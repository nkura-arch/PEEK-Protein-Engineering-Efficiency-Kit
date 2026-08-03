import os
import logging
import pandas as pd
from flask import Flask, request, render_template, jsonify, flash, redirect, url_for

logging.basicConfig(level=logging.DEBUG)

# Import the working algorithm
from silent_mutation_engine import (
    find_silent_sites, 
    filter_enzymes_not_in_sequence, 
    rank_results_by_preference,
    save_results_to_csv,
    get_preferred_enzymes,
    score_results_with_codon_optimization,
    export_scored_results_to_csv,
    get_top_codon_optimized_variants
)
from motif_library import get_motif_library
from Bio.Seq import Seq
from Bio.Restriction import CommOnly
from Bio.SeqUtils import MeltingTemp as mt
import math
import statistics

app = Flask(__name__)
app.config["TEMPLATES_AUTO_RELOAD"] = True
app.secret_key = os.environ.get("SESSION_SECRET", "dev-secret-key-for-peek")

# Let the browser cache static assets aggressively so switching tabs doesn't
# re-fetch CSS/JS every time. We bust the cache with a ?v=<mtime> query whenever
# a file changes (see asset_version below), so edits still show up immediately.
app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0  # no cache during dev


@app.context_processor
def inject_asset_version():
    """Expose asset_version(path) to templates for cache-busting static URLs."""
    def asset_version(filename):
        try:
            full = os.path.join(app.static_folder, filename)
            return str(int(os.path.getmtime(full)))
        except OSError:
            return "0"
    return {"asset_version": asset_version}

# ---------- CONFIG ----------
NEB_CORR = 18.0          # Empirical offset to match NEB
FORWARD_CORR = 9.0       # Additional correction for forward primers (non-SDM)
MIN_PRIMER_LEN = 15
MAX_COMBINED = 60
MAX_TM_DIFF = 5.0
NEB_ANNEAL_OFFSET = 3    # NEB-style annealing = avg_tm + 3

Q5_NA = 50               # mM monovalent
Q5_MG = 2                # mM Mg++
Q5_dNTP = 0.2            # mM dNTP
PRIMER_CONC = 500e-9     # 500 nM primer

# ---------- TM CONVERSION (PEEK <-> NEB) ----------
# Trendline: NEB_Tm = 0.585 * PEEK_Tm + 34.147
TM_SLOPE = 0.585
TM_INTERCEPT = 34.147

# NEB Tm targets (what user sees/wants)
NEB_TM_TARGET = 65.0     # Target NEB Tm (center of 60-70 range)
NEB_TM_MIN_SOFT = 60.0   # Soft minimum for optimal range
NEB_TM_MAX_SOFT = 70.0   # Soft maximum for optimal range  
NEB_TM_MIN_HARD = 55.0   # Hard minimum border
NEB_TM_MAX_HARD = 72.0   # Hard maximum border

def peek_to_neb_tm(peek_tm):
    """Convert PEEK Tm to NEB Tm using trendline equation"""
    return TM_SLOPE * peek_tm + TM_INTERCEPT

def neb_to_peek_tm(neb_tm):
    """Convert NEB Tm to PEEK Tm (inverse of trendline)"""
    return (neb_tm - TM_INTERCEPT) / TM_SLOPE

# Pre-calculated PEEK Tm targets from NEB targets
PEEK_TM_TARGET = neb_to_peek_tm(NEB_TM_TARGET)      # ~52.7°C
PEEK_TM_MIN_SOFT = neb_to_peek_tm(NEB_TM_MIN_SOFT)  # ~44.2°C
PEEK_TM_MAX_SOFT = neb_to_peek_tm(NEB_TM_MAX_SOFT)  # ~61.3°C
PEEK_TM_MIN_HARD = neb_to_peek_tm(NEB_TM_MIN_HARD)  # ~35.6°C
PEEK_TM_MAX_HARD = neb_to_peek_tm(NEB_TM_MAX_HARD)  # ~64.7°C

TARGET_TM = PEEK_TM_TARGET  # Use converted target for internal calculations


# ---------- UTILITIES ----------
def reverse_complement(seq):
    """Get reverse complement of DNA sequence"""
    return str(Seq(seq).reverse_complement())


def complement(seq):
    """Get complement of DNA sequence (no reversal)"""
    return str(Seq(seq).complement())


def validate_dna_sequence(sequence):
    """Validate DNA sequence input"""
    if not sequence:
        return False, "Please enter a DNA sequence"
    
    valid_bases = set('ATCG')
    sequence_bases = set(sequence.upper())
    
    if not sequence_bases.issubset(valid_bases):
        invalid_bases = sequence_bases - valid_bases
        return False, f"Invalid characters found: {', '.join(invalid_bases)}. Only A, T, C, G are allowed."
    
    return True, "Valid DNA sequence"


def tm_nn_peek(primer):
    """Calculate raw PEEK Tm (internal use for optimization)"""
    if not primer:
        return 0.0, 0.0
    try:
        tm_raw = mt.Tm_NN(
            primer,
            nn_table=mt.DNA_NN4,
            Na=Q5_NA,
            Mg=Q5_MG,
            dNTPs=Q5_dNTP,
            dnac1=PRIMER_CONC,
            dnac2=PRIMER_CONC
        )
    except Exception:
        tm_raw = 0.0
    
    tm_peek = tm_raw + NEB_CORR
    gc = 100.0 * (primer.count("G") + primer.count("C")) / len(primer)
    return round(tm_peek, 1), round(gc, 1)


def tm_nn_neb(primer):
    """Calculate NEB Tm using PEEK-to-NEB conversion trendline.
    
    Returns the Tm that should match NEB's calculator output.
    Uses: NEB_Tm = 0.585 * PEEK_Tm + 34.147
    """
    if not primer:
        return 0.0, 0.0
    
    peek_tm, gc = tm_nn_peek(primer)
    neb_tm = peek_to_neb_tm(peek_tm)
    return round(neb_tm, 1), round(gc, 1)


def neb_q5_anneal_temp(tm_f, tm_r):
    """NEB rule for Q5: Ta = lower Tm + 1°C, rounded up to whole number, max 72°C"""
    anneal = int(math.ceil(min(tm_f, tm_r) + 1.0))
    return min(anneal, 72)


def neb_q5_annealing_temp(tm_fwd, tm_rev):
    """Alias for backwards compatibility"""
    return neb_q5_anneal_temp(tm_fwd, tm_rev)


def calculate_tm(seq):
    """Calculate Tm using the new NEB-corrected method"""
    return tm_nn_neb(seq)[0]


def calculate_gc_content(seq):
    """Calculate GC content percentage"""
    gc_count = seq.count('G') + seq.count('C')
    return round((gc_count / len(seq)) * 100, 1) if seq else 0


def design_primer_pair(avail_forward, avail_reverse_rc, site_prefix=""):
    """
    Primer pair optimizer targeting NEB Tm 60-70°C with 55-72°C hard borders.
    
    Prioritizes:
    1. Both primers within optimal 60-70°C range
    2. Small Tm difference between primers
    3. Average Tm close to 65°C target
    4. Shorter combined length
    """
    best_optimal = None  # Both primers in 60-70°C
    best_acceptable = None  # Both primers in 55-72°C
    max_f = min(len(avail_forward)+len(site_prefix), MAX_COMBINED - MIN_PRIMER_LEN)
    
    for f_len in range(MIN_PRIMER_LEN, max_f+1):
        f_take = max(0, f_len - len(site_prefix))
        f_primer = (site_prefix + avail_forward[:f_take])[:f_len]
        f_tm, f_gc = tm_nn_neb(f_primer)
        
        # Skip if forward primer outside hard borders
        if f_tm < NEB_TM_MIN_HARD or f_tm > NEB_TM_MAX_HARD:
            continue
        
        max_r = min(len(avail_reverse_rc), MAX_COMBINED - len(f_primer))
        if max_r < MIN_PRIMER_LEN:
            continue
            
        for r_len in range(MIN_PRIMER_LEN, max_r+1):
            r_primer = avail_reverse_rc[:r_len]
            r_tm, r_gc = tm_nn_neb(r_primer)
            
            # Skip if reverse primer outside hard borders
            if r_tm < NEB_TM_MIN_HARD or r_tm > NEB_TM_MAX_HARD:
                continue
            
            tm_diff = abs(f_tm - r_tm)
            if tm_diff > MAX_TM_DIFF:
                continue
                
            avg_tm = (f_tm + r_tm) / 2
            combined_len = len(f_primer) + len(r_primer)
            
            # Score: prioritize small Tm diff, close to target, shorter length
            score = (tm_diff, abs(avg_tm - NEB_TM_TARGET), combined_len)
            
            # Check if both primers in optimal range
            f_in_optimal = NEB_TM_MIN_SOFT <= f_tm <= NEB_TM_MAX_SOFT
            r_in_optimal = NEB_TM_MIN_SOFT <= r_tm <= NEB_TM_MAX_SOFT
            
            result = {
                "f_primer": f_primer, "f_tm": f_tm, "f_gc": f_gc,
                "r_primer": r_primer, "r_tm": r_tm, "r_gc": r_gc,
                "avg_tm": round(avg_tm, 1), "combined_len": combined_len,
                "tm_diff": round(tm_diff, 1),
                "recommended_anneal": neb_q5_anneal_temp(f_tm, r_tm),
                "score": score
            }
            
            if f_in_optimal and r_in_optimal:
                if best_optimal is None or score < best_optimal["score"]:
                    best_optimal = result
            else:
                if best_acceptable is None or score < best_acceptable["score"]:
                    best_acceptable = result
    
    # Prefer optimal range, fall back to acceptable
    return best_optimal if best_optimal else best_acceptable


def design_optimal_primer_pair(avail_forward, avail_reverse_rc,
                               min_len=MIN_PRIMER_LEN, max_total=MAX_COMBINED,
                               target_tm=TARGET_TM, max_diff=MAX_TM_DIFF,
                               enzyme_site_prefix=""):
    """Wrapper for compatibility with existing code"""
    return design_primer_pair(avail_forward, avail_reverse_rc, site_prefix=enzyme_site_prefix)


def get_enzyme_site_dict(dna_sequence=None):
    """Get dictionary of enzyme names to recognition sites, optionally filtered by DNA sequence"""
    enzyme_dict = {}
    for enzyme in CommOnly:
        site = str(enzyme.site)
        if 'N' not in site and len(site) <= 8:  # Filter out ambiguous and very long sites
            # If DNA sequence provided, only include enzymes whose sites don't already exist
            if dna_sequence is None or site not in dna_sequence:
                enzyme_dict[enzyme.__name__] = site
    return enzyme_dict


def aa_to_dna(aa_seq):
    """Convert amino acid sequence to DNA using simplified codon table"""
    table = {'A':'GCT','R':'CGT','N':'AAT','D':'GAT','C':'TGT',
             'Q':'CAA','E':'GAA','G':'GGT','H':'CAT','I':'ATT',
             'L':'TTA','K':'AAA','M':'ATG','F':'TTT','P':'CCT',
             'S':'TCT','T':'ACT','W':'TGG','Y':'TAT','V':'GTT','*':'TAA'}
    dna = ""
    for aa in aa_seq:
        if aa not in table:
            raise ValueError(f"Invalid AA: {aa}")
        dna += table[aa]
    return dna


def extend_primer_to_tm(seq, min_tm=None, max_tm=None, target_tm=None, min_len=18):
    """
    Extend primer to optimal NEB Tm range (60-70°C target, 55-72°C hard borders).

    Prioritizes finding a primer length that gives Tm in the optimal 60-70°C range,
    but will accept primers within the hard borders if optimal range isn't achievable.
    min_len: minimum primer length to consider (default 18).
    """
    if min_tm is None:
        min_tm = NEB_TM_MIN_HARD  # 55°C
    if max_tm is None:
        max_tm = NEB_TM_MAX_HARD  # 72°C
    if target_tm is None:
        target_tm = NEB_TM_TARGET  # 65°C

    primer_len = max(18, min_len)
    max_len = 40
    
    best_optimal = None  # Best primer in 60-70°C range
    best_acceptable = None  # Best primer in 55-72°C range
    
    while primer_len <= min(max_len, len(seq)):
        candidate = seq[:primer_len]
        tm, gc = tm_nn_neb(candidate)
        
        # Check if within hard borders
        if min_tm <= tm <= max_tm:
            score = abs(tm - target_tm)
            
            # Check if within optimal range (60-70°C)
            if NEB_TM_MIN_SOFT <= tm <= NEB_TM_MAX_SOFT:
                if best_optimal is None or score < best_optimal[3]:
                    best_optimal = (candidate, tm, gc, score)
            
            # Track best acceptable regardless
            if best_acceptable is None or score < best_acceptable[3]:
                best_acceptable = (candidate, tm, gc, score)
        
        primer_len += 1
    
    # Prefer optimal range, fall back to acceptable
    if best_optimal:
        return best_optimal[0], best_optimal[1], best_optimal[2]
    if best_acceptable:
        return best_acceptable[0], best_acceptable[1], best_acceptable[2]
    
    # Fallback: return max length primer if nothing in range
    tm, gc = tm_nn_neb(seq[:min(max_len, len(seq))])
    return seq[:min(max_len, len(seq))], tm, gc


def design_enzyme_insertion_primers(dna, position, enzyme_site, min_tm=55, max_tm=72):
    """
    Design primers for inserting an enzyme site at a cut site.
    
    Forward primer: 5' enzyme tail + top strand starting at cut site, extend right for Tm.
    Reverse primer: bottom strand starting at cut site, extend left (5'→3') for Tm.
    """
    pos = position - 1  # 0-based index
    arm_len = 20

    # Homology arms for reference
    left_arm = dna[max(0, pos - arm_len):pos]
    right_arm = dna[pos:pos + arm_len]

    # Forward primer: enzyme tail + top strand sequence right of cut site.
    # Enforce minimum 15 nt genomic binding tail (NEB Q5 requires adequate template contact).
    # For high-GC enzyme sites (like AscI/SrfI) the combined Tm may exceed the 72°C hard cap;
    # in that case limit the binding tail to 15 nt and flag it with f_binding_tm for display.
    site_len = len(enzyme_site)
    min_binding = 15
    f_seq = enzyme_site + dna[pos:pos + 40]
    f_primer, f_tm, f_gc = extend_primer_to_tm(f_seq, min_tm, max_tm,
                                                 min_len=site_len + min_binding)
    # If the primer still exceeds max_tm (high-GC site pushed it over), cap the binding tail.
    if f_tm > max_tm:
        f_primer = f_seq[:site_len + min_binding]
        f_tm, f_gc = tm_nn_neb(f_primer)
    f_binding_tail = f_primer[site_len:]
    f_binding_tm, _ = tm_nn_neb(f_binding_tail) if len(f_binding_tail) >= 8 else (None, None)

    # Reverse primer: bottom strand sequence left of cut site, 5'→3'
    bottom_left_seq = dna[max(0, pos - 40):pos]
    r_primer_seq = reverse_complement(bottom_left_seq)  # now 5'→3'
    r_primer, r_tm, r_gc = extend_primer_to_tm(r_primer_seq, min_tm, max_tm)

    avg_tm = round((f_tm + r_tm) / 2, 1)
    tm_diff = round(abs(f_tm - r_tm), 1)
    anneal_temp = neb_q5_annealing_temp(f_tm, r_tm)

    return {
        'left_arm': left_arm,
        'enzyme_site': enzyme_site,
        'insert_seq': enzyme_site,
        'right_arm': right_arm,
        'f_primer': f_primer,
        'r_primer': r_primer,
        'f_tm': round(f_tm, 1),
        'r_tm': round(r_tm, 1),
        'f_gc': round(f_gc, 1),
        'r_gc': round(r_gc, 1),
        'f_binding_tail': f_binding_tail,
        'f_binding_tm': round(f_binding_tm, 1) if f_binding_tm is not None else None,
        'avg_tm': avg_tm,
        'tm_diff': tm_diff,
        'annealing_temp_Q5': anneal_temp,
        'recommended_anneal': anneal_temp,
        'cut_site_position': position
    }


def design_enzyme_primers(dna, position, enzyme_name, enzyme_site):
    """Design primers for enzyme insertion"""
    result = design_enzyme_insertion_primers(dna, position, enzyme_site)
    result['enzyme_name'] = enzyme_name
    result['insert_nt'] = enzyme_site
    return result


def design_insertion_primers(dna, position, insert_nt, arm_len=20, max_primer_len=25):
    """
    Design primers that bind to the homology arms flanking an insertion site.
    No insert included in the primers - they bind to the homology arms only.

    insert_nt is the exact nucleotide sequence to insert. Callers that have an
    amino acid sequence instead should convert it with aa_to_dna() first.
    """
    pos = position - 1

    left_arm = dna[max(0, pos - arm_len):pos]
    right_arm = dna[pos:min(len(dna), pos + arm_len)]
    homology_full = left_arm + insert_nt + right_arm

    f_primer_seq = left_arm[-max_primer_len:]
    tm_f, gc_f = tm_nn_neb(f_primer_seq)

    r_primer_seq = reverse_complement(right_arm[:max_primer_len])
    tm_r, gc_r = tm_nn_neb(r_primer_seq)

    avg_tm = round((tm_f + tm_r) / 2, 1)
    tm_diff = round(abs(tm_f - tm_r), 1)
    anneal_temp = neb_q5_annealing_temp(tm_f, tm_r)

    return {
        'left_arm': left_arm,
        'right_arm': right_arm,
        'insert_nt': insert_nt,
        'insert_seq': insert_nt,
        'homology_full': homology_full,
        'f_primer': f_primer_seq,
        'r_primer': r_primer_seq,
        'f_tm': tm_f,
        'r_tm': tm_r,
        'f_gc': gc_f,
        'r_gc': gc_r,
        'avg_tm': avg_tm,
        'tm_diff': tm_diff,
        'annealing_temp_Q5': anneal_temp,
        'recommended_anneal': anneal_temp,
        'insertion_point': f"Position {position} in sequence"
    }


def design_indel_primers(dna, insert_seq, operation='amplify', tolerance_mode='edge', tolerance_bp=30):
    """
    Design PCR primers for an indel operation around a target region.

    operation='amplify': primers point INWARD, flanking (edge mode) or
      overlapping (tolerant mode) the target region -- for amplifying it.
        Forward: anchored at the region's 5' boundary (or `tolerance_bp`
          upstream of it in tolerant mode), extends rightward INTO the region.
        Reverse: anchored at the region's 3' boundary (or `tolerance_bp'
          downstream of it in tolerant mode), extends leftward INTO the
          region on the bottom strand.

    operation='delete': primers point OUTWARD, back-to-back around the
      target region -- for cutting it out (inverse-PCR / deletion primers).
        Forward: anchored just after the region's 3' boundary (or further
          downstream in tolerant mode), extends rightward AWAY from the
          region into the remaining sequence.
        Reverse: anchored just before the region's 5' boundary (or further
          upstream in tolerant mode), extends leftward AWAY from the region
          into the remaining sequence, on the bottom strand.

    Tolerance gives leniency for edge/junction confirmation in amplify mode,
    or a clean margin removed around the target region in delete mode.
    """
    dna = dna.upper()
    insert_seq = insert_seq.upper()

    insert_start = dna.find(insert_seq)
    if insert_start == -1:
        raise ValueError("Target region not found in the provided DNA sequence")

    insert_len = len(insert_seq)
    insert_end = insert_start + insert_len  # 0-based, exclusive

    is_tolerant = tolerance_mode == 'tolerant'
    is_delete = operation == 'delete'
    tol = tolerance_bp if is_tolerant else 0

    if is_delete:
        fwd_anchor = min(len(dna), insert_end + tol)
        rev_anchor = max(0, insert_start - tol)
    else:
        fwd_anchor = max(0, insert_start - tol)
        rev_anchor = min(len(dna), insert_end + tol)

    # ---------- FORWARD PRIMER ----------
    f_seq = dna[fwd_anchor:fwd_anchor + 40]
    if not f_seq:
        raise ValueError("Not enough sequence to design a forward primer at this boundary -- try a smaller tolerance or provide more flanking sequence.")
    f_primer, f_tm, f_gc = extend_primer_to_tm(f_seq, min_len=18)
    fwd_start_0 = fwd_anchor
    fwd_end_0 = fwd_anchor + len(f_primer)

    # ---------- REVERSE PRIMER ----------
    r_window_start = max(0, rev_anchor - 40)
    r_template = dna[r_window_start:rev_anchor]
    if not r_template:
        raise ValueError("Not enough sequence to design a reverse primer at this boundary -- try a smaller tolerance or provide more flanking sequence.")
    r_seq = reverse_complement(r_template)  # 5'->3' on bottom strand, starting at rev_anchor
    r_primer, r_tm, r_gc = extend_primer_to_tm(r_seq, min_len=18)
    rev_end_0 = rev_anchor
    rev_start_0 = rev_anchor - len(r_primer)

    avg_tm = round((f_tm + r_tm) / 2, 1)
    tm_diff = round(abs(f_tm - r_tm), 1)
    anneal_temp = neb_q5_annealing_temp(f_tm, r_tm)

    amplicon_len = None
    deletion_size = None
    remaining_len = None
    if is_delete:
        deletion_size = fwd_anchor - rev_anchor
        remaining_len = len(dna) - deletion_size
    else:
        amplicon_len = rev_anchor - fwd_anchor

    return {
        'operation': operation,
        'insert_seq': insert_seq,
        'insert_len': insert_len,
        'insert_start': insert_start + 1,   # 1-based for display
        'insert_end': insert_end,           # 1-based inclusive == 0-based exclusive
        'insert_start_0': insert_start,
        'insert_end_0': insert_end,
        'tolerance_mode': tolerance_mode,
        'tolerance_bp': tol,
        'f_primer': f_primer,
        'r_primer': r_primer,
        'f_tm': f_tm,
        'r_tm': r_tm,
        'f_gc': f_gc,
        'r_gc': r_gc,
        'avg_tm': avg_tm,
        'tm_diff': tm_diff,
        'annealing_temp_Q5': anneal_temp,
        'recommended_anneal': anneal_temp,
        'amplicon_len': amplicon_len,
        'deletion_size': deletion_size,
        'remaining_len': remaining_len,
        'fwd_anchor': fwd_anchor + 1,        # 1-based for display
        'rev_anchor': rev_anchor,            # 1-based inclusive display == 0-based exclusive
        'fwd_start_0': fwd_start_0,
        'fwd_end_0': fwd_end_0,
        'rev_start_0': rev_start_0,
        'rev_end_0': rev_end_0,
    }


def annealing_tm(seq1, seq2=None, salt_conc=50):
    """
    Calculate the annealing Tm for two complementary DNA strands.
    
    Parameters:
    seq1       : Top strand sequence (5'->3')
    seq2       : Bottom strand sequence (5'->3'). If None, assumes perfect complement of seq1
    salt_conc  : Monovalent salt concentration in mM (default 50 mM)
    
    Returns:
    Tm in Celsius
    """
    import math
    
    if seq2 is None:
        seq2 = reverse_complement(seq1)
    
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be same length to anneal")
    
    A = seq1.count('A') + seq2.count('A')
    T = seq1.count('T') + seq2.count('T')
    G = seq1.count('G') + seq2.count('G')
    C = seq1.count('C') + seq2.count('C')
    N = len(seq1)
    
    tm = 64.9 + 41 * (G + C - 16.4)/N
    tm += 16.6 * math.log10(salt_conc / 1000)
    
    return round(tm, 1)


def design_enzyme_and_insert(dna, position, enzyme_name, enzyme_site, insert_seq, cut_offset=1, overhang_len=4):
    """
    Design SDM primers to insert a restriction enzyme site and create a double-stranded insert
    for ligation of a user-defined sequence.
    
    Parameters:
    dna         : Original plasmid sequence
    position    : 1-based insertion position for the enzyme site
    enzyme_name : Name of enzyme (for tracking)
    enzyme_site : Recognition sequence of the enzyme
    insert_seq  : Nucleotide sequence to insert via ligation
    cut_offset  : Offset of cut from 5' end of enzyme site (for sticky ends)
    overhang_len: Number of bases in sticky end
    
    Returns:
    Dictionary with SDM primers info and double-stranded insert sequence.
    """
    
    # --- SDM Primer Design ---
    pos = position - 1  # convert to 0-based
    insert_len = len(enzyme_site)
    flank_len = max(10, (30 - insert_len) // 2)
    
    start = max(0, pos - flank_len)
    end = min(len(dna), pos + flank_len)
    left_flank = dna[start:pos]
    right_flank = dna[pos:end]
    
    f_primer = left_flank + enzyme_site + right_flank
    r_primer = reverse_complement(f_primer)
    
    tm_f, gc_f = tm_nn_neb(f_primer)
    tm_r, gc_r = tm_nn_neb(r_primer)
    avg_tm = (tm_f + tm_r) / 2
    recommended_anneal = int(round(avg_tm - 3))
    
    # --- Homology arms (optional for SDM synthesis) ---
    arm_len = 20
    left_arm = dna[max(0, pos - arm_len):pos]
    right_arm = dna[pos:min(len(dna), pos + arm_len)]
    homology_full = left_arm + enzyme_site + right_arm
    
    # --- Double-stranded insert for ligation ---
    # Overhangs determined from cut_offset and overhang_len
    top_overhang = insert_seq[:overhang_len] if len(insert_seq) >= overhang_len else insert_seq
    bottom_overhang = reverse_complement(top_overhang)
    
    top_strand = top_overhang + insert_seq + bottom_overhang
    bottom_strand = reverse_complement(top_strand)
    
    # Calculate annealing Tm for ligation oligos
    oligo_anneal_tm = annealing_tm(top_strand, bottom_strand)
    
    return {
        # SDM primers
        'enzyme_name': enzyme_name,
        'enzyme_site': enzyme_site,
        'insert_nt': enzyme_site,
        'left_arm': left_arm,
        'right_arm': right_arm,
        'homology_full': homology_full,
        'f_primer': f_primer,
        'r_primer': r_primer,
        'f_tm': round(tm_f, 1),
        'r_tm': round(tm_r, 1),
        'f_gc': gc_f,
        'r_gc': gc_r,
        'avg_tm': round(avg_tm, 1),
        'recommended_anneal': recommended_anneal,
        'insertion_point': f"Position {position} in sequence",
        # Ligation insert
        'ligation_insert_seq': insert_seq,
        'insert_top': top_strand,
        'insert_bottom': bottom_strand,
        'overhang_len': overhang_len,
        'oligo_anneal_tm': oligo_anneal_tm
    }


def design_sdm_primers(plasmid_sequence, original_window, mutated_window, position=None, enzyme_name=None):
    """
    Benchling-correct Q5 SDM primer design:
    - Forward: top strand, rightward
    - Reverse: bottom strand, leftward
    - Reverse is NOT reverse-complement of forward
    
    The forward primer uses the user's mutated_window directly as the core,
    with flanking sequences taken from the original plasmid.
    """
    dna = plasmid_sequence.upper()
    orig = original_window.upper().replace(" ", "")
    mut = mutated_window.upper().replace(" ", "")

    if len(orig) != len(mut):
        raise ValueError("Original and mutated windows must be the same length")

    window_start = dna.find(orig)
    if window_start == -1:
        raise ValueError("Original window not found in plasmid")
    
    window_end = window_start + len(orig)

    # The silent-mutation window (orig/mut) is itself always codon-aligned to
    # the gene's real reading frame (see get_codon_aligned_window), so the
    # frame is simply window_start's residue mod 3 -- no need to redetect it.
    frame_offset = window_start % 3

    # Identify mutation span for reporting
    diff_indices = [i for i in range(len(orig)) if orig[i] != mut[i]]
    if not diff_indices:
        raise ValueError("No mutation detected")

    first_mut = diff_indices[0]
    last_mut = diff_indices[-1] + 1
    mut_start_abs = window_start + first_mut
    mut_end_abs = window_start + last_mut
    orig_segment = orig[first_mut:last_mut]
    mut_segment = mut[first_mut:last_mut]

    # ---------- FORWARD PRIMER ----------
    # Build primer: [left extension from plasmid] + [mutated window core] + [right extension from plasmid]
    # Can trim the mutated window OR extend beyond it using plasmid flanking sequence

    fwd_primer = None
    tm_fwd = gc_fwd = None
    best_score = None
    best_start = None

    min_flank = 5  # Minimum bp on each side of the mutation
    max_extend = 15  # Maximum bp to extend beyond the user's window

    # Try combinations: trim from window OR extend beyond window
    for left_trim in range(-max_extend, len(mut) - min_flank):
        for right_trim in range(-max_extend, len(mut) - max(0, left_trim) - min_flank):
            # Calculate plasmid positions
            primer_start_in_plasmid = window_start + left_trim
            primer_end_in_plasmid = window_end - right_trim

            # Bounds check
            if primer_start_in_plasmid < 0 or primer_end_in_plasmid > len(dna):
                continue

            # Codon-boundary check: both ends of the primer must fall exactly
            # at a codon boundary in the gene's reading frame, never mid-codon.
            # Otherwise the primer's 5'/3' edge could land inside a mutated
            # codon, splitting an amino acid across the primer edge.
            if (primer_start_in_plasmid - frame_offset) % 3 != 0:
                continue
            if (primer_end_in_plasmid - frame_offset) % 3 != 0:
                continue


            # Build the primer:
            # - Left extension from plasmid (if left_trim < 0)
            # - Core from mutated window
            # - Right extension from plasmid (if right_trim < 0)
            
            left_ext = ""
            if left_trim < 0:
                left_ext = dna[primer_start_in_plasmid:window_start]
                core_left = 0
            else:
                core_left = left_trim
            
            right_ext = ""
            if right_trim < 0:
                right_ext = dna[window_end:primer_end_in_plasmid]
                core_right = len(mut)
            else:
                core_right = len(mut) - right_trim
            
            core = mut[core_left:core_right]
            candidate = left_ext + core + right_ext
            
            primer_len = len(candidate)
            if primer_len < 15 or primer_len > 30:
                continue
            
            # Verify the mutation is still included with minimum flanks
            # Mutation position relative to primer start
            mut_pos_in_primer = first_mut - left_trim if left_trim >= 0 else first_mut + len(left_ext)
            right_flank_len = len(candidate) - mut_pos_in_primer - len(mut_segment)
            
            if mut_pos_in_primer < min_flank or right_flank_len < min_flank:
                continue
            
            tm, gc = tm_nn_neb(candidate)
            
            # Score: prefer optimal range, then closeness to 65°C target
            in_optimal = NEB_TM_MIN_SOFT <= tm <= NEB_TM_MAX_SOFT
            in_acceptable = NEB_TM_MIN_HARD <= tm <= NEB_TM_MAX_HARD
            
            if not in_acceptable:
                continue
            
            # Score: (not in optimal range, distance from target, length)
            score = (0 if in_optimal else 1, abs(tm - NEB_TM_TARGET), primer_len)
            
            if best_score is None or score < best_score:
                best_score = score
                fwd_primer = candidate
                tm_fwd = tm
                gc_fwd = gc
                best_start = primer_start_in_plasmid
    
    if fwd_primer is None:
        raise ValueError("Could not design suitable forward primer - Tm outside 55-72°C range")
    
    start = best_start  # Use for reverse primer calculation

    # ---------- REVERSE PRIMER (Benchling-correct) ----------
    # Same genomic start as forward (3' end), goes LEFT on bottom strand.
    # The reverse primer never overlaps the mutated codons -- it sits
    # entirely in the untouched flanking plasmid sequence -- so its length
    # is independently searched for the best Tm, in whole-codon (3 bp) steps
    # the same way the forward primer is, and lengthened or shortened as
    # needed to land in range. It previously just copied the forward
    # primer's length with no Tm check at all, so it could end up outside
    # the acceptable 55-72 C window even when the forward primer was fine.
    rev_primer = None
    tm_rev = gc_rev = None
    best_rev_score = None

    min_rev_len = 15   # multiple of 3
    max_rev_len = 30   # multiple of 3

    for rev_len in range(min_rev_len, max_rev_len + 1, 3):
        if start - rev_len < 0:
            continue
        rev_template_top = dna[start - rev_len:start]
        candidate_rev = complement(rev_template_top)[::-1]
        tm_r, gc_r = tm_nn_neb(candidate_rev)

        r_in_optimal = NEB_TM_MIN_SOFT <= tm_r <= NEB_TM_MAX_SOFT
        r_in_acceptable = NEB_TM_MIN_HARD <= tm_r <= NEB_TM_MAX_HARD
        if not r_in_acceptable:
            continue

        # Prefer the optimal range, then closeness to target, then
        # closeness to the forward primer's Tm (keeps annealing well matched).
        r_score = (0 if r_in_optimal else 1, abs(tm_r - NEB_TM_TARGET), abs(tm_r - tm_fwd))
        if best_rev_score is None or r_score < best_rev_score:
            best_rev_score = r_score
            rev_primer = candidate_rev
            tm_rev = tm_r
            gc_rev = gc_r

    if rev_primer is None:
        raise ValueError("Could not design suitable reverse primer - Tm outside 55-72°C range")

    # ---------- NEB Q5 ANNEALING ----------
    anneal_temp = neb_q5_annealing_temp(tm_fwd, tm_rev)
    avg_tm = round((tm_fwd + tm_rev) / 2, 1)
    tm_diff = round(abs(tm_fwd - tm_rev), 1)

    label = f"SDM for {enzyme_name}" if enzyme_name else "Silent Mutation (Q5)"

    return {
        'label': label,
        'enzyme_name': label,
        'forward_primer_5to3': fwd_primer,
        'reverse_primer_5to3': rev_primer,
        'f_primer': fwd_primer,
        'r_primer': rev_primer,
        'f_tm': tm_fwd,
        'r_tm': tm_rev,
        'f_gc': gc_fwd,
        'r_gc': gc_rev,
        'avg_tm': avg_tm,
        'tm_diff': tm_diff,
        'annealing_temp_Q5': anneal_temp,
        'recommended_anneal': anneal_temp,
        'mutation_region': f"Position {mut_start_abs + 1}–{mut_end_abs}: {orig_segment} → {mut_segment}",
        'mutation_start': mut_start_abs + 1,
        'mutation_end': mut_end_abs,
        'original_seq': orig_segment,
        'mutated_seq': mut_segment,
        'primer_length': len(fwd_primer),
        'tm_difference': tm_diff,
        'fwd_primer_start_0': best_start,
    }


# ---------- SEQUENCE VISUALIZER SIDEBAR DATA ----------
# Curated, click-to-insert items shown in the Sequence Visualizer sidebar.
_VIZ_MOTIFS = get_motif_library()
_VIZ_ENZYMES_CACHE = None


def get_common_enzyme_sites():
    """Sorted (name, site) list of common restriction enzymes for the sidebar (cached).

    Only unambiguous ATCG-only recognition sites are included, since clicking an
    item inserts the site as literal DNA (IUPAC ambiguity codes like W/N would be
    rejected by DNA validation).
    """
    global _VIZ_ENZYMES_CACHE
    if _VIZ_ENZYMES_CACHE is None:
        site_dict = get_enzyme_site_dict(None)
        unambiguous = {
            name: site for name, site in site_dict.items()
            if set(site.upper()) <= {'A', 'T', 'C', 'G'}
        }
        _VIZ_ENZYMES_CACHE = sorted(unambiguous.items(), key=lambda kv: kv[0])
    return _VIZ_ENZYMES_CACHE


@app.context_processor
def inject_viz_sidebar_data():
    """Make the motif + enzyme sidebar lists available to all templates."""
    return {
        'viz_motifs': _VIZ_MOTIFS,
        'viz_enzymes': get_common_enzyme_sites(),
    }


@app.route('/')
def index():
    return redirect('/silent-mutations')

@app.route('/home')
def home():
    return redirect('/silent-mutations')


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/silent-mutations', methods=['GET', 'POST'])
def silent_mutations():
    if request.method == 'GET':
        return render_template('silent_mutations.html')

    try:
        # Get form data
        dna = request.form.get('dna', '').strip().upper().replace('\n', '').replace(' ', '')
        cut_pos = int(request.form.get('cut_pos', 0))
        tolerance = int(request.form.get('tolerance', 3))
        window_codons = int(request.form.get('window_codons', 10))
        max_changes = int(request.form.get('max_changes', 3))
        host_organism = request.form.get('host_organism', 'ecoli')
        # An unchecked checkbox is omitted from the POST, so missing == disabled.
        enable_optimization = request.form.get('enable_optimization', '') == 'on'

        # Validate inputs
        is_valid, message = validate_dna_sequence(dna)
        if not is_valid:
            flash(f"Invalid DNA sequence: {message}", "error")
            return render_template('silent_mutations.html', error=message,
                                 dna=dna, cut_pos=cut_pos, tolerance=tolerance,
                                 window_codons=window_codons, max_changes=max_changes)

        if cut_pos <= 0 or cut_pos > len(dna):
            flash(f"Cut position ({cut_pos}) must be between 1 and {len(dna)}", "error")
            return render_template('silent_mutations.html', error="Invalid cut position",
                                 dna=dna, cut_pos=cut_pos, tolerance=tolerance,
                                 window_codons=window_codons, max_changes=max_changes)

        # Convert to 0-based indexing
        cut_pos_0based = cut_pos - 1

        # Filter enzymes not in sequence
        available_enzymes = filter_enzymes_not_in_sequence(dna)

        # Find silent sites. coding_start_pos is left unset so the reading
        # frame at cut_pos is auto-detected -- this works whether `dna` is
        # just the CDS or a whole plasmid pasted with cut_pos in absolute
        # coordinates.
        results = find_silent_sites(
            dna,
            cut_pos_0based,
            tolerance,
            window_codons,
            available_enzymes,
            max_changes=max_changes
        )

        if not results:
            flash("No valid silent mutations found with the given parameters.", "warning")
            return render_template('silent_mutations.html',
                                 error="No results found. Try increasing tolerance or window size.",
                                 dna=dna, cut_pos=cut_pos, tolerance=tolerance,
                                 window_codons=window_codons, max_changes=max_changes)

        # Score or rank results
        if enable_optimization:
            try:
                results_df = score_results_with_codon_optimization(results, dna, host_organism=host_organism)
                ranked_results = results_df.to_dict('records')
                top_variants = get_top_codon_optimized_variants(results_df, top_n=10)
            except Exception as opt_error:
                print(f"Codon optimization error: {opt_error}")
                results_df, matched_preferred = rank_results_by_preference(results)
                ranked_results = results_df.to_dict('records')
                top_variants = ranked_results[:10]
        else:
            results_df, matched_preferred = rank_results_by_preference(results)
            ranked_results = results_df.to_dict('records')
            top_variants = ranked_results[:10]

        # Get summary stats
        summary = {
            'total_results': len(results),
            'optimization_enabled': enable_optimization,
            'host_organism': host_organism,
            'unique_enzymes': results_df['enzyme'].astype(str).nunique() if not results_df.empty else 0,
            'avg_changes': results_df['removed_bases'].str.len().mean() if not results_df.empty else 0,
        }

        if not enable_optimization and 'matched_preferred' in locals():
            summary['preferred_enzymes'] = len(matched_preferred)
            summary['matched_preferred_list'] = matched_preferred

        return render_template('silent_mutations.html',
                             results=ranked_results,
                             top_variants=top_variants,
                             summary=summary,
                             success=True,
                             dna=dna,
                             cut_pos=cut_pos,
                             tolerance=tolerance,
                             window_codons=window_codons,
                             max_changes=max_changes,
                             host_organism=host_organism,
                             enable_optimization=enable_optimization)

    except ValueError as e:
        flash(f"Invalid input: {str(e)}", "error")
        return render_template('silent_mutations.html', error="Invalid numerical input")
    except Exception as e:
        flash(f"An error occurred: {str(e)}", "error")
        return render_template('silent_mutations.html', error="Processing error occurred")
@app.route('/download-csv')
def download_csv():
    """Download results as CSV"""
    try:
        # This would need to be called after results are generated
        # For now, return a simple response
        return "CSV download feature - results need to be stored in session first"
    except Exception as e:
        flash(f"Error downloading CSV: {str(e)}", "error")
        return redirect(url_for('silent_mutations'))


@app.route('/design-sdm-primers', methods=['POST'])
def design_sdm_primers_route():
    """Route for SDM primer design from silent mutations page"""
    try:
        dna = request.form.get('dna', '').strip().upper().replace('\n', '').replace(' ', '')
        window_start_str = request.form.get('window_start', '0')
        window_start = int(window_start_str) if window_start_str else 0
        original_nt = request.form.get('original_nt', '')
        mutated_nt = request.form.get('mutated_nt', '')
        enzyme_name = request.form.get('enzyme_name', '').strip()
        
        # Preserve form state for re-rendering
        cut_pos = request.form.get('cut_pos', '')
        tolerance = request.form.get('tolerance', '30')
        window_codons = request.form.get('window_codons', '10')
        max_changes = request.form.get('max_changes', '3')
        host_organism = request.form.get('host_organism', 'ecoli')
        enable_optimization = request.form.get('enable_optimization', '') == 'on'
        
        # Clean up sequences if they exist
        if original_nt:
            original_nt = original_nt.strip().upper().replace(' ', '')
        if mutated_nt:
            mutated_nt = mutated_nt.strip().upper().replace(' ', '')
        
        # Validate inputs
        if not dna:
            flash("DNA sequence is required", "error")
            return redirect(url_for('silent_mutations'))
        
        if not original_nt or not mutated_nt:
            flash(f"Both original and mutated sequences are required. Got orig='{original_nt}' mut='{mutated_nt}'", "error")
            return redirect(url_for('silent_mutations'))
        
        # Check for None before len comparison
        orig_len = len(original_nt) if original_nt else 0
        mut_len = len(mutated_nt) if mutated_nt else 0
        if orig_len != mut_len:
            flash(f"Original ({orig_len}) and mutated ({mut_len}) sequences must be the same length", "error")
            return redirect(url_for('silent_mutations'))
        
        # Validate DNA sequences
        is_valid_orig, msg_orig = validate_dna_sequence(original_nt)
        is_valid_mut, msg_mut = validate_dna_sequence(mutated_nt)
        if not is_valid_orig or not is_valid_mut:
            flash(f"Invalid nucleotide sequence", "error")
            return redirect(url_for('silent_mutations'))
        
        # Design SDM primers
        logging.debug(f"Calling design_sdm_primers with dna len={len(dna)}, orig={original_nt[:20]}..., mut={mutated_nt[:20]}...")
        try:
            sdm_result = design_sdm_primers(dna, original_nt, mutated_nt, window_start, enzyme_name=enzyme_name if enzyme_name else None)
            logging.debug(f"SDM result: {sdm_result}")
        except Exception as sdm_error:
            logging.error(f"SDM primer design error: {sdm_error}")
            flash(f"SDM primer design error: {sdm_error}", "error")
            return redirect(url_for('silent_mutations'))
        
        if sdm_result is None:
            flash("Could not design SDM primers. Check that the original sequence matches the DNA at the specified position.", "error")
            return redirect(url_for('silent_mutations'))
        
        # Re-run silent mutation analysis to preserve results
        ranked_results = None
        top_variants = None
        summary = None
        
        try:
            cut_pos_int = int(cut_pos) if cut_pos else 0
            tolerance_int = int(tolerance) if tolerance else 30
            window_codons_int = int(window_codons) if window_codons else 10
            max_changes_int = int(max_changes) if max_changes else 3
            
            if cut_pos_int > 0 and cut_pos_int <= len(dna):
                cut_pos_0based = cut_pos_int - 1
                available_enzymes = filter_enzymes_not_in_sequence(dna)
                
                results = find_silent_sites(
                    dna,
                    cut_pos_0based,
                    tolerance_int,
                    window_codons_int,
                    available_enzymes,
                    max_changes=max_changes_int
                )
                
                if results:
                    if enable_optimization:
                        try:
                            results_df = score_results_with_codon_optimization(results, dna, host_organism=host_organism)
                            ranked_results = results_df.to_dict('records')
                            top_variants = get_top_codon_optimized_variants(results_df, top_n=10)
                        except:
                            results_df, matched_preferred = rank_results_by_preference(results)
                            ranked_results = results_df.to_dict('records')
                            top_variants = ranked_results[:10]
                    else:
                        results_df, matched_preferred = rank_results_by_preference(results)
                        ranked_results = results_df.to_dict('records')
                        top_variants = ranked_results[:10]
                    
                    summary = {
                        'total_results': len(results),
                        'optimization_enabled': enable_optimization,
                        'host_organism': host_organism,
                        'unique_enzymes': results_df['enzyme'].astype(str).nunique() if not results_df.empty else 0,
                    }
        except:
            pass  # If re-running fails, just show SDM result without mutation results
        
        # Return with SDM result and preserved form state + mutation results
        return render_template('silent_mutations.html',
                             sdm_result=sdm_result,
                             dna=dna,
                             cut_pos=cut_pos,
                             tolerance=tolerance,
                             window_codons=window_codons,
                             max_changes=max_changes,
                             host_organism=host_organism,
                             enable_optimization=enable_optimization,
                             results=ranked_results,
                             top_variants=top_variants,
                             summary=summary,
                             success=True if ranked_results else False)
        
    except ValueError as e:
        flash(f"Error: {str(e)}", "error")
        return redirect(url_for('silent_mutations'))
    except Exception as e:
        flash(f"An error occurred: {str(e)}", "error")
        return redirect(url_for('silent_mutations'))


@app.route('/primer-design', methods=['GET', 'POST'])
def primer_design():
    if request.method == 'GET':
        # Get available enzymes for the dropdown
        available_enzymes = get_enzyme_site_dict()
        motif_library = get_motif_library()
        return render_template('primer_design.html',
                             available_enzymes=available_enzymes,
                             motif_library=motif_library,
                             enzyme_result=None,
                             insert_result=None,
                             insert_enzyme_result=None,
                             ligation_result=None,
                             amplify_result=None)

    try:
        # Get form data
        dna = request.form.get('dna', '').strip().upper().replace('\n', '').replace(' ', '')
        pos = int(request.form.get('pos', 0) or 0)
        mode = request.form.get('mode', 'enzyme')

        # Validate DNA sequence
        is_valid, message = validate_dna_sequence(dna)
        if not is_valid:
            flash(f"Invalid DNA sequence: {message}", "error")
            available_enzymes = get_enzyme_site_dict()
            motif_library = get_motif_library()
            return render_template('primer_design.html',
                                 error=message,
                                 available_enzymes=available_enzymes,
                                 motif_library=motif_library,
                                 dna=dna, pos=pos, mode=mode)

        # Validate position (not used in amplify mode, which locates the insert by sequence match)
        if mode != 'amplify' and (pos <= 0 or pos > len(dna)):
            flash(f"Position ({pos}) must be between 1 and {len(dna)}", "error")
            available_enzymes = get_enzyme_site_dict()
            motif_library = get_motif_library()
            return render_template('primer_design.html',
                                 error="Invalid position",
                                 available_enzymes=available_enzymes,
                                 motif_library=motif_library,
                                 dna=dna, pos=pos, mode=mode)

        # Get available enzymes for the dropdown, filtered to exclude those already in sequence
        available_enzymes = get_enzyme_site_dict(dna)
        motif_library = get_motif_library()

        enzyme_result = None
        insert_result = None
        insert_enzyme_result = None
        ligation_result = None
        amplify_result = None

        if mode == 'enzyme':
            # Enzyme insertion mode
            enzyme_name = request.form.get('enzyme_name')
            if not enzyme_name:
                flash("Please select an enzyme", "error")
                return render_template('primer_design.html',
                                         error="No enzyme selected",
                                         available_enzymes=available_enzymes,
                                         motif_library=motif_library,
                                         dna=dna, pos=pos, mode=mode)
            
            # Check if enzyme exists in our list
            if enzyme_name not in available_enzymes:
                flash("Selected enzyme not found", "error")
                return render_template('primer_design.html',
                                         error="Invalid enzyme selection",
                                         available_enzymes=available_enzymes,
                                         motif_library=motif_library,
                                         dna=dna, pos=pos, mode=mode)
            
            enzyme_site = available_enzymes[enzyme_name]
            
            # Design enzyme insertion primers
            enzyme_result = design_enzyme_primers(dna, pos, enzyme_name, enzyme_site)
            
            # Check if primer design failed
            if enzyme_result is None:
                flash("Could not design primer pair meeting constraints (min len, max combined, ΔTm). Consider providing more flanking sequence.", "error")
                return render_template('primer_design.html',
                                         error="Primer design failed",
                                         available_enzymes=available_enzymes,
                                         motif_library=motif_library,
                                         dna=dna, pos=pos, mode=mode,
                                         enzyme_name=enzyme_name)
        
        elif mode == 'insert':
            # Amino acid OR nucleotide insertion mode
            aa_seq = request.form.get('aa_seq', '').strip().upper()
            insert_type = request.form.get('insert_type', 'aa').strip().lower()
            if insert_type not in ('aa', 'dna'):
                insert_type = 'aa'
            enable_ligation = request.form.get('enable_ligation') == 'on'
            ligation_enzyme = request.form.get('ligation_enzyme', '')
            overhang_len = int(request.form.get('overhang_len', 4))
            insert_enzyme_name = request.form.get('insert_enzyme_name', '').strip()
            insert_owned_enzymes = request.form.get('insert_owned_enzymes', '').strip()

            seq_label = "DNA" if insert_type == 'dna' else "amino acid"

            if not aa_seq:
                flash(f"Please enter a {seq_label} sequence", "error")
                return render_template('primer_design.html',
                                     error=f"No {seq_label} sequence provided",
                                     available_enzymes=available_enzymes,
                                     motif_library=motif_library,
                                     dna=dna, pos=pos, mode=mode, insert_type=insert_type,
                                     enable_ligation=enable_ligation,
                                     ligation_enzyme=ligation_enzyme, overhang_len=overhang_len,
                                     insert_enzyme_name=insert_enzyme_name, insert_owned_enzymes=insert_owned_enzymes)

            # Validate and resolve to a nucleotide insert
            if insert_type == 'dna':
                is_valid_insert, msg_insert = validate_dna_sequence(aa_seq)
                if not is_valid_insert:
                    flash(f"Invalid DNA sequence: {msg_insert}", "error")
                    return render_template('primer_design.html',
                                         error="Invalid DNA sequence",
                                         available_enzymes=available_enzymes,
                                         motif_library=motif_library,
                                         dna=dna, pos=pos, mode=mode, aa_seq=aa_seq, insert_type=insert_type,
                                         enable_ligation=enable_ligation, ligation_enzyme=ligation_enzyme,
                                         overhang_len=overhang_len,
                                         insert_enzyme_name=insert_enzyme_name, insert_owned_enzymes=insert_owned_enzymes)
                resolved_insert_nt = aa_seq
            else:
                valid_aas = set('ARNDCQEGHILKMFPSTWYV')
                if not set(aa_seq).issubset(valid_aas):
                    invalid_aas = set(aa_seq) - valid_aas
                    flash(f"Invalid amino acids: {', '.join(invalid_aas)}", "error")
                    return render_template('primer_design.html',
                                         error="Invalid amino acid sequence",
                                         available_enzymes=available_enzymes,
                                         motif_library=motif_library,
                                         dna=dna, pos=pos, mode=mode, aa_seq=aa_seq, insert_type=insert_type,
                                         enable_ligation=enable_ligation, ligation_enzyme=ligation_enzyme,
                                         overhang_len=overhang_len,
                                         insert_enzyme_name=insert_enzyme_name, insert_owned_enzymes=insert_owned_enzymes)
                resolved_insert_nt = aa_to_dna(aa_seq)

            # Design insertion primers using the resolved nucleotide insert
            insert_result = design_insertion_primers(dna, pos, resolved_insert_nt)

            # Also design enzyme insertion primers if enzyme was selected
            insert_enzyme_result = None
            if insert_enzyme_name and insert_enzyme_name in available_enzymes:
                enzyme_site = available_enzymes[insert_enzyme_name]
                insert_enzyme_result = design_enzyme_primers(dna, pos, insert_enzyme_name, enzyme_site)

            # Check if primer design failed
            if insert_result is None:
                flash("Could not design primer pair meeting constraints (min len, max combined, ΔTm). Consider providing more flanking sequence.", "error")
                return render_template('primer_design.html',
                                         error="Primer design failed",
                                         available_enzymes=available_enzymes,
                                         motif_library=motif_library,
                                         dna=dna, pos=pos, mode=mode, aa_seq=aa_seq, insert_type=insert_type,
                                         enable_ligation=enable_ligation, ligation_enzyme=ligation_enzyme,
                                         overhang_len=overhang_len)

            # Check if homology sequence is too long
            if len(insert_result['homology_full']) > 60:
                flash(f"Warning: Homology sequence is {len(insert_result['homology_full'])} nt (>60 nt may be difficult to synthesize)", "warning")

            # If restriction-ligation mode is enabled, also generate ligation oligos
            if enable_ligation and ligation_enzyme:
                if ligation_enzyme in available_enzymes:
                    enzyme_site = available_enzymes[ligation_enzyme]
                    ligation_result = design_enzyme_and_insert(dna, pos, ligation_enzyme, enzyme_site, resolved_insert_nt, overhang_len=overhang_len)
                    ligation_result['aa_seq'] = aa_seq
                else:
                    flash(f"Selected enzyme '{ligation_enzyme}' not found for ligation", "warning")

        elif mode == 'amplify':
            # Indel primers: amplify a target region, or cut it out (delete)
            insert_amplify_seq = request.form.get('insert_amplify_seq', '').strip().upper().replace('\n', '').replace(' ', '')
            operation = request.form.get('operation', 'amplify')
            if operation not in ('amplify', 'delete'):
                operation = 'amplify'
            tolerance_mode = request.form.get('tolerance_mode', 'edge')
            tolerance_bp_str = request.form.get('tolerance_bp', '30').strip()

            if not insert_amplify_seq:
                flash("Please enter the target region sequence", "error")
                return render_template('primer_design.html',
                                     error="No target sequence provided",
                                     available_enzymes=available_enzymes,
                                     motif_library=motif_library,
                                     dna=dna, pos=pos, mode=mode,
                                     insert_amplify_seq=insert_amplify_seq,
                                     operation=operation,
                                     tolerance_mode=tolerance_mode,
                                     tolerance_bp=tolerance_bp_str)

            is_valid_ins, msg_ins = validate_dna_sequence(insert_amplify_seq)
            if not is_valid_ins:
                flash(f"Invalid target sequence: {msg_ins}", "error")
                return render_template('primer_design.html',
                                     error="Invalid target sequence",
                                     available_enzymes=available_enzymes,
                                     motif_library=motif_library,
                                     dna=dna, pos=pos, mode=mode,
                                     insert_amplify_seq=insert_amplify_seq,
                                     operation=operation,
                                     tolerance_mode=tolerance_mode,
                                     tolerance_bp=tolerance_bp_str)

            try:
                tolerance_bp = int(tolerance_bp_str) if tolerance_bp_str else 30
            except ValueError:
                tolerance_bp = 30
            tolerance_bp = max(0, tolerance_bp)

            try:
                amplify_result = design_indel_primers(
                    dna, insert_amplify_seq, operation=operation,
                    tolerance_mode=tolerance_mode, tolerance_bp=tolerance_bp)
            except ValueError as amp_error:
                flash(str(amp_error), "error")
                return render_template('primer_design.html',
                                     error=str(amp_error),
                                     available_enzymes=available_enzymes,
                                     motif_library=motif_library,
                                     dna=dna, pos=pos, mode=mode,
                                     insert_amplify_seq=insert_amplify_seq,
                                     operation=operation,
                                     tolerance_mode=tolerance_mode,
                                     tolerance_bp=tolerance_bp_str)

        return render_template('primer_design.html',
                                 available_enzymes=available_enzymes,
                                 motif_library=motif_library,
                                 enzyme_result=enzyme_result,
                                 insert_result=insert_result,
                                 insert_enzyme_result=insert_enzyme_result,
                                 ligation_result=ligation_result,
                                 amplify_result=amplify_result,
                                 dna=dna,
                                 pos=pos,
                                 mode=mode,
                                 aa_seq=request.form.get('aa_seq', '') if mode == 'insert' else '',
                                 insert_type=request.form.get('insert_type', 'aa') if mode == 'insert' else 'aa',
                                 enzyme_name=request.form.get('enzyme_name', '') if mode == 'enzyme' else '',
                                 insert_enzyme_name=request.form.get('insert_enzyme_name', '') if mode == 'insert' else '',
                                 insert_owned_enzymes=request.form.get('insert_owned_enzymes', '') if mode == 'insert' else '',
                                 enable_ligation=request.form.get('enable_ligation') == 'on' if mode == 'insert' else False,
                                 ligation_enzyme=request.form.get('ligation_enzyme', '') if mode == 'insert' else '',
                                 overhang_len=request.form.get('overhang_len', '4') if mode == 'insert' else '4',
                                 insert_amplify_seq=request.form.get('insert_amplify_seq', '') if mode == 'amplify' else '',
                                 operation=request.form.get('operation', 'amplify') if mode == 'amplify' else 'amplify',
                                 tolerance_mode=request.form.get('tolerance_mode', 'edge') if mode == 'amplify' else 'edge',
                                 tolerance_bp=request.form.get('tolerance_bp', '30') if mode == 'amplify' else '30')
        
    except ValueError as e:
        flash(f"Invalid input: {str(e)}", "error")
        available_enzymes = get_enzyme_site_dict()
        motif_library = get_motif_library()
        owned_enzymes = request.form.get('owned_enzymes', '').strip()
        return render_template('primer_design.html', 
                             error="Invalid numerical input",
                             available_enzymes=available_enzymes,
                             motif_library=motif_library,
                             owned_enzymes=owned_enzymes)
    except Exception as e:
        flash(f"An error occurred: {str(e)}", "error")
        available_enzymes = get_enzyme_site_dict()
        motif_library = get_motif_library()
        owned_enzymes = request.form.get('owned_enzymes', '').strip()
        return render_template('primer_design.html', 
                             error="Processing error occurred",
                             available_enzymes=available_enzymes,
                             motif_library=motif_library,
                             owned_enzymes=owned_enzymes)


def build_sequence_viz(dna, mutated_dna, mutation_start_0, mutation_end_0,
                        fwd_start_0, fwd_end_0, rev_start_0, rev_end_0,
                        context_pad=90, row_size=60, center_0=None):
    """
    Build a context window of the sequence around the mutation for visualization.
    Returns rows of annotated nucleotides for rendering.
    center_0 overrides context centering when there is no mutation (e.g. primers-only mode).
    """
    seq_len = len(dna)
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    # Determine context window
    if mutation_start_0 is not None:
        ctx_start = max(0, mutation_start_0 - context_pad)
        ctx_end = min(seq_len, (mutation_end_0 if mutation_end_0 else mutation_start_0 + 1) + context_pad)
    elif center_0 is not None:
        ctx_start = max(0, center_0 - context_pad)
        ctx_end = min(seq_len, center_0 + context_pad)
    else:
        ctx_start = 0
        ctx_end = min(seq_len, 180)

    # Snap to row boundaries
    ctx_start = (ctx_start // row_size) * row_size
    ctx_end = min(seq_len, ((ctx_end // row_size) + 1) * row_size)

    rows = []
    for row_start in range(ctx_start, ctx_end, row_size):
        row_end = min(row_start + row_size, ctx_end)
        nts = []
        comp_nts = []

        # Pixel constants: must match CSS exactly
        NT_PX = 11    # .nt { width: 11px }
        SPACER_PX = 44  # .strand-label(28) + .strand-dir(16)

        # Forward primer bar for this row — pixel-precise positioning
        fwd_bar = None
        if fwd_start_0 is not None and fwd_end_0 is not None:
            bar_s = max(fwd_start_0, row_start)
            bar_e = min(fwd_end_0, row_end)
            if bar_s < bar_e:
                fwd_bar = {
                    'left_px': SPACER_PX + (bar_s - row_start) * NT_PX,
                    'width_px': (bar_e - bar_s) * NT_PX,
                    'is_start': fwd_start_0 >= row_start,
                    'is_end': fwd_end_0 <= row_end,
                }

        # Reverse primer bar — pixel-precise positioning
        rev_bar = None
        if rev_start_0 is not None and rev_end_0 is not None:
            bar_s = max(rev_start_0, row_start)
            bar_e = min(rev_end_0, row_end)
            if bar_s < bar_e:
                rev_bar = {
                    'left_px': SPACER_PX + (bar_s - row_start) * NT_PX,
                    'width_px': (bar_e - bar_s) * NT_PX,
                    'is_start': rev_start_0 >= row_start,
                    'is_end': rev_end_0 <= row_end,
                }

        for i in range(row_start, row_end):
            base = dna[i] if i < len(dna) else 'N'
            mut_base = mutated_dna[i] if mutated_dna and i < len(mutated_dna) else base
            comp_base = complement_map.get(base, 'N')

            in_mut = mutation_start_0 is not None and mutation_start_0 <= i < mutation_end_0
            in_fwd = fwd_start_0 is not None and fwd_start_0 <= i < fwd_end_0
            in_rev = rev_start_0 is not None and rev_start_0 <= i < rev_end_0

            nts.append({
                'pos': i + 1,
                'base': base,
                'mut_base': mut_base,
                'changed': base != mut_base,
                'in_mut': in_mut,
                'in_fwd': in_fwd,
                'in_rev': in_rev,
            })
            comp_nts.append({
                'pos': i + 1,
                'base': comp_base,
                'in_mut': in_mut,
                'in_rev': in_rev,
            })

        # Ruler ticks: pixel-positioned labels every 10 bases
        ruler_ticks = []
        for i in range(row_start, row_end):
            pos_1 = i + 1
            if pos_1 % 10 == 0 or i == row_start:
                ruler_ticks.append({
                    'px': SPACER_PX + (i - row_start) * NT_PX,
                    'pos': pos_1,
                })

        rows.append({
            'start': row_start + 1,
            'end': row_end,
            'nts': nts,
            'comp_nts': comp_nts,
            'fwd_bar': fwd_bar,
            'rev_bar': rev_bar,
            'ruler_ticks': ruler_ticks,
        })

    return rows


@app.route('/sequence-visualizer', methods=['GET', 'POST'])
def sequence_visualizer():
    if request.method == 'GET':
        return render_template('sequence_visualizer.html')

    try:
        dna = request.form.get('dna', '').strip().upper().replace('\n', '').replace('\r', '').replace(' ', '')
        mode = request.form.get('mode', 'sdm')

        is_valid, message = validate_dna_sequence(dna)
        if not is_valid:
            flash(f"Invalid DNA sequence: {message}", "error")
            return render_template('sequence_visualizer.html', dna=dna, mode=mode)

        primer_result = None
        mutated_dna = dna
        mutation_start_0 = None
        mutation_end_0 = None
        fwd_start_0 = None
        fwd_end_0 = None
        rev_start_0 = None
        rev_end_0 = None
        insert_label = None
        center_0 = None

        if mode == 'sdm':
            original_window = request.form.get('original_window', '').strip().upper().replace(' ', '')
            mutated_window = request.form.get('mutated_window', '').strip().upper().replace(' ', '')

            if not original_window or not mutated_window:
                flash("Please enter both the original and mutated sequence windows.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                       original_window=original_window, mutated_window=mutated_window)

            if len(original_window) != len(mutated_window):
                flash("Original and mutated windows must be the same length.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                       original_window=original_window, mutated_window=mutated_window)

            is_valid_o, _ = validate_dna_sequence(original_window)
            is_valid_m, _ = validate_dna_sequence(mutated_window)
            if not is_valid_o or not is_valid_m:
                flash("Windows must contain only A, T, C, G.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                       original_window=original_window, mutated_window=mutated_window)

            if dna.find(original_window) == -1:
                flash("Original window not found in the plasmid sequence.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                       original_window=original_window, mutated_window=mutated_window)

            primer_result = design_sdm_primers(dna, original_window, mutated_window)

            mut_start_1 = primer_result['mutation_start']
            mut_end_1 = primer_result['mutation_end']
            mutation_start_0 = mut_start_1 - 1
            mutation_end_0 = mut_end_1

            window_start_0 = dna.find(original_window)
            mutated_dna = dna[:window_start_0] + mutated_window + dna[window_start_0 + len(original_window):]

            fwd_start_0_raw = primer_result.get('fwd_primer_start_0')
            if fwd_start_0_raw is not None:
                fwd_start_0 = fwd_start_0_raw
                fwd_end_0 = fwd_start_0 + len(primer_result['f_primer'])
                rev_end_0 = fwd_start_0
                rev_start_0 = max(0, rev_end_0 - len(primer_result['r_primer']))
            else:
                fwd_start_0 = max(0, mutation_start_0 - len(primer_result['f_primer']) // 2)
                fwd_end_0 = fwd_start_0 + len(primer_result['f_primer'])
                rev_end_0 = fwd_start_0
                rev_start_0 = max(0, rev_end_0 - len(primer_result['r_primer']))

        elif mode == 'insertion':
            position_str = request.form.get('position', '1').strip()
            insert_seq = request.form.get('insert_seq', '').strip().upper().replace(' ', '')
            insert_type = request.form.get('insert_type', 'dna')

            try:
                position = int(position_str)
            except ValueError:
                flash("Position must be a valid integer.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode)

            if position < 1 or position > len(dna):
                flash(f"Position must be between 1 and {len(dna)}.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode, position=position_str)

            if not insert_seq:
                flash("Please enter an insert sequence.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode, position=position_str)

            if insert_type == 'aa':
                valid_aas = set('ARNDCQEGHILKMFPSTWYV')
                if not set(insert_seq).issubset(valid_aas):
                    flash("Invalid amino acid characters in insert.", "error")
                    return render_template('sequence_visualizer.html', dna=dna, mode=mode, position=position_str, insert_seq=insert_seq)
                insert_dna = aa_to_dna(insert_seq)
                insert_label = f"AA insert: {insert_seq} → {insert_dna}"
            else:
                is_valid_ins, msg_ins = validate_dna_sequence(insert_seq)
                if not is_valid_ins:
                    flash(f"Invalid insert DNA: {msg_ins}", "error")
                    return render_template('sequence_visualizer.html', dna=dna, mode=mode, position=position_str, insert_seq=insert_seq)
                insert_dna = insert_seq
                insert_label = f"DNA insert: {insert_dna}"

            pos0 = position - 1
            mutated_dna = dna[:pos0] + insert_dna + dna[pos0:]
            mutation_start_0 = pos0
            mutation_end_0 = pos0 + len(insert_dna)

            primer_result = design_enzyme_insertion_primers(dna, position, insert_dna)

            fwd_start_0 = pos0
            fwd_end_0 = min(len(dna), pos0 + len(primer_result['f_primer']) - len(insert_dna))
            rev_end_0 = pos0
            rev_start_0 = max(0, pos0 - len(primer_result['r_primer']))
            insert_label = insert_label

        elif mode == 'amplify':
            insert_amplify_seq = request.form.get('insert_amplify_seq', '').strip().upper().replace('\n', '').replace(' ', '')
            operation = request.form.get('operation', 'amplify')
            if operation not in ('amplify', 'delete'):
                operation = 'amplify'
            tolerance_mode = request.form.get('tolerance_mode', 'edge')
            tolerance_bp_str = request.form.get('tolerance_bp', '30').strip()

            if not insert_amplify_seq:
                flash("Please enter the target region sequence.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                       operation=operation,
                                       tolerance_mode=tolerance_mode, tolerance_bp=tolerance_bp_str)

            is_valid_ins, msg_ins = validate_dna_sequence(insert_amplify_seq)
            if not is_valid_ins:
                flash(f"Invalid target sequence: {msg_ins}", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                       insert_amplify_seq=insert_amplify_seq,
                                       operation=operation,
                                       tolerance_mode=tolerance_mode, tolerance_bp=tolerance_bp_str)

            try:
                tolerance_bp = int(tolerance_bp_str) if tolerance_bp_str else 30
            except ValueError:
                tolerance_bp = 30
            tolerance_bp = max(0, tolerance_bp)

            try:
                primer_result = design_indel_primers(
                    dna, insert_amplify_seq, operation=operation,
                    tolerance_mode=tolerance_mode, tolerance_bp=tolerance_bp)
            except ValueError as amp_error:
                flash(str(amp_error), "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                       insert_amplify_seq=insert_amplify_seq,
                                       operation=operation,
                                       tolerance_mode=tolerance_mode, tolerance_bp=tolerance_bp_str)

            mutation_start_0 = primer_result['insert_start_0']
            mutation_end_0 = primer_result['insert_end_0']
            mutated_dna = dna

            fwd_start_0 = primer_result['fwd_start_0']
            fwd_end_0 = primer_result['fwd_end_0']
            rev_start_0 = primer_result['rev_start_0']
            rev_end_0 = primer_result['rev_end_0']

            insert_preview = insert_amplify_seq if len(insert_amplify_seq) <= 40 else insert_amplify_seq[:40] + '…'
            if operation == 'delete':
                insert_label = f"Region to remove: {insert_preview} ({primer_result['insert_len']} nt)"
            else:
                insert_label = f"Insert to amplify: {insert_preview} ({primer_result['insert_len']} nt)"

        elif mode == 'primers':
            fwd_input = request.form.get('fwd_primer_input', '').strip().upper().replace(' ', '')
            rev_input = request.form.get('rev_primer_input', '').strip().upper().replace(' ', '')

            if not fwd_input and not rev_input:
                flash("Please enter at least one primer sequence.", "error")
                return render_template('sequence_visualizer.html', dna=dna, mode=mode)

            # Validate
            if fwd_input:
                is_valid_f, msg_f = validate_dna_sequence(fwd_input)
                if not is_valid_f:
                    flash(f"Invalid forward primer: {msg_f}", "error")
                    return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                           fwd_primer_input=fwd_input, rev_primer_input=rev_input)
            if rev_input:
                is_valid_r, msg_r = validate_dna_sequence(rev_input)
                if not is_valid_r:
                    flash(f"Invalid reverse primer: {msg_r}", "error")
                    return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                           fwd_primer_input=fwd_input, rev_primer_input=rev_input)

            # Find forward primer binding site in top strand
            if fwd_input:
                fwd_pos = dna.find(fwd_input)
                if fwd_pos == -1:
                    flash(f"Forward primer not found in the plasmid sequence. Make sure it matches a substring of the top strand (5'→3').", "error")
                    return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                           fwd_primer_input=fwd_input, rev_primer_input=rev_input)
                fwd_start_0 = fwd_pos
                fwd_end_0 = fwd_pos + len(fwd_input)
            else:
                fwd_start_0 = fwd_end_0 = None

            # Find reverse primer binding site: search for its reverse complement in the top strand
            if rev_input:
                rev_rc = reverse_complement(rev_input)
                rev_pos = dna.find(rev_rc)
                if rev_pos == -1:
                    flash(f"Reverse primer binding site not found. The reverse complement of your reverse primer should appear in the top strand.", "error")
                    return render_template('sequence_visualizer.html', dna=dna, mode=mode,
                                           fwd_primer_input=fwd_input, rev_primer_input=rev_input)
                rev_start_0 = rev_pos
                rev_end_0 = rev_pos + len(rev_input)
            else:
                rev_start_0 = rev_end_0 = None

            # No mutation overlay in this mode
            mutation_start_0 = None
            mutation_end_0 = None
            mutated_dna = dna

            # Center view on the amplicon midpoint
            if fwd_start_0 is not None and rev_end_0 is not None:
                center_0 = (fwd_start_0 + rev_end_0) // 2
            elif fwd_start_0 is not None:
                center_0 = (fwd_start_0 + fwd_end_0) // 2
            elif rev_start_0 is not None:
                center_0 = (rev_start_0 + rev_end_0) // 2
            else:
                center_0 = None

            # Build primer stats
            f_tm, f_gc = tm_nn_neb(fwd_input) if fwd_input else (0.0, 0.0)
            r_tm, r_gc = tm_nn_neb(rev_input) if rev_input else (0.0, 0.0)
            avg_tm = round((f_tm + r_tm) / 2, 1) if (fwd_input and rev_input) else (f_tm or r_tm)
            tm_diff = round(abs(f_tm - r_tm), 1) if (fwd_input and rev_input) else 0.0
            anneal_temp = neb_q5_anneal_temp(f_tm, r_tm) if (fwd_input and rev_input) else None

            amplicon_len = None
            if fwd_end_0 is not None and rev_start_0 is not None and rev_start_0 > fwd_end_0:
                amplicon_len = rev_end_0 - fwd_start_0

            primer_result = {
                'f_primer': fwd_input,
                'r_primer': rev_input,
                'f_tm': f_tm,
                'r_tm': r_tm,
                'f_gc': f_gc,
                'r_gc': r_gc,
                'avg_tm': avg_tm,
                'tm_diff': tm_diff,
                'recommended_anneal': anneal_temp,
                'amplicon_len': amplicon_len,
            }

        # Determine center for view (used when no mutation)
        viz_center = center_0 if mode == 'primers' else None

        viz_rows = build_sequence_viz(
            dna, mutated_dna,
            mutation_start_0, mutation_end_0,
            fwd_start_0, fwd_end_0,
            rev_start_0, rev_end_0,
            center_0=viz_center
        )

        return render_template('sequence_visualizer.html',
                               dna=dna,
                               mutated_dna=mutated_dna,
                               mode=mode,
                               primer_result=primer_result,
                               viz_rows=viz_rows,
                               mutation_start_0=mutation_start_0,
                               mutation_end_0=mutation_end_0,
                               fwd_start_0=fwd_start_0,
                               fwd_end_0=fwd_end_0,
                               rev_start_0=rev_start_0,
                               rev_end_0=rev_end_0,
                               insert_label=insert_label,
                               original_window=request.form.get('original_window', ''),
                               mutated_window=request.form.get('mutated_window', ''),
                               insert_seq=request.form.get('insert_seq', ''),
                               position=request.form.get('position', ''),
                               insert_type=request.form.get('insert_type', 'dna'),
                               fwd_primer_input=request.form.get('fwd_primer_input', ''),
                               rev_primer_input=request.form.get('rev_primer_input', ''),
                               insert_amplify_seq=request.form.get('insert_amplify_seq', ''),
                               operation=request.form.get('operation', 'amplify'),
                               tolerance_mode=request.form.get('tolerance_mode', 'edge'),
                               tolerance_bp=request.form.get('tolerance_bp', '30'),
                               success=True)

    except ValueError as e:
        flash(f"Error: {str(e)}", "error")
        return render_template('sequence_visualizer.html',
                               dna=request.form.get('dna', ''),
                               mode=request.form.get('mode', 'sdm'))
    except Exception as e:
        logging.exception("Sequence visualizer error")
        flash(f"An error occurred: {str(e)}", "error")
        return render_template('sequence_visualizer.html',
                               dna=request.form.get('dna', ''),
                               mode=request.form.get('mode', 'sdm'))


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)