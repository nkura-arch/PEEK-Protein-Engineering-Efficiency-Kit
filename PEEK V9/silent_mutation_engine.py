# Silent Mutation Detection Engine
# Extracted from the original bioinformatics algorithm

from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Restriction import CommOnly
import itertools
import csv
import pandas as pd
from codon_optimizer import CodonOptimizer


# Setup codon tables
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
SYNONYMOUS = {}
for codon, aa in standard_table.forward_table.items():
    SYNONYMOUS.setdefault(aa, []).append(codon)
CODON_TO_AA = standard_table.forward_table


def translate(seq):
    """Translate DNA sequence to amino acids"""
    return str(Seq(seq).translate())


def get_synonymous(codon):
    """Get synonymous codons for a given codon"""
    aa = CODON_TO_AA.get(codon, None)
    if not aa:
        return [codon]
    return [syn for syn in SYNONYMOUS[aa] if syn != codon]


def site_in_window(seq, enzyme, cut_pos, tolerance):
    """Check if enzyme recognition site is within tolerance of cut position"""
    rec = str(enzyme.site)
    rec_len = len(rec)
    cut_index = enzyme.fst5 if enzyme.fst5 is not None else 0
    for i in range(len(seq) - rec_len + 1):
        if seq[i:i + rec_len] == rec:
            cut_site = i + cut_index
            if abs(cut_site - cut_pos) <= tolerance:
                return cut_site
    return None


def is_valid_window(seq):
    """Check if window is valid: non-empty, a whole number of codons, and free
    of internal stop codons."""
    if not seq or len(seq) % 3 != 0:
        return False
    return "*" not in translate(seq)


def _round_to_frame(value, frame_offset, direction):
    """Round `value` to the nearest position congruent to frame_offset (mod 3).

    direction='up' rounds up (never decreases below `value`), used for lower
    bounds so we never start before the window we were asked for.
    direction='down' rounds down (never increases above `value`), used for
    upper bounds so we never extend past the window/sequence we were asked
    for. Either way, the result always lands exactly on a codon boundary.
    """
    remainder = (value - frame_offset) % 3
    if remainder == 0:
        return value
    if direction == 'up':
        return value + (3 - remainder)
    return value - remainder


STOP_CODONS = {"TAA", "TAG", "TGA"}


def _stop_free_run_length(dna, position, frame_offset):
    """Length of the contiguous run of non-stop codons, in the given frame,
    that covers `position`. Used to auto-detect which of the 3 possible
    reading frames is the real coding frame at a given position: only the
    true frame will typically run for many uninterrupted codons, while the
    other two frames hit a stop codon almost immediately."""
    n = len(dna)
    codon_start = _round_to_frame(position, frame_offset, 'down')

    left = codon_start
    while left - 3 >= 0 and dna[left - 3:left] not in STOP_CODONS:
        left -= 3

    right = codon_start
    while right + 3 <= n and dna[right:right + 3] not in STOP_CODONS:
        right += 3

    return right - left


def detect_coding_frame(dna, cut_pos):
    """Auto-detect the reading frame (as a coding_start_pos-equivalent 0, 1,
    or 2) that applies at cut_pos, without requiring the caller to know
    where the CDS starts.

    Works whether `dna` is just a CDS or a much larger sequence (e.g. a
    whole plasmid) -- it picks whichever of the 3 possible frames has the
    longest stop-codon-free run through cut_pos, which is virtually always
    the real coding frame for any non-trivial open reading frame.
    """
    best_offset = 0
    best_length = -1
    for frame_offset in (0, 1, 2):
        length = _stop_free_run_length(dna, cut_pos, frame_offset)
        if length > best_length:
            best_length = length
            best_offset = frame_offset
    return best_offset


def get_codon_aligned_window(dna_len, cut_pos, window_codons, coding_start_pos=0):
    """Get a codon-aligned window around cut position.

    Guarantees (checked, not just assumed):
      - start and end are both exactly on a codon boundary relative to
        coding_start_pos (never mid-codon)
      - 0 <= start <= end <= dna_len
      - (end - start) is a whole number of codons

    If no valid in-frame window exists within [0, dna_len] (e.g. dna_len is
    too short to contain even one full codon in this frame), an empty
    zero-length window at the nearest valid boundary is returned instead of
    ever returning a misaligned or out-of-bounds pair.
    """
    frame_offset = coding_start_pos % 3

    # The in-frame codon boundary at or before cut_pos.
    codon_start = _round_to_frame(cut_pos, frame_offset, 'down')

    half_window = window_codons // 2
    raw_start = codon_start - half_window * 3
    raw_end = codon_start + (window_codons - half_window) * 3

    # The valid range for any boundary is [0, dna_len]; snap the requested
    # raw start/end into that range FIRST, then round to frame -- rounding
    # must always move inward (start up, end down) so the result can never
    # slip outside [0, dna_len].
    start = _round_to_frame(max(0, raw_start), frame_offset, 'up')
    end = _round_to_frame(min(dna_len, raw_end), frame_offset, 'down')

    # Degenerate case: no in-frame codon boundary fits in [0, dna_len] at all
    # (only possible for pathologically short sequences, e.g. dna_len smaller
    # than the frame offset itself). Collapse to an empty window rather than
    # ever hand back start > end or a mid-codon boundary.
    if start > end:
        start = end = 0

    # Hard guarantee -- catches any future change to the math above that
    # would otherwise silently reintroduce a mid-codon boundary. An empty
    # window carries no codons, so alignment only needs to hold when
    # start != end.
    if start != end:
        assert (start - frame_offset) % 3 == 0, "window start is not codon-aligned"
        assert (end - frame_offset) % 3 == 0, "window end is not codon-aligned"
    assert 0 <= start <= end <= dna_len, "window is out of bounds"

    return start, end, codon_start


def find_silent_sites(dna, cut_pos, tolerance, window_codons, enzymes, coding_start_pos=None, max_changes=3):
    """Find silent mutations that introduce restriction enzyme sites.

    coding_start_pos is optional. If not given, the reading frame at
    cut_pos is auto-detected (see detect_coding_frame) -- callers don't need
    to know or supply where the CDS starts, even if `dna` is a whole plasmid
    rather than just the gene.
    """
    results = []
    if coding_start_pos is None:
        coding_start_pos = detect_coding_frame(dna, cut_pos)
    start, end, codon_start = get_codon_aligned_window(len(dna), cut_pos, window_codons, coding_start_pos)

    # Defense in depth: re-verify the boundaries returned above are actually
    # on codon boundaries before slicing. If this ever fails, a mutation
    # could split a codon, so we bail out instead of returning bad results.
    frame_offset = coding_start_pos % 3
    if (start - frame_offset) % 3 != 0 or (end - frame_offset) % 3 != 0 or start > end:
        return []

    window = dna[start:end]

    if not is_valid_window(window):
        return []
    
    codons = [window[i:i + 3] for i in range(0, len(window), 3)]
    codon_indices = [i for i, codon in enumerate(codons) if get_synonymous(codon)]
    orig_aa = translate(window)
    
    # Prioritize preferred enzymes for faster results
    preferred_enzymes = get_preferred_enzymes()
    preferred_set = set(preferred_enzymes)
    priority_enzymes = [e for e in enzymes if str(e) in preferred_set]
    other_enzymes = [e for e in enzymes if str(e) not in preferred_set]
    
    # Process preferred enzymes first, then others
    enzyme_list = priority_enzymes + other_enzymes
    
    # Limit search space for performance
    max_combinations = min(1000, len(enzyme_list) * 10)
    combinations_tested = 0
    
    for k in range(1, min(max_changes + 1, len(codon_indices) + 1)):
        for combo in itertools.combinations(codon_indices, k):
            if combinations_tested >= max_combinations:
                break
                
            syn_options = [get_synonymous(codons[i]) for i in combo]
            for replacements in itertools.product(*syn_options):
                if combinations_tested >= max_combinations:
                    break
                    
                combinations_tested += 1
                new_codons = codons[:]
                for idx, new_codon in zip(combo, replacements):
                    new_codons[idx] = new_codon
                mutated = ''.join(new_codons)
                if translate(mutated) != orig_aa:
                    continue
                    
                for enzyme in enzyme_list:
                    cut = site_in_window(mutated, enzyme, cut_pos - start, tolerance)
                    if cut is not None:
                        # Determine removed sequence (difference from original)
                        removed_bases = ''.join([o for o, m in zip(window, mutated) if o != m])
                        side = "5'" if mutated.find(str(enzyme.site)) < len(mutated) // 2 else "3'"
                        
                        results.append({
                            'enzyme': enzyme,
                            'new_cut': start + cut,
                            'orig_window': window,
                            'mutated': mutated,
                            'aa': translate(mutated),
                            'window_start': start,
                            'window_end': end,
                            'codon_start': codon_start,
                            'removed_bases': removed_bases,
                            'removed_side': side
                        })
                        
                        # Stop early if we found enough preferred results
                        if len(results) >= 50 and str(enzyme) in preferred_set:
                            return results
            
            if combinations_tested >= max_combinations:
                break
    
    return results


def save_results_to_csv(results, filename="silent_mutation_results.csv"):
    """Save results to CSV file"""
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([
            "Result #", "Enzyme", "Recognition Site", "New Cut Position (0-based)",
            "Original Window Start", "Original Window End", "Codon Start in Window",
            "Original Window", "Mutated Window", "Amino Acids",
            "Removed Bases", "Removed Side (5'/3')"
        ])
        for i, result in enumerate(results, 1):
            writer.writerow([
                i,
                str(result['enzyme']),
                result['enzyme'].site,
                result['new_cut'],
                result['window_start'],
                result['window_end'],
                result['codon_start'],
                result['orig_window'],
                result['mutated'],
                result['aa'],
                result['removed_bases'],
                result['removed_side']
            ])
    return filename


def get_preferred_enzymes():
    """Get list of preferred enzymes"""
    return ["SalI", "EagI", "RsrII", "KpnI", "MluI", "SbfI", "NcoI", "DpnI"]


def filter_enzymes_not_in_sequence(dna):
    """Filter enzymes that are not already present in the sequence"""
    return [e for e in CommOnly if str(e.site) not in dna]


def rank_results_by_preference(results):
    """Rank results by enzyme preference"""
    preferred = get_preferred_enzymes()
    results_df = pd.DataFrame(results)
    
    if results_df.empty:
        return results_df, []
    
    # Add preference score (lower is better)
    results_df['preference_score'] = results_df['enzyme'].astype(str).apply(
        lambda x: preferred.index(x) if x in preferred else len(preferred)
    )
    
    # Sort by preference, then by number of changes (fewer is better)
    results_df = results_df.sort_values(['preference_score', 'removed_bases'])
    
    # Get matched preferred enzymes  
    matched_series = results_df[results_df['enzyme'].astype(str).isin(preferred)]['enzyme'].astype(str)
    matched_enzymes = matched_series.drop_duplicates().tolist()
    
    return results_df, matched_enzymes


def score_results_with_codon_optimization(results, original_seq, 
                                         host_organism="ecoli",
                                         weights=None):
    """
    Score silent mutation results using algorithmic codon optimization
    
    Args:
        results: List of results from find_silent_sites
        original_seq: Original DNA sequence
        host_organism: 'ecoli', 'human', or custom codon usage dict
        weights: Custom scoring weights
    
    Returns: DataFrame with scored and ranked results
    """
    if not results:
        return pd.DataFrame()
    
    optimizer = CodonOptimizer(host_organism=host_organism)
    
    scored_results = []
    for result in results:
        mutated = result.get('mutated')
        orig_window = result.get('orig_window') or result.get('original_window')

        # Use the codon-aligned original window when scoring (ensure equal lengths)
        if not orig_window or not mutated or len(orig_window) != len(mutated):
            # Skip inconsistent entries
            continue

        # Check if this mutation introduces a restriction site (by design true)
        has_restriction_site = True

        # Score the variant using the orig_window as the original sequence
        scores = optimizer.score_variant(
            orig_window,
            mutated,
            has_restriction_site=has_restriction_site,
            weights=weights
        )

        # Merge scores into result and include both sequences explicitly
        result_with_scores = result.copy()
        result_with_scores.update(scores)
        result_with_scores['original_seq'] = orig_window
        result_with_scores['variant_seq'] = mutated
        scored_results.append(result_with_scores)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(scored_results)
    
    # Sort by total score (highest first)
    if 'total_score' in results_df.columns:
        results_df = results_df.sort_values('total_score', ascending=False)
        # Assign rank (1 = best)
        results_df = results_df.reset_index(drop=True)
        results_df['rank'] = results_df.index + 1
    
    return results_df


def export_scored_results_to_csv(results_df, filename="codon_optimized_results.csv"):
    """Export scored results to CSV with all scoring metrics"""
    if results_df.empty:
        return filename
    
    # Select relevant columns for export
    export_cols = [
        'enzyme', 'mutated', 'aa', 'codon_changes',
        'cai', 'rare_penalty', 'structure_penalty', 
        'mutation_penalty', 'restriction_bonus', 'total_score'
    ]
    
    # Only include columns that exist
    available_cols = [col for col in export_cols if col in results_df.columns]
    
    results_df[available_cols].to_csv(filename, index=False)
    return filename


def get_top_codon_optimized_variants(results_df, top_n=5):
    """Get top N variants by codon optimization score"""
    if results_df.empty or 'total_score' not in results_df.columns:
        return []
    
    top_variants = results_df.nlargest(top_n, 'total_score')
    return top_variants.to_dict('records')
