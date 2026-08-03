# Algorithmic Codon Optimization with Explicit Scoring
# Deterministic, transparent, and scientifically defensible

from Bio.Seq import Seq
from Bio.Data import CodonTable
import RNA  # ViennaRNA for secondary structure prediction
import math
from typing import Dict, List, Tuple
import json


# Standard codon usage table for E. coli (commonly used model organism)
CODON_USAGE_ECOLI = {
    'TTT': 26, 'TTC': 26, 'TTA': 4, 'TTG': 12,
    'TCT': 15, 'TCC': 14, 'TCA': 5, 'TCG': 6,
    'TAT': 15, 'TAC': 15, 'TAA': 1, 'TAG': 0,
    'TGT': 10, 'TGC': 10, 'TGA': 1, 'TGG': 10,
    'CTT': 12, 'CTC': 10, 'CTA': 3, 'CTG': 47,
    'CCT': 13, 'CCC': 4, 'CCA': 6, 'CCG': 27,
    'CAT': 14, 'CAC': 14, 'CAA': 4, 'CAG': 37,
    'CGT': 36, 'CGC': 36, 'CGA': 7, 'CGG': 11,
    'ATT': 38, 'ATC': 29, 'ATA': 7, 'ATG': 26,
    'ACT': 18, 'ACC': 27, 'ACA': 8, 'ACG': 16,
    'AAT': 17, 'AAC': 17, 'AAA': 24, 'AAG': 24,
    'AGT': 7, 'AGC': 16, 'AGA': 5, 'AGG': 4,
    'GTT': 27, 'GTC': 20, 'GTA': 11, 'GTG': 37,
    'GCT': 26, 'GCC': 30, 'GCA': 16, 'GCG': 34,
    'GAT': 34, 'GAC': 34, 'GAA': 29, 'GAG': 16,
    'GGT': 35, 'GGC': 37, 'GGA': 12, 'GGG': 12,
}

# Alternative codon usage tables
CODON_USAGE_HUMAN = {
    'TTT': 45, 'TTC': 55, 'TTA': 27, 'TTG': 13,
    'TCT': 31, 'TCC': 28, 'TCA': 24, 'TCG': 6,
    'TAT': 44, 'TAC': 56, 'TAA': 1, 'TAG': 0,
    'TGT': 45, 'TGC': 55, 'TGA': 1, 'TGG': 100,
    'CTT': 32, 'CTC': 21, 'CTA': 13, 'CTG': 39,
    'CCT': 33, 'CCC': 20, 'CCA': 27, 'CCG': 11,
    'CAT': 42, 'CAC': 58, 'CAA': 25, 'CAG': 75,
    'CGT': 18, 'CGC': 36, 'CGA': 10, 'CGG': 19,
    'ATT': 46, 'ATC': 51, 'ATA': 24, 'ATG': 100,
    'ACT': 30, 'ACC': 46, 'ACA': 26, 'ACG': 11,
    'AAT': 46, 'AAC': 54, 'AAA': 43, 'AAG': 57,
    'AGT': 35, 'AGC': 70, 'AGA': 21, 'AGG': 20,
    'GTT': 37, 'GTC': 25, 'GTA': 11, 'GTG': 54,
    'GCT': 40, 'GCC': 41, 'GCA': 26, 'GCG': 15,
    'GAT': 46, 'GAC': 54, 'GAA': 42, 'GAG': 58,
    'GGT': 37, 'GGC': 49, 'GGA': 25, 'GGG': 16,
}

# Amino acid to codons mapping
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
AA_TO_CODONS = {}
for codon, aa in standard_table.forward_table.items():
    if aa not in AA_TO_CODONS:
        AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)


class CodonOptimizer:
    """Algorithmic codon optimization with explicit scoring"""
    
    def __init__(self, host_organism: str = "ecoli"):
        """
        Initialize optimizer with specified host organism
        
        Args:
            host_organism: 'ecoli', 'human', or custom codon usage dict
        """
        if host_organism == "ecoli":
            self.codon_usage = CODON_USAGE_ECOLI
        elif host_organism == "human":
            self.codon_usage = CODON_USAGE_HUMAN
        else:
            self.codon_usage = host_organism
        
        # Pre-compute frequency maximums per amino acid
        self.max_frequencies = self._compute_max_frequencies()
        
    def _compute_max_frequencies(self) -> Dict[str, float]:
        """Compute max codon frequency for each amino acid"""
        max_freq = {}
        for aa, codons in AA_TO_CODONS.items():
            max_freq[aa] = max(self.codon_usage.get(codon, 1) for codon in codons)
        return max_freq
    
    def _translate_to_aa(self, seq: str) -> str:
        """Translate DNA sequence to amino acids"""
        return str(Seq(seq).translate())
    
    def _compute_cai(self, seq: str) -> float:
        """
        Compute Codon Adaptation Index (CAI)
        
        CAI = exp(sum(ln(w_i)) / n) where w_i = codon_freq / max_freq_for_aa
        Range: [0, 1], where 1 is optimal
        """
        if len(seq) % 3 != 0:
            return 0.0
        
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        
        # Get amino acids
        aa_seq = self._translate_to_aa(seq)
        
        # Skip stop codons
        if '*' in aa_seq:
            return 0.0
        
        weights = []
        for codon, aa in zip(codons, aa_seq):
            if aa == '*':
                continue
            max_freq = self.max_frequencies.get(aa, 1)
            # If codon usage data missing for this codon, use a small default
            codon_freq = self.codon_usage.get(codon, None)
            if codon_freq is None or max_freq <= 0:
                # assign a small weight to penalize unknown codons
                weight = 1e-3
            else:
                weight = float(codon_freq) / float(max_freq)
                # ensure weight is in (0, 1]
                if weight <= 0:
                    weight = 1e-6
                elif weight > 1.0:
                    weight = 1.0
            weights.append(weight)
        
        if not weights:
            return 0.0
        
        # Geometric mean with safe log for very small weights
        eps = 1e-6
        log_sum = sum(math.log(w if w > 0 else eps) for w in weights)
        cai = math.exp(log_sum / len(weights))
        # Sanity constrain
        if math.isnan(cai) or cai < 0:
            return 0.0
        return max(0.0, min(cai, 1.0))  # Cap between 0 and 1
    
    def _penalize_rare_codon_clusters(self, seq: str, window_size: int = 6, 
                                     rarity_threshold: float = 0.5) -> float:
        """
        Penalize sequences with clusters of rare codons
        
        A rare codon is one with frequency < rarity_threshold * max_freq_for_aa
        Clusters detected in sliding windows cause penalties
        
        Returns: penalty value (higher = worse)
        """
        if len(seq) % 3 != 0:
            return 1.0
        
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        aa_seq = self._translate_to_aa(seq)
        
        if '*' in aa_seq or len(codons) < window_size:
            return 0.0
        
        # Identify rare codons
        rare_codons = []
        for codon, aa in zip(codons, aa_seq):
            max_freq = self.max_frequencies.get(aa, 1)
            codon_freq = self.codon_usage.get(codon, 1)
            if codon_freq < rarity_threshold * max_freq:
                rare_codons.append(True)
            else:
                rare_codons.append(False)
        
        # Count clusters in sliding window
        penalty = 0.0
        cluster_weight = 0.1  # Each cluster costs 0.1
        
        for i in range(len(rare_codons) - window_size + 1):
            window = rare_codons[i:i+window_size]
            if sum(window) >= 2:  # 2 or more rare codons = cluster
                penalty += cluster_weight
        
        return penalty
    
    def _compute_rrna_structure_penalty(self, seq: str, window_size: int = 30) -> float:
        """
        Compute mRNA secondary structure penalty using minimum free energy (ΔG)
        
        Stable secondary structures (negative ΔG) slow translation
        Penalty = |min(ΔG)| * weight for overly stable structures
        
        Returns: penalty value
        """
        if len(seq) < window_size:
            return 0.0
        
        max_penalty = 0.0
        window_stride = window_size // 2
        
        for i in range(0, len(seq) - window_size + 1, window_stride):
            window = seq[i:i+window_size]
            try:
                # Compute minimum free energy
                structure, energy = RNA.fold(window)
                # Penalize overly stable structures (energy < -2 kcal/mol)
                if energy < -2.0:
                    penalty = abs(energy) * 0.05  # Weight factor
                    max_penalty = max(max_penalty, penalty)
            except:
                continue
        
        return max_penalty
    
    def _count_codon_changes(self, original: str, mutated: str) -> int:
        """Count number of codon changes between sequences"""
        if len(original) % 3 != 0 or len(mutated) % 3 != 0:
            return 0
        
        changes = 0
        for i in range(0, len(original), 3):
            if original[i:i+3] != mutated[i:i+3]:
                changes += 1
        return changes
    
    def _compute_mutation_penalty(self, original: str, mutated: str, 
                                 penalty_weight: float = 0.05) -> float:
        """
        Prefer conservative edits with fewer codon changes
        
        Returns: penalty value
        """
        changes = self._count_codon_changes(original, mutated)
        return changes * penalty_weight
    
    def score_variant(self, original_seq: str, variant_seq: str, 
                     has_restriction_site: bool = False,
                     weights: Dict[str, float] = None) -> Dict[str, float]:
        """
        Score a silent mutation variant using the full scoring algorithm
        
        Args:
            original_seq: Original DNA sequence
            variant_seq: Mutated DNA sequence
            has_restriction_site: Whether variant contains desired restriction site
            weights: Custom scoring weights {w1, w2, w3, w4, w5}
                    Default: w1=1.0, w2=0.5, w3=0.3, w4=0.1, w5=1.0
        
        Returns: Dict with breakdown of all scoring components
        """
        # Default weights
        if weights is None:
            weights = {
                'w_cai': 1.0,           # CAI contribution
                'w_rare': 0.5,          # Rare codon penalty
                'w_structure': 0.3,     # mRNA structure penalty
                'w_mutation': 0.1,      # Minimal change preference
                'w_restriction': 1.0    # Restriction site bonus
            }
        
        # Verify sequences are synonymous (same amino acids)
        orig_aa = self._translate_to_aa(original_seq)
        var_aa = self._translate_to_aa(variant_seq)
        
        if orig_aa != var_aa or '*' in orig_aa:
            return {
                'cai': 0.0,
                'rare_penalty': 10.0,
                'structure_penalty': 10.0,
                'mutation_penalty': 10.0,
                'restriction_bonus': 0.0,
                'total_score': -100.0,
                'valid': False
            }
        
        # Step 2: CAI score
        cai = self._compute_cai(variant_seq)
        
        # Step 3: Rare codon cluster penalty
        rare_penalty = self._penalize_rare_codon_clusters(variant_seq)
        
        # Step 4: mRNA secondary structure penalty
        structure_penalty = self._compute_rrna_structure_penalty(variant_seq)
        
        # Step 5: Mutation penalty
        mutation_penalty = self._compute_mutation_penalty(original_seq, variant_seq)
        
        # Step 6: Restriction site bonus
        restriction_bonus = 1.0 if has_restriction_site else 0.0
        
        # Final score (higher is better)
        total_score = (
            weights['w_cai'] * cai
            - weights['w_rare'] * rare_penalty
            - weights['w_structure'] * structure_penalty
            - weights['w_mutation'] * mutation_penalty
            + weights['w_restriction'] * restriction_bonus
        )
        
        return {
            'cai': round(cai, 4),
            'rare_penalty': round(rare_penalty, 4),
            'structure_penalty': round(structure_penalty, 4),
            'mutation_penalty': round(mutation_penalty, 4),
            'restriction_bonus': round(restriction_bonus, 4),
            'total_score': round(total_score, 4),
            'valid': True,
            'codon_changes': self._count_codon_changes(original_seq, variant_seq)
        }
    
    def rank_variants(self, original_seq: str, variants: List[str], 
                     restriction_sites: List[bool] = None,
                     weights: Dict[str, float] = None) -> List[Dict]:
        """
        Score and rank multiple variants
        
        Args:
            original_seq: Original sequence
            variants: List of variant sequences
            restriction_sites: List of booleans for restriction site presence
            weights: Custom scoring weights
        
        Returns: Ranked list of scored variants
        """
        if restriction_sites is None:
            restriction_sites = [False] * len(variants)
        
        scored_variants = []
        for i, (variant, has_site) in enumerate(zip(variants, restriction_sites)):
            score = self.score_variant(original_seq, variant, has_site, weights)
            score['variant_index'] = i
            score['sequence'] = variant
            scored_variants.append(score)
        
        # Sort by total_score (highest first)
        ranked = sorted(scored_variants, key=lambda x: x['total_score'], reverse=True)
        
        # Add rank
        for rank, variant in enumerate(ranked, 1):
            variant['rank'] = rank
        
        return ranked
