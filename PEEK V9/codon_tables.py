# Host-specific codon usage tables for CAI calculation
# Data sourced from Kazusa Codon Usage Database
# Values represent frequency per 1000 codons

# E. coli K12 codon usage (highly expressed genes)
ECOLI_CODON_USAGE = {
    'TTT': 22.0, 'TTC': 16.5, 'TTA': 14.0, 'TTG': 13.4,
    'CTT': 12.2, 'CTC': 10.5, 'CTA': 4.0, 'CTG': 50.7,
    'ATT': 30.3, 'ATC': 24.1, 'ATA': 5.0, 'ATG': 27.1,
    'GTT': 26.1, 'GTC': 15.0, 'GTA': 11.0, 'GTG': 25.6,
    'TAT': 16.5, 'TAC': 12.2, 'TAA': 2.0, 'TAG': 0.3,
    'CAT': 12.8, 'CAC': 9.4, 'CAA': 15.1, 'CAG': 28.8,
    'AAT': 18.3, 'AAC': 21.4, 'AAA': 34.6, 'AAG': 11.5,
    'GAT': 32.1, 'GAC': 19.2, 'GAA': 39.1, 'GAG': 18.3,
    'TCT': 8.6, 'TCC': 8.8, 'TCA': 7.5, 'TCG': 8.7,
    'CCT': 7.2, 'CCC': 5.6, 'CCA': 8.5, 'CCG': 22.6,
    'ACT': 9.0, 'ACC': 22.9, 'ACA': 7.6, 'ACG': 14.4,
    'GCT': 15.4, 'GCC': 25.3, 'GCA': 20.5, 'GCG': 32.4,
    'TGT': 5.2, 'TGC': 6.3, 'TGA': 1.0, 'TGG': 15.2,
    'CGT': 20.7, 'CGC': 21.4, 'CGA': 3.7, 'CGG': 5.7,
    'AGT': 9.1, 'AGC': 15.8, 'AGA': 2.4, 'AGG': 1.5,
    'GGT': 24.7, 'GGC': 28.1, 'GGA': 8.4, 'GGG': 11.2,
}

# Saccharomyces cerevisiae (Yeast) codon usage
YEAST_CODON_USAGE = {
    'TTT': 26.1, 'TTC': 18.4, 'TTA': 26.2, 'TTG': 27.2,
    'CTT': 12.3, 'CTC': 5.4, 'CTA': 13.4, 'CTG': 10.5,
    'ATT': 30.3, 'ATC': 17.2, 'ATA': 17.8, 'ATG': 20.9,
    'GTT': 22.1, 'GTC': 11.8, 'GTA': 11.8, 'GTG': 10.8,
    'TAT': 18.8, 'TAC': 14.8, 'TAA': 1.1, 'TAG': 0.5,
    'CAT': 13.7, 'CAC': 7.8, 'CAA': 27.3, 'CAG': 12.2,
    'AAT': 35.7, 'AAC': 24.8, 'AAA': 41.9, 'AAG': 30.8,
    'GAT': 37.6, 'GAC': 20.2, 'GAA': 45.6, 'GAG': 19.2,
    'TCT': 23.5, 'TCC': 14.2, 'TCA': 18.7, 'TCG': 8.6,
    'CCT': 13.5, 'CCC': 6.8, 'CCA': 18.3, 'CCG': 5.3,
    'ACT': 20.3, 'ACC': 12.7, 'ACA': 17.8, 'ACG': 8.0,
    'GCT': 21.2, 'GCC': 12.6, 'GCA': 16.2, 'GCG': 6.2,
    'TGT': 8.1, 'TGC': 4.8, 'TGA': 0.7, 'TGG': 10.4,
    'CGT': 6.4, 'CGC': 2.6, 'CGA': 3.0, 'CGG': 1.7,
    'AGT': 14.2, 'AGC': 9.8, 'AGA': 21.3, 'AGG': 9.2,
    'GGT': 23.9, 'GGC': 9.8, 'GGA': 10.9, 'GGG': 6.0,
}

# Human (Homo sapiens) codon usage
HUMAN_CODON_USAGE = {
    'TTT': 17.6, 'TTC': 20.3, 'TTA': 7.7, 'TTG': 12.9,
    'CTT': 13.2, 'CTC': 19.6, 'CTA': 7.2, 'CTG': 39.6,
    'ATT': 16.0, 'ATC': 20.8, 'ATA': 7.5, 'ATG': 22.0,
    'GTT': 11.0, 'GTC': 14.5, 'GTA': 7.1, 'GTG': 28.1,
    'TAT': 12.2, 'TAC': 15.3, 'TAA': 1.0, 'TAG': 0.8,
    'CAT': 10.9, 'CAC': 15.1, 'CAA': 12.3, 'CAG': 34.2,
    'AAT': 17.0, 'AAC': 19.1, 'AAA': 24.4, 'AAG': 31.9,
    'GAT': 21.8, 'GAC': 25.1, 'GAA': 29.0, 'GAG': 39.6,
    'TCT': 15.2, 'TCC': 17.7, 'TCA': 12.2, 'TCG': 4.4,
    'CCT': 17.5, 'CCC': 19.8, 'CCA': 16.9, 'CCG': 6.9,
    'ACT': 13.1, 'ACC': 18.9, 'ACA': 15.1, 'ACG': 6.1,
    'GCT': 18.4, 'GCC': 27.7, 'GCA': 15.8, 'GCG': 7.4,
    'TGT': 10.6, 'TGC': 12.6, 'TGA': 1.6, 'TGG': 13.2,
    'CGT': 4.5, 'CGC': 10.4, 'CGA': 6.2, 'CGG': 11.4,
    'AGT': 12.1, 'AGC': 19.5, 'AGA': 12.2, 'AGG': 12.0,
    'GGT': 10.8, 'GGC': 22.2, 'GGA': 16.5, 'GGG': 16.5,
}

# CHO (Chinese Hamster Ovary) cells - similar to human but adjusted
CHO_CODON_USAGE = {
    'TTT': 16.9, 'TTC': 21.8, 'TTA': 6.7, 'TTG': 12.0,
    'CTT': 12.5, 'CTC': 20.5, 'CTA': 6.8, 'CTG': 41.2,
    'ATT': 15.2, 'ATC': 22.4, 'ATA': 6.8, 'ATG': 22.5,
    'GTT': 10.5, 'GTC': 15.2, 'GTA': 6.8, 'GTG': 29.5,
    'TAT': 11.5, 'TAC': 16.2, 'TAA': 0.9, 'TAG': 0.7,
    'CAT': 10.2, 'CAC': 16.0, 'CAA': 11.8, 'CAG': 35.8,
    'AAT': 16.2, 'AAC': 20.5, 'AAA': 23.5, 'AAG': 33.2,
    'GAT': 20.8, 'GAC': 26.5, 'GAA': 28.2, 'GAG': 41.2,
    'TCT': 14.5, 'TCC': 18.8, 'TCA': 11.5, 'TCG': 4.8,
    'CCT': 16.8, 'CCC': 20.8, 'CCA': 16.2, 'CCG': 7.2,
    'ACT': 12.5, 'ACC': 20.2, 'ACA': 14.5, 'ACG': 6.5,
    'GCT': 17.8, 'GCC': 29.2, 'GCA': 15.2, 'GCG': 7.8,
    'TGT': 10.2, 'TGC': 13.2, 'TGA': 1.5, 'TGG': 13.5,
    'CGT': 4.2, 'CGC': 11.2, 'CGA': 6.0, 'CGG': 12.0,
    'AGT': 11.5, 'AGC': 20.2, 'AGA': 11.8, 'AGG': 11.5,
    'GGT': 10.2, 'GGC': 23.5, 'GGA': 16.2, 'GGG': 17.0,
}

# Insect cells (Sf9/Sf21) codon usage
INSECT_CODON_USAGE = {
    'TTT': 20.5, 'TTC': 22.8, 'TTA': 11.2, 'TTG': 18.5,
    'CTT': 14.8, 'CTC': 15.2, 'CTA': 8.5, 'CTG': 28.2,
    'ATT': 22.5, 'ATC': 21.8, 'ATA': 10.2, 'ATG': 22.8,
    'GTT': 15.2, 'GTC': 16.5, 'GTA': 9.8, 'GTG': 24.5,
    'TAT': 14.5, 'TAC': 16.8, 'TAA': 1.2, 'TAG': 0.8,
    'CAT': 12.2, 'CAC': 14.5, 'CAA': 18.5, 'CAG': 25.8,
    'AAT': 22.5, 'AAC': 21.2, 'AAA': 32.5, 'AAG': 28.8,
    'GAT': 26.5, 'GAC': 24.2, 'GAA': 34.5, 'GAG': 32.2,
    'TCT': 16.8, 'TCC': 16.2, 'TCA': 14.5, 'TCG': 6.2,
    'CCT': 15.8, 'CCC': 14.2, 'CCA': 18.2, 'CCG': 7.5,
    'ACT': 14.8, 'ACC': 18.2, 'ACA': 16.5, 'ACG': 7.2,
    'GCT': 18.5, 'GCC': 24.5, 'GCA': 17.2, 'GCG': 8.5,
    'TGT': 9.8, 'TGC': 11.2, 'TGA': 1.2, 'TGG': 12.5,
    'CGT': 6.2, 'CGC': 9.5, 'CGA': 7.2, 'CGG': 9.8,
    'AGT': 13.5, 'AGC': 17.2, 'AGA': 14.5, 'AGG': 10.8,
    'GGT': 14.2, 'GGC': 20.5, 'GGA': 18.2, 'GGG': 14.5,
}

# Map host names to codon tables
CODON_TABLES = {
    'ecoli': ECOLI_CODON_USAGE,
    'yeast': YEAST_CODON_USAGE,
    'human': HUMAN_CODON_USAGE,
    'cho': CHO_CODON_USAGE,
    'insect': INSECT_CODON_USAGE,
}

# Host display names
HOST_NAMES = {
    'ecoli': 'E. coli K12',
    'yeast': 'S. cerevisiae (Yeast)',
    'human': 'Human (H. sapiens)',
    'cho': 'CHO Cells',
    'insect': 'Insect Cells (Sf9)',
}

# Codon to amino acid mapping
CODON_TO_AA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Amino acid to codons mapping
AA_TO_CODONS = {}
for codon, aa in CODON_TO_AA.items():
    if aa not in AA_TO_CODONS:
        AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)


def get_codon_table(host):
    """Get codon usage table for specified host"""
    return CODON_TABLES.get(host, ECOLI_CODON_USAGE)


def get_host_names():
    """Get dictionary of host IDs to display names"""
    return HOST_NAMES


def get_synonymous_codons(codon):
    """Get all synonymous codons for a given codon"""
    aa = CODON_TO_AA.get(codon)
    if not aa or aa == '*':
        return []
    return [c for c in AA_TO_CODONS[aa] if c != codon]
