"""Curated motif library for quick insertion in primer design (amino acid sequences).

Provides MOTIF_LIBRARY dict and a helper get_motif_library().
"""

MOTIF_LIBRARY = {
    # Purification tags
    'his6': {
        'name': '6xHis',
        'sequence': 'HHHHHH',
        'category': 'purification',
        'description': 'Hexahistidine tag for IMAC purification (His6).',
        'preferred_position': 'C-terminal',
        'notes': 'Small, widely used; may affect metal-binding proteins.',
        'risk_internal': False,
        'risk_terminal': False
    },
    'flag': {
        'name': 'FLAG',
        'sequence': 'DYKDDDDK',
        'category': 'purification',
        'description': 'FLAG epitope tag for affinity detection and purification.',
        'preferred_position': 'N-terminal',
        'notes': 'Acidic; sometimes removed by proteases.',
        'risk_internal': True,
        'risk_terminal': False
    },
    'ha': {
        'name': 'HA',
        'sequence': 'YPYDVPDYA',
        'category': 'purification',
        'description': 'Hemagglutinin (HA) epitope tag.',
        'preferred_position': 'C-terminal',
        'notes': 'Small and commonly used for detection.',
        'risk_internal': True,
        'risk_terminal': False
    },
    'myc': {
        'name': 'Myc',
        'sequence': 'EQKLISEEDL',
        'category': 'purification',
        'description': 'c-Myc epitope tag for detection.',
        'preferred_position': 'C-terminal',
        'notes': '',
        'risk_internal': True,
        'risk_terminal': False
    },
    'strep': {
        'name': 'Strep',
        'sequence': 'WSHPQFEK',
        'category': 'purification',
        'description': 'Strep-tag II for gentle affinity purification.',
        'preferred_position': 'N-terminal',
        'notes': 'Very specific binding to Strep-Tactin.',
        'risk_internal': False,
        'risk_terminal': False
    },

    # Cleavage motifs
    'tev': {
        'name': 'TEV',
        'sequence': 'ENLYFQG',
        'category': 'cleavage',
        'description': 'Tobacco Etch Virus protease recognition site (TEV).',
        'preferred_position': 'internal',
        'notes': 'Efficient and specific; cleavage occurs between Q and G.',
        'risk_internal': False,
        'risk_terminal': False
    },
    'hrv3c': {
        'name': 'HRV 3C',
        'sequence': 'LEVLFQGP',
        'category': 'cleavage',
        'description': 'Human rhinovirus 3C protease site (HRV 3C).',
        'preferred_position': 'internal',
        'notes': 'Often used as an alternative to TEV.',
        'risk_internal': False,
        'risk_terminal': False
    },
    'thrombin': {
        'name': 'Thrombin',
        'sequence': 'LVPRGS',
        'category': 'cleavage',
        'description': 'Thrombin protease recognition site.',
        'preferred_position': 'internal',
        'notes': 'Can sometimes cleave at cryptic sites.',
        'risk_internal': True,
        'risk_terminal': False
    },
    'factor_xa': {
        'name': 'Factor Xa',
        'sequence': 'IEGR',
        'category': 'cleavage',
        'description': 'Factor Xa protease recognition sequence.',
        'preferred_position': 'internal',
        'notes': 'Short sequence; may appear elsewhere by chance.',
        'risk_internal': True,
        'risk_terminal': False
    },
    'enterokinase': {
        'name': 'Enterokinase',
        'sequence': 'DDDDK',
        'category': 'cleavage',
        'description': 'Enterokinase cleavage site (Asp-Asp-Asp-Asp-Lys).',
        'preferred_position': 'internal',
        'notes': 'Four Asp residues - acidic patch may influence local structure.',
        'risk_internal': True,
        'risk_terminal': False
    },

    # Linkers
    'gs_short': {
        'name': 'GS (short)',
        'sequence': 'GGGGS',
        'category': 'linker',
        'description': 'Flexible glycine-serine linker (short).',
        'preferred_position': 'internal',
        'notes': 'Useful to separate domains with flexibility.',
        'risk_internal': False,
        'risk_terminal': False
    },
    'gs_medium': {
        'name': 'GS (medium)',
        'sequence': 'GGGGSGGGGS',
        'category': 'linker',
        'description': 'Flexible glycine-serine linker (medium).',
        'preferred_position': 'internal',
        'notes': '',
        'risk_internal': False,
        'risk_terminal': False
    },
    'gs_long': {
        'name': 'GS (long)',
        'sequence': 'GGGGSGGGGSGGGGS',
        'category': 'linker',
        'description': 'Flexible glycine-serine linker (long).',
        'preferred_position': 'internal',
        'notes': 'Long linkers increase separation but may add disorder.',
        'risk_internal': False,
        'risk_terminal': False
    },
    'eaak': {
        'name': 'EAAAK',
        'sequence': 'EAAAK',
        'category': 'linker',
        'description': 'Rigid helical linker motif (EAAAK).',
        'preferred_position': 'internal',
        'notes': 'Imposes rigidity between domains.',
        'risk_internal': False,
        'risk_terminal': False
    },
    'pro_rich': {
        'name': 'Proline-rich',
        'sequence': 'PPPP',
        'category': 'linker',
        'description': 'Proline-rich short motif for rigidity.',
        'preferred_position': 'internal',
        'notes': 'May introduce kinks and structural constraints.',
        'risk_internal': True,
        'risk_terminal': False
    },

    # Combination motifs
    'his_tev': {
        'name': 'His6-TEV',
        'sequence': 'HHHHHHENLYFQG',
        'category': 'combination',
        'description': '6xHis followed by TEV site for removable His-tag.',
        'preferred_position': 'N-terminal',
        'notes': 'Commonly used for purification then tag removal.',
        'risk_internal': False,
        'risk_terminal': False
    },
    'flag_tev': {
        'name': 'FLAG-TEV',
        'sequence': 'DYKDDDDKENLYFQG',
        'category': 'combination',
        'description': 'FLAG tag followed by TEV site.',
        'preferred_position': 'N-terminal',
        'notes': '',
        'risk_internal': True,
        'risk_terminal': False
    }
}


def get_motif_library():
    """Return the motif library dictionary.

    The consumer (template) may iterate categories and entries as needed.
    """
    return MOTIF_LIBRARY
