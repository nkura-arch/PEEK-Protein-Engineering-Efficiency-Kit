# Overview

PEEK (Protein Engineering Efficiency Kit) v1.2 is a Flask-based web application designed for molecular biology and protein engineering workflows. The application provides tools for:
1. **Silent Mutations + CAI** - Finding codon-preserving mutations to introduce restriction enzyme sites with integrated algorithmic codon optimization scoring
2. **Primer Design** - Designing primers for enzyme insertion, amino acid insertion, or Site-Directed Mutagenesis (SDM) with nearest-neighbor Tm calculations

It leverages BioPython for sequence manipulation, ViennaRNA for mRNA secondary structure analysis, and restriction enzyme analysis, making it a specialized tool for researchers working with DNA sequences and protein engineering.

# User Preferences

Preferred communication style: Simple, everyday language.
Design approach: Scientific rigor, reproducibility - NO AI/LLM-based predictions, only algorithmic approaches.
Feature philosophy: All features must be removable if requested - treated as experimental additions.

# System Architecture

## Frontend Architecture
- **Framework**: Traditional server-rendered Flask templates with Bootstrap 5 dark theme
- **Styling**: Bootstrap CSS framework with custom CSS for molecular biology-specific UI components
- **JavaScript**: Vanilla JavaScript for client-side interactions (copy-to-clipboard, form validation, DNA sequence validation, motif library selection)
- **Template Engine**: Jinja2 templates with a base template providing consistent navigation and layout

## Backend Architecture
- **Framework**: Flask web framework with modular routing structure
- **Core Logic**: BioPython-based sequence analysis and manipulation
- **Session Management**: Flask sessions with configurable secret key
- **Request Handling**: Standard Flask request/response pattern with form processing

## Key Functional Components
- **Codon Table Management**: Host-specific codon usage tables from Kazusa database (E. coli, Yeast, Human, CHO, Insect cells)
- **Restriction Enzyme Database**: Integration with BioPython's restriction enzyme collection (CommOnly enzymes)
- **Silent Mutation Engine**: Algorithm for generating codon-preserving mutations within specified tolerance windows
- **Codon Optimizer Engine**: CAI-based optimization with rare codon cluster detection, GC content analysis, homopolymer scoring
- **Primer Design Engine**: Tools for calculating melting temperatures and optimizing primer sequences for PCR
- **Motif Library**: Curated collection of purification tags, cleavage sites, and linkers for functional insertions

## Codon Optimization Scoring Formula (ViennaRNA-based)
```
Score = 1.0·CAI - 0.5·RarePenalty - 0.3·StructurePenalty - 0.1·MutationPenalty + 1.0·RestrictionBonus
```
- **CAI**: Codon Adaptation Index (0-1, geometric mean of codon frequencies)
- **RarePenalty**: Clusters of codons <50% max frequency in sliding windows
- **StructurePenalty**: mRNA secondary structure penalty using ViennaRNA (ΔG < -2 kcal/mol)
- **MutationPenalty**: Minimal change preference (0.05 per codon change)
- **RestrictionBonus**: +1.0 if target restriction site is present

## Data Processing Pipeline
- **Sequence Validation**: DNA sequence input validation and uppercase normalization
- **Translation Services**: DNA to amino acid translation with stop codon detection
- **Window-based Analysis**: Codon-aligned sequence windows for mutation analysis
- **Enzyme Site Detection**: Pattern matching for restriction enzyme recognition sites
- **CAI Calculation**: Geometric mean of relative adaptiveness values per host

## Application Structure
- **app.py**: Main Flask application with routes
- **silent_mutation_engine.py**: Silent mutation detection algorithm with scoring integration
- **codon_optimizer.py**: ViennaRNA-based codon optimization with CAI, structure, and penalty scoring
- **motif_library.py**: Curated functional motif sequences
- **templates/**: Jinja2 HTML templates
- **static/**: CSS and JavaScript assets

# External Dependencies

## Core Scientific Libraries
- **BioPython**: Primary dependency for sequence analysis, restriction enzymes, and genetic code tables
- **ViennaRNA**: mRNA secondary structure prediction using minimum free energy (ΔG) calculations
- **Pandas**: Data manipulation and CSV processing for results export

## Web Framework Stack
- **Flask**: Web application framework and routing
- **Bootstrap 5**: Frontend CSS framework with dark theme variant
- **Font Awesome**: Icon library for UI elements

## Development Dependencies
- **Python Standard Library**: CSV, IO, itertools, math for data processing
- **Environment Variables**: Session secret key configuration via environment

## Client-side Resources
- **CDN Dependencies**: Bootstrap CSS and Font Awesome served from CDNs
- **No Database**: Application operates without persistent data storage, processing sequences in memory

# Recent Changes (January 2026)

- **PEEK v1.2 integrated** - Complete ViennaRNA-based codon optimization:
  - Select host organism (E. coli or Human) directly in Silent Mutations form
  - Enable/disable algorithmic scoring with checkbox
  - Results show CAI scores, rare codon penalties, mRNA structure penalties, and total P-scores
  - Top 20 preview with toggle to full results table
  - CSV export with all scoring metrics
- Removed standalone Codon Optimizer page (functionality fully integrated into Silent Mutations + CAI)
- Updated navigation to "Silent Mutations + CAI"
- Added Motif Library with 15+ curated sequences (His-tags, FLAG, TEV, linkers, etc.)
- Enhanced Primer Design page with motif library dropdown for quick insertions
- **Added Site-Directed Mutagenesis (SDM) mode** to Primer Design page:
  - New third mode option for Q5 mutagenesis-style silent mutation primer design
  - Takes original and mutated nucleotide sequences from Silent Mutations results
  - Designs forward primer containing mutation + downstream sequence
  - Designs reverse primer as reverse-complement upstream sequence
  - Returns primer lengths, GC%, ΔTm, and recommended annealing temperature
  - Compatible with NEB Q5 Site-Directed Mutagenesis Kit (E0554)
  - Feature is experimental and can be removed if requested
