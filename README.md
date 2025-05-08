# Yeast Transposable Element Analyzer

This tool is undergoing development.

## Features

- Identifies all Ty1-Ty5 transposable elements in the yeast genome
- Determines whether each TE is active or inactive
- For inactive TEs, identifies the likely cause of inactivation:
  - Truncation
  - Critical mutations
  - Epigenetic silencing
  - Recombination
- Analyzes nearby genomic features
- Generates comprehensive reports

## Prerequisites

On DRAC:
- StdEnv/2020 gcc/9.3.0 scipy-stack python blast+ muscle


## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/yeast-te-analyzer.git
cd yeast-te-analyzer
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Download TE consensus sequences:
```bash
python download_consensus.py
```

## Usage

1. Download the S. cerevisiae reference genome from NCBI:
   - Genome: GCF_000146045.2_R64_genomic.fna.gz

2. Run the analysis:
```bash
python yeast_te_analyzer.py \
  --genome GCF_000146045.2_R64_genomic.fna.gz \
  --transcriptome GCF_000146045.2_R64_rna.fna.gz \
  --output results
```

## Output

The script generates a tab-separated file (`results_te_analysis.tsv`) containing:
- TE family (Ty1-Ty5)
- Chromosome location
- Start and end positions
- Strand
- Activity status
- Inactivation reason (if applicable)
- Length
- Nearby genomic features
