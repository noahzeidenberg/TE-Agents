import os
import argparse
import gzip
from dataclasses import dataclass
from enum import Enum, auto
from typing import List, Dict, Optional, Set, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import requests
from collections import defaultdict
import re
import subprocess
import tempfile
from Bio import Entrez
import shutil
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import json
from pathlib import Path
from itertools import product

Entrez.email = "nzeidenb@uoguelph.ca"

class TEFamily(Enum):
    """Yeast TE families"""
    TY1 = "Ty1"
    TY2 = "Ty2"
    TY3 = "Ty3"
    TY4 = "Ty4"
    TY5 = "Ty5"

class TEStatus(Enum):
    """TE activity status"""
    ACTIVE = "active"
    TRUNCATED = "truncated"
    MUTATED = "mutated"
    SILENCED = "silenced"
    UNKNOWN = "unknown"

class InactivationReason(Enum):
    TRUNCATION = "Truncation"
    MUTATION = "Mutation"
    SILENCING = "Silencing"
    NONE = "None"

@dataclass
class TEAnnotation:
    family: TEFamily
    chromosome: str
    start: int
    end: int
    strand: str
    sequence: str
    _status: TEStatus = TEStatus.UNKNOWN
    _inactivation_reason: InactivationReason = InactivationReason.NONE
    
    def __hash__(self):
        # Hash only the immutable identifying attributes
        return hash((self.family, self.chromosome, self.start, self.end, self.strand))
    
    @property
    def status(self) -> TEStatus:
        """Get the TE status"""
        return self._status
    
    @status.setter
    def status(self, value: TEStatus):
        """Set the TE status"""
        object.__setattr__(self, '_status', value)
    
    @property
    def inactivation_reason(self) -> InactivationReason:
        return self._inactivation_reason
    
    @inactivation_reason.setter
    def inactivation_reason(self, value: InactivationReason):
        object.__setattr__(self, '_inactivation_reason', value)
    
    def get_id(self) -> str:
        """Generate a unique identifier for this TE"""
        return f"{self.chromosome}_{self.start}_{self.end}_{self.family.value}"
    
    def __eq__(self, other):
        if not isinstance(other, TEAnnotation):
            return NotImplemented
        return (self.family == other.family and
                self.chromosome == other.chromosome and
                self.start == other.start and
                self.end == other.end and
                self.strand == other.strand)

class SpacedSeed:
    def __init__(self, pattern: str):
        """
        Initialize a spaced seed with a pattern string.
        '1' represents match position, '0' represents don't care position
        Example: '11011' requires matches at positions 0,1,3,4
        """
        self.pattern = pattern
        self.weight = pattern.count('1')  # Number of required matches
        self.length = len(pattern)
        self.match_positions = [i for i, c in enumerate(pattern) if c == '1']

class YeastTEAnalyzer:
    def __init__(self, genome_file: str, gff_file: str):
        """Initialize the analyzer with genome and GFF annotation files"""
        self.genome_file = genome_file
        self.gff_file = gff_file
        self.output_dir = "repeatmasker_output"
        self.te_annotations = []
        self.gff_features = []
        self.nearby_features = {}
        
        # Create output directory
        Path(self.output_dir).mkdir(exist_ok=True)
        
        # Initialize consensus sequences first
        self._load_te_consensus()
        
    def _setup_blast_db(self):
        """Set up BLAST database for the genome"""
        print("Setting up BLAST database...")
        subprocess.run([
            "makeblastdb",
            "-in", self.genome_file,
            "-dbtype", "nucl",
            "-out", self.blast_db
        ])

    def _load_te_consensus(self):
        """Load consensus sequences from local files"""
        print("Loading TE consensus sequences from local files...")
        
        # Define local file paths for each TE family
        te_files = {
            TEFamily.TY1: "consensus/S288C_YALWdelta1_YALWdelta1_genomic.fsa",
            TEFamily.TY2: "consensus/S288C_YGRWTy2-2_genomic.fsa",
            TEFamily.TY3: "consensus/S288C_YGRWTy3-1_genomic.fsa",
            TEFamily.TY4: "consensus/S288C_YHLWTy4-1_genomic.fsa",
            TEFamily.TY5: "consensus/S288C_YCLWTy5-1_genomic.fsa"
        }
        
        for family, file_path in te_files.items():
            try:
                # Verify file exists
                if not os.path.exists(file_path):
                    print(f"Warning: Consensus file not found for {family.value}: {file_path}")
                    continue
                    
                # Read and validate the consensus sequence
                with open(file_path) as f:
                    sequence_data = f.read()
                    
                # Store the file path directly since it's already in FASTA format
                self.te_consensus[family] = file_path
                
                # Debug: Print sequence info
                seq = ''.join(line.strip() for line in sequence_data.split('\n')[1:])
                print(f"Loaded consensus for {family.value}")
                print(f"Sequence length: {len(seq)}")
                print(f"First 50 bp: {seq[:50]}")
                
            except Exception as e:
                print(f"Error loading consensus for {family.value}: {e}")
                continue

    def _load_genome(self):
        """Load the genome sequences"""
        print("Loading genome sequences...")
        with open(self.genome_file, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                # Print the first few characters of each header to help debug
                print(f"Found sequence with ID: {record.id}")
                self.genome_sequences[record.id] = str(record.seq)

    def _load_gff_annotations(self):
        """Load gene annotations from GFF file"""
        print("Loading GFF annotations...")
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                if fields[2] in ['gene', 'tRNA_gene']:
                    chromosome = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                    gene_id = attributes.get('ID', 'unknown')
                    
                    if chromosome not in self.gene_annotations:
                        self.gene_annotations[chromosome] = []
                    self.gene_annotations[chromosome].append({
                        'id': gene_id,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'type': fields[2]
                    })

    def _generate_spaced_seeds(self) -> List[SpacedSeed]:
        """Generate optimized spaced seeds for TE detection"""
        # These patterns are optimized for TE detection based on literature
        # Multiple seeds increase sensitivity
        seed_patterns = [
            '1111100001111',  # Optimized for ~70-80% similarity
            '11011011000011',  # Good for finding TEs with internal deletions
            '111010010100111'  # Effective for finding degraded LTRs
        ]
        return [SpacedSeed(pattern) for pattern in seed_patterns]

    def _generate_seed_hits(self, sequence: str, seed: SpacedSeed) -> Set[str]:
        """Generate all seed hits from a sequence using a spaced seed pattern"""
        hits = set()
        for i in range(len(sequence) - seed.length + 1):
            # Extract characters at match positions only
            hit = ''.join(sequence[i + pos] for pos in seed.match_positions)
            hits.add(hit)
        return hits

    def _run_repeatmasker(self):
        """Run RepeatMasker on the genome"""
        print("Running RepeatMasker analysis...")
        
        # Create a custom library from our consensus sequences
        library_file = "consensus/yeast_te_library.fa"
        with open(library_file, 'w') as lib:
            for te_file in Path("consensus").glob("S288C_*.fsa"):
                with open(te_file) as f:
                    lib.write(f.read() + "\n")
        
        # Run RepeatMasker with custom settings
        cmd = [
            "RepeatMasker",
            "-lib", library_file,          # Use custom library
            "-species", "Saccharomyces",   # Specify species
            "-pa", "4",                    # Use 4 parallel processes
            "-gff",                        # Output GFF format
            "-nolow",                      # Don't mask low complexity
            "-no_is",                      # Don't mask bacterial insertion elements
            "-dir", self.output_dir,       # Output directory
            self.genome_file               # Input genome
        ]
        
        print("Running command:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        
        # Parse RepeatMasker output
        self._parse_repeatmasker_output()

    def _parse_repeatmasker_output(self):
        """Parse RepeatMasker output files"""
        # Parse the .out file which contains detailed annotations
        out_file = Path(self.output_dir) / f"{Path(self.genome_file).name}.out"
        
        print(f"Parsing RepeatMasker output: {out_file}")
        
        with open(out_file) as f:
            # Skip header lines
            for _ in range(3):
                next(f)
            
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 15:
                    te = TEAnnotation(
                        family=self._get_te_family(fields[10]),
                        chromosome=fields[4],
                        start=int(fields[5]),
                        end=int(fields[6]),
                        strand=fields[8],
                        sequence="",  # We'll fill this in later
                        _status=self._determine_status(fields),
                        _inactivation_reason=self._determine_inactivation(fields)
                    )
                    self.te_annotations.append(te)

    def _get_te_family(self, name: str) -> TEFamily:
        """Map RepeatMasker names to TE families"""
        name = name.upper()
        if "TY1" in name:
            return TEFamily.TY1
        elif "TY2" in name:
            return TEFamily.TY2
        elif "TY3" in name:
            return TEFamily.TY3
        elif "TY4" in name:
            return TEFamily.TY4
        elif "TY5" in name:
            return TEFamily.TY5
        else:
            return None

    def _determine_status(self, fields) -> TEStatus:
        """Determine TE status from RepeatMasker output"""
        div_pct = float(fields[1])  # Divergence percentage
        del_pct = float(fields[2])  # Deletion percentage
        ins_pct = float(fields[3])  # Insertion percentage
        
        if div_pct < 5 and del_pct < 5 and ins_pct < 5:
            return TEStatus.ACTIVE
        elif del_pct > 20:
            return TEStatus.TRUNCATED
        elif div_pct > 20:
            return TEStatus.MUTATED
        else:
            return TEStatus.UNKNOWN

    def _determine_inactivation(self, fields) -> InactivationReason:
        """Determine inactivation reason from RepeatMasker output"""
        div_pct = float(fields[1])
        del_pct = float(fields[2])
        
        if del_pct > 20:
            return InactivationReason.TRUNCATION
        elif div_pct > 20:
            return InactivationReason.MUTATION
        else:
            return InactivationReason.NONE

    def analyze(self):
        """Run the complete TE analysis"""
        # Run RepeatMasker
        self._run_repeatmasker()
        
        # Load GFF annotations for feature analysis
        self._load_gff_annotations()
        
        # Analyze nearby features and silencing
        for te in self.te_annotations:
            if self._check_silencing_marks(te):
                te.status = TEStatus.SILENCED
                te.inactivation_reason = InactivationReason.SILENCING
        
        # Print summary
        self._print_summary()

    def _print_summary(self):
        """Print analysis summary"""
        total = len(self.te_annotations)
        status_counts = {}
        family_counts = {}
        inactivation_counts = {}
        
        for te in self.te_annotations:
            status_counts[te.status.value] = status_counts.get(te.status.value, 0) + 1
            family_counts[te.family.value] = family_counts.get(te.family.value, 0) + 1
            inactivation_counts[te.inactivation_reason.value] = inactivation_counts.get(te.inactivation_reason.value, 0) + 1
        
        print("\nAnalysis Summary:")
        print(f"Total TEs found: {total}")
        
        print("\nTE Family Distribution:")
        for family, count in family_counts.items():
            print(f"{family}: {count}")
        
        print("\nTE Status Distribution:")
        for status, count in status_counts.items():
            print(f"{status}: {count}")
        
        print("\nInactivation Reason Distribution:")
        for reason, count in inactivation_counts.items():
            print(f"{reason}: {count}")

    def _cleanup(self):
        """Clean up temporary files"""
        # Remove BLAST database files
        for ext in [".nhr", ".nin", ".nsq"]:
            try:
                os.remove(self.blast_db + ext)
            except OSError:
                pass
        
        # Remove temporary consensus files
        try:
            shutil.rmtree("temp_consensus")
        except OSError:
            pass

    def __del__(self):
        """Destructor to ensure cleanup"""
        self._cleanup()

    def generate_report(self) -> pd.DataFrame:
        """Generate a report of TE analysis results"""
        if not self.te_annotations:
            print("No TEs found in the analysis.")
            return pd.DataFrame(columns=[
                'Family', 'Chromosome', 'Start', 'End', 'Strand',
                'Status', 'Inactivation_Reason', 'Length', 'Nearby_Features'
            ])
        
        report_data = []
        for te in self.te_annotations:
            report_data.append({
                'Family': te.family.value,
                'Chromosome': te.chromosome,
                'Start': te.start,
                'End': te.end,
                'Strand': te.strand,
                'Status': te.status.value,
                'Inactivation_Reason': te.inactivation_reason.value,
                'Length': len(te.sequence),
                'Nearby_Features': len(self.nearby_features.get(te.get_id(), []))
            })
        
        return pd.DataFrame(report_data)

    def save_results(self, output_prefix: str):
        """Save analysis results to files"""
        # Generate and save report
        report = self.generate_report()
        report.to_csv(f"{output_prefix}_report.csv", index=False)
        
        # Save summary statistics
        summary = {
            'total_tes': len(self.te_annotations),
            'active_tes': len([te for te in self.te_annotations if te.status == TEStatus.ACTIVE]),
            'truncated_tes': len([te for te in self.te_annotations if te.inactivation_reason == InactivationReason.TRUNCATION]),
            'mutated_tes': len([te for te in self.te_annotations if te.inactivation_reason == InactivationReason.MUTATION]),
            'silenced_tes': len([te for te in self.te_annotations if te.inactivation_reason == InactivationReason.SILENCING])
        }
        
        with open(f"{output_prefix}_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)

def main():
    parser = argparse.ArgumentParser(description='Analyze TEs in yeast genome')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('--gff', required=True, help='Path to GFF annotation file')
    parser.add_argument('--output', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    analyzer = YeastTEAnalyzer(args.genome, args.gff)
    analyzer.analyze()
    
    # Generate and save report
    report = analyzer.generate_report()
    if not report.empty:
        report.to_csv(f"{args.output}_report.csv", index=False)
        status_counts = report['Status'].value_counts()
        print("\nAnalysis Summary:")
        print(f"Total TEs found: {len(report)}")
        print("\nTE Status Distribution:")
        for status, count in status_counts.items():
            print(f"{status}: {count}")
    else:
        print("\nNo TEs were found in the analysis.")

if __name__ == "__main__":
    main() 