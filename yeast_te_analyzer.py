import os
import argparse
import gzip
from dataclasses import dataclass
from enum import Enum, auto
from typing import List, Dict, Optional, Set, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
import requests
from collections import defaultdict
import re
import subprocess
import tempfile
from Bio import Entrez
import shutil
import json
from pathlib import Path
from itertools import product
import sys
import logging

Entrez.email = "nzeidenb@uoguelph.ca"

# At the start of the script:
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('yeast_te_analyzer.log'),
        logging.StreamHandler()
    ]
)

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
        self.output_dir = "edta_output"
        self.te_annotations = []
        self.gff_features = []
        self.nearby_features = {}
        self.gene_annotations = defaultdict(list)
        
        # Create output directory
        Path(self.output_dir).mkdir(exist_ok=True)
        
        # Prepare curated library
        self._prepare_curated_library()

    def _prepare_curated_library(self):
        """Prepare curated TE library from consensus sequences"""
        print("Preparing curated TE library...")
        
        # Create consensus directory if it doesn't exist
        Path("consensus").mkdir(exist_ok=True)
        
        # Check if consensus files exist
        consensus_files = list(Path("consensus").glob("S288C_*.fsa"))
        if not consensus_files:
            print("Warning: No consensus sequences found in consensus directory")
            print("Please ensure S288C_*.fsa files are present in the consensus directory")
            raise FileNotFoundError("Missing consensus sequence files")
        
        library_file = "consensus/yeast_te_library.fasta"
        print(f"Creating library file: {library_file}")
        with open(library_file, 'w') as lib:
            for te_file in consensus_files:
                print(f"Adding sequences from: {te_file}")
                with open(te_file) as f:
                    content = f.read()
                    if not content.startswith('>'):
                        print(f"Warning: {te_file} may not be in FASTA format")
                    lib.write(content + "\n")
        
        print(f"Verifying library file...")
        with open(library_file) as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                raise ValueError(f"Library file {library_file} is not in proper FASTA format")
        
        self.curated_library = library_file

    def _run_edta(self):
        """Run RepeatMasker analysis"""
        print("Running RepeatMasker analysis...")
        
        # Get absolute paths and prepare files
        abs_genome = os.path.abspath(self.genome_file)
        abs_lib = os.path.abspath(self.curated_library)
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Get available resources
        threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
        
        # Use full path to RepeatMasker
        repeatmasker_path = os.path.expanduser("~/scratch/TE-Agents/tools_fresh/RepeatMasker/RepeatMasker")
        
        # Run RepeatMasker using only our curated library
        rm_cmd = [
            repeatmasker_path,
            "-pa", str(threads),
            "-lib", abs_lib,  # Use our curated yeast TE library
            "-gff",  # Output in GFF format
            "-dir", self.output_dir,
            "-nolow",  # Skip low complexity repeats
            "-no_is",  # Skip bacterial insertion elements
            "-engine", "crossmatch",  # Use crossmatch instead of ncbi/rmblast
            abs_genome
        ]
        
        print("Running RepeatMasker command:", " ".join(rm_cmd))
        try:
            subprocess.run(rm_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"RepeatMasker analysis failed: {e}")
            raise

    def _parse_edta_output(self):
        """Parse RepeatMasker output files"""
        # Update to parse RepeatMasker output instead of EDTA
        anno_file = os.path.join(self.output_dir, os.path.basename(self.genome_file) + ".out.gff")
        
        print(f"Parsing RepeatMasker annotations: {anno_file}")
        
        if not os.path.exists(anno_file):
            print(f"Warning: Annotation file not found at {anno_file}")
            return
        
        with open(anno_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                # Parse RepeatMasker GFF format
                attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                name = attributes.get('Target', '').split()[0]
                
                te_family = self._get_te_family(name)
                if te_family is None:
                    continue
                
                te = TEAnnotation(
                    family=te_family,
                    chromosome=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    sequence="",
                    _status=self._determine_status(attributes),
                    _inactivation_reason=self._determine_inactivation(attributes)
                )
                self.te_annotations.append(te)

    def _determine_status(self, attributes: Dict[str, str]) -> TEStatus:
        """Determine TE status from EDTA attributes"""
        # EDTA provides detailed classification
        classification = attributes.get('Classification', '').lower()
        identity = float(attributes.get('Identity', '0'))
        
        if 'intact' in classification and identity >= 90:
            return TEStatus.ACTIVE
        elif 'truncated' in classification:
            return TEStatus.TRUNCATED
        elif 'degraded' in classification or identity < 80:
            return TEStatus.MUTATED
        else:
            return TEStatus.UNKNOWN

    def analyze(self):
        """Run the complete TE analysis"""
        # Run EDTA
        self._run_edta()
        
        # Parse EDTA output
        self._parse_edta_output()
        
        # Load GFF annotations for feature analysis
        self._load_gff_annotations()
        
        # Check for silencing marks
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
        # Remove temporary consensus files
        try:
            shutil.rmtree("temp_consensus")
        except OSError:
            pass
        
        # Clean up EDTA temporary files
        try:
            for temp_file in Path(self.output_dir).glob("*.temp*"):
                temp_file.unlink()
        except OSError:
            pass
        
        # Clean up copied input files
        try:
            if hasattr(self, 'genome_file'):
                os.remove(os.path.basename(self.genome_file))
            if hasattr(self, 'gff_file'):
                os.remove(os.path.basename(self.gff_file))
            if hasattr(self, 'curated_library'):
                os.remove(os.path.basename(self.curated_library))
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

    def _get_te_family(self, name: str) -> TEFamily:
        """Map EDTA names to TE families"""
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

    def _determine_inactivation(self, attributes: Dict[str, str]) -> InactivationReason:
        """Determine inactivation reason from EDTA attributes"""
        # EDTA provides detailed classification
        classification = attributes.get('Classification', '').lower()
        
        if 'truncated' in classification:
            return InactivationReason.TRUNCATION
        elif 'degraded' in classification:
            return InactivationReason.MUTATION
        else:
            return InactivationReason.NONE

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

    def _check_silencing_marks(self, te: TEAnnotation) -> bool:
        """Check if a TE is silenced by nearby features"""
        # Implementation of the method to check for silencing marks
        # This is a placeholder and should be implemented based on your specific requirements
        return False  # Placeholder return, actual implementation needed

def main():
    parser = argparse.ArgumentParser(description='Analyze TEs in yeast genome')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('--gff', required=True, help='Path to GFF annotation file')
    parser.add_argument('--output', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    try:
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
    except Exception as e:
        print(f"Error during analysis: {e}")
        raise

if __name__ == "__main__":
    main() 