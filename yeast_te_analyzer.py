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
        self.blast_db = "temp_genome_db"
        self.genome_sequences = {}
        self.te_annotations = []
        self.te_consensus = {}
        self.gff_features = []  # Store GFF features
        self.nearby_features = {}  # Will use TE ID as key instead of TE object
        self.gene_annotations = {}  # Store gene locations from GFF
        
        # Initialize consensus sequences first
        self._load_te_consensus()
        
        self._setup_blast_db()
        self._load_gff_annotations()
        
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

    def _identify_te_locations(self):
        """Identify TE locations using BLAST"""
        print("Identifying TE locations...")
        
        # First, get the actual chromosome IDs from the genome
        valid_chromosomes = set(self.genome_sequences.keys())
        print(f"Valid chromosome IDs: {', '.join(valid_chromosomes)}")
        
        for family in TEFamily:
            if family not in self.te_consensus:
                print(f"Skipping {family.value}: no consensus sequence available")
                continue
            
            consensus_file = self.te_consensus[family]
            
            # Create a BLAST database with parameters optimized for TEs
            output_file = f"temp_blast_{family.value}.txt"
            blastn_cline = NcbiblastnCommandline(
                query=consensus_file,
                db=self.blast_db,
                outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                out=output_file,
                word_size=11,
                evalue=1e-5,        # More permissive e-value for short sequences
                dust='no',
                soft_masking='false',
                task='blastn',
                gapopen=5,
                gapextend=2,
                reward=2,
                penalty=-3,
                max_target_seqs=1000
            )
            
            print(f"\nRunning BLAST search for {family.value}...")
            stdout, stderr = blastn_cline()
            
            # Parse tab-delimited BLAST results
            hit_count = 0
            if os.path.exists(output_file):
                with open(output_file) as f:
                    content = f.read()
                    if not content.strip():
                        print(f"No BLAST hits found for {family.value}")
                        print("Debug: Checking BLAST parameters and consensus sequence...")
                        continue
                        
                    for line in content.strip().split('\n'):
                        fields = line.split('\t')
                        if len(fields) >= 12:
                            pident = float(fields[2])
                            length = int(fields[3])
                            
                            # Adjust filtering based on consensus length
                            min_length = 50  # Minimum length for delta elements
                            min_identity = 80  # Minimum identity percentage
                            
                            if length >= min_length and pident >= min_identity:
                                chrom_id = fields[1]
                                if chrom_id in valid_chromosomes:
                                    start = min(int(fields[8]), int(fields[9]))
                                    end = max(int(fields[8]), int(fields[9]))
                                    
                                    # Get the actual sequence from the genome
                                    sequence = self.genome_sequences[chrom_id][start-1:end]
                                    
                                    te = TEAnnotation(
                                        family=family,
                                        chromosome=chrom_id,
                                        start=start,
                                        end=end,
                                        strand='+' if int(fields[8]) < int(fields[9]) else '-',
                                        sequence=sequence
                                    )
                                    self.te_annotations.append(te)
                                    hit_count += 1
            
            print(f"Found {hit_count} hits for {family.value}")
            os.remove(output_file)

    def _merge_overlapping_hits(self):
        """Merge overlapping TE hits to avoid redundant annotations"""
        print("Merging overlapping hits...")
        merged = []
        
        # Sort by chromosome and start position
        sorted_tes = sorted(
            self.te_annotations,
            key=lambda x: (x.chromosome, x.start)
        )
        
        if not sorted_tes:
            return
        
        current = sorted_tes[0]
        
        for next_te in sorted_tes[1:]:
            if (current.chromosome == next_te.chromosome and 
                next_te.start <= current.end + 100):  # Allow small gaps
                # Merge the TEs
                current = TEAnnotation(
                    family=current.family,
                    chromosome=current.chromosome,
                    start=min(current.start, next_te.start),
                    end=max(current.end, next_te.end),
                    strand=current.strand,
                    sequence=current.sequence  # Keep the longer sequence
                    if len(current.sequence) > len(next_te.sequence)
                    else next_te.sequence
                )
            else:
                merged.append(current)
                current = next_te
        
        merged.append(current)
        self.te_annotations = merged

    def _analyze_nearby_features(self, te: TEAnnotation) -> List[str]:
        """Analyze genomic features near the TE"""
        nearby = []
        window = 1000  # Look 1kb upstream and downstream
        
        for feature in self.gff_features:
            if feature['chromosome'] != te.chromosome:
                continue
                
            # Check if feature is within window of TE
            if (abs(feature['start'] - te.start) <= window or 
                abs(feature['end'] - te.end) <= window):
                nearby.append(feature)
        
        # Store nearby features using TE ID as key
        self.nearby_features[te.get_id()] = nearby
        return nearby

    def _check_silencing_marks(self, te: TEAnnotation) -> bool:
        """Check for silencing marks near the TE"""
        te_id = te.get_id()
        if te_id not in self.nearby_features:
            self._analyze_nearby_features(te)
            
        nearby = self.nearby_features.get(te_id, [])
        
        # Check for tRNA genes or other silencing-associated features within 1kb
        for feature in nearby:
            feature_type = feature.get('type', '').lower()
            if any(mark in feature_type for mark in ['trna', 'telomer', 'centromer', 'silenc']):
                return True
        
        return False

    def _analyze_te_status(self, te: TEAnnotation):
        """Analyze the status of a TE"""
        # First check truncation as it's the most common inactivation mechanism
        if self._check_truncation(te):
            te.status = TEStatus.TRUNCATED
            te.inactivation_reason = InactivationReason.TRUNCATION
            return
        
        # Then check for mutations
        if self._check_mutations(te):
            te.status = TEStatus.MUTATED
            te.inactivation_reason = InactivationReason.MUTATION
            return
        
        # Finally check for silencing marks
        if self._check_silencing_marks(te):
            te.status = TEStatus.SILENCED
            te.inactivation_reason = InactivationReason.SILENCING
            return
        
        # If none of the above, consider it potentially active
        te.status = TEStatus.ACTIVE
        te.inactivation_reason = InactivationReason.NONE

    def _check_truncation(self, te: TEAnnotation) -> bool:
        """Check if TE is truncated based on expected lengths"""
        # Define expected lengths for each TE family (in base pairs)
        expected_lengths = {
            TEFamily.TY1: 5900,  # ~5.9 kb
            TEFamily.TY2: 5900,  # ~5.9 kb
            TEFamily.TY3: 5400,  # ~5.4 kb
            TEFamily.TY4: 6200,  # ~6.2 kb
            TEFamily.TY5: 5400   # ~5.4 kb
        }
        
        # Allow for some variation (e.g., 80% of expected length)
        min_length = expected_lengths[te.family] * 0.8
        actual_length = te.end - te.start + 1
        
        return actual_length < min_length

    def _run_muscle_alignment(self, sequence1_file: str, sequence2_file: str, output_file: str) -> None:
        """Run MUSCLE alignment using subprocess"""
        try:
            result = subprocess.run(
                [
                    'muscle',
                    '-in', sequence1_file,
                    '-in2', sequence2_file,
                    '-out', output_file,
                    '-quiet',
                    '-maxiters', '1'
                ],
                capture_output=True,
                text=True,
                check=True
            )
            if result.stderr:
                print(f"MUSCLE stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"MUSCLE alignment failed: {e}")
            raise

    def _check_mutations(self, te: TEAnnotation) -> bool:
        """Check for inactivating mutations in the TE sequence"""
        consensus_file = self.te_consensus[te.family]
        temp_sequence_file = f"temp_te_{te.get_id()}.fasta"
        output_file = f"temp_alignment_{te.get_id()}.fasta"
        
        try:
            # Write TE sequence to temporary file
            with open(temp_sequence_file, 'w') as f:
                f.write(f">TE_{te.get_id()}\n{te.sequence}\n")
            
            # Run MUSCLE alignment
            self._run_muscle_alignment(consensus_file, temp_sequence_file, output_file)
            
            # Analyze alignment for mutations
            with open(output_file) as f:
                alignments = list(SeqIO.parse(f, "fasta"))
                if len(alignments) == 2:
                    consensus_seq = str(alignments[0].seq)
                    te_seq = str(alignments[1].seq)
                    
                    # Calculate sequence identity
                    matches = sum(1 for a, b in zip(consensus_seq, te_seq) if a == b and a != '-' and b != '-')
                    aligned_positions = sum(1 for a, b in zip(consensus_seq, te_seq) if a != '-' and b != '-')
                    
                    if aligned_positions == 0:
                        return True
                    
                    identity = (matches / aligned_positions) * 100
                    
                    # Check for frameshift mutations (gaps of non-3 length)
                    gaps = [len(gap) for gap in te_seq.split('-') if gap]
                    has_frameshift = any(gap % 3 != 0 for gap in gaps)
                    
                    # Check for premature stop codons
                    has_premature_stop = False
                    for i in range(0, len(te_seq) - 2, 3):
                        codon = te_seq[i:i+3].upper()
                        if codon in ['TAA', 'TAG', 'TGA']:
                            has_premature_stop = True
                            break
                    
                    return identity < 90 or has_frameshift or has_premature_stop
            
            return False
        
        finally:
            # Clean up temporary files
            for file in [temp_sequence_file, output_file]:
                try:
                    Path(file).unlink(missing_ok=True)
                except Exception as e:
                    print(f"Warning: Could not remove temporary file {file}: {e}")

    def _check_recombination(self, te: TEAnnotation) -> bool:
        """Check for recombination"""
        if len(te.sequence) < 300:  # Typical LTR length
            te.inactivation_reason = InactivationReason.SILENCING
            return True
        return False

    def _has_critical_mutations(self, mutations: List[str]) -> bool:
        """Check if mutations affect critical regions"""
        critical_regions = {
            TEFamily.TY1: [
                (292, 1850),  # GAG region
                (1850, 5561)  # POL region
            ],
            TEFamily.TY2: [
                (292, 1850),
                (1850, 5561)
            ],
            TEFamily.TY3: [
                (238, 1628),  # GAG
                (1628, 4952)  # POL
            ],
            TEFamily.TY4: [
                (271, 1674),
                (1674, 5376)
            ],
            TEFamily.TY5: [
                (280, 1713),
                (1713, 5369)
            ]
        }
        
        for mutation in mutations:
            match = re.match(r'[ATCG](\d+)[ATCG]', mutation)
            if match:
                pos = int(match.group(1))
                for family, regions in critical_regions.items():
                    for start, end in regions:
                        if start <= pos <= end:
                            return True
        return False

    def analyze(self):
        """Run the complete TE analysis"""
        # Load required data
        print("Loading GFF annotations...")
        self._load_gff_annotations()
        
        print("Loading genome sequences...")
        self._load_genome()
        
        # Configure BLAST database with optimized parameters
        print("Setting up BLAST database...")
        makeblastdb_cline = NcbimakeblastdbCommandline(
            dbtype="nucl",
            input_file=self.genome_file,
            out=self.blast_db,
            parse_seqids=True
        )
        makeblastdb_cline()
        
        # Identify TEs using spaced seeds
        self._identify_te_locations()
        
        # Merge overlapping hits
        self._merge_overlapping_hits()
        
        # Continue with status analysis
        print("\nAnalyzing TE status...")
        total = len(self.te_annotations)
        for i, te in enumerate(self.te_annotations, 1):
            self._analyze_te_status(te)
            if i % 10 == 0:  # Progress update every 10 TEs
                print(f"Processed {i}/{total} TEs...")
        
        # Print summary statistics
        status_counts = {}
        for te in self.te_annotations:
            status_counts[te.status.value] = status_counts.get(te.status.value, 0) + 1
        
        print("\nAnalysis Summary:")
        print(f"Total TEs found: {total}")
        print("\nTE Status Distribution:")
        for status, count in status_counts.items():
            print(f"{status}: {count}")

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