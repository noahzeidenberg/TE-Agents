import os
import argparse
import gzip
from dataclasses import dataclass
from enum import Enum
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
    INACTIVE_TRUNCATED = "inactive_truncated"
    INACTIVE_MUTATED = "inactive_mutated"
    INACTIVE_SILENCED = "inactive_silenced"
    INACTIVE_RECOMBINED = "inactive_recombined"
    UNKNOWN = "unknown"

@dataclass
class TEAnnotation:
    """Stores information about a TE instance"""
    family: TEFamily
    chromosome: str
    start: int
    end: int
    strand: str
    status: TEStatus
    inactivation_reason: Optional[str] = None
    sequence: Optional[str] = None
    mutations: Optional[List[str]] = None
    nearby_features: Optional[Dict[str, str]] = None

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
        self.nearby_features = {}  # Initialize nearby_features as empty dict
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
        """Load consensus sequences for each TE family from NCBI"""
        print("Fetching TE consensus sequences from NCBI...")
        
        # Define NCBI accession numbers for each TE family
        te_accessions = {
            TEFamily.TY1: "M18706.1",
            TEFamily.TY2: "X03840.1",
            TEFamily.TY3: "M34549.1",
            TEFamily.TY4: "X67284.1",
            TEFamily.TY5: "U19263.1"
        }
        
        # Create a temporary directory for consensus files
        os.makedirs("temp_consensus", exist_ok=True)
        
        for family, accession in te_accessions.items():
            # Fetch sequence from NCBI
            consensus_file = f"temp_consensus/{family.value}.fasta"
            
            try:
                # Use Entrez to fetch the sequence
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype="fasta",
                    retmode="text"
                )
                
                # Save to temporary file
                with open(consensus_file, 'w') as f:
                    f.write(handle.read())
                
                self.te_consensus[family] = consensus_file
                print(f"Retrieved consensus for {family.value}: {accession}")
                
            except Exception as e:
                print(f"Error fetching consensus for {family.value}: {e}")
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

    def _identify_te_locations(self):
        """Identify TE locations using BLAST"""
        print("Identifying TE locations...")
        
        # First, get the actual chromosome IDs from the genome
        valid_chromosomes = set(self.genome_sequences.keys())
        print(f"Valid chromosome IDs: {', '.join(valid_chromosomes)}")
        
        for family in TEFamily:
            consensus_file = self.te_consensus[family]
            
            # Run BLAST
            output_file = f"temp_blast_{family.value}.xml"
            blastn_cline = NcbiblastnCommandline(
                query=consensus_file,
                db=self.blast_db,
                outfmt=5,
                out=output_file,
                word_size=11,
                evalue=1e-10,
                dust="no"
            )
            stdout, stderr = blastn_cline()
            
            # Parse BLAST results
            with open(output_file) as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            # Filter hits by identity and length
                            identity = (hsp.identities / hsp.align_length) * 100
                            if identity >= 80 and hsp.align_length >= 100:
                                # Extract the actual chromosome ID from the BLAST title
                                # The title format might be like "gnl|BL_ORD_ID|3 NC_001135.5"
                                blast_title = alignment.title.split()
                                
                                # Try to find a matching chromosome ID
                                chrom_id = None
                                for title_part in blast_title:
                                    if title_part in valid_chromosomes:
                                        chrom_id = title_part
                                        break
                                
                                if chrom_id is None:
                                    print(f"Warning: Could not find valid chromosome ID in BLAST title: {alignment.title}")
                                    continue
                                
                                te = TEAnnotation(
                                    family=family,
                                    chromosome=chrom_id,
                                    start=hsp.sbjct_start,
                                    end=hsp.sbjct_end,
                                    strand='+' if hsp.sbjct_start < hsp.sbjct_end else '-',
                                    status=TEStatus.UNKNOWN,
                                    sequence=hsp.sbjct
                                )
                                self.te_annotations.append(te)
            
            os.remove(output_file)

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
        
        # Store nearby features for this TE
        self.nearby_features[te] = nearby
        return nearby

    def _check_silencing_marks(self, te: TEAnnotation) -> bool:
        """Check for silencing marks near the TE"""
        if te not in self.nearby_features:
            self._analyze_nearby_features(te)
            
        nearby = self.nearby_features.get(te, [])
        
        # Check for tRNA genes or other silencing-associated features
        for feature in nearby:
            if 'tRNA' in feature.get('type', ''):
                return True
            # Add other silencing mark checks here
        
        return False

    def _analyze_te_status(self, te: TEAnnotation):
        """Analyze the status of a TE"""
        # First analyze nearby features
        self._analyze_nearby_features(te)
        
        # Then check various status conditions
        if self._check_truncation(te):
            te.status = TEStatus.INACTIVE_TRUNCATED
        elif self._check_mutations(te):
            te.status = TEStatus.INACTIVE_MUTATED
        elif self._check_silencing_marks(te):
            te.status = TEStatus.INACTIVE_SILENCED
        else:
            te.status = TEStatus.ACTIVE

    def _check_truncation(self, te: TEAnnotation) -> bool:
        """Check for truncation"""
        sequence = self.genome_sequences[te.chromosome][te.start:te.end]
        if len(sequence) < 0.9 * len(self.te_consensus[te.family]):
            te.inactivation_reason = f"Truncated: {len(sequence)}bp vs consensus {len(self.te_consensus[te.family])}bp"
            return True
        return False

    def _check_mutations(self, te: TEAnnotation) -> bool:
        """Check for mutations"""
        sequence = self.genome_sequences[te.chromosome][te.start:te.end]
        consensus = self.te_consensus[te.family]
        mutations = self._find_mutations(sequence, consensus)
        if self._has_critical_mutations(mutations):
            te.inactivation_reason = "Critical mutations in functional regions"
            te.mutations = mutations
            return True
        return False

    def _check_recombination(self, te: TEAnnotation) -> bool:
        """Check for recombination"""
        if len(te.sequence) < 300:  # Typical LTR length
            te.inactivation_reason = "Inactivated by recombination"
            return True
        return False

    def _find_mutations(self, sequence: str, consensus: str) -> List[str]:
        """Find mutations compared to consensus sequence"""
        mutations = []
        
        # Align sequences
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(f">query\n{sequence}\n>consensus\n{consensus}\n")
            temp_file = f.name
            
        # Run muscle alignment
        output_file = f"{temp_file}.aln"
        subprocess.run([
            "muscle",
            "-in", temp_file,
            "-out", output_file
        ])
        
        # Parse alignment and find mutations
        aligned_seqs = {}
        for record in SeqIO.parse(output_file, "fasta"):
            aligned_seqs[record.id] = str(record.seq)
            
        query_seq = aligned_seqs["query"]
        consensus_seq = aligned_seqs["consensus"]
        
        for i, (q, c) in enumerate(zip(query_seq, consensus_seq)):
            if q != c and q != '-' and c != '-':
                mutations.append(f"{c}{i+1}{q}")
                
        os.remove(temp_file)
        os.remove(output_file)
        return mutations

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
        
        # Set up BLAST database
        print("Setting up BLAST database...")
        makeblastdb_cline = NcbiblastdbCommandline(
            dbtype="nucl",
            input_file=self.genome_file,
            out=self.blast_db
        )
        makeblastdb_cline()
        
        # Identify and analyze TEs
        self._identify_te_locations()
        
        # Analyze each TE
        for te in self.te_annotations:
            self._analyze_te_status(te)
        
        # Clean up temporary files
        self._cleanup()

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
        """Generate a detailed report of findings"""
        data = []
        for te in self.te_annotations:
            data.append({
                'Family': te.family.value,
                'Chromosome': te.chromosome,
                'Start': te.start,
                'End': te.end,
                'Strand': te.strand,
                'Status': te.status.value,
                'Inactivation_Reason': te.inactivation_reason,
                'Length': te.end - te.start,
                'Nearby_Features': str(te.nearby_features)
            })
        return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Analyze TEs in S. cerevisiae genome')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('--gff', required=True, help='Path to GFF annotation file')
    parser.add_argument('--output', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    analyzer = YeastTEAnalyzer(args.genome, args.gff)
    analyzer.analyze()
    
    # Generate and save report
    report = analyzer.generate_report()
    report.to_csv(f"{args.output}_te_analysis.tsv", sep='\t', index=False)
    
    # Print summary statistics
    print("\nAnalysis Summary:")
    print(f"Total TEs found: {len(analyzer.te_annotations)}")
    status_counts = report['Status'].value_counts()
    print("\nTE Status Distribution:")
    for status, count in status_counts.items():
        print(f"{status}: {count}")

if __name__ == "__main__":
    main() 