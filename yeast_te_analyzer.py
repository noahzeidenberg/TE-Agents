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
        self.genome_sequences = {}
        self.te_annotations = []
        self.gene_annotations = {}  # Store gene locations from GFF
        self.te_consensus = self._load_te_consensus()
        self._setup_blast_db()
        self._load_gff_annotations()
        
    def _setup_blast_db(self):
        """Set up BLAST database for the genome"""
        print("Setting up BLAST database...")
        self.blast_db = "temp_genome_db"
        subprocess.run([
            "makeblastdb",
            "-in", self.genome_file,
            "-dbtype", "nucl",
            "-out", self.blast_db
        ])

    def _load_te_consensus(self) -> Dict[TEFamily, str]:
        """Load consensus sequences for each TE family from NCBI"""
        print("Fetching TE consensus sequences from NCBI...")
        
        # NCBI accession numbers for S. cerevisiae TE consensus sequences
        consensus_accessions = {
            TEFamily.TY1: "M18706.1",  # Ty1 consensus
            TEFamily.TY2: "X03840.1",  # Ty2 consensus
            TEFamily.TY3: "M34549.1",  # Ty3 consensus
            TEFamily.TY4: "X67284.1",  # Ty4 consensus
            TEFamily.TY5: "U19263.1"   # Ty5 consensus
        }
        
        consensus = {}
        for family, accession in consensus_accessions.items():
            try:
                # Fetch sequence from NCBI
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
                response = requests.get(url)
                response.raise_for_status()
                
                # Create temporary file for the consensus sequence
                with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
                    f.write(response.text)
                    consensus[family] = f.name
                    print(f"Retrieved consensus for {family.value}: {accession}")
                    
            except requests.exceptions.RequestException as e:
                print(f"Error fetching consensus for {family.value}: {e}")
                continue
        
        if not consensus:
            raise ValueError("Failed to retrieve any consensus sequences from NCBI")
        
        return consensus

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

    def _analyze_te_status(self, te: TEAnnotation):
        """Determine if a TE is active or inactive and why"""
        sequence = self.genome_sequences[te.chromosome][te.start:te.end]
        
        # Check for truncation
        if len(sequence) < 0.9 * len(self.te_consensus[te.family]):
            te.status = TEStatus.INACTIVE_TRUNCATED
            te.inactivation_reason = f"Truncated: {len(sequence)}bp vs consensus {len(self.te_consensus[te.family])}bp"
            return

        # Check for mutations in key regions
        mutations = self._find_mutations(sequence, self.te_consensus[te.family])
        if self._has_critical_mutations(mutations):
            te.status = TEStatus.INACTIVE_MUTATED
            te.inactivation_reason = "Critical mutations in functional regions"
            te.mutations = mutations
            return

        # Check for silencing marks
        if self._check_silencing_marks(te):
            te.status = TEStatus.INACTIVE_SILENCED
            te.inactivation_reason = "Epigenetically silenced"
            return

        # Check for recombination
        if self._check_recombination(te):
            te.status = TEStatus.INACTIVE_RECOMBINED
            te.inactivation_reason = "Inactivated by recombination"
            return

        te.status = TEStatus.ACTIVE

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

    def _check_silencing_marks(self, te: TEAnnotation) -> bool:
        """Check for epigenetic silencing marks"""
        # Look for known silencing patterns in flanking regions
        window = 500
        chromosome = self.genome_sequences[te.chromosome]
        upstream = chromosome[max(0, te.start - window):te.start]
        downstream = chromosome[te.end:min(len(chromosome), te.end + window)]
        
        # Check for tRNA genes nearby (known to attract silencing)
        if 'tRNA' in str(self.nearby_features):
            return True
            
        # Check for telomeric regions
        if te.start < 10000 or te.end > len(chromosome) - 10000:
            return True
            
        return False

    def _check_recombination(self, te: TEAnnotation) -> bool:
        """Check for recombination events"""
        # Look for solo LTRs or truncated elements
        if len(te.sequence) < 300:  # Typical LTR length
            return True
            
        # Check for hybrid elements
        other_tes = [t for t in self.te_annotations 
                    if t.chromosome == te.chromosome and 
                    abs(t.start - te.end) < 1000]
        
        if other_tes:
            return True
            
        return False

    def _analyze_nearby_features(self, te: TEAnnotation):
        """Analyze genomic features near the TE using GFF annotations"""
        window = 1000  # 1kb window
        features = {}
        
        if te.chromosome in self.gene_annotations:
            nearby_genes = []
            for gene in self.gene_annotations[te.chromosome]:
                # Check if gene is within window of TE
                if (abs(gene['start'] - te.end) < window or 
                    abs(gene['end'] - te.start) < window):
                    nearby_genes.append({
                        'id': gene['id'],
                        'type': gene['type'],
                        'distance': min(
                            abs(gene['start'] - te.end),
                            abs(gene['end'] - te.start)
                        )
                    })
            if nearby_genes:
                features['nearby_genes'] = nearby_genes
        
        te.nearby_features = features

    def analyze(self):
        """Main analysis function"""
        self._load_genome()
        self._identify_te_locations()
        
        for te in self.te_annotations:
            self._analyze_te_status(te)
            self._analyze_nearby_features(te)

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

    def __del__(self):
        """Cleanup temporary files"""
        # Remove BLAST database files
        for ext in ['.nhr', '.nin', '.nsq']:
            try:
                os.remove(self.blast_db + ext)
            except (OSError, AttributeError):
                pass
            
        # Remove temporary consensus files
        if hasattr(self, 'te_consensus'):
            for temp_file in self.te_consensus.values():
                try:
                    os.remove(temp_file)
                except OSError:
                    pass

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