import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Set
from enum import Enum
import random
import multiprocessing as mp
from functools import partial
import concurrent.futures
import torch
import pickle
import gzip
from Bio import SeqIO

class TEType(Enum):
    LTR_RETROTRANSPOSON = 0
    SINE = 1
    LINE = 2
    DNA_TRANSPOSON = 3
    
    @property
    def name_str(self) -> str:
        return {
            0: "LTR_retrotransposon",
            1: "SINE",
            2: "LINE",
            3: "DNA_transposon"
        }[self.value]

@dataclass
class TEProperties:
    """Properties of a transposable element"""
    te_type: TEType
    transposition_rate: float  # Probability of transposition per time step
    copy_mechanism: str  # "cut_and_paste" or "copy_and_paste"
    target_site_bias: float  # Bias towards certain genomic features (0-1)
    silencing_sensitivity: float  # How sensitive to host silencing (0-1)
    fitness_impact: float  # Impact on host fitness (-1 to 1)
    # New properties from reference code
    length: int  # Length of the TE in base pairs
    death_rate: float  # Rate at which TE becomes inactive/dead
    excision_rate: float  # Rate at which TE excises itself (for DNA transposons)
    progeny_rate: float  # Number of copies produced per transposition
    # Inactivation probabilities
    initial_active_rate: float = 0.1  # Proportion of TEs that start active (typically 5-15%)
    truncation_rate: float = 0.4  # Probability of being truncated
    mutation_rate: float = 0.3  # Probability of being mutated
    recombination_rate: float = 0.2  # Probability of being inactivated by recombination
    silencing_rate: float = 0.1  # Probability of being epigenetically silenced
    preferred_insertion_sites: List[Tuple[int, int]] = None  # Preferred genomic regions for insertion
    epigenetic_sensitivity: float = 0.5  # How sensitive to epigenetic modifications
    host_tissue_specificity: Dict[str, float] = None  # Tissue-specific activity levels
    horizontal_transfer_rate: float = 0.0  # Rate of horizontal transfer between hosts
    regulatory_impact: float = 0.0  # Impact on nearby gene regulation
    structural_impact: float = 0.0  # Impact on genome structure (e.g., chromosomal rearrangements)
    age_dependent_activity: bool = False  # Whether activity changes with age
    host_species_specificity: Dict[str, float] = None  # Species-specific activity levels

class InactivationMechanism(Enum):
    """Types of TE inactivation mechanisms"""
    ACTIVE = "active"
    TRUNCATED = "truncated"
    MUTATED = "mutated"
    RECOMBINED = "recombined"
    SILENCED = "silenced"

class TransposableElement:
    """Represents a single transposable element"""
    def __init__(self, 
                 te_type: TEType,
                 location: int,
                 properties: TEProperties,
                 is_silenced: bool = False,
                 is_dead: bool = False,
                 inactivation_mechanism: InactivationMechanism = InactivationMechanism.ACTIVE,
                 truncation_length: Optional[int] = None):
        self.te_type = te_type
        self.location = location
        self.properties = properties
        self.is_silenced = is_silenced
        self.is_dead = is_dead
        self.age = 0  # Time steps since creation
        self.copies_made = 0  # Number of copies this TE has produced
        self.last_transposition = 0  # Time step of last transposition
        self.epigenetic_state = {}  # Current epigenetic state
        self.inactivation_mechanism = inactivation_mechanism
        self.truncation_length = truncation_length  # For truncated TEs, stores the remaining length
        
    def is_active(self) -> bool:
        """Check if the TE is currently active"""
        return (not self.is_dead and 
                not self.is_silenced and 
                self.inactivation_mechanism == InactivationMechanism.ACTIVE)

    def get_inactivation_details(self) -> Dict:
        """Get details about the TE's inactivation state"""
        return {
            "mechanism": self.inactivation_mechanism.value,
            "is_dead": self.is_dead,
            "is_silenced": self.is_silenced,
            "truncation_length": self.truncation_length if self.truncation_length else self.properties.length
        }
        
    def attempt_transposition(self, genome_size: int) -> Optional[int]:
        """Attempt to transpose to a new location"""
        if self.is_silenced or self.is_dead:
            return None
            
        # Check death rate
        if random.random() < self.properties.death_rate:
            self.is_dead = True
            return None
            
        # Check excision rate for DNA transposons
        if (self.properties.copy_mechanism == "cut_and_paste" and 
            random.random() < self.properties.excision_rate):
            return None
            
        # Age-dependent activity
        if self.properties.age_dependent_activity:
            activity_factor = max(0.0, 1.0 - (self.age / 1000))
        else:
            activity_factor = 1.0
            
        # Check transposition rate with modifiers
        if random.random() > self.properties.transposition_rate * activity_factor:
            return None
            
        # Determine new location with bias
        if self.properties.preferred_insertion_sites:
            # Choose from preferred sites with bias
            preferred_sites = self.properties.preferred_insertion_sites
            weights = [1.0 - self.properties.target_site_bias for _ in preferred_sites]
            chosen_site = random.choices(preferred_sites, weights=weights)[0]
            new_location = random.randint(chosen_site[0], chosen_site[1])
        else:
            # Random insertion with slight bias
            bias = self.properties.target_site_bias
            if random.random() < bias:
                # Bias towards certain regions (e.g., gene-rich regions)
                new_location = int(random.gauss(genome_size/2, genome_size/4))
            else:
                new_location = random.randint(0, genome_size - 1)
                
        # Ensure location is within bounds
        new_location = max(0, min(genome_size - 1, new_location))
        
        # Update transposition history
        self.last_transposition = self.age
        self.copies_made += 1
        
        return new_location
        
    def update_epigenetic_state(self, genome: 'GenomeEnvironment'):
        """Update TE's epigenetic state based on genome context"""
        # Check nearby epigenetic marks
        for offset in range(-100, 101):
            pos = self.location + offset
            if pos in genome.epigenetic_marks:
                self.epigenetic_state[pos] = genome.epigenetic_marks[pos]
                
    def get_activity_level(self, tissue: str = None, species: str = None) -> float:
        """Get current activity level based on context"""
        activity = 1.0
        
        # Apply tissue-specific activity
        if tissue and self.properties.host_tissue_specificity:
            activity *= self.properties.host_tissue_specificity.get(tissue, 1.0)
            
        # Apply species-specific activity
        if species and self.properties.host_species_specificity:
            activity *= self.properties.host_species_specificity.get(species, 1.0)
            
        # Apply epigenetic effects
        if self.epigenetic_state:
            active_marks = sum(1 for mark in self.epigenetic_state.values() 
                             if mark == "active_te")
            activity *= (1 + 0.1 * active_marks)
            
        return max(0.0, min(1.0, activity))

class GenomeEnvironment:
    """Represents the host genome environment"""
    def __init__(self, 
                 size: int,
                 gene_locations: List[Tuple[int, int]] = None,
                 regulatory_regions: List[Tuple[int, int]] = None,
                 epigenetic_marks: Dict[int, str] = None,
                 chromosome_structure: Dict[str, List[Tuple[int, int]]] = None,
                 pathway_network: 'PathwayNetwork' = None):
        self.size = size
        self.gene_locations = gene_locations or []
        self.regulatory_regions = regulatory_regions or []
        self.te_copies: List[TransposableElement] = []
        self.silencing_machinery = 0.5  # Host's TE silencing capability (0-1)
        # New properties from reference code
        self.epigenetic_marks = epigenetic_marks or {}  # Position -> mark type
        self.chromosome_structure = chromosome_structure or {}  # Chromosome -> regions
        self.recombination_rate = 0.01  # Rate of recombination events
        self.mutation_rate = 0.0001  # Rate of point mutations
        self.gene_density = len(self.gene_locations) / size  # Genes per base pair
        self.te_density = 0.0  # Updated in update_density()
        self.host_fitness = 1.0  # Current host fitness
        self.genome_stability = 1.0  # Measure of genome stability
        self.pathway_network = pathway_network
        self.gene_to_location = self._create_gene_location_map()
        
    def _create_gene_location_map(self) -> Dict[str, Tuple[int, int]]:
        """Create mapping of gene IDs to their locations"""
        return {f"gene_{i}": loc for i, loc in enumerate(self.gene_locations)}
        
    def get_affected_genes(self) -> Set[str]:
        """Get set of genes affected by TE insertions"""
        affected_genes = set()
        
        for te in self.te_copies:
            if te.is_active():
                # Check for direct gene disruption
                for gene_id, (start, end) in self.gene_to_location.items():
                    if start <= te.location <= end:
                        affected_genes.add(gene_id)
                        
                # Check for regulatory region effects
                for i, (start, end) in enumerate(self.regulatory_regions):
                    if start <= te.location <= end:
                        # Add genes within 10kb of regulatory region
                        for gene_id, (g_start, g_end) in self.gene_to_location.items():
                            if abs(g_start - end) < 10000 or abs(g_end - start) < 10000:
                                affected_genes.add(gene_id)
                                
        return affected_genes
        
    def add_te(self, te: TransposableElement):
        """Add a TE to the genome"""
        self.te_copies.append(te)
        
    def remove_te(self, te: TransposableElement):
        """Remove a TE from the genome"""
        if te in self.te_copies:
            self.te_copies.remove(te)
            
    def update_silencing(self):
        """Update silencing status of TEs based on host machinery"""
        for te in self.te_copies:
            if random.random() < self.silencing_machinery * te.properties.silencing_sensitivity:
                te.is_silenced = True
            else:
                te.is_silenced = False
                
    def calculate_fitness_impact(self) -> float:
        """Calculate overall fitness impact of TEs"""
        total_impact = 0
        for te in self.te_copies:
            if not te.is_silenced:
                total_impact += te.properties.fitness_impact
        return total_impact

    def update_density(self):
        """Update TE density in the genome"""
        self.te_density = len(self.te_copies) / self.size
        
    def update_epigenetic_marks(self):
        """Update epigenetic marks based on TE activity"""
        for te in self.te_copies:
            if not te.is_silenced:
                # Add epigenetic marks around active TEs
                for offset in range(-100, 101):
                    pos = te.location + offset
                    if 0 <= pos < self.size:
                        self.epigenetic_marks[pos] = "active_te"
                        
    def calculate_genome_stability(self) -> float:
        """Calculate genome stability based on TE distribution and activity"""
        stability = 1.0
        
        # Penalize high TE density
        if self.te_density > 0.5:
            stability *= (1 - (self.te_density - 0.5))
            
        # Penalize clustered TEs
        te_locations = sorted([te.location for te in self.te_copies])
        for i in range(len(te_locations) - 1):
            if te_locations[i + 1] - te_locations[i] < 1000:
                stability *= 0.95
                
        # Penalize active TEs in gene regions
        active_tes = [te for te in self.te_copies if not te.is_silenced]
        for te in active_tes:
            for gene_start, gene_end in self.gene_locations:
                if gene_start <= te.location <= gene_end:
                    stability *= 0.9
                    
        return max(0.0, min(1.0, stability))
        
    def update_fitness(self):
        """Update host fitness based on TE activity, genome stability, and pathway disruption"""
        # Base fitness from TE impact
        self.host_fitness = 1.0 + self.calculate_fitness_impact()
        
        # Adjust for genome stability
        self.genome_stability = self.calculate_genome_stability()
        self.host_fitness *= self.genome_stability
        
        # Adjust for pathway disruption if pathway network exists
        if self.pathway_network:
            affected_genes = self.get_affected_genes()
            self.pathway_network.calculate_pathway_disruption(affected_genes)
            pathway_impact = self.pathway_network.calculate_fitness_impact()
            self.host_fitness *= (1.0 - pathway_impact)
        
        # Ensure fitness stays in reasonable range
        self.host_fitness = max(0.0, min(1.0, self.host_fitness))
        
    def update(self):
        """Update all genome properties"""
        self.update_density()
        self.update_epigenetic_marks()
        self.update_silencing()
        self.update_fitness()

def _process_te_batch(te_batch, genome_size, tissue_type, species):
    """Process a batch of TEs in parallel"""
    results = []
    for te in te_batch:
        activity = te.get_activity_level(tissue_type, species)
        new_location = te.attempt_transposition(genome_size)
        results.append((te, new_location, activity))
    return results

class TESimulation:
    """Main simulation class with parallel processing and GPU acceleration"""
    def __init__(self, 
                 genome_size: int,
                 initial_te_count: int,
                 te_properties: Dict[TEType, TEProperties],
                 tissue_type: str = None,
                 species: str = None,
                 n_processes: int = None,
                 use_gpu: bool = False):
        self.genome = GenomeEnvironment(genome_size)
        self.te_properties = te_properties
        self.time_step = 0
        self.tissue_type = tissue_type
        self.species = species
        self.history = []
        self.n_processes = n_processes or mp.cpu_count()
        self.use_gpu = use_gpu and torch.cuda.is_available()
        self.device = torch.device("cuda" if self.use_gpu else "cpu")
        self.initialize_tes(initial_te_count)
        
    def save_checkpoint(self, filename: str):
        """Save simulation state to file"""
        state = {
            'time_step': self.time_step,
            'genome': self.genome,
            'history': self.history,
            'te_properties': self.te_properties
        }
        with open(filename, 'wb') as f:
            pickle.dump(state, f)
            
    def load_checkpoint(self, filename: str):
        """Load simulation state from file"""
        with open(filename, 'rb') as f:
            state = pickle.load(f)
        self.time_step = state['time_step']
        self.genome = state['genome']
        self.history = state['history']
        self.te_properties = state['te_properties']
        
    def initialize_tes(self, count: int):
        """Initialize TEs using vectorized operations with realistic inactive proportions"""
        te_types = np.random.choice(list(self.te_properties.keys()), size=count)
        locations = np.random.randint(0, self.genome.size, size=count)
        
        self.genome.te_copies = []
        
        for te_type, location in zip(te_types, locations):
            props = self.te_properties[te_type]
            
            # Determine if TE starts active based on initial_active_rate
            if random.random() < props.initial_active_rate:
                # Create active TE
                te = TransposableElement(
                    te_type=te_type,
                    location=location,
                    properties=props,
                    inactivation_mechanism=InactivationMechanism.ACTIVE
                )
            else:
                # Choose inactivation mechanism based on relative probabilities
                total_prob = (props.truncation_rate + props.mutation_rate + 
                            props.recombination_rate + props.silencing_rate)
                
                rand = random.random() * total_prob
                cumulative = 0
                
                if rand < (cumulative := props.truncation_rate):
                    # Create truncated TE
                    truncation_length = random.randint(1, props.length - 1)
                    te = TransposableElement(
                        te_type=te_type,
                        location=location,
                        properties=props,
                        is_dead=True,
                        inactivation_mechanism=InactivationMechanism.TRUNCATED,
                        truncation_length=truncation_length
                    )
                elif rand < (cumulative := cumulative + props.mutation_rate):
                    # Create mutated TE
                    te = TransposableElement(
                        te_type=te_type,
                        location=location,
                        properties=props,
                        is_dead=True,
                        inactivation_mechanism=InactivationMechanism.MUTATED
                    )
                elif rand < (cumulative := cumulative + props.recombination_rate):
                    # Create recombined TE
                    te = TransposableElement(
                        te_type=te_type,
                        location=location,
                        properties=props,
                        is_dead=True,
                        inactivation_mechanism=InactivationMechanism.RECOMBINED
                    )
                else:
                    # Create silenced TE
                    te = TransposableElement(
                        te_type=te_type,
                        location=location,
                        properties=props,
                        is_silenced=True,
                        inactivation_mechanism=InactivationMechanism.SILENCED
                    )
            
            self.genome.te_copies.append(te)
            
    def _parallel_te_processing(self, te_copies):
        """Process TEs in parallel using multiple CPU cores"""
        # Split TEs into batches for parallel processing
        batch_size = max(1, len(te_copies) // (self.n_processes * 4))
        te_batches = [te_copies[i:i + batch_size] 
                     for i in range(0, len(te_copies), batch_size)]
        
        # Process batches in parallel
        process_batch = partial(_process_te_batch, 
                              genome_size=self.genome.size,
                              tissue_type=self.tissue_type,
                              species=self.species)
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            results = list(executor.map(process_batch, te_batches))
            
        # Flatten results
        return [item for batch_result in results for item in batch_result]

    def _parallel_te_processing_gpu(self, te_copies):
        """Process TEs in parallel using GPU"""
        if not te_copies:
            return []
            
        # Prepare data for GPU
        locations = torch.tensor([te.location for te in te_copies], device=self.device)
        types = torch.tensor([te.te_type.value for te in te_copies], dtype=torch.long, device=self.device)
        
        # Generate random numbers in batch
        rand_vals = torch.rand(len(te_copies), 4, device=self.device)
        
        # Vectorized operations on GPU
        is_active = torch.tensor(
            [not te.is_silenced and not te.is_dead for te in te_copies],
            dtype=torch.bool,
            device=self.device
        )
        
        # Calculate transposition probabilities
        trans_rates = torch.tensor(
            [te.properties.transposition_rate for te in te_copies],
            dtype=torch.float32,
            device=self.device
        )
        
        # Determine which TEs will transpose
        will_transpose = (rand_vals[:, 0] < trans_rates) & is_active
        
        # Calculate new locations for transposing TEs
        new_locations = torch.empty_like(locations)
        bias_mask = rand_vals[:, 1] < torch.tensor(
            [te.properties.target_site_bias for te in te_copies],
            dtype=torch.float32,
            device=self.device
        )
        
        # Biased locations
        new_locations[bias_mask] = torch.normal(
            mean=torch.tensor(self.genome.size/2, dtype=torch.float32),
            std=torch.tensor(self.genome.size/4, dtype=torch.float32),
            size=(bias_mask.sum(),),
            device=self.device
        ).long()
        
        # Random locations
        new_locations[~bias_mask] = torch.randint(
            0, self.genome.size,
            size=((~bias_mask).sum(),),
            device=self.device
        )
        
        # Ensure locations are within bounds
        new_locations = torch.clamp(new_locations, 0, self.genome.size - 1)
        
        # Process results
        results = []
        for i, te in enumerate(te_copies):
            if will_transpose[i].item():
                activity = te.get_activity_level(self.tissue_type, self.species)
                results.append((te, new_locations[i].item(), activity))
            else:
                results.append((te, None, 0.0))
                
        return results
        
    def step(self):
        """Perform one simulation step with GPU acceleration if available"""
        # Update genome properties
        self.genome.update()
        
        # Process TEs
        if len(self.genome.te_copies) > 100:
            if self.use_gpu:
                te_results = self._parallel_te_processing_gpu(self.genome.te_copies)
            else:
                te_results = self._parallel_te_processing(self.genome.te_copies)
        else:
            te_results = [(te, te.attempt_transposition(self.genome.size),
                          te.get_activity_level(self.tissue_type, self.species))
                         for te in self.genome.te_copies]
        
        # Process results and create new TEs
        new_tes = []
        for te, new_location, activity in te_results:
            if new_location is not None:
                if te.properties.copy_mechanism == "copy_and_paste":
                    num_copies = int(te.properties.progeny_rate * activity)
                    if num_copies > 0:
                        # Vectorized TE creation
                        new_locations = torch.randint(
                            0, self.genome.size,
                            size=(num_copies,),
                            device=self.device
                        ).tolist()
                        
                        new_tes.extend([
                            TransposableElement(
                                te_type=te.te_type,
                                location=loc,
                                properties=te.properties
                            ) for loc in new_locations
                        ])
                else:  # cut_and_paste
                    te.location = new_location
        
        # Bulk add new TEs
        if new_tes:
            self.genome.te_copies.extend(new_tes)
        
        # Update ages (vectorized)
        if self.use_gpu and self.genome.te_copies:
            ages = torch.tensor([te.age for te in self.genome.te_copies], device=self.device)
            ages += 1
            for te, age in zip(self.genome.te_copies, ages.tolist()):
                te.age = age
        else:
            for te in self.genome.te_copies:
                te.age += 1
        
        # Record statistics
        self.history.append(self.get_statistics())
        self.time_step += 1

    def get_statistics(self) -> Dict:
        """Get current simulation statistics using vectorized operations"""
        te_copies = np.array(self.genome.te_copies)
        if len(te_copies) == 0:
            return self._empty_statistics()

        # Vectorized calculations
        active_mask = np.array([te.is_active() for te in te_copies])
        dead_mask = np.array([te.is_dead for te in te_copies])
        ages = np.array([te.age for te in te_copies])
        
        # Count inactivation mechanisms
        inactivation_counts = {
            mechanism: sum(1 for te in te_copies if te.inactivation_mechanism == mechanism)
            for mechanism in InactivationMechanism
        }
        
        stats = {
            "time_step": self.time_step,
            "total_tes": len(te_copies),
            "active_tes": np.sum(active_mask),
            "dead_tes": np.sum(dead_mask),
            "fitness_impact": self.genome.calculate_fitness_impact(),
            "genome_stability": self.genome.genome_stability,
            "host_fitness": self.genome.host_fitness,
            "te_density": self.genome.te_density,
            "te_types": {
                te_type: np.sum([te.te_type == te_type for te in te_copies])
                for te_type in TEType
            },
            "inactivation_mechanisms": {
                mechanism.value: count
                for mechanism, count in inactivation_counts.items()
            },
            "epigenetic_marks": len(self.genome.epigenetic_marks),
            "average_te_age": np.mean(ages) if len(ages) > 0 else 0
        }
        
        return stats

    def _empty_statistics(self) -> Dict:
        """Return empty statistics when no TEs exist"""
        return {
            "time_step": self.time_step,
            "total_tes": 0,
            "active_tes": 0,
            "dead_tes": 0,
            "fitness_impact": 0,
            "genome_stability": self.genome.genome_stability,
            "host_fitness": self.genome.host_fitness,
            "te_density": 0,
            "te_types": {te_type: 0 for te_type in TEType},
            "inactivation_mechanisms": {mechanism.value: 0 for mechanism in InactivationMechanism},
            "epigenetic_marks": 0,
            "average_te_age": 0
        }

    def _load_genome(self):
        """Load the genome sequences"""
        print("Loading genome sequences...")
        try:
            # First try opening as gzip
            with gzip.open(self.genome_file, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    self.genome_sequences[record.id] = str(record.seq)
        except gzip.BadGzipFile:
            # If not gzipped, open as regular file
            with open(self.genome_file, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    self.genome_sequences[record.id] = str(record.seq) 