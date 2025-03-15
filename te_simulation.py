import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from enum import Enum
import random
import multiprocessing as mp
from functools import partial
import concurrent.futures

class TEType(Enum):
    LTR_RETROTRANSPOSON = "LTR_retrotransposon"
    SINE = "SINE"
    LINE = "LINE"
    DNA_TRANSPOSON = "DNA_transposon"

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
    preferred_insertion_sites: List[Tuple[int, int]] = None  # Preferred genomic regions for insertion
    epigenetic_sensitivity: float = 0.5  # How sensitive to epigenetic modifications
    host_tissue_specificity: Dict[str, float] = None  # Tissue-specific activity levels
    horizontal_transfer_rate: float = 0.0  # Rate of horizontal transfer between hosts
    regulatory_impact: float = 0.0  # Impact on nearby gene regulation
    structural_impact: float = 0.0  # Impact on genome structure (e.g., chromosomal rearrangements)
    age_dependent_activity: bool = False  # Whether activity changes with age
    host_species_specificity: Dict[str, float] = None  # Species-specific activity levels

class TransposableElement:
    """Represents a single transposable element"""
    def __init__(self, 
                 te_type: TEType,
                 location: int,
                 properties: TEProperties,
                 is_silenced: bool = False,
                 is_dead: bool = False):
        self.te_type = te_type
        self.location = location
        self.properties = properties
        self.is_silenced = is_silenced
        self.is_dead = is_dead
        self.age = 0  # Time steps since creation
        self.copies_made = 0  # Number of copies this TE has produced
        self.last_transposition = 0  # Time step of last transposition
        self.epigenetic_state = {}  # Current epigenetic state
        
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
                 chromosome_structure: Dict[str, List[Tuple[int, int]]] = None):
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
        """Update host fitness based on TE activity and genome stability"""
        # Base fitness from TE impact
        self.host_fitness = 1.0 + self.calculate_fitness_impact()
        
        # Adjust for genome stability
        self.genome_stability = self.calculate_genome_stability()
        self.host_fitness *= self.genome_stability
        
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
    """Main simulation class with parallel processing"""
    def __init__(self, 
                 genome_size: int,
                 initial_te_count: int,
                 te_properties: Dict[TEType, TEProperties],
                 tissue_type: str = None,
                 species: str = None,
                 n_processes: int = None):
        self.genome = GenomeEnvironment(genome_size)
        self.te_properties = te_properties
        self.time_step = 0
        self.tissue_type = tissue_type
        self.species = species
        self.history = []
        self.n_processes = n_processes or mp.cpu_count()
        self.initialize_tes(initial_te_count)
        
    def initialize_tes(self, count: int):
        """Initialize TEs using vectorized operations"""
        te_types = np.random.choice(list(self.te_properties.keys()), size=count)
        locations = np.random.randint(0, self.genome.size, size=count)
        
        # Vectorized creation of TEs
        self.genome.te_copies = [
            TransposableElement(
                te_type=te_type,
                location=loc,
                properties=self.te_properties[te_type]
            ) for te_type, loc in zip(te_types, locations)
        ]

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

    def step(self):
        """Perform one simulation step with parallel processing"""
        # Update genome properties (vectorized where possible)
        self.genome.update()
        
        # Update epigenetic states in parallel
        active_tes = [te for te in self.genome.te_copies if not te.is_silenced]
        
        # Process TEs in parallel
        if len(self.genome.te_copies) > 100:  # Only parallelize for large numbers of TEs
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
                    # Vectorized creation of new copies
                    num_copies = int(te.properties.progeny_rate * activity)
                    if num_copies > 0:
                        new_tes.extend([
                            TransposableElement(
                                te_type=te.te_type,
                                location=new_location,
                                properties=te.properties
                            ) for _ in range(num_copies)
                        ])
                else:  # cut_and_paste
                    te.location = new_location

        # Bulk add new TEs
        if new_tes:
            self.genome.te_copies.extend(new_tes)

        # Vectorized age update
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
        active_mask = np.array([not te.is_silenced for te in te_copies])
        dead_mask = np.array([te.is_dead for te in te_copies])
        ages = np.array([te.age for te in te_copies])
        
        return {
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
            "epigenetic_marks": len(self.genome.epigenetic_marks),
            "average_te_age": np.mean(ages) if len(ages) > 0 else 0
        }

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
            "epigenetic_marks": 0,
            "average_te_age": 0
        } 