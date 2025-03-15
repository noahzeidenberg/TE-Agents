# Methods: Agent-Based Model of Transposable Element Dynamics

## Model Overview

We developed an agent-based model to simulate the dynamics of transposable elements (TEs) within a host genome. The model incorporates both ecological and molecular perspectives, allowing for the study of TE-host interactions at multiple scales. The implementation is written in Python and uses NumPy for numerical computations.

## Core Components

### 1. Transposable Element Properties

Each TE is characterized by a set of biological properties:

- **Basic Properties**
  - Type (LTR retrotransposon, SINE, LINE, DNA transposon)
  - Length
  - Transposition rate
  - Copy mechanism (cut-and-paste or copy-and-paste)
  - Target site bias
  - Silencing sensitivity
  - Impact on host

- **Lifecycle Properties**
  - Death rate
  - Excision rate
  - Progeny rate
  - Age-dependent activity

- **Environmental Interactions**
  - Epigenetic sensitivity
  - Tissue-specific activity
  - Species-specific activity
  - Horizontal transfer rate
  - Regulatory impact
  - Structural impact

### 2. Genome Environment

The host genome is modeled as a continuous environment with the following features:

- **Structure**
  - Total size (base pairs)
  - Gene locations
  - Regulatory regions
  - Chromosome structure
  - Epigenetic marks

- **Properties**
  - Recombination rate
  - Mutation rate
  - Gene density
  - TE density
  - Genome stability
  - Host fitness

### 3. Simulation Dynamics

The simulation proceeds in discrete time steps, with the following processes occurring each step:

1. **Genome Updates**
   - Update TE density
   - Update epigenetic marks
   - Update silencing status
   - Calculate genome stability
   - Update host fitness

2. **TE Updates**
   - Update epigenetic state
   - Calculate activity level
   - Attempt transposition
   - Create new copies
   - Update age

3. **Statistics Recording**
   - Total TE count
   - Active vs. silenced TEs
   - TE type distribution
   - Fitness impact
   - Genome stability
   - Epigenetic mark distribution

## Implementation Details

### Data Structures

1. **TEProperties (dataclass)**
   ```python
   @dataclass
   class TEProperties:
       te_type: TEType
       transposition_rate: float
       copy_mechanism: str
       target_site_bias: float
       silencing_sensitivity: float
       fitness_impact: float
       length: int
       death_rate: float
       excision_rate: float
       progeny_rate: float
       # ... additional properties
   ```

2. **TransposableElement (class)**
   - Location tracking
   - State management (silenced/dead)
   - Age tracking
   - Transposition history
   - Epigenetic state

3. **GenomeEnvironment (class)**
   - TE collection
   - Genome features
   - Environmental properties
   - Update mechanisms

### Key Algorithms

1. **Transposition Process**
   ```python
   def attempt_transposition(self, genome_size: int) -> Optional[int]:
       # Check viability
       if self.is_silenced or self.is_dead:
           return None
           
       # Apply death rate
       if random.random() < self.properties.death_rate:
           self.is_dead = True
           return None
           
       # Calculate activity level
       activity = self.get_activity_level()
       
       # Determine new location
       new_location = self._determine_insertion_site(genome_size)
       
       return new_location
   ```

2. **Genome Stability Calculation**
   ```python 
   def calculate_genome_stability(self) -> float:  # currently penalizing if TEs are densely packed, but ecosystem engineering implies that could be a helpful mechanism
       stability = 1.0
       
       # Density penalty
       if self.te_density > 0.5:
           stability *= (1 - (self.te_density - 0.5))
           
       # Clustering penalty
       for i in range(len(te_locations) - 1):
           if te_locations[i + 1] - te_locations[i] < 1000:
               stability *= 0.95
               
       return max(0.0, min(1.0, stability))
   ```

## Visualization Methods

The model includes comprehensive visualization capabilities:

1. **TE Distribution Plot**
   - Histogram of TE locations
   - Active vs. silenced TE distribution
   - Genome position mapping

2. **TE Type Distribution**
   - Bar plot of TE types
   - Population composition
   - Type-specific statistics

3. **Simulation History**
   - TE population dynamics
   - Fitness impact over time
   - Active vs. total TE ratio

4. **TE Density Heatmap**
   - Window-based density calculation
   - 2D heatmap representation
   - Spatial distribution analysis

## Parameter Settings

Default parameter values used in simulations:

```python
# TE Properties
te_properties = {
    TEType.LTR_RETROTRANSPOSON: TEProperties(
        transposition_rate=0.1,
        target_site_bias=0.3,
        silencing_sensitivity=0.7,
        fitness_impact=-0.2
    ),
    # ... other TE types
}

# Genome Parameters
genome_size = 1000000  # 1 million base pairs
initial_te_count = 100
simulation_steps = 100
```

## Usage

The simulation can be run using the following command:

```python
simulation = TESimulation(genome_size, initial_te_count, te_properties)
for _ in range(n_steps):
    simulation.step()
```

## Data Collection and Analysis

Statistics are collected at each time step and include:

1. **Population Metrics**
   - Total TE count
   - Active TE count
   - Dead TE count
   - TE type distribution

2. **Environmental Metrics**
   - Genome stability
   - Host fitness
   - TE density
   - Epigenetic mark distribution

3. **Temporal Metrics**
   - Average TE age
   - Transposition frequency
   - Fitness impact over time

## Limitations and Assumptions

1. **Model Limitations**
   - Discrete time steps
   - Simplified genome structure
   - Homogeneous environment
   - Independent TE interactions

2. **Biological Assumptions**
   - Constant mutation rates
   - Simplified epigenetic effects
   - Linear fitness impacts
   - Independent TE families

## Future Improvements

1. **Model Enhancements**
   - Continuous time simulation
   - Complex genome structure
   - Heterogeneous environment
   - TE-TE interactions

2. **Biological Realism**
   - Variable mutation rates
   - Complex epigenetic networks
   - Non-linear fitness effects
   - TE family interactions 