from te_simulation import TESimulation, TEType, TEProperties
from visualization import SimulationVisualizer
from tqdm import tqdm
import time
import multiprocessing as mp
import os
import torch

def main():
    # Get SLURM environment variables
    n_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', mp.cpu_count()))
    has_gpu = torch.cuda.is_available()
    
    print(f"Available CPU cores: {n_cpus}")
    print(f"GPU available: {has_gpu}")
    if has_gpu:
        print(f"GPU device: {torch.cuda.get_device_name(0)}")
    
    start_time = time.time()

    # Define TE properties for different types
    te_properties = {
        TEType.LTR_RETROTRANSPOSON: TEProperties(
            te_type=TEType.LTR_RETROTRANSPOSON,
            transposition_rate=0.1,
            copy_mechanism="copy_and_paste",
            target_site_bias=0.3,
            silencing_sensitivity=0.7,
            fitness_impact=-0.2,
            length=5000,  # Typical LTR retrotransposon length
            death_rate=0.01,
            excision_rate=0.001,
            progeny_rate=2.0,
            # Inactivation parameters
            initial_active_rate=0.15,  # Higher activity rate for LTRs
            truncation_rate=0.35,
            mutation_rate=0.25,
            recombination_rate=0.3,
            silencing_rate=0.1
        ),
        TEType.SINE: TEProperties(
            te_type=TEType.SINE,
            transposition_rate=0.05,
            copy_mechanism="copy_and_paste",
            target_site_bias=0.5,
            silencing_sensitivity=0.8,
            fitness_impact=-0.1,
            length=300,  # Typical SINE length
            death_rate=0.02,
            excision_rate=0.001,
            progeny_rate=3.0,
            # Inactivation parameters
            initial_active_rate=0.1,
            truncation_rate=0.45,  # Higher truncation rate due to small size
            mutation_rate=0.25,
            recombination_rate=0.2,
            silencing_rate=0.1
        ),
        TEType.LINE: TEProperties(
            te_type=TEType.LINE,
            transposition_rate=0.08,
            copy_mechanism="copy_and_paste",
            target_site_bias=0.4,
            silencing_sensitivity=0.6,
            fitness_impact=-0.3,
            length=6000,  # Typical LINE length
            death_rate=0.015,
            excision_rate=0.001,
            progeny_rate=1.5,
            # Inactivation parameters
            initial_active_rate=0.05,  # Lower activity rate due to length
            truncation_rate=0.5,  # Higher truncation rate due to length
            mutation_rate=0.3,
            recombination_rate=0.15,
            silencing_rate=0.05
        ),
        TEType.DNA_TRANSPOSON: TEProperties(
            te_type=TEType.DNA_TRANSPOSON,
            transposition_rate=0.15,
            copy_mechanism="cut_and_paste",
            target_site_bias=0.2,
            silencing_sensitivity=0.5,
            fitness_impact=-0.15,
            length=2000,  # Typical DNA transposon length
            death_rate=0.02,
            excision_rate=0.05,  # Higher excision rate for DNA transposons
            progeny_rate=1.0,
            # Inactivation parameters
            initial_active_rate=0.12,
            truncation_rate=0.3,
            mutation_rate=0.35,  # Higher mutation rate
            recombination_rate=0.25,
            silencing_rate=0.1
        )
    }
    
    # Initialize simulation with HPC parameters
    genome_size = 10000000  # Increased to 10 million base pairs for HPC
    initial_te_count = 1000  # Increased initial TE count
    
    print(f"\nInitializing simulation with {n_cpus} processes...")
    simulation = TESimulation(
        genome_size=genome_size,
        initial_te_count=initial_te_count,
        te_properties=te_properties,
        n_processes=n_cpus,
        use_gpu=has_gpu
    )
    
    # Run simulation with checkpointing
    n_steps = 1000  # Increased number of steps for HPC
    history = []
    checkpoint_interval = 100
    
    print("\nRunning simulation...")
    with tqdm(total=n_steps, desc="Simulation progress") as pbar:
        for step in range(n_steps):
            step_start = time.time()
            simulation.step()
            step_time = time.time() - step_start
            
            # Save state and history at checkpoints
            if (step + 1) % checkpoint_interval == 0:
                checkpoint_file = f"checkpoint_step_{step+1}.pkl"
                simulation.save_checkpoint(checkpoint_file)
                print(f"\nCheckpoint saved: {checkpoint_file}")
            
            history.append(simulation.get_statistics())
            pbar.set_postfix({
                "Step time": f"{step_time:.2f}s",
                "Active TEs": history[-1]['active_tes'],
                "Fitness": f"{history[-1]['host_fitness']:.3f}"
            })
            pbar.update(1)
    
    total_time = time.time() - start_time
    print(f"\nTotal simulation time: {total_time:.2f} seconds")
    print(f"Average time per step: {total_time/n_steps:.2f} seconds")
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    visualizer = SimulationVisualizer()
    
    visualizer.plot_te_distribution(genome_size, simulation.genome.te_copies)
    visualizer.plot_te_type_distribution(simulation.genome.te_copies)
    visualizer.plot_simulation_history(history)
    visualizer.plot_te_heatmap(genome_size, simulation.genome.te_copies)
    
    if simulation.genome.pathway_network:
        visualizer.plot_pathway_disruption(simulation.genome.pathway_network)
    
    # Print final statistics
    print("\nFinal Statistics:")
    final_stats = simulation.get_statistics()
    print(f"Total TEs: {final_stats['total_tes']}")
    print(f"Active TEs: {final_stats['active_tes']}")
    print(f"Fitness Impact: {final_stats['fitness_impact']:.3f}")
    print("\nTE Type Distribution:")
    for te_type, count in final_stats['te_types'].items():
        print(f"{te_type.value}: {count}")

if __name__ == "__main__":
    main() 