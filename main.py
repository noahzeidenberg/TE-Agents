from te_simulation import TESimulation, TEType, TEProperties
from visualization import SimulationVisualizer
from tqdm import tqdm
import time
import multiprocessing as mp

def main():
    print(f"Available CPU cores: {mp.cpu_count()}")
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
            progeny_rate=2.0
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
            progeny_rate=3.0
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
            progeny_rate=1.5
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
            progeny_rate=1.0  # Usually 1 for cut-and-paste
        )
    }
    
    # Initialize simulation with parallel processing
    genome_size = 1000000  # 1 million base pairs
    initial_te_count = 100
    n_processes = mp.cpu_count()  # Use all available CPU cores
    
    print(f"\nInitializing simulation with {n_processes} processes...")
    simulation = TESimulation(
        genome_size=genome_size,
        initial_te_count=initial_te_count,
        te_properties=te_properties,
        n_processes=n_processes
    )
    
    # Run simulation
    n_steps = 100
    history = []
    
    print("\nRunning simulation...")
    with tqdm(total=n_steps, desc="Simulation progress") as pbar:
        for step in range(n_steps):
            step_start = time.time()
            simulation.step()
            step_time = time.time() - step_start
            history.append(simulation.get_statistics())
            pbar.set_postfix({"Step time": f"{step_time:.2f}s"})
            pbar.update(1)
    
    total_time = time.time() - start_time
    print(f"\nTotal simulation time: {total_time:.2f} seconds")
    print(f"Average time per step: {total_time/n_steps:.2f} seconds")
    
    # Visualize results
    print("\nGenerating visualizations...")
    visualizer = SimulationVisualizer()
    
    # Plot TE distribution
    visualizer.plot_te_distribution(genome_size, simulation.genome.te_copies)
    
    # Plot TE type distribution
    visualizer.plot_te_type_distribution(simulation.genome.te_copies)
    
    # Plot simulation history
    visualizer.plot_simulation_history(history)
    
    # Plot TE density heatmap
    visualizer.plot_te_heatmap(genome_size, simulation.genome.te_copies)
    
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