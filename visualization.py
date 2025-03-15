import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import List, Dict
from te_simulation import TEType
import os

class SimulationVisualizer:
    """Handles visualization of TE simulation results"""
    
    def __init__(self):
        # Create output directory if it doesn't exist
        self.output_dir = "simulation_plots"
        os.makedirs(self.output_dir, exist_ok=True)
    
    def plot_te_distribution(self, genome_size: int, te_copies: List['TransposableElement']):
        """Plot the distribution of TEs along the genome"""
        plt.figure(figsize=(12, 4))
        
        # Create histogram of TE locations
        locations = [te.location for te in te_copies]
        plt.hist(locations, bins=50, alpha=0.6, label='All TEs')
        
        # Plot active vs silenced TEs
        active_locations = [te.location for te in te_copies if not te.is_silenced]
        silenced_locations = [te.location for te in te_copies if te.is_silenced]
        
        plt.hist(active_locations, bins=50, alpha=0.4, label='Active TEs', color='green')
        plt.hist(silenced_locations, bins=50, alpha=0.4, label='Silenced TEs', color='red')
        
        plt.xlabel('Genome Position')
        plt.ylabel('Number of TEs')
        plt.title('Distribution of TEs Along the Genome')
        plt.legend()
        plt.savefig(os.path.join(self.output_dir, 'te_distribution.png'))
        plt.close()
        
    def plot_te_type_distribution(self, te_copies: List['TransposableElement']):
        """Plot the distribution of TE types"""
        plt.figure(figsize=(8, 6))
        
        # Count TEs by type
        type_counts = {te_type: 0 for te_type in TEType}
        for te in te_copies:
            type_counts[te.te_type] += 1
            
        # Create bar plot
        plt.bar([t.value for t in type_counts.keys()], 
                list(type_counts.values()))
        
        plt.xlabel('TE Type')
        plt.ylabel('Number of Copies')
        plt.title('Distribution of TE Types')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'te_types.png'))
        plt.close()
        
    def plot_simulation_history(self, history: List[Dict]):
        """Plot simulation statistics over time"""
        # Extract data
        time_steps = [h['time_step'] for h in history]
        total_tes = [h['total_tes'] for h in history]
        active_tes = [h['active_tes'] for h in history]
        fitness_impact = [h['fitness_impact'] for h in history]
        
        # Create subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Plot TE counts
        ax1.plot(time_steps, total_tes, label='Total TEs', color='blue')
        ax1.plot(time_steps, active_tes, label='Active TEs', color='green')
        ax1.set_xlabel('Time Step')
        ax1.set_ylabel('Number of TEs')
        ax1.set_title('TE Population Dynamics')
        ax1.legend()
        
        # Plot fitness impact
        ax2.plot(time_steps, fitness_impact, color='red')
        ax2.set_xlabel('Time Step')
        ax2.set_ylabel('Fitness Impact')
        ax2.set_title('Host Fitness Impact Over Time')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'simulation_history.png'))
        plt.close()
        
    def plot_te_heatmap(self, genome_size: int, te_copies: List['TransposableElement'], 
                       window_size: int = 1000):
        """Create a heatmap of TE density along the genome"""
        # Create density array
        n_windows = genome_size // window_size
        density = np.zeros(n_windows)
        
        # Calculate TE density in windows
        for te in te_copies:
            window_idx = te.location // window_size
            if window_idx < len(density):
                density[window_idx] += 1
                
        # Calculate grid dimensions for square-ish layout
        n_rows = int(np.floor(np.sqrt(n_windows)))
        n_cols = n_rows  # Make it square
        
        # Trim density array to fit square grid
        density = density[:n_rows * n_cols]
        heatmap = density.reshape(n_rows, n_cols)
        
        # Plot heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap, cmap='YlOrRd')
        plt.xlabel('Genome Position (Windows)')
        plt.ylabel('Genome Position (Windows)')
        plt.title('TE Density Heatmap')
        plt.savefig(os.path.join(self.output_dir, 'te_heatmap.png'))
        plt.close() 