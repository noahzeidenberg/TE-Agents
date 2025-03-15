import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import List, Dict
from te_simulation import TEType, InactivationMechanism
import os
import networkx as nx

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
        """Plot the distribution of TE types and their inactivation states"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Count TEs by type
        type_counts = {te_type: 0 for te_type in TEType}
        for te in te_copies:
            type_counts[te.te_type] += 1
            
        # Plot TE types
        ax1.bar([t.name_str for t in type_counts.keys()], 
                list(type_counts.values()))
        ax1.set_xlabel('TE Type')
        ax1.set_ylabel('Number of Copies')
        ax1.set_title('Distribution of TE Types')
        ax1.tick_params(axis='x', rotation=45)
        
        # Count inactivation mechanisms
        inactivation_counts = {
            mechanism: sum(1 for te in te_copies if te.inactivation_mechanism == mechanism)
            for mechanism in InactivationMechanism
        }
        
        # Plot inactivation mechanisms
        colors = {'active': 'green', 'truncated': 'red', 'mutated': 'orange',
                 'recombined': 'purple', 'silenced': 'blue'}
        ax2.bar([m.value for m in inactivation_counts.keys()],
                list(inactivation_counts.values()),
                color=[colors[m.value] for m in inactivation_counts.keys()])
        ax2.set_xlabel('Inactivation Mechanism')
        ax2.set_ylabel('Number of TEs')
        ax2.set_title('Distribution of Inactivation Mechanisms')
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'te_types.png'))
        plt.close()
        
    def plot_simulation_history(self, history: List[Dict]):
        """Plot simulation statistics over time"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Extract data
        time_steps = [h['time_step'] for h in history]
        total_tes = [h['total_tes'] for h in history]
        active_tes = [h['active_tes'] for h in history]
        fitness_impact = [h['fitness_impact'] for h in history]
        
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
        
        # Plot inactivation mechanisms over time
        mechanisms = list(history[0]['inactivation_mechanisms'].keys())
        for mechanism in mechanisms:
            counts = [h['inactivation_mechanisms'][mechanism] for h in history]
            ax3.plot(time_steps, counts, label=mechanism)
        ax3.set_xlabel('Time Step')
        ax3.set_ylabel('Number of TEs')
        ax3.set_title('Inactivation Mechanisms Over Time')
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Plot proportions of active vs inactive
        inactive_tes = [h['total_tes'] - h['active_tes'] for h in history]
        ax4.stackplot(time_steps, [active_tes, inactive_tes],
                     labels=['Active', 'Inactive'],
                     colors=['green', 'red'])
        ax4.set_xlabel('Time Step')
        ax4.set_ylabel('Number of TEs')
        ax4.set_title('Active vs Inactive TEs')
        ax4.legend()
        
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
        
    def plot_pathway_disruption(self, pathway_network: 'PathwayNetwork'):
        """Plot pathway disruption patterns"""
        plt.figure(figsize=(15, 10))
        
        # Create subplots
        gs = plt.GridSpec(2, 2)
        ax1 = plt.subplot(gs[0, 0])  # Pathway type disruption
        ax2 = plt.subplot(gs[0, 1])  # Most affected pathways
        ax3 = plt.subplot(gs[1, :])  # Pathway network graph
        
        # Plot disruption by pathway type
        type_disruption = {}
        for pathway_id, score in pathway_network.disruption_scores.items():
            pathway_type = pathway_network.pathways[pathway_id].type
            if pathway_type not in type_disruption:
                type_disruption[pathway_type] = []
            type_disruption[pathway_type].append(score)
            
        types = list(type_disruption.keys())
        means = [np.mean(scores) for scores in type_disruption.values()]
        sems = [np.std(scores)/np.sqrt(len(scores)) for scores in type_disruption.values()]
        
        ax1.bar([t.value for t in types], means, yerr=sems)
        ax1.set_xlabel('Pathway Type')
        ax1.set_ylabel('Average Disruption Score')
        ax1.set_title('Disruption by Pathway Type')
        ax1.tick_params(axis='x', rotation=45)
        
        # Plot most affected pathways
        top_pathways = sorted(
            pathway_network.disruption_scores.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]
        
        names = [pathway_network.pathways[pid].name for pid, _ in top_pathways]
        scores = [score for _, score in top_pathways]
        
        ax2.barh(range(len(names)), scores)
        ax2.set_yticks(range(len(names)))
        ax2.set_yticklabels(names)
        ax2.set_xlabel('Disruption Score')
        ax2.set_title('Most Affected Pathways')
        
        # Plot pathway network graph
        # Add nodes
        for pid, pathway in pathway_network.pathways.items():
            G.add_node(pid, 
                      type=pathway.type.value,
                      disruption=pathway_network.disruption_scores[pid])
            
        # Add edges based on shared genes
        for pid1 in pathway_network.pathways:
            genes1 = pathway_network.pathways[pid1].genes
            for pid2 in pathway_network.pathways:
                if pid1 < pid2:  # Avoid duplicate edges
                    genes2 = pathway_network.pathways[pid2].genes
                    shared = len(genes1.intersection(genes2))
                    if shared > 0:
                        G.add_edge(pid1, pid2, weight=shared)
        
        # Draw network
        pos = nx.spring_layout(G)
        
        # Node colors based on disruption
        node_colors = [G.nodes[node]['disruption'] for node in G.nodes()]
        
        # Node sizes based on number of genes
        node_sizes = [len(pathway_network.pathways[pid].genes) for pid in G.nodes()]
        
        nx.draw_networkx(
            G, pos,
            ax=ax3,
            node_color=node_colors,
            node_size=[s*100 for s in node_sizes],
            cmap=plt.cm.YlOrRd,
            with_labels=False,
            alpha=0.7
        )
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd)
        sm.set_array([])
        plt.colorbar(sm, ax=ax3, label='Disruption Score')
        
        ax3.set_title('Pathway Interaction Network')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'pathway_analysis.png'))
        plt.close() 