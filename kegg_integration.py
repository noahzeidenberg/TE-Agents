import numpy as np
from enum import Enum
from typing import Dict, List, Set, Tuple
from dataclasses import dataclass

class PathwayType(Enum):
    METABOLIC = "metabolic"
    SIGNALING = "signaling"
    REGULATORY = "regulatory"

@dataclass
class KEGGPathway:
    """Represents a KEGG pathway"""
    id: str
    name: str
    type: PathwayType
    genes: Set[str]  # Gene IDs in pathway
    interactions: List[Tuple[str, str]]  # Gene-gene interactions
    robustness: float  # Pathway's robustness to perturbation (0-1)
    
class PathwayNetwork:
    """Manages KEGG pathway interactions"""
    def __init__(self):
        self.pathways: Dict[str, KEGGPathway] = {}
        self.gene_to_pathways: Dict[str, Set[str]] = {}  # Maps genes to pathway IDs
        self.disruption_scores: Dict[str, float] = {}  # Pathway disruption scores
        
    def add_pathway(self, pathway: KEGGPathway):
        """Add a pathway to the network"""
        self.pathways[pathway.id] = pathway
        for gene in pathway.genes:
            if gene not in self.gene_to_pathways:
                self.gene_to_pathways[gene] = set()
            self.gene_to_pathways[gene].add(pathway.id)
        self.disruption_scores[pathway.id] = 0.0
        
    def calculate_pathway_disruption(self, affected_genes: Set[str]) -> Dict[str, float]:
        """Calculate disruption scores for all pathways based on affected genes"""
        disruption_scores = {}
        for pathway_id, pathway in self.pathways.items():
            affected_in_pathway = affected_genes.intersection(pathway.genes)
            if affected_in_pathway:
                # Calculate disruption based on:
                # 1. Proportion of pathway genes affected
                # 2. Pathway robustness
                # 3. Network centrality of affected genes
                proportion_affected = len(affected_in_pathway) / len(pathway.genes)
                centrality_factor = self._calculate_gene_centrality(affected_in_pathway, pathway)
                
                disruption = proportion_affected * (1 - pathway.robustness) * centrality_factor
                disruption_scores[pathway_id] = min(1.0, disruption)
            else:
                disruption_scores[pathway_id] = 0.0
                
        self.disruption_scores = disruption_scores
        return disruption_scores
        
    def _calculate_gene_centrality(self, genes: Set[str], pathway: KEGGPathway) -> float:
        """Calculate the average centrality of affected genes in the pathway"""
        if not genes or not pathway.interactions:
            return 1.0
            
        # Create adjacency matrix
        gene_list = list(pathway.genes)
        n = len(gene_list)
        adj_matrix = np.zeros((n, n))
        
        # Fill adjacency matrix
        gene_to_idx = {gene: idx for idx, gene in enumerate(gene_list)}
        for g1, g2 in pathway.interactions:
            if g1 in gene_to_idx and g2 in gene_to_idx:
                i, j = gene_to_idx[g1], gene_to_idx[g2]
                adj_matrix[i,j] = adj_matrix[j,i] = 1
                
        # Calculate degree centrality for affected genes
        affected_indices = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        if not affected_indices:
            return 1.0
            
        degrees = np.sum(adj_matrix, axis=0)
        max_degree = max(degrees) if max(degrees) > 0 else 1
        centralities = [degrees[i]/max_degree for i in affected_indices]
        
        return np.mean(centralities)
        
    def get_cascade_effects(self, initial_pathways: Set[str], threshold: float = 0.3) -> Set[str]:
        """Identify pathways affected through cascade effects"""
        affected = initial_pathways.copy()
        to_check = initial_pathways.copy()
        
        while to_check:
            pathway_id = to_check.pop()
            pathway = self.pathways[pathway_id]
            
            # Find connected pathways through shared genes
            for gene in pathway.genes:
                connected_pathways = self.gene_to_pathways.get(gene, set())
                for connected_id in connected_pathways:
                    if connected_id not in affected:
                        # Calculate influence score
                        influence = (self.disruption_scores[pathway_id] * 
                                  (1 - self.pathways[connected_id].robustness))
                        if influence > threshold:
                            affected.add(connected_id)
                            to_check.add(connected_id)
                            
        return affected
        
    def calculate_fitness_impact(self) -> float:
        """Calculate overall fitness impact from pathway disruptions"""
        if not self.disruption_scores:
            return 0.0
            
        # Weight disruption scores by pathway importance
        weighted_scores = []
        for pathway_id, disruption in self.disruption_scores.items():
            pathway = self.pathways[pathway_id]
            
            # Weight based on pathway type
            type_weight = {
                PathwayType.METABOLIC: 1.0,
                PathwayType.SIGNALING: 0.8,
                PathwayType.REGULATORY: 0.6
            }[pathway.type]
            
            # Consider cascade effects
            cascade_effect = len(self.get_cascade_effects({pathway_id}))
            cascade_weight = 1.0 + (0.1 * cascade_effect)
            
            weighted_scores.append(disruption * type_weight * cascade_weight)
            
        # Combine scores non-linearly
        return 1.0 - (1.0 / (1.0 + np.sum(weighted_scores))) 