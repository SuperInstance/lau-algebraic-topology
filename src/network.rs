//! Agent network topology analysis using algebraic topology.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use crate::chain_complex::{ChainComplex, HomologyGroup};
use crate::euler::euler_from_betti;
use nalgebra::DMatrix;

/// An agent in the network.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Agent {
    /// Agent identifier.
    pub id: String,
    /// Metadata.
    pub metadata: HashMap<String, String>,
}

/// A connection between agents.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Connection {
    /// Source agent.
    pub from: String,
    /// Target agent.
    pub to: String,
    /// Weight (optional).
    pub weight: f64,
}

/// An agent network topology.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AgentNetwork {
    /// Agents (vertices).
    pub agents: Vec<Agent>,
    /// Connections (edges).
    pub connections: Vec<Connection>,
    /// Whether connections are bidirectional.
    pub bidirectional: bool,
}

impl AgentNetwork {
    /// Create a new empty network.
    pub fn new() -> Self {
        AgentNetwork {
            agents: Vec::new(),
            connections: Vec::new(),
            bidirectional: true,
        }
    }

    /// Add an agent.
    pub fn add_agent(&mut self, id: &str) {
        self.agents.push(Agent {
            id: id.to_string(),
            metadata: HashMap::new(),
        });
    }

    /// Add a connection.
    pub fn connect(&mut self, from: &str, to: &str) {
        self.connections.push(Connection {
            from: from.to_string(),
            to: to.to_string(),
            weight: 1.0,
        });
    }

    /// Build the adjacency matrix.
    pub fn adjacency_matrix(&self) -> DMatrix<f64> {
        let n = self.agents.len();
        let mut mat = DMatrix::zeros(n, n);
        let idx: HashMap<&str, usize> = self.agents.iter()
            .enumerate()
            .map(|(i, a)| (a.id.as_str(), i))
            .collect();

        for conn in &self.connections {
            if let (Some(&i), Some(&j)) = (idx.get(conn.from.as_str()), idx.get(conn.to.as_str())) {
                mat[(i, j)] += conn.weight;
                if self.bidirectional {
                    mat[(j, i)] += conn.weight;
                }
            }
        }
        mat
    }

    /// Compute connected components.
    pub fn connected_components(&self) -> Vec<Vec<String>> {
        let n = self.agents.len();
        let mut parent: Vec<usize> = (0..n).collect();
        let idx: HashMap<&str, usize> = self.agents.iter()
            .enumerate()
            .map(|(i, a)| (a.id.as_str(), i))
            .collect();

        fn find(parent: &mut Vec<usize>, i: usize) -> usize {
            if parent[i] != i {
                parent[i] = find(parent, parent[i]);
            }
            parent[i]
        }

        for conn in &self.connections {
            if let (Some(&i), Some(&j)) = (idx.get(conn.from.as_str()), idx.get(conn.to.as_str())) {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
                if self.bidirectional {
                    let ri2 = find(&mut parent, i);
                    let rj2 = find(&mut parent, j);
                    if ri2 != rj2 {
                        parent[ri2] = rj2;
                    }
                }
            }
        }

        let mut components: HashMap<usize, Vec<String>> = HashMap::new();
        for (i, agent) in self.agents.iter().enumerate() {
            let root = find(&mut parent, i);
            components.entry(root).or_default().push(agent.id.clone());
        }
        components.into_values().collect()
    }

    /// Number of connected components.
    pub fn num_components(&self) -> usize {
        self.connected_components().len()
    }

    /// Check if network is fully connected.
    pub fn is_connected_when_complete(&self) -> bool {
        self.num_components() == 1
    }

    /// Compute the homology of the network (viewed as a simplicial complex).
    /// We build a Vietoris-Rips complex with parameter r.
    pub fn network_homology(&self, max_dim: usize) -> Vec<HomologyGroup> {
        let n = self.agents.len();
        if n == 0 {
            return vec![];
        }

        // Build adjacency as a simplicial complex
        // 0-simplices: agents
        // 1-simplices: connections
        // Higher simplices: cliques in the connection graph

        let idx: HashMap<&str, usize> = self.agents.iter()
            .enumerate()
            .map(|(i, a)| (a.id.as_str(), i))
            .collect();

        let mut adj: HashSet<(usize, usize)> = HashSet::new();
        for conn in &self.connections {
            if let (Some(&i), Some(&j)) = (idx.get(conn.from.as_str()), idx.get(conn.to.as_str())) {
                adj.insert((i.min(j), i.max(j)));
            }
        }

        // Build simplices by dimension
        let mut simplices: Vec<HashSet<Vec<usize>>> = Vec::new();

        // 0-simplices
        let mut s0 = HashSet::new();
        for i in 0..n {
            s0.insert(vec![i]);
        }
        simplices.push(s0);

        // 1-simplices
        let mut s1 = HashSet::new();
        for &(i, j) in &adj {
            s1.insert(vec![i, j]);
        }
        simplices.push(s1);

        // Higher simplices: find cliques
        for dim in 2..=max_dim {
            let mut sd = HashSet::new();
            // A (dim+1)-clique exists iff all pairs are adjacent
            if let Some(prev) = simplices.get(dim - 1) {
                for simplex in prev {
                    // Try extending by each vertex
                    for v in 0..n {
                        if simplex.contains(&v) {
                            continue;
                        }
                        // Check if v is adjacent to all vertices in simplex
                        let all_connected = simplex.iter().all(|&u| {
                            adj.contains(&(u.min(v), u.max(v)))
                        });
                        if all_connected {
                            let mut new_simplex = simplex.clone();
                            new_simplex.push(v);
                            new_simplex.sort();
                            if new_simplex.len() == dim + 1 {
                                sd.insert(new_simplex);
                            }
                        }
                    }
                }
            }
            simplices.push(sd);
        }

        // Build chain complex
        let ranks: Vec<usize> = simplices.iter().map(|s| s.len()).collect();
        let mut boundary_maps = Vec::new();

        for d in 0..simplices.len().saturating_sub(1) {
            let sources: Vec<_> = simplices.get(d + 1).cloned().unwrap_or_default().into_iter().collect();
            let targets: Vec<_> = simplices.get(d).cloned().unwrap_or_default().into_iter().collect();

            let target_idx: HashMap<Vec<usize>, usize> = targets.iter()
                .enumerate()
                .map(|(i, s)| (s.clone(), i))
                .collect();

            let mut mat = DMatrix::zeros(targets.len(), sources.len());

            for (j, simplex) in sources.iter().enumerate() {
                for (i, _) in simplex.iter().enumerate() {
                    let mut face: Vec<usize> = simplex.clone();
                    face.remove(i);
                    face.sort();
                    if let Some(&row) = target_idx.get(&face) {
                        mat[(row, j)] += if i % 2 == 0 { 1 } else { -1 };
                    }
                }
            }

            boundary_maps.push(mat);
        }

        let cc = ChainComplex::from_ranks_and_maps(ranks, boundary_maps);
        cc.all_homology(max_dim)
    }

    /// Analyze the network topology.
    pub fn analyze(&self) -> NetworkAnalysis {
        let components = self.connected_components();
        let num_components = components.len();
        let homology = self.network_homology(3);
        let betti: Vec<usize> = homology.iter().map(|h| h.rank()).collect();
        let euler = euler_from_betti(&betti);

        NetworkAnalysis {
            num_agents: self.agents.len(),
            num_connections: self.connections.len(),
            num_components,
            betti_numbers: betti.clone(),
            euler_characteristic: euler,
            homology_groups: homology,
            components,
            is_connected: num_components == 1,
            holes_1d: betti.get(1).copied().unwrap_or(0),
            voids_2d: betti.get(2).copied().unwrap_or(0),
            voids_3d: betti.get(3).copied().unwrap_or(0),
        }
    }
}

/// Analysis of a network topology.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NetworkAnalysis {
    /// Number of agents.
    pub num_agents: usize,
    /// Number of connections.
    pub num_connections: usize,
    /// Number of connected components.
    pub num_components: usize,
    /// Betti numbers.
    pub betti_numbers: Vec<usize>,
    /// Euler characteristic.
    pub euler_characteristic: i64,
    /// Homology groups.
    pub homology_groups: Vec<HomologyGroup>,
    /// Connected components.
    pub components: Vec<Vec<String>>,
    /// Whether the network is connected.
    pub is_connected: bool,
    /// Number of 1-dimensional holes (loops).
    pub holes_1d: usize,
    /// Number of 2-dimensional voids.
    pub voids_2d: usize,
    /// Number of 3-dimensional voids.
    pub voids_3d: usize,
}

impl std::fmt::Display for NetworkAnalysis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Network Analysis:")?;
        writeln!(f, "  Agents: {}", self.num_agents)?;
        writeln!(f, "  Connections: {}", self.num_connections)?;
        writeln!(f, "  Connected: {}", if self.is_connected { "yes" } else { "no" })?;
        writeln!(f, "  Components: {}", self.num_components)?;
        writeln!(f, "  Betti numbers: {:?}", self.betti_numbers)?;
        writeln!(f, "  Euler characteristic: {}", self.euler_characteristic)?;
        writeln!(f, "  1D holes: {}", self.holes_1d)?;
        writeln!(f, "  2D voids: {}", self.voids_2d)?;
        writeln!(f, "  3D voids: {}", self.voids_3d)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_network() {
        let net = AgentNetwork::new();
        assert_eq!(net.agents.len(), 0);
    }

    #[test]
    fn test_connected_network() {
        let mut net = AgentNetwork::new();
        net.add_agent("A");
        net.add_agent("B");
        net.add_agent("C");
        net.connect("A", "B");
        net.connect("B", "C");
        assert!(net.is_connected_when_complete());
    }

    #[test]
    fn test_disconnected_network() {
        let mut net = AgentNetwork::new();
        net.add_agent("A");
        net.add_agent("B");
        net.add_agent("C");
        net.connect("A", "B");
        // C is disconnected
        assert_eq!(net.num_components(), 2);
    }

    #[test]
    fn test_triangle_network() {
        let mut net = AgentNetwork::new();
        net.add_agent("A");
        net.add_agent("B");
        net.add_agent("C");
        net.connect("A", "B");
        net.connect("B", "C");
        net.connect("A", "C");
        let analysis = net.analyze();
        assert!(analysis.is_connected);
        // Triangle fills the hole, so β_1 = 0
        assert_eq!(analysis.holes_1d, 0);
    }

    #[test]
    fn test_square_network_hole() {
        let mut net = AgentNetwork::new();
        net.add_agent("A");
        net.add_agent("B");
        net.add_agent("C");
        net.add_agent("D");
        net.connect("A", "B");
        net.connect("B", "C");
        net.connect("C", "D");
        net.connect("D", "A");
        // Square has a 1D hole (no diagonal)
        let analysis = net.analyze();
        assert!(analysis.is_connected);
        assert_eq!(analysis.holes_1d, 1);
    }

    #[test]
    fn test_two_triangles() {
        let mut net = AgentNetwork::new();
        for c in "ABCDEF".chars() {
            net.add_agent(&c.to_string());
        }
        // Two triangles sharing no vertices
        net.connect("A", "B");
        net.connect("B", "C");
        net.connect("A", "C");
        net.connect("D", "E");
        net.connect("E", "F");
        net.connect("D", "F");
        let analysis = net.analyze();
        assert_eq!(analysis.num_components, 2);
    }

    #[test]
    fn test_tetrahedron_network() {
        let mut net = AgentNetwork::new();
        net.add_agent("A");
        net.add_agent("B");
        net.add_agent("C");
        net.add_agent("D");
        // Complete graph K4
        net.connect("A", "B");
        net.connect("A", "C");
        net.connect("A", "D");
        net.connect("B", "C");
        net.connect("B", "D");
        net.connect("C", "D");
        let analysis = net.analyze();
        assert!(analysis.is_connected);
        // K4 is a tetrahedron = S^2 boundary filled, so β_1 = 0
        assert_eq!(analysis.holes_1d, 0);
    }

    #[test]
    fn test_network_analysis_display() {
        let mut net = AgentNetwork::new();
        net.add_agent("A");
        net.add_agent("B");
        net.connect("A", "B");
        let analysis = net.analyze();
        let display = format!("{}", analysis);
        assert!(display.contains("Agents: 2"));
    }

    #[test]
    fn test_adjacency_matrix() {
        let mut net = AgentNetwork::new();
        net.add_agent("A");
        net.add_agent("B");
        net.connect("A", "B");
        let adj = net.adjacency_matrix();
        assert_eq!(adj[(0, 1)], 1.0);
        assert_eq!(adj[(1, 0)], 1.0); // bidirectional
    }
}
