//! Simplicial complexes and simplicial homology.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap, HashSet};
use crate::chain_complex::{ChainComplex, HomologyGroup};

/// A vertex index.
pub type Vertex = usize;

/// A simplex represented as a sorted set of vertices.
pub type Simplex = BTreeSet<Vertex>;

/// A simplicial complex.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SimplicialComplex {
    /// All simplices, indexed by dimension.
    simplices: Vec<HashSet<BTreeSet<Vertex>>>,
    /// All vertices.
    vertices: HashSet<Vertex>,
    /// Name for display.
    pub name: String,
}

impl SimplicialComplex {
    /// Create a new empty simplicial complex.
    pub fn new(name: &str) -> Self {
        SimplicialComplex {
            simplices: Vec::new(),
            vertices: HashSet::new(),
            name: name.to_string(),
        }
    }

    /// Add a simplex and all its faces.
    pub fn add_simplex(&mut self, vertices: &[Vertex]) {
        let simplex: BTreeSet<Vertex> = vertices.iter().copied().collect();
        self.vertices.extend(vertices.iter().copied());

        let dim = simplex.len().saturating_sub(1);

        // Ensure we have enough dimensions
        while self.simplices.len() <= dim {
            self.simplices.push(HashSet::new());
        }

        // Add the simplex and all its faces
        self.add_simplex_recursive(&simplex);
    }

    fn add_simplex_recursive(&mut self, simplex: &BTreeSet<Vertex>) {
        let dim = simplex.len().saturating_sub(1);
        while self.simplices.len() <= dim {
            self.simplices.push(HashSet::new());
        }
        if !self.simplices[dim].contains(simplex) {
            self.simplices[dim].insert(simplex.clone());
            if simplex.len() > 1 {
                // Add all (dim-1)-faces
                for v in simplex.iter() {
                    let mut face = simplex.clone();
                    face.remove(v);
                    self.add_simplex_recursive(&face);
                }
            }
        }
    }

    /// Get number of simplices of given dimension.
    pub fn num_simplices(&self, dim: usize) -> usize {
        self.simplices.get(dim).map(|s| s.len()).unwrap_or(0)
    }

    /// Get the dimension (highest dimension of any simplex).
    pub fn dimension(&self) -> usize {
        self.simplices.len().saturating_sub(1)
    }

    /// Get all simplices of given dimension, sorted.
    pub fn simplices_of_dim(&self, dim: usize) -> Vec<BTreeSet<Vertex>> {
        let mut result: Vec<_> = self.simplices.get(dim)
            .map(|s| s.iter().cloned().collect())
            .unwrap_or_default();
        result.sort();
        result
    }

    /// Build the boundary matrix for dimension n+1 -> n.
    /// boundary_maps[n] represents d_{n+1}: C_{n+1} -> C_n.
    fn boundary_matrix(&self, dim: usize) -> DMatrix<i64> {
        // Boundary map d_{dim+1}: C_{dim+1} -> C_{dim}
        let targets = self.simplices_of_dim(dim);
        let sources = self.simplices_of_dim(dim + 1);

        if sources.is_empty() || targets.is_empty() {
            return DMatrix::zeros(targets.len(), sources.len());
        }

        // Build index maps
        let target_idx: HashMap<BTreeSet<Vertex>, usize> = targets.iter()
            .enumerate()
            .map(|(i, s)| (s.clone(), i))
            .collect();

        let mut mat = DMatrix::zeros(targets.len(), sources.len());

        for (j, simplex) in sources.iter().enumerate() {
            let vertices: Vec<_> = simplex.iter().copied().collect();
            // Boundary = sum over faces with sign (-1)^i
            for (i, _) in vertices.iter().enumerate() {
                let mut face = simplex.clone();
                face.remove(&vertices[i]);
                if let Some(&row) = target_idx.get(&face) {
                    mat[(row, j)] += if i % 2 == 0 { 1 } else { -1 };
                }
            }
        }

        mat
    }

    /// Build the chain complex for this simplicial complex.
    pub fn chain_complex(&self) -> ChainComplex {
        let max_dim = self.dimension();
        let mut boundary_maps = Vec::new();

        let mut ranks = Vec::new();
        for d in 0..=max_dim {
            ranks.push(self.num_simplices(d));
        }

        for d in 0..max_dim {
            boundary_maps.push(self.boundary_matrix(d));
        }

        ChainComplex::from_ranks_and_maps(ranks, boundary_maps)
    }

    /// Compute homology group H_n.
    pub fn homology(&self, n: usize) -> HomologyGroup {
        self.chain_complex().homology(n)
    }

    /// Compute all homology groups up to the dimension of the complex.
    pub fn all_homology(&self) -> Vec<HomologyGroup> {
        self.chain_complex().all_homology(self.dimension())
    }

    /// Compute Betti numbers.
    pub fn betti_numbers(&self) -> Vec<usize> {
        self.chain_complex().betti_numbers(self.dimension())
    }

    /// Euler characteristic: χ = Σ (-1)^n * (number of n-simplices).
    pub fn euler_characteristic(&self) -> i64 {
        let mut chi = 0i64;
        for d in 0..=self.dimension() {
            chi += if d % 2 == 0 { 1 } else { -1 } * self.num_simplices(d) as i64;
        }
        chi
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.num_simplices(1)
    }
}

// ============ Construction helpers for well-known spaces ============

/// Build a simplicial complex homeomorphic to S^1 (circle).
pub fn circle() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("S^1");
    // Triangle: 3 vertices, 3 edges
    k.add_simplex(&[0, 1]);
    k.add_simplex(&[1, 2]);
    k.add_simplex(&[0, 2]);
    k
}

/// Build S^2 (2-sphere) as the boundary of a tetrahedron.
pub fn sphere_2() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("S^2");
    // 4 vertices, 4 triangles (boundary of tetrahedron)
    k.add_simplex(&[0, 1, 2]);
    k.add_simplex(&[0, 1, 3]);
    k.add_simplex(&[0, 2, 3]);
    k.add_simplex(&[1, 2, 3]);
    k
}

/// Build S^n for arbitrary n using the boundary of an (n+1)-simplex.
pub fn sphere_n(n: usize) -> SimplicialComplex {
    let mut k = SimplicialComplex::new(&format!("S^{}", n));
    // The boundary of an (n+1)-simplex is homeomorphic to S^n
    // We need n+2 vertices
    let vertices: Vec<usize> = (0..=n + 1).collect();
    // Add all (n)-simplices (i.e., all faces of dimension n of the (n+1)-simplex)
    // Each is obtained by removing one vertex
    for i in 0..=n + 1 {
        let face: Vec<usize> = vertices.iter().copied().filter(|&v| v != i).collect();
        k.add_simplex(&face);
    }
    k
}

/// Build a torus T^2 as a simplicial complex (minimal triangulation).
pub fn torus() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("T^2");
    // 7-vertex triangulation of the torus
    // Using 7 vertices with specific triangles
    let triangles = [
        [0, 1, 3], [1, 3, 4], [1, 2, 4], [2, 4, 5], [0, 2, 5],
        [0, 3, 5], [3, 4, 6], [3, 5, 6], [4, 5, 6], [0, 1, 6],
        [1, 2, 6], [0, 2, 6],
    ];
    for tri in &triangles {
        k.add_simplex(tri);
    }
    k
}

/// Build the real projective plane RP².
pub fn rp2() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("RP²");
    // 6-vertex triangulation of RP²
    let triangles = [
        [0, 1, 2], [0, 1, 3], [0, 2, 4], [0, 3, 4],
        [1, 2, 5], [1, 3, 5], [2, 4, 5], [3, 4, 5],
    ];
    for tri in &triangles {
        k.add_simplex(tri);
    }
    k
}

/// Build the Klein bottle.
pub fn klein_bottle() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("Klein bottle");
    // Triangulation of the Klein bottle (8 vertices)
    let triangles = [
        [0, 1, 5], [1, 5, 6], [1, 2, 6], [2, 6, 7],
        [0, 2, 7], [0, 4, 5], [4, 5, 6], [3, 4, 6],
        [3, 6, 7], [0, 3, 7], [0, 1, 3], [1, 3, 4],
        [1, 4, 5], [2, 5, 6], [2, 6, 7], [0, 2, 7],
        [0, 2, 3], [2, 3, 4], [0, 4, 7], [4, 6, 7],
    ];
    for tri in &triangles {
        k.add_simplex(tri);
    }
    k
}

/// Build the Möbius strip.
pub fn mobius_strip() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("Möbius strip");
    // Möbius strip triangulation
    let triangles = [
        [0, 1, 3], [1, 3, 4], [1, 2, 4], [2, 4, 5], [0, 2, 5],
    ];
    for tri in &triangles {
        k.add_simplex(tri);
    }
    k
}

/// Build a single point.
pub fn point() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("pt");
    k.add_simplex(&[0]);
    k
}

/// Build a line segment (contractible).
pub fn interval() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("I");
    k.add_simplex(&[0, 1]);
    k
}

/// Build a solid triangle (disk, contractible).
pub fn disk() -> SimplicialComplex {
    let mut k = SimplicialComplex::new("D^2");
    k.add_simplex(&[0, 1, 2]);
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_homology() {
        let p = point();
        let h0 = p.homology(0);
        assert_eq!(h0.free_rank, 1);
        assert!(h0.torsion.is_empty());
        let h1 = p.homology(1);
        assert!(h1.is_trivial());
    }

    #[test]
    fn test_circle_homology() {
        let s1 = circle();
        let h0 = s1.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = s1.homology(1);
        assert_eq!(h1.free_rank, 1);
        let h2 = s1.homology(2);
        assert!(h2.is_trivial());
    }

    #[test]
    fn test_s2_homology() {
        let s2 = sphere_2();
        let h0 = s2.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = s2.homology(1);
        assert!(h1.is_trivial());
        let h2 = s2.homology(2);
        assert_eq!(h2.free_rank, 1);
    }

    #[test]
    fn test_sphere_n_homology() {
        // S^3
        let s3 = sphere_n(3);
        let h0 = s3.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = s3.homology(1);
        assert!(h1.is_trivial());
        let h2 = s3.homology(2);
        assert!(h2.is_trivial());
        let h3 = s3.homology(3);
        assert_eq!(h3.free_rank, 1);
    }

    #[test]
    fn test_disk_contractible() {
        let d = disk();
        let h0 = d.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = d.homology(1);
        assert!(h1.is_trivial());
    }

    #[test]
    fn test_circle_euler() {
        let s1 = circle();
        // 3 vertices, 3 edges: χ = 3 - 3 = 0
        assert_eq!(s1.euler_characteristic(), 0);
    }

    #[test]
    fn test_s2_euler() {
        let s2 = sphere_2();
        // 4 vertices, 6 edges, 4 triangles: χ = 4 - 6 + 4 = 2
        assert_eq!(s2.euler_characteristic(), 2);
    }

    #[test]
    fn test_interval_contractible() {
        let i = interval();
        let h0 = i.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = i.homology(1);
        assert!(h1.is_trivial());
    }

    #[test]
    fn test_mobius_strip_homology() {
        let m = mobius_strip();
        let h0 = m.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = m.homology(1);
        assert_eq!(h1.free_rank, 1); // Z for Möbius strip H_1
    }
}
