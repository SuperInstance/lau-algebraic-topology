//! Singular homology: singular simplices and chain maps.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use crate::chain_complex::{ChainComplex, HomologyGroup};

/// A singular n-simplex is a continuous map σ: Δ^n → X.
/// We represent it abstractly by its label and the space it maps into.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SingularSimplex {
    /// Dimension of the simplex.
    pub dimension: usize,
    /// Label/identifier.
    pub label: String,
    /// Vertices (as indices into the space).
    pub vertices: Vec<usize>,
}

impl SingularSimplex {
    /// Create a new singular simplex.
    pub fn new(dimension: usize, label: &str, vertices: Vec<usize>) -> Self {
        SingularSimplex {
            dimension,
            label: label.to_string(),
            vertices,
        }
    }

    /// The i-th face map: returns the face with vertex i removed.
    pub fn face(&self, i: usize) -> SingularSimplex {
        let mut verts = self.vertices.clone();
        verts.remove(i);
        SingularSimplex::new(
            self.dimension.saturating_sub(1),
            &format!("∂_{}({})", i, self.label),
            verts,
        )
    }
}

/// A singular chain: a formal linear combination of singular simplices.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SingularChain {
    /// Dimension.
    pub dimension: usize,
    /// Coefficients indexed by simplex index.
    pub coefficients: Vec<i64>,
}

impl SingularChain {
    /// Create a zero chain of given dimension with n generators.
    pub fn zero(dimension: usize, n: usize) -> Self {
        SingularChain {
            dimension,
            coefficients: vec![0; n],
        }
    }

    /// Add two chains.
    pub fn add(&self, other: &SingularChain) -> SingularChain {
        let coeffs: Vec<i64> = self.coefficients.iter()
            .zip(other.coefficients.iter())
            .map(|(&a, &b)| a + b)
            .collect();
        SingularChain {
            dimension: self.dimension,
            coefficients: coeffs,
        }
    }

    /// Scale by integer.
    pub fn scale(&self, n: i64) -> SingularChain {
        SingularChain {
            dimension: self.dimension,
            coefficients: self.coefficients.iter().map(|&c| c * n).collect(),
        }
    }
}

/// A singular chain complex for a space represented abstractly.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SingularChainComplex {
    /// The underlying chain complex.
    pub complex: ChainComplex,
    /// Labels for simplices at each dimension.
    pub labels: Vec<Vec<String>>,
}

impl SingularChainComplex {
    /// Create from a pre-built chain complex.
    pub fn from_chain_complex(complex: ChainComplex, labels: Vec<Vec<String>>) -> Self {
        SingularChainComplex { complex, labels }
    }

    /// Compute H_n.
    pub fn homology(&self, n: usize) -> HomologyGroup {
        self.complex.homology(n)
    }

    /// Compute all homology groups.
    pub fn all_homology(&self, max_n: usize) -> Vec<HomologyGroup> {
        self.complex.all_homology(max_n)
    }
}

/// A chain map between singular chain complexes.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ChainMap {
    /// Maps at each dimension: matrix representing the linear map.
    pub maps: Vec<DMatrix<i64>>,
}

impl ChainMap {
    /// Create a new chain map.
    pub fn new(maps: Vec<DMatrix<i64>>) -> Self {
        ChainMap { maps }
    }

    /// The induced map on homology H_n.
    /// For simplicity, we compute this as a matrix acting on representatives.
    pub fn induced_homology_matrix(&self, _n: usize) -> Option<DMatrix<i64>> {
        // This is a simplified version
        self.maps.first().cloned()
    }
}

/// Build a singular chain complex equivalent to a simplicial one.
/// For testing, we just delegate to the simplicial computation.
pub fn singular_homology_of_space(complex: &ChainComplex, n: usize) -> HomologyGroup {
    complex.homology(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplicial;

    #[test]
    fn test_singular_simplex_face() {
        let sigma = SingularSimplex::new(2, "σ", vec![0, 1, 2]);
        let f0 = sigma.face(0);
        assert_eq!(f0.vertices, vec![1, 2]);
        let f1 = sigma.face(1);
        assert_eq!(f1.vertices, vec![0, 2]);
    }

    #[test]
    fn test_singular_chain_add() {
        let a = SingularChain { dimension: 1, coefficients: vec![1, 2, 3] };
        let b = SingularChain { dimension: 1, coefficients: vec![4, 5, 6] };
        let c = a.add(&b);
        assert_eq!(c.coefficients, vec![5, 7, 9]);
    }

    #[test]
    fn test_singular_homology_circle() {
        let s1 = simplicial::circle();
        let cc = s1.chain_complex();
        let h0 = singular_homology_of_space(&cc, 0);
        assert_eq!(h0.free_rank, 1);
        let h1 = singular_homology_of_space(&cc, 1);
        assert_eq!(h1.free_rank, 1);
    }

    #[test]
    fn test_singular_chain_complex() {
        let s2 = simplicial::sphere_2();
        let cc = s2.chain_complex();
        let scc = SingularChainComplex::from_chain_complex(
            cc,
            vec![vec!["v0".into(), "v1".into(), "v2".into(), "v3".into()], vec![], vec![]],
        );
        let h0 = scc.homology(0);
        assert_eq!(h0.free_rank, 1);
    }
}
