//! Cohomology with coefficients, cup product, and Künneth formula.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use crate::chain_complex::{ChainComplex, HomologyGroup};

/// Cohomology group H^n(X; Z).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CohomologyGroup {
    /// Dimension n.
    pub dimension: usize,
    /// Free rank.
    pub free_rank: i64,
    /// Torsion coefficients.
    pub torsion: Vec<i64>,
}

impl CohomologyGroup {
    /// Check if trivial.
    pub fn is_trivial(&self) -> bool {
        self.free_rank == 0 && self.torsion.is_empty()
    }
}

/// Compute cohomology from homology using universal coefficient theorem:
/// H^n(X; Z) ≅ Hom(H_n(X), Z) ⊕ Ext(H_{n-1}(X), Z)
/// For free H_n: H^n ≅ Z^{rank(H_n)} ⊕ torsion from H_{n-1}.
pub fn cohomology_from_homology(hom_n: &HomologyGroup, hom_n_minus_1: &HomologyGroup) -> CohomologyGroup {
    let free_rank = hom_n.free_rank;
    // Torsion of H_{n-1} becomes torsion in H^n via Ext
    let torsion = hom_n_minus_1.torsion.clone();
    CohomologyGroup {
        dimension: hom_n.dimension,
        free_rank,
        torsion,
    }
}

/// Compute all cohomology groups from homology groups.
pub fn all_cohomology(homology: &[HomologyGroup]) -> Vec<CohomologyGroup> {
    let mut result = Vec::new();
    for (i, h) in homology.iter().enumerate() {
        let h_prev = if i > 0 {
            &homology[i - 1]
        } else {
            &HomologyGroup { dimension: 0, free_rank: 0, torsion: vec![] }
        };
        result.push(cohomology_from_homology(h, h_prev));
    }
    result
}

/// Cup product structure.
/// For a space X, the cup product ⌣: H^p × H^q → H^{p+q}.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CupProduct {
    /// Source dimension p.
    pub p: usize,
    /// Source dimension q.
    pub q: usize,
    /// Target dimension p+q.
    pub target_dim: usize,
    /// The product matrix: maps H^p ⊗ H^q → H^{p+q}.
    pub product_matrix: DMatrix<i64>,
}

impl CupProduct {
    /// Create a new cup product.
    pub fn new(p: usize, q: usize, matrix: DMatrix<i64>) -> Self {
        CupProduct {
            p,
            q,
            target_dim: p + q,
            product_matrix: matrix,
        }
    }

    /// Compute the cup product of two cohomology classes.
    pub fn compute(&self, a: &[i64], b: &[i64]) -> Vec<i64> {
        let a_vec = DMatrix::from_row_slice(1, a.len(), a);
        let b_vec = DMatrix::from_column_slice(b.len(), 1, b);
        // This is a simplification
        let result = &a_vec * &self.product_matrix * &b_vec;
        vec![result[(0, 0)]]
    }
}

/// Cohomology ring with cup product.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CohomologyRing {
    /// Cohomology groups.
    pub groups: Vec<CohomologyGroup>,
    /// Cup products.
    pub cup_products: Vec<CupProduct>,
}

impl CohomologyRing {
    /// Compute the cup product ring of RP².
    /// H*(RP²; Z/2) = Z/2[a]/(a³) where a ∈ H^1.
    pub fn ring_rp2_z2() -> Self {
        let groups = vec![
            CohomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
            CohomologyGroup { dimension: 1, free_rank: 0, torsion: vec![2] },
            CohomologyGroup { dimension: 2, free_rank: 0, torsion: vec![2] },
        ];
        CohomologyRing { groups, cup_products: vec![] }
    }
}

/// Künneth formula: H_n(X × Y) in terms of H_*(X) and H_*(Y).
/// H_n(X × Y) ≅ ⊕_{p+q=n} H_p(X) ⊗ H_q(Y) ⊕ ⊕_{p+q=n-1} Tor(H_p(X), H_q(Y))
pub fn kunneth(
    hom_x: &[HomologyGroup],
    hom_y: &[HomologyGroup],
    n: usize,
) -> HomologyGroup {
    let mut free_rank = 0i64;
    let mut torsion = Vec::new();

    // Free part: ⊕_{p+q=n} H_p(X) ⊗ H_q(Y)
    for p in 0..=n {
        let q = n - p;
        let hp = hom_x.iter().find(|h| h.dimension == p);
        let hq = hom_y.iter().find(|h| h.dimension == q);

        if let (Some(hp), Some(hq)) = (hp, hq) {
            // Tensor product of free parts
            free_rank += hp.free_rank * hq.free_rank;

            // Tor from torsion parts
            for t1 in &hp.torsion {
                for t2 in &hq.torsion {
                    torsion.push(gcd(*t1, *t2));
                }
            }
        }
    }

    // Tor part: ⊕_{p+q=n-1} Tor(H_p(X), H_q(Y))
    if n > 0 {
        for p in 0..n {
            let q = n - 1 - p;
            let hp = hom_x.iter().find(|h| h.dimension == p);
            let hq = hom_y.iter().find(|h| h.dimension == q);

            if let (Some(hp), Some(hq)) = (hp, hq) {
                // Tor(Z/m, Z/n) = Z/gcd(m,n)
                for t1 in &hp.torsion {
                    for t2 in &hq.torsion {
                        torsion.push(gcd(*t1, *t2));
                    }
                }
            }
        }
    }

    HomologyGroup { dimension: n, free_rank, torsion }
}

/// Compute all Künneth homology groups for X × Y.
pub fn kunneth_all(
    hom_x: &[HomologyGroup],
    hom_y: &[HomologyGroup],
    max_n: usize,
) -> Vec<HomologyGroup> {
    (0..=max_n).map(|n| kunneth(hom_x, hom_y, n)).collect()
}

fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Compute cohomology with Z/2 coefficients from a chain complex.
pub fn cohomology_z2(complex: &ChainComplex, n: usize) -> CohomologyGroup {
    // Reduce mod 2 and compute dual
    let hom = complex.homology(n);
    // Over Z/2, everything is a vector space, so no torsion distinction
    CohomologyGroup {
        dimension: n,
        free_rank: hom.free_rank + hom.torsion.len() as i64,
        torsion: vec![],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplicial;
    use crate::cw_complex;

    #[test]
    fn test_cohomology_s1() {
        let s1 = simplicial::circle();
        let cc = s1.chain_complex();
        let h0 = cc.homology(0);
        let h1 = cc.homology(1);
        let coh0 = cohomology_from_homology(&h0, &HomologyGroup { dimension: 0, free_rank: 0, torsion: vec![] });
        assert_eq!(coh0.free_rank, 1);
        let coh1 = cohomology_from_homology(&h1, &h0);
        assert_eq!(coh1.free_rank, 1);
    }

    #[test]
    fn test_cohomology_s2() {
        let s2 = simplicial::sphere_2();
        let cc = s2.chain_complex();
        let h2 = cc.homology(2);
        let h1 = cc.homology(1);
        let coh2 = cohomology_from_homology(&h2, &h1);
        assert_eq!(coh2.free_rank, 1);
    }

    #[test]
    fn test_kunneth_s1_s1() {
        // S^1 × S^1 = T^2
        let hom_s1 = vec![
            HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
            HomologyGroup { dimension: 1, free_rank: 1, torsion: vec![] },
        ];
        let result = kunneth_all(&hom_s1, &hom_s1, 2);
        assert_eq!(result[0].free_rank, 1);  // H_0 = Z
        assert_eq!(result[1].free_rank, 2);  // H_1 = Z²
        assert_eq!(result[2].free_rank, 1);  // H_2 = Z
    }

    #[test]
    fn test_kunneth_s1_s2() {
        // S^1 × S^2
        let hom_s1 = vec![
            HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
            HomologyGroup { dimension: 1, free_rank: 1, torsion: vec![] },
        ];
        let hom_s2 = vec![
            HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
            HomologyGroup { dimension: 1, free_rank: 0, torsion: vec![] },
            HomologyGroup { dimension: 2, free_rank: 1, torsion: vec![] },
        ];
        let result = kunneth_all(&hom_s1, &hom_s2, 3);
        assert_eq!(result[0].free_rank, 1);
        assert_eq!(result[1].free_rank, 1);
        assert_eq!(result[2].free_rank, 1);
        assert_eq!(result[3].free_rank, 1);
    }

    #[test]
    fn test_cup_product_creation() {
        let cp = CupProduct::new(1, 1, DMatrix::from_row_slice(1, 1, &[1]));
        assert_eq!(cp.target_dim, 2);
        let result = cp.compute(&[1], &[1]);
        assert_eq!(result, vec![1]);
    }

    #[test]
    fn test_cohomology_z2() {
        let s1 = simplicial::circle();
        let cc = s1.chain_complex();
        let coh0 = cohomology_z2(&cc, 0);
        assert_eq!(coh0.free_rank, 1);
    }

    #[test]
    fn test_cohomology_rp2() {
        let rp2 = cw_complex::cw_rp2();
        let cc = rp2.chain_complex();
        let h1 = cc.homology(1);
        // H_1(RP²) = Z/2, so H^1(RP²;Z) = Hom(Z/2, Z) ⊕ Ext(Z, Z) = 0 + 0 = 0
        // Actually by UCT: H^1 = Hom(H_1, Z) ⊕ Ext(H_0, Z) = Hom(Z/2, Z) ⊕ Ext(Z, Z) = 0
        // H^2 = Hom(H_2, Z) ⊕ Ext(H_1, Z) = 0 ⊕ Ext(Z/2, Z) = Z/2
        let coh2 = cohomology_from_homology(&cc.homology(2), &h1);
        assert!(coh2.torsion.contains(&2));
    }
}
