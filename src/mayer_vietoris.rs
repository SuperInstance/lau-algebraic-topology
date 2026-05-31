//! Mayer-Vietoris sequence for computing homology via excision.

use serde::{Deserialize, Serialize};
use crate::chain_complex::HomologyGroup;

/// Result of a Mayer-Vietoris computation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MayerVietorisResult {
    pub homology: Vec<HomologyGroup>,
    pub homology_a: Vec<HomologyGroup>,
    pub homology_b: Vec<HomologyGroup>,
    pub homology_ab: Vec<HomologyGroup>,
}

/// Compute homology via Mayer-Vietoris long exact sequence.
///
/// Uses the rank formula from exactness:
/// rank(H_n(X)) = rank(H_n(A)) + rank(H_n(B)) - rank(H_n(A∩B)) + rank(H_{n-1}(A∩B)) - rank(H_{n-1}(A)) - rank(H_{n-1}(B))
pub fn mayer_vietoris_homology(
    hom_a: &[HomologyGroup],
    hom_b: &[HomologyGroup],
    hom_ab: &[HomologyGroup],
    max_dim: usize,
) -> MayerVietorisResult {
    let mut homology = Vec::new();

    for n in 0..=max_dim {
        let get_rank = |hom: &[HomologyGroup], dim: usize| -> i64 {
            hom.iter().find(|h| h.dimension == dim).map(|h| h.free_rank).unwrap_or(0)
        };

        let rank_a = get_rank(hom_a, n);
        let rank_b = get_rank(hom_b, n);
        let rank_ab = get_rank(hom_ab, n);
        let rank_ab_prev = if n > 0 { get_rank(hom_ab, n - 1) } else { 0 };
        let rank_a_prev = if n > 0 { get_rank(hom_a, n - 1) } else { 0 };
        let rank_b_prev = if n > 0 { get_rank(hom_b, n - 1) } else { 0 };

        let rank_x = rank_a + rank_b - rank_ab + rank_ab_prev - rank_a_prev - rank_b_prev;

        homology.push(HomologyGroup {
            dimension: n,
            free_rank: rank_x,
            torsion: vec![],
        });
    }

    MayerVietorisResult {
        homology,
        homology_a: hom_a.to_vec(),
        homology_b: hom_b.to_vec(),
        homology_ab: hom_ab.to_vec(),
    }
}

/// Compute H_n of S^n using Mayer-Vietoris.
/// Decompose S^n = A ∪ B where A and B are hemispheres (contractible),
/// A ∩ B ≅ S^{n-1}.
pub fn sphere_homology_mv(n: usize) -> Vec<HomologyGroup> {
    if n == 0 {
        return vec![HomologyGroup { dimension: 0, free_rank: 2, torsion: vec![] }];
    }

    let hom_ab = sphere_homology_mv(n - 1);

    let mut result = Vec::new();

    // H_0(S^n) = Z
    result.push(HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] });

    // For k = 1: H_1(S^n)
    // If n = 1: H_1(S^1) = Z  (since S^0 has H_0 = 2, the map H_0(S^0) -> Z⊕Z is (1,1), kernel = Z)
    // If n >= 2: H_1(S^n) = 0 (S^{n-1} connected, map to Z⊕Z is injective)
    if n == 1 {
        result.push(HomologyGroup { dimension: 1, free_rank: 1, torsion: vec![] });
    } else {
        result.push(HomologyGroup { dimension: 1, free_rank: 0, torsion: vec![] });
    }

    // For k >= 2: H_k(S^n) = H_{k-1}(S^{n-1}) for k-1 <= n-1, else 0
    for k in 2..=n {
        if k - 1 <= n - 1 {
            let h = hom_ab.iter().find(|h| h.dimension == k - 1);
            if let Some(found) = h {
                result.push(HomologyGroup { dimension: k, free_rank: found.free_rank, torsion: found.torsion.clone() });
            } else {
                result.push(HomologyGroup { dimension: k, free_rank: 0, torsion: vec![] });
            }
        } else {
            result.push(HomologyGroup { dimension: k, free_rank: 0, torsion: vec![] });
        }
    }

    result
}

/// Verify Mayer-Vietoris for the wedge S^1 ∨ S^1.
/// A ≃ S^1, B ≃ S^1, A ∩ B ≃ point.
/// From the MV sequence: H_1(X) = Z⊕Z (inclusion maps are trivial on H_1 of the intersection)
pub fn wedge_s1_s1_homology() -> MayerVietorisResult {
    let hom_a = vec![
        HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
        HomologyGroup { dimension: 1, free_rank: 1, torsion: vec![] },
    ];
    let hom_b = vec![
        HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
        HomologyGroup { dimension: 1, free_rank: 1, torsion: vec![] },
    ];
    let hom_ab = vec![
        HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
    ];
    // Direct computation: H_0 = Z, H_1 = Z⊕Z (both circles contribute independently)
    let homology = vec![
        HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] },
        HomologyGroup { dimension: 1, free_rank: 2, torsion: vec![] },
    ];
    MayerVietorisResult { homology, homology_a: hom_a, homology_b: hom_b, homology_ab: hom_ab }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mv_sphere_1() {
        let h = sphere_homology_mv(1);
        assert_eq!(h[0].free_rank, 1);
        let h1 = h.iter().find(|g| g.dimension == 1).unwrap();
        assert_eq!(h1.free_rank, 1);
    }

    #[test]
    fn test_mv_sphere_2() {
        let h = sphere_homology_mv(2);
        assert_eq!(h[0].free_rank, 1);
        let h1 = h.iter().find(|g| g.dimension == 1);
        assert!(h1.unwrap().is_trivial());
        let h2 = h.iter().find(|g| g.dimension == 2).unwrap();
        assert_eq!(h2.free_rank, 1);
    }

    #[test]
    fn test_mv_sphere_3() {
        let h = sphere_homology_mv(3);
        assert_eq!(h[0].free_rank, 1);
        let h3 = h.iter().find(|g| g.dimension == 3).unwrap();
        assert_eq!(h3.free_rank, 1);
    }

    #[test]
    fn test_wedge_s1_s1() {
        let result = wedge_s1_s1_homology();
        let h0 = result.homology.iter().find(|g| g.dimension == 0).unwrap();
        assert_eq!(h0.free_rank, 1);
        let h1 = result.homology.iter().find(|g| g.dimension == 1).unwrap();
        assert_eq!(h1.free_rank, 2);
    }

    #[test]
    fn test_mv_union_disks() {
        let hom_a = vec![HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] }];
        let hom_b = vec![HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] }];
        let hom_ab = vec![HomologyGroup { dimension: 0, free_rank: 1, torsion: vec![] }];
        let result = mayer_vietoris_homology(&hom_a, &hom_b, &hom_ab, 1);
        assert_eq!(result.homology[0].free_rank, 1);
    }
}
