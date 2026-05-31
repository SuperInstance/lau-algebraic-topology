//! Euler characteristic computation from Betti numbers.

use serde::{Deserialize, Serialize};

/// Compute Euler characteristic from Betti numbers.
/// χ = Σ_{n=0}^{∞} (-1)^n β_n
pub fn euler_from_betti(betti: &[usize]) -> i64 {
    betti.iter()
        .enumerate()
        .map(|(n, &b)| if n % 2 == 0 { b as i64 } else { -(b as i64) })
        .sum()
}

/// Compute Euler characteristic from homology ranks.
pub fn euler_from_homology_ranks(ranks: &[i64]) -> i64 {
    ranks.iter()
        .enumerate()
        .map(|(n, &r)| if n % 2 == 0 { r } else { -r })
        .sum()
}

/// Euler characteristic info for a space.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EulerInfo {
    /// Name of the space.
    pub space: String,
    /// Euler characteristic.
    pub euler: i64,
    /// Betti numbers.
    pub betti: Vec<usize>,
}

impl EulerInfo {
    /// Euler characteristic of S^n.
    pub fn sphere(n: usize) -> Self {
        let mut betti = vec![0usize; n + 1];
        if n == 0 {
            betti[0] = 2; // S^0 = two points
        } else {
            betti[0] = 1;
            betti[n] = 1;
        }
        EulerInfo {
            space: format!("S^{}", n),
            euler: euler_from_betti(&betti),
            betti,
        }
    }

    /// Euler characteristic of T^n (n-torus).
    pub fn torus(n: usize) -> Self {
        // β_k(T^n) = C(n, k)
        let betti: Vec<usize> = (0..=n)
            .map(|k| binomial(n, k))
            .collect();
        EulerInfo {
            space: format!("T^{}", n),
            euler: euler_from_betti(&betti),
            betti,
        }
    }

    /// Euler characteristic of RP².
    pub fn rp2() -> Self {
        EulerInfo {
            space: "RP²".into(),
            euler: 1,
            betti: vec![1, 0, 0], // Over Z: β_0=1, β_1=0, β_2=0
        }
    }

    /// Euler characteristic of the Klein bottle.
    pub fn klein_bottle() -> Self {
        EulerInfo {
            space: "Klein bottle".into(),
            euler: 0,
            betti: vec![1, 1, 0],
        }
    }

    /// Euler characteristic of a point.
    pub fn point() -> Self {
        EulerInfo {
            space: "pt".into(),
            euler: 1,
            betti: vec![1],
        }
    }

    /// Euler characteristic of CP^n.
    pub fn cpn(n: usize) -> Self {
        let betti: Vec<usize> = (0..=2 * n)
            .map(|k| if k % 2 == 0 && k <= 2 * n { 1 } else { 0 })
            .collect();
        EulerInfo {
            space: format!("CP^{}", n),
            euler: n as i64 + 1,
            betti,
        }
    }
}

fn binomial(n: usize, k: usize) -> usize {
    if k > n { return 0; }
    let k = k.min(n - k);
    let mut result = 1usize;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_sphere_0() {
        let info = EulerInfo::sphere(0);
        assert_eq!(info.euler, 2); // Two points
    }

    #[test]
    fn test_euler_sphere_1() {
        let info = EulerInfo::sphere(1);
        assert_eq!(info.euler, 0);
    }

    #[test]
    fn test_euler_sphere_2() {
        let info = EulerInfo::sphere(2);
        assert_eq!(info.euler, 2);
    }

    #[test]
    fn test_euler_sphere_even() {
        for n in (0..10).step_by(2) {
            let info = EulerInfo::sphere(n);
            assert_eq!(info.euler, 2, "S^{} should have χ=2", n);
        }
    }

    #[test]
    fn test_euler_sphere_odd() {
        for n in (1..10).step_by(2) {
            let info = EulerInfo::sphere(n);
            assert_eq!(info.euler, 0, "S^{} should have χ=0", n);
        }
    }

    #[test]
    fn test_euler_torus_2() {
        let info = EulerInfo::torus(2);
        assert_eq!(info.euler, 0); // 1 - 2 + 1 = 0
    }

    #[test]
    fn test_euler_torus_3() {
        let info = EulerInfo::torus(3);
        // β = (1, 3, 3, 1) → χ = 1 - 3 + 3 - 1 = 0
        assert_eq!(info.euler, 0);
    }

    #[test]
    fn test_euler_rp2() {
        let info = EulerInfo::rp2();
        assert_eq!(info.euler, 1);
    }

    #[test]
    fn test_euler_klein() {
        let info = EulerInfo::klein_bottle();
        assert_eq!(info.euler, 0);
    }

    #[test]
    fn test_euler_cp2() {
        let info = EulerInfo::cpn(2);
        assert_eq!(info.euler, 3); // 1 + 1 + 1 = 3
    }

    #[test]
    fn test_euler_from_betti() {
        assert_eq!(euler_from_betti(&[1, 0, 1]), 2);
        assert_eq!(euler_from_betti(&[1, 2, 1]), 0);
        assert_eq!(euler_from_betti(&[1, 1]), 0);
    }
}
