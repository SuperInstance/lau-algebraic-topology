//! Chain complexes, boundary operators, and homology computation.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Element of Z (integer ring).
pub type Z = i64;

/// Type alias for dynamic matrix over Z.
pub type ZMatrix = DMatrix<Z>;

/// A chain complex over Z.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ChainComplex {
    /// Boundary matrices: boundary_maps[k] is d_{k+1}: C_{k+1} -> C_k.
    pub boundary_maps: Vec<ZMatrix>,
    /// Rank (number of generators) at each dimension.
    pub ranks: Vec<usize>,
}

impl ChainComplex {
    /// Create a chain complex with given ranks and boundary maps.
    pub fn from_ranks_and_maps(ranks: Vec<usize>, boundary_maps: Vec<ZMatrix>) -> Self {
        ChainComplex { boundary_maps, ranks }
    }

    /// Create a trivial chain complex.
    pub fn trivial() -> Self {
        ChainComplex { boundary_maps: vec![], ranks: vec![] }
    }

    #[allow(dead_code)]
    fn gcd(a: Z, b: Z) -> Z {
        let (mut a, mut b) = (a.abs(), b.abs());
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }

    /// Compute the Smith normal form of a matrix over Z.
    pub fn smith_normal_form(mat: &ZMatrix) -> (ZMatrix, ZMatrix, ZMatrix) {
        let m = mat.nrows();
        let n = mat.ncols();
        let mut s = mat.clone();
        let mut u: ZMatrix = DMatrix::identity(m, m);
        let mut v: ZMatrix = DMatrix::identity(n, n);

        let min_dim = m.min(n);

        for k in 0..min_dim {
            // Repeatedly reduce until position (k,k) divides everything in the submatrix
            loop {
                // Find the smallest nonzero entry in the submatrix s[k..m, k..n]
                let mut pivot_row = k;
                let mut pivot_col = k;
                let mut pivot_val = 0i64;

                for i in k..m {
                    for j in k..n {
                        let val = s[(i, j)];
                        if val != 0 && (pivot_val == 0 || val.abs() < pivot_val.abs()) {
                            pivot_val = val;
                            pivot_row = i;
                            pivot_col = j;
                        }
                    }
                }

                if pivot_val == 0 {
                    break; // All zeros in submatrix
                }

                // Move pivot to (k,k)
                if pivot_row != k {
                    for j in 0..n {
                        let tmp = s[(k, j)];
                        s[(k, j)] = s[(pivot_row, j)];
                        s[(pivot_row, j)] = tmp;
                    }
                    for j in 0..m {
                        let tmp = u[(k, j)];
                        u[(k, j)] = u[(pivot_row, j)];
                        u[(pivot_row, j)] = tmp;
                    }
                }
                if pivot_col != k {
                    for i in 0..m {
                        let tmp = s[(i, k)];
                        s[(i, k)] = s[(i, pivot_col)];
                        s[(i, pivot_col)] = tmp;
                    }
                    for i in 0..n {
                        let tmp = v[(i, k)];
                        v[(i, k)] = v[(i, pivot_col)];
                        v[(i, pivot_col)] = tmp;
                    }
                }

                // Make pivot positive
                if s[(k, k)] < 0 {
                    for j in 0..n {
                        s[(k, j)] = -s[(k, j)];
                    }
                    for j in 0..m {
                        u[(k, j)] = -u[(k, j)];
                    }
                }

                // Try to eliminate all entries in row k and column k
                let mut all_divide = true;

                // Eliminate column k entries below pivot
                for i in (k + 1)..m {
                    if s[(i, k)] != 0 {
                        let d = s[(k, k)];
                        if d == 0 { continue; }
                        let q = s[(i, k)] / d;
                        for j in 0..n {
                            s[(i, j)] -= q * s[(k, j)];
                        }
                        for j in 0..m {
                            u[(i, j)] -= q * u[(k, j)];
                        }
                        if s[(i, k)] != 0 {
                            all_divide = false;
                        }
                    }
                }

                // Eliminate row k entries right of pivot
                for j in (k + 1)..n {
                    if s[(k, j)] != 0 {
                        let d = s[(k, k)];
                        if d == 0 { continue; }
                        let q = s[(k, j)] / d;
                        for i in 0..m {
                            s[(i, j)] -= q * s[(i, k)];
                        }
                        for i in 0..n {
                            v[(i, j)] -= q * v[(i, k)];
                        }
                        if s[(k, j)] != 0 {
                            all_divide = false;
                        }
                    }
                }

                if all_divide {
                    // Check that s[(k,k)] divides every entry in the remaining submatrix
                    let d = s[(k, k)];
                    let mut divides_all = true;
                    for i in (k + 1)..m {
                        for j in (k + 1)..n {
                            if s[(i, j)] != 0 && s[(i, j)] % d != 0 {
                                // Add row i to row k to bring a non-divisible entry into row k
                                for jj in 0..n {
                                    s[(k, jj)] += s[(i, jj)];
                                }
                                for jj in 0..m {
                                    u[(k, jj)] += u[(i, jj)];
                                }
                                divides_all = false;
                                break;
                            }
                        }
                        if !divides_all { break; }
                    }
                    if divides_all {
                        break; // Done with this position
                    }
                }
                // Otherwise loop again to reduce further
            }
        }

        // Normalize: make all diagonal entries positive
        for k in 0..min_dim {
            if s[(k, k)] < 0 {
                for j in 0..n {
                    s[(k, j)] = -s[(k, j)];
                }
                for j in 0..m {
                    u[(k, j)] = -u[(k, j)];
                }
            }
        }

        (u, s, v)
    }

    /// Get diagonal entries from a possibly non-square matrix.
    fn diagonal_entries(mat: &ZMatrix) -> Vec<Z> {
        let min_dim = mat.nrows().min(mat.ncols());
        (0..min_dim).map(|i| mat[(i, i)]).collect()
    }

    /// Compute dimension of kernel of a matrix over Z.
    pub fn kernel_dimension(mat: &ZMatrix) -> usize {
        let (_, s, _) = Self::smith_normal_form(mat);
        let diag = Self::diagonal_entries(&s);
        let rank = diag.iter().filter(|&&x| x != 0).count();
        mat.ncols() - rank
    }

    /// Compute rank of image of a matrix.
    pub fn image_rank(mat: &ZMatrix) -> usize {
        let (_, s, _) = Self::smith_normal_form(mat);
        let diag = Self::diagonal_entries(&s);
        diag.iter().filter(|&&x| x != 0).count()
    }

    /// Compute image rank and torsion coefficients.
    fn image_rank_and_torsion(mat: &ZMatrix) -> (usize, Vec<Z>) {
        let (_, s, _) = Self::smith_normal_form(mat);
        let diag = Self::diagonal_entries(&s);
        let mut torsion = Vec::new();
        let mut rank = 0;
        for &d in &diag {
            if d == 0 {
                break;
            }
            if d == 1 {
                rank += 1;
            } else {
                rank += 1;
                torsion.push(d);
            }
        }
        (rank, torsion)
    }

    /// Compute homology group H_n.
    /// boundary_maps[k] = d_{k+1}: C_{k+1} -> C_k
    /// H_n = ker(d_n) / im(d_{n+1})
    /// d_n = boundary_maps[n-1] for n >= 1, d_0 = 0
    pub fn homology(&self, n: usize) -> HomologyGroup {
        let rank_cn = self.ranks.get(n).copied().unwrap_or(0);
        if rank_cn == 0 {
            return HomologyGroup { dimension: n, free_rank: 0, torsion: vec![] };
        }

        // ker(d_n)
        let ker_dim = if n == 0 {
            rank_cn // d_0 = 0
        } else if n - 1 < self.boundary_maps.len() {
            Self::kernel_dimension(&self.boundary_maps[n - 1])
        } else {
            rank_cn
        };

        // im(d_{n+1}) and torsion
        let (img_dim, torsion) = if n < self.boundary_maps.len() {
            Self::image_rank_and_torsion(&self.boundary_maps[n])
        } else {
            (0, vec![])
        };

        HomologyGroup {
            dimension: n,
            free_rank: ker_dim as i64 - img_dim as i64,
            torsion,
        }
    }

    /// Compute all homology groups up to max_n.
    pub fn all_homology(&self, max_n: usize) -> Vec<HomologyGroup> {
        (0..=max_n).map(|n| self.homology(n)).collect()
    }

    /// Compute Betti numbers.
    pub fn betti_numbers(&self, max_n: usize) -> Vec<usize> {
        self.all_homology(max_n)
            .iter()
            .map(|h| h.free_rank.max(0) as usize)
            .collect()
    }

    /// Check validity: d_k ∘ d_{k+1} = 0.
    pub fn is_valid(&self) -> bool {
        for k in 1..self.boundary_maps.len() {
            let product = &self.boundary_maps[k] * &self.boundary_maps[k - 1];
            if product.iter().any(|&x| x != 0) {
                return false;
            }
        }
        true
    }
}

/// A homology group H_n.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HomologyGroup {
    pub dimension: usize,
    pub free_rank: i64,
    pub torsion: Vec<Z>,
}

impl HomologyGroup {
    pub fn is_trivial(&self) -> bool {
        self.free_rank == 0 && self.torsion.is_empty()
    }

    pub fn rank(&self) -> usize {
        self.free_rank.max(0) as usize
    }
}

impl fmt::Display for HomologyGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_trivial() {
            write!(f, "H_{} = 0", self.dimension)
        } else {
            write!(f, "H_{} = Z^{}", self.dimension, self.free_rank)?;
            for t in &self.torsion {
                write!(f, " ⊕ Z/{}", t)?;
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trivial_complex() {
        let cc = ChainComplex::trivial();
        assert!(cc.homology(0).is_trivial());
    }

    #[test]
    fn test_smith_normal_form_identity() {
        let mat: ZMatrix = DMatrix::from_row_slice(2, 2, &[1, 0, 0, 1]);
        let (_, s, _) = ChainComplex::smith_normal_form(&mat);
        let d = ChainComplex::diagonal_entries(&s);
        assert_eq!(d, vec![1, 1]);
    }

    #[test]
    fn test_smith_normal_form_diagonal() {
        let mat: ZMatrix = DMatrix::from_row_slice(2, 2, &[2, 0, 0, 3]);
        let (_, s, _) = ChainComplex::smith_normal_form(&mat);
        let d = ChainComplex::diagonal_entries(&s);
        assert_eq!(d, vec![1, 6]);
    }

    #[test]
    fn test_kernel_dimension() {
        let mat: ZMatrix = DMatrix::from_row_slice(1, 2, &[1, -1]);
        assert_eq!(ChainComplex::kernel_dimension(&mat), 1);
    }

    #[test]
    fn test_image_rank() {
        let mat: ZMatrix = DMatrix::from_row_slice(1, 2, &[1, -1]);
        assert_eq!(ChainComplex::image_rank(&mat), 1);
    }
}
