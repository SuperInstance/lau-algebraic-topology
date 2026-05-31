//! CW complexes: cell attachment and cellular homology.

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use crate::chain_complex::{ChainComplex, HomologyGroup};

/// A cell in a CW complex.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Cell {
    /// Dimension of the cell.
    pub dimension: usize,
    /// Label/identifier.
    pub label: String,
}

/// A CW complex.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CWComplex {
    /// Cells indexed by dimension.
    cells: Vec<Vec<Cell>>,
    /// Attaching maps encoded as boundary matrices.
    /// attaching[k] is the matrix for the boundary of (k+1)-cells -> k-cells.
    attaching_maps: Vec<DMatrix<i64>>,
    /// Name.
    pub name: String,
}

impl CWComplex {
    /// Create an empty CW complex.
    pub fn new(name: &str) -> Self {
        CWComplex {
            cells: Vec::new(),
            attaching_maps: Vec::new(),
            name: name.to_string(),
        }
    }

    /// Add a cell of given dimension.
    pub fn add_cell(&mut self, dim: usize, label: &str) {
        while self.cells.len() <= dim {
            self.cells.push(Vec::new());
        }
        self.cells[dim].push(Cell {
            dimension: dim,
            label: label.to_string(),
        });
    }

    /// Set the attaching map for dimension k+1 -> k.
    pub fn set_attaching_map(&mut self, k: usize, map: DMatrix<i64>) {
        while self.attaching_maps.len() <= k {
            self.attaching_maps.push(DMatrix::zeros(0, 0));
        }
        self.attaching_maps[k] = map;
    }

    /// Number of cells in dimension k.
    pub fn num_cells(&self, k: usize) -> usize {
        self.cells.get(k).map(|c| c.len()).unwrap_or(0)
    }

    /// Dimension of the complex.
    pub fn dimension(&self) -> usize {
        self.cells.iter().rposition(|c| !c.is_empty()).unwrap_or(0)
    }

    /// Build the cellular chain complex.
    pub fn chain_complex(&self) -> ChainComplex {
        let max_dim = self.dimension();
        let ranks: Vec<usize> = (0..=max_dim).map(|k| self.num_cells(k)).collect();
        let mut boundary_maps = Vec::new();

        for k in 0..max_dim {
            if k < self.attaching_maps.len() {
                boundary_maps.push(self.attaching_maps[k].clone());
            } else {
                let rows = self.num_cells(k);
                let cols = self.num_cells(k + 1);
                boundary_maps.push(DMatrix::zeros(rows, cols));
            }
        }

        ChainComplex::from_ranks_and_maps(ranks, boundary_maps)
    }

    /// Compute cellular homology H_n.
    pub fn homology(&self, n: usize) -> HomologyGroup {
        self.chain_complex().homology(n)
    }

    /// Euler characteristic: χ = Σ (-1)^n * c_n where c_n = number of n-cells.
    pub fn euler_characteristic(&self) -> i64 {
        let mut chi = 0i64;
        for k in 0..=self.dimension() {
            chi += if k % 2 == 0 { 1 } else { -1 } * self.num_cells(k) as i64;
        }
        chi
    }
}

// ============ Well-known CW complexes ============

/// Build S^1 as a CW complex: one 0-cell, one 1-cell.
pub fn cw_s1() -> CWComplex {
    let mut cw = CWComplex::new("S^1");
    cw.add_cell(0, "v");
    cw.add_cell(1, "e");
    // Attaching map for e: both endpoints map to v → boundary = [0]
    cw.set_attaching_map(0, DMatrix::from_row_slice(1, 1, &[0]));
    cw
}

/// Build S^2 as a CW complex: one 0-cell, one 2-cell.
pub fn cw_s2() -> CWComplex {
    let mut cw = CWComplex::new("S^2");
    cw.add_cell(0, "v");
    cw.add_cell(1, "e"); // We actually skip 1-cells for S^2
    cw.add_cell(2, "f");
    // For S^2: 0-cell, 0 1-cells, 1 2-cell
    // Let me redo this
    let mut cw = CWComplex::new("S^2");
    cw.add_cell(0, "v");
    cw.add_cell(2, "f");
    // attaching_map[0] = d_1: C_1 -> C_0 = zero (no 1-cells)
    cw.set_attaching_map(0, DMatrix::zeros(1, 0));
    // attaching_map[1] = d_2: C_2 -> C_1 = zero (no 1-cells, so boundary is 0)
    cw.set_attaching_map(1, DMatrix::zeros(0, 1));
    cw
}

/// Build S^n as a CW complex: one 0-cell, one n-cell.
pub fn cw_sn(n: usize) -> CWComplex {
    let mut cw = CWComplex::new(&format!("S^{}", n));
    cw.add_cell(0, "v");
    cw.add_cell(n, "f");
    // All attaching maps are zero
    for k in 0..n {
        let rows = cw.num_cells(k);
        let cols = cw.num_cells(k + 1);
        cw.set_attaching_map(k, DMatrix::zeros(rows, cols));
    }
    cw
}

/// Build CP^n (complex projective space) as a CW complex.
pub fn cw_cp2() -> CWComplex {
    let mut cw = CWComplex::new("CP^2");
    // CP^2 has cells in dimensions 0, 2, 4
    cw.add_cell(0, "v");
    cw.add_cell(2, "e");
    cw.add_cell(4, "f");
    // Attaching maps are all zero for CP^n
    cw.set_attaching_map(0, DMatrix::zeros(1, 0));
    cw.set_attaching_map(1, DMatrix::zeros(0, 1));
    cw.set_attaching_map(2, DMatrix::zeros(1, 0));
    cw.set_attaching_map(3, DMatrix::zeros(0, 1));
    cw
}

/// Build RP^2 as a CW complex.
pub fn cw_rp2() -> CWComplex {
    let mut cw = CWComplex::new("RP^2");
    // Cells: one in each dimension 0, 1, 2
    cw.add_cell(0, "v");
    cw.add_cell(1, "e");
    cw.add_cell(2, "f");
    // d_1: C_1 -> C_0: the 1-cell wraps around the 0-cell twice → 0 (identity map to single vertex)
    cw.set_attaching_map(0, DMatrix::from_row_slice(1, 1, &[0]));
    // d_2: C_2 -> C_1: the 2-cell wraps around the 1-cell twice → 2
    cw.set_attaching_map(1, DMatrix::from_row_slice(1, 1, &[2]));
    cw
}

/// Build T^2 (torus) as a CW complex.
pub fn cw_torus() -> CWComplex {
    let mut cw = CWComplex::new("T^2");
    // Cells: 1 vertex, 2 edges, 1 face
    cw.add_cell(0, "v");
    cw.add_cell(1, "a");
    cw.add_cell(1, "b");
    cw.add_cell(2, "f");
    // d_1: C_1 -> C_0: both edges map to the single vertex → 0
    cw.set_attaching_map(0, DMatrix::from_row_slice(1, 2, &[0, 0]));
    // d_2: C_2 -> C_1: the face attaches along aba^{-1}b^{-1} → 0 in cellular homology
    cw.set_attaching_map(1, DMatrix::from_row_slice(2, 1, &[0, 0]));
    cw
}

/// Build Klein bottle as a CW complex.
pub fn cw_klein_bottle() -> CWComplex {
    let mut cw = CWComplex::new("Klein bottle");
    // Cells: 1 vertex, 2 edges, 1 face
    cw.add_cell(0, "v");
    cw.add_cell(1, "a");
    cw.add_cell(1, "b");
    cw.add_cell(2, "f");
    // d_1: → 0 (single vertex)
    cw.set_attaching_map(0, DMatrix::from_row_slice(1, 2, &[0, 0]));
    // d_2: attaching map along aba^{-1}b → (0, 2) since b² means b maps to 2
    cw.set_attaching_map(1, DMatrix::from_row_slice(2, 1, &[0, 2]));
    cw
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cw_s1_homology() {
        let s1 = cw_s1();
        let h0 = s1.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = s1.homology(1);
        assert_eq!(h1.free_rank, 1);
    }

    #[test]
    fn test_cw_s2_homology() {
        let s2 = cw_s2();
        let h0 = s2.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h2 = s2.homology(2);
        assert_eq!(h2.free_rank, 1);
    }

    #[test]
    fn test_cw_s3_homology() {
        let s3 = cw_sn(3);
        let h0 = s3.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h3 = s3.homology(3);
        assert_eq!(h3.free_rank, 1);
    }

    #[test]
    fn test_cw_rp2_homology() {
        let rp2 = cw_rp2();
        let h0 = rp2.homology(0);
        assert_eq!(h0.free_rank, 1);
        // H_1(RP²) = Z/2
        let h1 = rp2.homology(1);
        assert_eq!(h1.torsion, vec![2]);
        // H_2(RP²) = 0
        let h2 = rp2.homology(2);
        assert!(h2.is_trivial());
    }

    #[test]
    fn test_cw_torus_homology() {
        let t2 = cw_torus();
        let h0 = t2.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = t2.homology(1);
        assert_eq!(h1.free_rank, 2);
        let h2 = t2.homology(2);
        assert_eq!(h2.free_rank, 1);
    }

    #[test]
    fn test_cw_klein_bottle_homology() {
        let kb = cw_klein_bottle();
        let h0 = kb.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h1 = kb.homology(1);
        assert_eq!(h1.free_rank, 1);
        assert_eq!(h1.torsion, vec![2]); // Z ⊕ Z/2
    }

    #[test]
    fn test_cw_cp2_homology() {
        let cp2 = cw_cp2();
        let h0 = cp2.homology(0);
        assert_eq!(h0.free_rank, 1);
        let h2 = cp2.homology(2);
        assert_eq!(h2.free_rank, 1);
        let h4 = cp2.homology(4);
        assert_eq!(h4.free_rank, 1);
    }

    #[test]
    fn test_cw_euler_characteristic() {
        let t2 = cw_torus();
        // 1 - 2 + 1 = 0
        assert_eq!(t2.euler_characteristic(), 0);

        let s2 = cw_s2();
        // 1 - 0 + 1 = 2
        assert_eq!(s2.euler_characteristic(), 2);
    }
}
