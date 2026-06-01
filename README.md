# lau-algebraic-topology

> Fundamental algebraic topology: simplicial/singular homology, cohomology, Mayer-Vietoris, Poincaré duality, homotopy groups, CW complexes, and agent network topology analysis

## What This Does

Fundamental algebraic topology: simplicial/singular homology, cohomology, Mayer-Vietoris, Poincaré duality, homotopy groups, CW complexes, and agent network topology analysis. Part of the PLATO/LAU ecosystem — a mathematically rigorous framework for building educational agents that learn, teach, and evolve.

## The Key Idea

This crate implements the core abstractions needed for its domain, with a focus on correctness, composability, and conservation guarantees. Every public type is serializable (serde), every algorithm is tested, and every invariant is verified.

## Install

```bash
cargo add lau-algebraic-topology
```

## Quick Start

See the API Reference below for complete usage. Key entry points:

```rust
use lau_algebraic_topology::*;
// See types and methods below for complete usage
```

## API Reference

```rust
pub struct Cell 
pub struct CWComplex 
    pub fn new(name: &str) -> Self 
    pub fn add_cell(&mut self, dim: usize, label: &str) 
    pub fn set_attaching_map(&mut self, k: usize, map: DMatrix<i64>) 
    pub fn num_cells(&self, k: usize) -> usize 
    pub fn dimension(&self) -> usize 
    pub fn chain_complex(&self) -> ChainComplex 
    pub fn homology(&self, n: usize) -> HomologyGroup 
    pub fn euler_characteristic(&self) -> i64 
pub fn cw_s1() -> CWComplex 
pub fn cw_s2() -> CWComplex 
pub fn cw_sn(n: usize) -> CWComplex 
pub fn cw_cp2() -> CWComplex 
pub fn cw_rp2() -> CWComplex 
pub fn cw_torus() -> CWComplex 
pub fn cw_klein_bottle() -> CWComplex 
pub fn euler_from_betti(betti: &[usize]) -> i64 
pub fn euler_from_homology_ranks(ranks: &[i64]) -> i64 
pub struct EulerInfo 
    pub fn sphere(n: usize) -> Self 
    pub fn torus(n: usize) -> Self 
    pub fn rp2() -> Self 
    pub fn klein_bottle() -> Self 
    pub fn point() -> Self 
    pub fn cpn(n: usize) -> Self 
pub struct ChainComplex 
    pub fn from_ranks_and_maps(ranks: Vec<usize>, boundary_maps: Vec<ZMatrix>) -> Self 
    pub fn trivial() -> Self 
    pub fn smith_normal_form(mat: &ZMatrix) -> (ZMatrix, ZMatrix, ZMatrix) 
    pub fn kernel_dimension(mat: &ZMatrix) -> usize 
    pub fn image_rank(mat: &ZMatrix) -> usize 
    pub fn homology(&self, n: usize) -> HomologyGroup 
    pub fn all_homology(&self, max_n: usize) -> Vec<HomologyGroup> 
    pub fn betti_numbers(&self, max_n: usize) -> Vec<usize> 
    pub fn is_valid(&self) -> bool 
pub struct HomologyGroup 
    pub fn is_trivial(&self) -> bool 
    pub fn rank(&self) -> usize 
pub struct SingularSimplex 
    pub fn new(dimension: usize, label: &str, vertices: Vec<usize>) -> Self 
    pub fn face(&self, i: usize) -> SingularSimplex 
pub struct SingularChain 
    pub fn zero(dimension: usize, n: usize) -> Self 
    pub fn add(&self, other: &SingularChain) -> SingularChain 
    pub fn scale(&self, n: i64) -> SingularChain 
pub struct SingularChainComplex 
    pub fn from_chain_complex(complex: ChainComplex, labels: Vec<Vec<String>>) -> Self 
    pub fn homology(&self, n: usize) -> HomologyGroup 
    pub fn all_homology(&self, max_n: usize) -> Vec<HomologyGroup> 
pub struct ChainMap 
    pub fn new(maps: Vec<DMatrix<i64>>) -> Self 
    pub fn induced_homology_matrix(&self, _n: usize) -> Option<DMatrix<i64>> 
pub fn singular_homology_of_space(complex: &ChainComplex, n: usize) -> HomologyGroup 
pub struct Agent 
pub struct Connection 
pub struct AgentNetwork 
    pub fn new() -> Self 
    pub fn add_agent(&mut self, id: &str) 
    pub fn connect(&mut self, from: &str, to: &str) 
```

## How It Works

Read the source in `src/` for full implementation details. All algorithms are documented with inline comments explaining the mathematical foundations.

## The Math

This crate implements formal mathematical constructs. See the source documentation for theorem statements and proofs of correctness.

## Testing

**76 tests** covering construction, serialization, correctness properties, edge cases, and composability with other lau-* crates.

## License

MIT
