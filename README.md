# lau-algebraic-topology

> **Algebraic topology in pure Rust** — simplicial complexes, homology, cohomology, Mayer-Vietoris, Poincaré duality, homotopy groups, CW complexes, and agent network topology analysis.

Part of the [PLATO/LAU](https://github.com/SuperInstance) ecosystem.

---

## What This Does

`lau-algebraic-topology` implements the fundamental tools of algebraic topology:

- **Simplicial complexes** — build, compute homology, Betti numbers, Euler characteristic
- **Chain complexes** — boundary operators, Smith normal form over ℤ, homology computation
- **Singular homology** — singular simplices, chains, chain maps
- **CW complexes** — cell attachment, cellular homology
- **Cohomology** — universal coefficient theorem, cup products, Künneth formula
- **Mayer-Vietoris sequence** — compute homology by decomposition
- **Poincaré duality** — verify duality for closed orientable manifolds
- **Homotopy groups** — π₁ and higher homotopy groups for standard spaces
- **Euler characteristic** — from Betti numbers or cell counts
- **Agent network analysis** — model agent networks as simplicial complexes, detect holes and voids

This is the topological backbone that lets the PLATO platform reason about connectivity, holes in agent interaction graphs, and conservation of topological invariants.

---

## The Key Idea

> **Topology sees shape through algebra.** Homology groups count holes; cohomology rings capture how they interact.

Instead of treating agent networks as flat graphs, `lau-algebraic-topology` lifts them to simplicial complexes and computes their homology:

1. **Agents** become vertices (0-simplices)
2. **Connections** become edges (1-simplices)
3. **Cliques** become higher simplices (triangles, tetrahedra, ...)
4. **Holes** are detected by homology groups (β₁ = number of loops, β₂ = number of voids)
5. **Euler characteristic** gives a single-number summary of the network's shape

This lets you detect when a network has structural gaps, redundant loops, or is simply connected — all algebraically.

---

## Install

```bash
cargo add lau-algebraic-topology
```

Or in `Cargo.toml`:

```toml
[dependencies]
lau-algebraic-topology = "0.1"
```

Requires **Rust 2021 edition**. Dependencies: `serde`, `nalgebra`.

---

## Quick Start

### Build a Simplicial Complex

```rust
use lau_algebraic_topology::simplicial::*;

let mut k = SimplicialComplex::new("my complex");
k.add_simplex(&[0, 1, 2]); // triangle with all faces
k.add_simplex(&[1, 2, 3]); // adjacent triangle

// How many simplices?
println!("Vertices: {}", k.num_vertices());  // 4
println!("Edges: {}", k.num_edges());        // 5
println!("Triangles: {}", k.num_simplices(2)); // 2
println!("Dimension: {}", k.dimension());     // 2
```

### Compute Homology

```rust
// Pre-built circle S¹
let s1 = circle();
let h0 = s1.homology(0); // H₀ = Z (connected)
let h1 = s1.homology(1); // H₁ = Z (one loop)
let h2 = s1.homology(2); // H₂ = 0

assert_eq!(h0.free_rank, 1);
assert_eq!(h1.free_rank, 1);
assert!(h2.is_trivial());

// Betti numbers: [1, 1]
let betti = s1.betti_numbers();
println!("Betti numbers: {:?}", betti);
```

### Analyze Agent Networks

```rust
use lau_algebraic_topology::network::*;

let mut net = AgentNetwork::new();
net.add_agent("Alice");
net.add_agent("Bob");
net.add_agent("Carol");
net.add_agent("Dave");

net.connect("Alice", "Bob");
net.connect("Bob", "Carol");
net.connect("Carol", "Dave");
net.connect("Dave", "Alice");
// Square: one 1D hole

let analysis = net.analyze();
println!("{}", analysis);
// Agents: 4, Connections: 4
// Connected: yes
// Betti numbers: [1, 1, 0, 0]
// 1D holes: 1 (the square's interior)
// Euler characteristic: 0
```

### Well-Known Spaces

```rust
use lau_algebraic_topology::simplicial::*;

let s2 = sphere_2();          // 2-sphere
let t2 = torus();             // Torus T²
let rp2 = rp2();              // Real projective plane
let kb = klein_bottle();      // Klein bottle
let mob = mobius_strip();     // Möbius strip
let s3 = sphere_n(3);         // 3-sphere
let pt = point();             // Single point
let d2 = disk();              // 2-disk (contractible)
```

### CW Complexes

```rust
use lau_algebraic_topology::cw_complex::*;

// S¹: one 0-cell, one 1-cell
let s1_cw = cw_s1();
assert_eq!(s1_cw.homology(1).free_rank, 1);

// RP²: one cell per dimension 0,1,2
let rp2_cw = cw_rp2();
let h1 = rp2_cw.homology(1);
assert_eq!(h1.torsion, vec![2]); // H₁(RP²) = Z/2

// Torus CW: 1 vertex, 2 edges, 1 face
let t2_cw = cw_torus();
assert_eq!(t2_cw.homology(1).free_rank, 2); // H₁(T²) = Z²
```

### Euler Characteristic

```rust
use lau_algebraic_topology::euler::*;

let info = EulerInfo::sphere(2);
assert_eq!(info.euler, 2);  // χ(S²) = 2

let info = EulerInfo::torus(2);
assert_eq!(info.euler, 0);  // χ(T²) = 0

let info = EulerInfo::cpn(2);
assert_eq!(info.euler, 3);  // χ(CP²) = 3

// Euler from Betti numbers: χ = β₀ - β₁ + β₂ - ...
assert_eq!(euler_from_betti(&[1, 0, 1]), 2); // S²
assert_eq!(euler_from_betti(&[1, 2, 1]), 0); // T²
```

### Homotopy Groups

```rust
use lau_algebraic_topology::homotopy::*;

// π₁(S¹) = Z
let g = pi1_s1();
assert_eq!(g.description, "Z (free group on 1 generator)");

// π₁(RP²) = Z/2
let g = pi1_space("RP^2");
assert!(g.description.contains("Z/2"));

// π₁(T²) = Z × Z
let g = pi1_space("T^2");
assert!(g.description.contains("Z × Z"));

// π₃(S²) = Z (Hopf fibration!)
let g = pik_sn(3, 2);
assert_eq!(g.description, "Z");

// Fundamental group presentations
let fg_torus = FundamentalGroup::of_torus();
// 2 generators, relation aba⁻¹b⁻¹, abelian
assert_eq!(fg_torus.num_generators, 2);
assert!(fg_torus.is_abelian);
```

---

## API Reference

### `simplicial` Module

#### `SimplicialComplex`
A simplicial complex: a set of simplices closed under taking faces.

| Method | Description |
|---|---|
| `new(name)` | Create empty complex |
| `add_simplex(&[v0, v1, ...])` | Add simplex and all faces |
| `num_simplices(dim)` | Count simplices of given dimension |
| `simplices_of_dim(dim)` | Get sorted list of simplices |
| `dimension()` | Highest dimension |
| `num_vertices()` / `num_edges()` | Convenience counts |
| `chain_complex()` | Build the chain complex |
| `homology(n)` | Compute Hₙ |
| `all_homology()` | All homology groups up to dimension |
| `betti_numbers()` | Betti numbers [β₀, β₁, ...] |
| `euler_characteristic()` | χ = Σ(-1)ⁿ · |simplicesₙ| |

**Pre-built spaces:** `circle()`, `sphere_2()`, `sphere_n(n)`, `torus()`, `rp2()`, `klein_bottle()`, `mobius_strip()`, `point()`, `interval()`, `disk()`.

### `chain_complex` Module

#### `ChainComplex`
A chain complex over ℤ with boundary matrices.

| Method | Description |
|---|---|
| `from_ranks_and_maps(ranks, maps)` | Construct from ranks and boundary matrices |
| `trivial()` | Zero complex |
| `smith_normal_form(&mat)` | SNF: returns (U, S, V) where S is diagonal |
| `kernel_dimension(&mat)` | dim ker via SNF |
| `image_rank(&mat)` | rank im via SNF |
| `homology(n)` | Hₙ = ker(dₙ) / im(dₙ₊₁) |
| `all_homology(max_n)` | All groups up to max_n |
| `betti_numbers(max_n)` | Free ranks |
| `is_valid()` | Check dₖ ∘ dₖ₊₁ = 0 |

#### `HomologyGroup`
Hₙ = ℤ^r ⊕ ℤ/t₁ ⊕ ℤ/t₂ ⊕ ...

| Field/Method | Description |
|---|---|
| `dimension` | n |
| `free_rank` | Rank of free part |
| `torsion` | Torsion coefficients |
| `is_trivial()` | free_rank = 0 and no torsion |
| `rank()` | max(free_rank, 0) as usize |

### `cw_complex` Module

#### `CWComplex`
A CW complex with explicit attaching maps.

| Method | Description |
|---|---|
| `new(name)` | Empty complex |
| `add_cell(dim, label)` | Add a cell |
| `set_attaching_map(k, matrix)` | Set boundary map for (k+1)-cells → k-cells |
| `num_cells(k)` | Count k-cells |
| `chain_complex()` | Cellular chain complex |
| `homology(n)` | Cellular homology |
| `euler_characteristic()` | χ = Σ(-1)ⁿ · cₙ |

**Pre-built:** `cw_s1()`, `cw_s2()`, `cw_sn(n)`, `cw_rp2()`, `cw_torus()`, `cw_klein_bottle()`, `cw_cp2()`.

### `cohomology` Module

| Function | Description |
|---|---|
| `cohomology_from_homology(&hₙ, &hₙ₋₁)` | Universal coefficient theorem: Hⁿ from Hₙ, Hₙ₋₁ |
| `all_cohomology(&homology)` | All cohomology groups |
| `kunneth(&hom_x, &hom_y, n)` | Künneth: Hₙ(X × Y) |
| `kunneth_all(&hom_x, &hom_y, max_n)` | All Künneth groups |
| `cohomology_z2(&complex, n)` | Hⁿ with ℤ/2 coefficients |

#### `CupProduct`
Cup product ⌣: Hᵖ × Hᵠ → Hᵖ⁺ᵠ with product matrix.

#### `CohomologyRing`
Cohomology ring with cup product structure.

### `singular` Module

#### `SingularSimplex`
A singular n-simplex σ: Δⁿ → X, represented by label and vertices.

| Method | Description |
|---|---|
| `new(dim, label, vertices)` | Create |
| `face(i)` | i-th face map (remove vertex i) |

#### `SingularChain`
Formal ℤ-linear combination of singular simplices.

| Method | Description |
|---|---|
| `zero(dim, n)` | Zero chain with n generators |
| `add(&other)` | Add chains |
| `scale(n)` | Multiply by integer |

#### `ChainMap`
A chain map between chain complexes, inducing maps on homology.

### `mayer_vietoris` Module

| Function | Description |
|---|---|
| `mayer_vietoris_homology(&hA, &hB, &hA∩B, max_dim)` | H(X) from H(A), H(B), H(A∩B) via exact sequence |
| `sphere_homology_mv(n)` | H(Sⁿ) via recursive Mayer-Vietoris |
| `wedge_s1_s1_homology()` | H(S¹ ∨ S¹) example |

### `poincare` Module

| Function | Description |
|---|---|
| `check_poincare_duality(&homology, n)` | Verify Hᵏ ≅ Hₙ₋ₖ for all k |
| `verify_poincare_s2()` / `_s3()` / `_sn(n)` / `_torus()` | Pre-built verifications |
| `intersection_form(betti_2k)` | Intersection form of a 4k-manifold |

### `homotopy` Module

| Function | Description |
|---|---|
| `pi1_s1()` | π₁(S¹) = ℤ |
| `pin_s1(n)` | πₙ(S¹) = 0 for n ≥ 2 |
| `pin_sn(n)` | πₙ(Sⁿ) = ℤ |
| `pik_sn(k, n)` | πₖ(Sⁿ) for known cases |
| `pi1_space(name)` | π₁ of named spaces |

#### `FundamentalGroup`
Presentation of π₁: generators, relations, abelianness.

| Constructor | Description |
|---|---|
| `of_s1()` | ℤ (1 generator, no relations) |
| `of_wedge_circles(k)` | Fₖ (free group on k generators) |
| `of_torus()` | ℤ × ℤ (2 generators, relation aba⁻¹b⁻¹) |
| `of_rp2()` | ℤ/2 (1 generator, relation a²) |

### `euler` Module

| Function/Type | Description |
|---|---|
| `euler_from_betti(&betti)` | χ = Σ(-1)ⁿβₙ |
| `euler_from_homology_ranks(&ranks)` | Same, from ranks |
| `EulerInfo::sphere(n)` | χ(Sⁿ) = 1 + (-1)ⁿ |
| `EulerInfo::torus(n)` | χ(Tⁿ) = 0 |
| `EulerInfo::rp2()` | χ(RP²) = 1 |
| `EulerInfo::klein_bottle()` | χ(K) = 0 |
| `EulerInfo::cpn(n)` | χ(CPⁿ) = n + 1 |
| `EulerInfo::point()` | χ(pt) = 1 |

### `network` Module

#### `AgentNetwork`
Model agent interactions as a simplicial complex.

| Method | Description |
|---|---|
| `new()` | Empty network |
| `add_agent(id)` | Add vertex |
| `connect(from, to)` | Add edge |
| `adjacency_matrix()` | n×n adjacency matrix |
| `connected_components()` | Union-find components |
| `num_components()` | Count |
| `network_homology(max_dim)` | Homology via clique complex |
| `analyze()` | Full `NetworkAnalysis` |

#### `NetworkAnalysis`
| Field | Description |
|---|---|
| `num_agents` / `num_connections` | Counts |
| `num_components` | Connected components |
| `is_connected` | Single component? |
| `betti_numbers` | [β₀, β₁, β₂, β₃] |
| `euler_characteristic` | χ |
| `holes_1d` | Loops (β₁) |
| `voids_2d` / `voids_3d` | Higher voids (β₂, β₃) |
| `components` | Agent IDs per component |

---

## How It Works

### Homology via Smith Normal Form

```
Chain complex: ... → Cₙ₊₁ --dₙ₊₁--> Cₙ --dₙ--> Cₙ₋₁ → ...
                              ↑              ↑
                          boundary map    boundary map

Hₙ = ker(dₙ) / im(dₙ₊₁)

Computation:
1. Build boundary matrices over ℤ
2. Compute Smith normal form: S = U · M · V
3. rank = number of nonzero diagonal entries of S
4. ker_dim = ncols - rank
5. free_rank = ker(dₙ) - im(dₙ₊₁)
6. torsion = diagonal entries > 1 of dₙ₊₁'s SNF
```

### Simplicial Boundary

The boundary of an n-simplex [v₀, v₁, ..., vₙ]:

```
∂[v₀, ..., vₙ] = Σᵢ₌₀ⁿ (-1)ⁱ [v₀, ..., v̂ᵢ, ..., vₙ]
```

where v̂ᵢ means "omit vertex i". This satisfies ∂² = 0.

### Network → Simplicial Complex (Clique Complex)

```
For a graph G = (V, E):
  0-simplices: all vertices
  1-simplices: all edges
  k-simplices: all (k+1)-cliques (complete subgraphs)

Then compute homology of the resulting complex.
```

A square (4 vertices, 4 edges, no diagonal) has β₁ = 1 (one hole).
Adding a diagonal fills the hole: β₁ = 0.

### Mayer-Vietoris

If X = A ∪ B, the long exact sequence:

```
... → Hₙ(A∩B) → Hₙ(A) ⊕ Hₙ(B) → Hₙ(X) → Hₙ₋₁(A∩B) → ...
```

From exactness, the rank formula:

```
rank(Hₙ(X)) = rank(Hₙ(A)) + rank(Hₙ(B)) - rank(Hₙ(A∩B))
             + rank(Hₙ₋₁(A∩B)) - rank(Hₙ₋₁(A)) - rank(Hₙ₋₁(B))
```

### Poincaré Duality

For a closed orientable n-manifold M:

```
Hᵏ(M; ℤ) ≅ Hₙ₋ₖ(M; ℤ)   for all k
```

Verified by computing cohomology via the universal coefficient theorem and checking rank equality.

---

## The Math

### Homology Groups

The n-th homology group of a space X:

$$H_n(X; \mathbb{Z}) = \ker(d_n) / \operatorname{im}(d_{n+1})$$

By the structure theorem for finitely generated abelian groups:

$$H_n \cong \mathbb{Z}^{\beta_n} \oplus \mathbb{Z}/t_1 \oplus \cdots \oplus \mathbb{Z}/t_k$$

where βₙ is the n-th **Betti number** (number of n-dimensional holes).

### Betti Numbers

| β₀ | Connected components |
| β₁ | Independent loops (1D holes) |
| β₂ | Voids (2D cavities) |
| βₙ | n-dimensional "holes" |

### Euler Characteristic

$$\chi(X) = \sum_{n=0}^{\infty} (-1)^n \beta_n = \sum_{n=0}^{\infty} (-1)^n c_n$$

where cₙ is the number of n-cells (simplices). This is a **topological invariant** — independent of the particular triangulation or CW decomposition.

### Well-Known Homology

| Space | H₀ | H₁ | H₂ | H₃ | χ |
|---|---|---|---|---|---|
| S⁰ | ℤ² | — | — | — | 2 |
| S¹ | ℤ | ℤ | — | — | 0 |
| S² | ℤ | 0 | ℤ | — | 2 |
| Sⁿ | ℤ | 0 | ... | 0, ℤ | 1+(−1)ⁿ |
| T² | ℤ | ℤ² | ℤ | — | 0 |
| RP² | ℤ | ℤ/2 | 0 | — | 1 |
| Klein bottle | ℤ | ℤ ⊕ ℤ/2 | 0 | — | 0 |
| CP² | ℤ | 0 | ℤ | 0 | 3 |
| pt | ℤ | — | — | — | 1 |

### Universal Coefficient Theorem for Cohomology

$$H^n(X; \mathbb{Z}) \cong \operatorname{Hom}(H_n(X), \mathbb{Z}) \oplus \operatorname{Ext}(H_{n-1}(X), \mathbb{Z})$$

Free parts of homology become free cohomology; torsion in Hₙ₋₁ becomes torsion in Hⁿ.

### Künneth Formula

$$H_n(X \times Y) \cong \bigoplus_{p+q=n} H_p(X) \otimes H_q(Y) \oplus \bigoplus_{p+q=n-1} \operatorname{Tor}(H_p(X), H_q(Y))$$

Example: H(T²) = H(S¹ × S¹) via Künneth with H(S¹) gives ℤ, ℤ², ℤ.

### Homotopy Groups of Spheres

| | n=1 | n=2 | n=3 | n=4 |
|---|---|---|---|---|
| **k=1** | ℤ | 0 | 0 | 0 |
| **k=2** | 0 | ℤ | 0 | 0 |
| **k=3** | 0 | ℤ | ℤ | 0 |
| **k=4** | 0 | ℤ/2 | ℤ/2 | ℤ |

---

## Testing

**76 tests** across 10 modules:

| Module | Tests | Coverage |
|---|---|---|
| `chain_complex` | 5 | SNF, kernel, image, trivial complex |
| `simplicial` | 9 | All pre-built spaces, homology, Euler |
| `singular` | 4 | Singular simplices, chains, face maps |
| `cw_complex` | 8 | CW S¹/S²/S³/RP²/T²/Klein/CP², Euler |
| `cohomology` | 7 | UCT, Künneth, cup product, ℤ/2 |
| `mayer_vietoris` | 5 | Spheres, wedge, union |
| `poincare` | 7 | S²/S³/S⁴/S⁵/T², failure case |
| `euler` | 11 | All spaces, from Betti, from ranks |
| `homotopy` | 11 | π₁, πₙ(Sⁿ), πₖ(Sⁿ), fundamental groups |
| `network` | 9 | Connected/disconnected, holes, adjacency |

Run: `cargo test`

---

## Caveats

- **Integer Smith normal form:** Exact computation over ℤ via iterative reduction. May be slow for large boundary matrices (O(n³) worst case).
- **Torsion detection:** Torsion is extracted from diagonal entries of the SNF. For complexes with many generators, this is accurate but potentially slow.
- **Simplicial homology only:** No persistent homology or Morse theory (yet).
- **Homotopy groups are lookup tables:** πₙ(Sᵐ) for arbitrary (n,m) is an open problem in mathematics. Only known cases are provided.
- **Network homology uses clique complex:** Higher simplices are all complete subgraphs. For large dense networks, this can be exponential.

---

## License

MIT
