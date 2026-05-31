# lau-algebraic-topology

Fundamental algebraic topology in Rust: simplicial/singular homology, cohomology, Mayer-Vietoris sequence, Poincaré duality, homotopy groups, CW complexes, and agent network topology analysis.

## Features

- **Simplicial homology**: Chain complexes, boundary operators, H_n computation via Smith normal form
- **Singular homology**: Singular simplices, chain maps
- **Mayer-Vietoris sequence**: Excision-based homology computation
- **Homotopy groups**: π_n computation for n=1,2,3
- **CW complexes**: Cell attachment, cellular homology
- **Cohomology**: Universal coefficient theorem, cup product, Künneth formula
- **Poincaré duality**: H^k ≅ H_{n-k} for closed orientable manifolds
- **Euler characteristic**: From alternating sum of Betti numbers
- **Network topology**: Agent network analysis — connectedness, holes, voids

## Usage

```rust
use lau_algebraic_topology::simplicial;

// Compute homology of S^2
let s2 = simplicial::sphere_2();
let h0 = s2.homology(0); // Z
let h1 = s2.homology(1); // 0
let h2 = s2.homology(2); // Z

// CW complex homology
use lau_algebraic_topology::cw_complex;
let rp2 = cw_complex::cw_rp2();
let h1 = rp2.homology(1); // Z/2

// Network topology
use lau_algebraic_topology::network::AgentNetwork;
let mut net = AgentNetwork::new();
net.add_agent("A");
net.add_agent("B");
net.add_agent("C");
net.add_agent("D");
net.connect("A", "B");
net.connect("B", "C");
net.connect("C", "D");
net.connect("D", "A");
let analysis = net.analyze();
println!("{}", analysis); // Shows 1D hole in square
```

## Spaces with Known Homology

| Space | H_0 | H_1 | H_2 | H_3 | χ |
|-------|-----|-----|-----|-----|---|
| S^1 | Z | Z | - | - | 0 |
| S^2 | Z | 0 | Z | - | 2 |
| S^n | Z | 0 | ... | Z | (-1)^n·2 or 0 |
| T^2 | Z | Z² | Z | - | 0 |
| RP² | Z | Z/2 | 0 | - | 1 |
| Klein bottle | Z | Z⊕Z/2 | 0 | - | 0 |

## License

MIT
