//! Homotopy groups π_n for n=1,2,3.

use serde::{Deserialize, Serialize};

/// A homotopy group result.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HomotopyGroup {
    /// Space name.
    pub space: String,
    /// Dimension n.
    pub dimension: usize,
    /// Whether the group is abelian (always true for n >= 2).
    pub is_abelian: bool,
    /// Description of the group.
    pub description: String,
}

/// Compute π_1(S^1) = Z.
pub fn pi1_s1() -> HomotopyGroup {
    HomotopyGroup {
        space: "S^1".into(),
        dimension: 1,
        is_abelian: false,
        description: "Z (free group on 1 generator)".into(),
    }
}

/// Compute π_n(S^1) = 0 for n >= 2.
pub fn pin_s1(n: usize) -> HomotopyGroup {
    HomotopyGroup {
        space: "S^1".into(),
        dimension: n,
        is_abelian: true,
        description: "0 (trivial)".into(),
    }
}

/// Compute π_n(S^n) = Z for all n >= 1.
pub fn pin_sn(n: usize) -> HomotopyGroup {
    HomotopyGroup {
        space: format!("S^{}", n),
        dimension: n,
        is_abelian: n >= 2,
        description: "Z".into(),
    }
}

/// π_k(S^n) for known cases.
pub fn pik_sn(k: usize, n: usize) -> HomotopyGroup {
    if k == n {
        return pin_sn(n);
    }
    if k < n {
        return HomotopyGroup {
            space: format!("S^{}", n),
            dimension: k,
            is_abelian: k >= 2,
            description: "0 (trivial)".into(),
        };
    }
    if n == 1 && k >= 2 {
        return pin_s1(k);
    }
    // Some known higher homotopy groups
    match (k, n) {
        (2, 2) => HomotopyGroup {
            space: "S^2".into(), dimension: 2, is_abelian: true,
            description: "Z".into(),
        },
        (3, 2) => HomotopyGroup {
            space: "S^2".into(), dimension: 3, is_abelian: true,
            description: "Z".into(),
        },
        (4, 2) => HomotopyGroup {
            space: "S^2".into(), dimension: 4, is_abelian: true,
            description: "Z/2".into(),
        },
        (3, 3) => HomotopyGroup {
            space: "S^3".into(), dimension: 3, is_abelian: true,
            description: "Z".into(),
        },
        (4, 3) => HomotopyGroup {
            space: "S^3".into(), dimension: 4, is_abelian: true,
            description: "Z/2".into(),
        },
        _ => HomotopyGroup {
            space: format!("S^{}", n), dimension: k, is_abelian: true,
            description: "Unknown / complex".into(),
        },
    }
}

/// π_1 of common spaces.
pub fn pi1_space(space: &str) -> HomotopyGroup {
    match space {
        "S^1" => pi1_s1(),
        "S^n" | "S^2" | "S^3" => HomotopyGroup {
            space: space.into(), dimension: 1, is_abelian: true,
            description: "0 (simply connected)".into(),
        },
        "T^2" | "torus" => HomotopyGroup {
            space: "T^2".into(), dimension: 1, is_abelian: false,
            description: "Z × Z (free abelian on 2 generators)".into(),
        },
        "RP^2" | "RP²" => HomotopyGroup {
            space: "RP²".into(), dimension: 1, is_abelian: false,
            description: "Z/2".into(),
        },
        "Klein" | "Klein bottle" => HomotopyGroup {
            space: "Klein bottle".into(), dimension: 1, is_abelian: false,
            description: "Z ⋊ Z (semidirect product)".into(),
        },
        "pt" | "point" => HomotopyGroup {
            space: "pt".into(), dimension: 1, is_abelian: true,
            description: "0 (trivial)".into(),
        },
        _ => HomotopyGroup {
            space: space.into(), dimension: 1, is_abelian: true,
            description: "Unknown".into(),
        },
    }
}

/// Compute the fundamental group π_1 from the 1-skeleton and 2-cells of a CW complex.
/// Returns a description of the group via generators and relations.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FundamentalGroup {
    /// Number of generators (from 1-cells).
    pub num_generators: usize,
    /// Relations from 2-cells.
    pub relations: Vec<String>,
    /// Whether the group is abelian.
    pub is_abelian: bool,
    /// Description.
    pub description: String,
}

impl FundamentalGroup {
    /// Compute π_1 of S^1: one generator, no relations.
    pub fn of_s1() -> Self {
        FundamentalGroup {
            num_generators: 1,
            relations: vec![],
            is_abelian: false,
            description: "Z (free cyclic)".into(),
        }
    }

    /// Compute π_1 of a wedge of k circles: k generators, no relations.
    pub fn of_wedge_circles(k: usize) -> Self {
        FundamentalGroup {
            num_generators: k,
            relations: vec![],
            is_abelian: false,
            description: format!("Free group F_{}", k),
        }
    }

    /// Compute π_1 of T^2: two generators, one relation aba^{-1}b^{-1}.
    pub fn of_torus() -> Self {
        FundamentalGroup {
            num_generators: 2,
            relations: vec!["aba^{-1}b^{-1}".into()],
            is_abelian: true,
            description: "Z × Z".into(),
        }
    }

    /// Compute π_1 of RP²: one generator, one relation a².
    pub fn of_rp2() -> Self {
        FundamentalGroup {
            num_generators: 1,
            relations: vec!["a^2".into()],
            is_abelian: false,
            description: "Z/2".into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pi1_s1() {
        let g = pi1_s1();
        assert_eq!(g.description, "Z (free group on 1 generator)");
    }

    #[test]
    fn test_pin_s1_trivial() {
        let g = pin_s1(2);
        assert_eq!(g.description, "0 (trivial)");
    }

    #[test]
    fn test_pi1_sn() {
        let g = pi1_space("S^2");
        assert!(g.description.contains("0"));
    }

    #[test]
    fn test_pi1_torus() {
        let g = pi1_space("T^2");
        assert!(g.description.contains("Z × Z"));
    }

    #[test]
    fn test_pi1_rp2() {
        let g = pi1_space("RP^2");
        assert!(g.description.contains("Z/2"));
    }

    #[test]
    fn test_fundamental_group_s1() {
        let fg = FundamentalGroup::of_s1();
        assert_eq!(fg.num_generators, 1);
        assert!(fg.relations.is_empty());
    }

    #[test]
    fn test_fundamental_group_torus() {
        let fg = FundamentalGroup::of_torus();
        assert_eq!(fg.num_generators, 2);
        assert!(fg.is_abelian);
    }

    #[test]
    fn test_fundamental_group_rp2() {
        let fg = FundamentalGroup::of_rp2();
        assert_eq!(fg.num_generators, 1);
        assert!(!fg.is_abelian);
    }

    #[test]
    fn test_pik_sn_identity() {
        let g = pik_sn(2, 2);
        assert_eq!(g.description, "Z");
    }

    #[test]
    fn test_pik_sn_below() {
        let g = pik_sn(1, 2);
        assert!(g.description.contains("0"));
    }

    #[test]
    fn test_pi1_klein() {
        let g = pi1_space("Klein bottle");
        assert!(g.description.contains("Z"));
    }
}
