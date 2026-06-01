//! Lebesgue decomposition: ν = ν_ac + ν_singular
//!
//! Any σ-finite measure ν can be uniquely decomposed as:
//! - ν_ac: absolutely continuous part (ν_ac ≪ μ)
//! - ν_s: singular part (ν_s ⊥ μ)
//! where ν_ac and ν_s are mutually singular.

use serde::{Serialize, Deserialize};
use crate::sigma_algebra::Subset;
use crate::measure::Measure;

/// Lebesgue decomposition of a measure relative to a reference.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LebesgueDecomposition {
    /// Absolutely continuous part: ν_ac ≪ μ.
    pub ac_part: Measure,
    /// Singular part: ν_s ⊥ μ.
    pub singular_part: Measure,
    /// The original measure ν.
    pub original_name: String,
    /// The reference measure μ.
    pub reference_name: String,
}

impl LebesgueDecomposition {
    /// Compute the Lebesgue decomposition of ν with respect to μ.
    ///
    /// For finite spaces:
    /// - ac_part: contributions where μ > 0 (ν_ac({x}) = ν({x}) if μ({x}) > 0)
    /// - singular_part: contributions where μ = 0
    pub fn compute(
        nu: &Measure,
        mu: &Measure,
        elements: &[String],
        sa: &crate::sigma_algebra::SigmaAlgebra,
    ) -> Self {
        let mut ac_values = Vec::new();
        let mut sing_values = Vec::new();

        for elem in elements {
            let singleton = Subset::from_slice(&[elem]);
            let nu_val = nu.measure_of(&singleton);
            let mu_val = mu.measure_of(&singleton);

            if mu_val > 0.0 {
                // Absolutely continuous part
                ac_values.push((singleton.clone(), nu_val));
            } else if nu_val > 0.0 {
                // Singular part: ν is positive where μ is zero
                sing_values.push((singleton.clone(), nu_val));
            }
        }

        let ac_part = Measure::new(&format!("{}_ac", nu.name()), sa, &ac_values).unwrap_or(
            Measure::new(&format!("{}_ac", nu.name()), sa, &[]).unwrap()
        );
        let singular_part = Measure::new(&format!("{}_sing", nu.name()), sa, &sing_values).unwrap_or(
            Measure::new(&format!("{}_sing", nu.name()), sa, &[]).unwrap()
        );

        LebesgueDecomposition {
            ac_part,
            singular_part,
            original_name: nu.name().to_string(),
            reference_name: mu.name().to_string(),
        }
    }

    /// Verify the decomposition: ν = ν_ac + ν_s.
    pub fn verify(&self, nu: &Measure, elements: &[String]) -> bool {
        for elem in elements {
            let singleton = Subset::from_slice(&[elem]);
            let nu_val = nu.measure_of(&singleton);
            let ac_val = self.ac_part.measure_of(&singleton);
            let sing_val = self.singular_part.measure_of(&singleton);
            if (nu_val - ac_val - sing_val).abs() > 1e-10 {
                return false;
            }
        }
        true
    }

    /// Is the decomposition pure (entirely ac or entirely singular)?
    pub fn is_pure(&self) -> PureType {
        let ac_has_mass = self.ac_part.total_mass() > 1e-10;
        let sing_has_mass = self.singular_part.total_mass() > 1e-10;
        match (ac_has_mass, sing_has_mass) {
            (true, false) => PureType::AbsolutelyContinuous,
            (false, true) => PureType::Singular,
            (true, true) => PureType::Mixed,
            (false, false) => PureType::Zero,
        }
    }
}

/// Classification of decomposition purity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PureType {
    /// Entirely absolutely continuous w.r.t. μ.
    AbsolutelyContinuous,
    /// Entirely singular w.r.t. μ.
    Singular,
    /// Both components present.
    Mixed,
    /// Zero measure.
    Zero,
}

/// Hahn decomposition: for a signed measure μ, partition Ω = P ∪ N
/// where μ is non-negative on P and non-positive on N.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HahnDecomposition {
    /// Positive set.
    pub positive: Subset,
    /// Negative set.
    pub negative: Subset,
}

impl HahnDecomposition {
    /// Compute the Hahn decomposition for a signed measure on a finite space.
    pub fn compute(
        values: &[(String, f64)],
    ) -> Self {
        let mut pos = std::collections::BTreeSet::new();
        let mut neg = std::collections::BTreeSet::new();

        for (elem, val) in values {
            if *val >= 0.0 {
                pos.insert(elem.clone());
            } else {
                neg.insert(elem.clone());
            }
        }

        HahnDecomposition {
            positive: Subset(pos),
            negative: Subset(neg),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sigma_algebra::SigmaAlgebra;

    fn setup_decomp() -> (SigmaAlgebra, Measure, Measure) {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c", "d"]));
        // μ is positive on {a, b}, zero on {c, d}
        let mu = Measure::new("mu", &sa, &[
            (Subset::from_slice(&["a"]), 0.3),
            (Subset::from_slice(&["b"]), 0.2),
            (Subset::from_slice(&["c"]), 0.0),
            (Subset::from_slice(&["d"]), 0.0),
        ]).unwrap();
        // ν has mass everywhere
        let nu = Measure::new("nu", &sa, &[
            (Subset::from_slice(&["a"]), 0.4),
            (Subset::from_slice(&["b"]), 0.3),
            (Subset::from_slice(&["c"]), 0.2),
            (Subset::from_slice(&["d"]), 0.1),
        ]).unwrap();
        (sa, mu, nu)
    }

    #[test]
    fn test_lebesgue_decomposition() {
        let (sa, mu, nu) = setup_decomp();
        let elems = vec!["a".into(), "b".into(), "c".into(), "d".into()];
        let decomp = LebesgueDecomposition::compute(&nu, &mu, &elems, &sa);
        
        // ac part: {a} → 0.4, {b} → 0.3 (where μ > 0)
        assert!((decomp.ac_part.measure_of(&Subset::from_slice(&["a"])) - 0.4).abs() < 1e-10);
        assert!((decomp.ac_part.measure_of(&Subset::from_slice(&["b"])) - 0.3).abs() < 1e-10);
        
        // singular part: {c} → 0.2, {d} → 0.1 (where μ = 0)
        assert!((decomp.singular_part.measure_of(&Subset::from_slice(&["c"])) - 0.2).abs() < 1e-10);
        assert!((decomp.singular_part.measure_of(&Subset::from_slice(&["d"])) - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_decomposition_verify() {
        let (sa, mu, nu) = setup_decomp();
        let elems = vec!["a".into(), "b".into(), "c".into(), "d".into()];
        let decomp = LebesgueDecomposition::compute(&nu, &mu, &elems, &sa);
        assert!(decomp.verify(&nu, &elems));
    }

    #[test]
    fn test_pure_ac() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let mu = Measure::uniform(&sa).unwrap();
        let nu = Measure::probability(&sa, &[
            ("a".into(), 0.6), ("b".into(), 0.4),
        ]).unwrap();
        let elems = vec!["a".into(), "b".into()];
        let decomp = LebesgueDecomposition::compute(&nu, &mu, &elems, &sa);
        assert_eq!(decomp.is_pure(), PureType::AbsolutelyContinuous);
    }

    #[test]
    fn test_pure_singular() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let mu = Measure::dirac(&sa, "a").unwrap();
        let nu = Measure::dirac(&sa, "b").unwrap();
        let elems = vec!["a".into(), "b".into()];
        let decomp = LebesgueDecomposition::compute(&nu, &mu, &elems, &sa);
        assert_eq!(decomp.is_pure(), PureType::Singular);
    }

    #[test]
    fn test_mixed_decomposition() {
        let (sa, mu, nu) = setup_decomp();
        let elems = vec!["a".into(), "b".into(), "c".into(), "d".into()];
        let decomp = LebesgueDecomposition::compute(&nu, &mu, &elems, &sa);
        assert_eq!(decomp.is_pure(), PureType::Mixed);
    }

    #[test]
    fn test_hahn_decomposition() {
        let values = vec![
            ("a".into(), 2.0),
            ("b".into(), -1.0),
            ("c".into(), 0.5),
            ("d".into(), -3.0),
        ];
        let hahn = HahnDecomposition::compute(&values);
        assert!(hahn.positive.contains("a"));
        assert!(hahn.positive.contains("c"));
        assert!(hahn.negative.contains("b"));
        assert!(hahn.negative.contains("d"));
    }

    #[test]
    fn test_hahn_all_positive() {
        let values = vec![("a".into(), 1.0), ("b".into(), 2.0)];
        let hahn = HahnDecomposition::compute(&values);
        assert!(hahn.negative.is_empty());
    }

    #[test]
    fn test_hahn_all_negative() {
        let values = vec![("a".into(), -1.0), ("b".into(), -2.0)];
        let hahn = HahnDecomposition::compute(&values);
        assert!(hahn.positive.is_empty());
    }

    #[test]
    fn test_decomposition_same_measure() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let mu = Measure::uniform(&sa).unwrap();
        let elems = vec!["a".into(), "b".into()];
        let decomp = LebesgueDecomposition::compute(&mu, &mu, &elems, &sa);
        assert_eq!(decomp.is_pure(), PureType::AbsolutelyContinuous);
    }

    #[test]
    fn test_decomposition_zero_reference() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a"]));
        let mu = Measure::new("zero", &sa, &[]).unwrap();
        let nu = Measure::dirac(&sa, "a").unwrap();
        let elems = vec!["a".into()];
        let decomp = LebesgueDecomposition::compute(&nu, &mu, &elems, &sa);
        assert_eq!(decomp.is_pure(), PureType::Singular);
    }
}
