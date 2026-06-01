//! Riesz representation theorem.
//!
//! Every positive linear functional L on C_c(X) (continuous functions with compact support)
//! can be represented as L(f) = ∫ f dμ for a unique regular Borel measure μ.
//!
//! For finite spaces, every positive linear functional on ℝ^n corresponds to integration
//! against a non-negative measure.

use serde::{Serialize, Deserialize};
use crate::sigma_algebra::Subset;
use crate::measure::Measure;

/// A linear functional on a finite-dimensional space of functions.
/// Represented by its values on basis functions (indicator functions of singletons).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearFunctional {
    name: String,
    /// Maps element → value of L(1_{element}).
    basis_values: std::collections::HashMap<String, f64>,
}

impl LinearFunctional {
    /// Create a linear functional from values on basis functions.
    pub fn new(name: &str, values: Vec<(&str, f64)>) -> Self {
        LinearFunctional {
            name: name.to_string(),
            basis_values: values.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
        }
    }

    /// Evaluate the functional on a function (given by values on elements).
    /// L(f) = Σ f(x_i) * L(1_{x_i})  — linearity.
    pub fn evaluate(&self, f_values: &std::collections::HashMap<String, f64>) -> f64 {
        let mut result = 0.0;
        for (elem, &f_val) in f_values {
            let basis_val = self.basis_values.get(elem).copied().unwrap_or(0.0);
            result += f_val * basis_val;
        }
        result
    }

    /// Is this a positive functional? L(f) ≥ 0 whenever f ≥ 0.
    pub fn is_positive(&self) -> bool {
        self.basis_values.values().all(|&v| v >= 0.0)
    }

    /// The norm of the functional.
    pub fn norm(&self) -> f64 {
        self.basis_values.values().map(|v| v.abs()).sum()
    }

    /// Scale the functional.
    pub fn scale(&self, c: f64) -> LinearFunctional {
        LinearFunctional {
            name: format!("{}*{}", c, self.name),
            basis_values: self.basis_values.iter().map(|(k, &v)| (k.clone(), v * c)).collect(),
        }
    }

    /// Add two functionals.
    pub fn add(&self, other: &LinearFunctional) -> LinearFunctional {
        let mut values = self.basis_values.clone();
        for (k, &v) in &other.basis_values {
            *values.entry(k.clone()).or_insert(0.0) += v;
        }
        LinearFunctional {
            name: format!("{}+{}", self.name, other.name),
            basis_values: values,
        }
    }

    /// Get the name.
    pub fn name(&self) -> &str {
        &self.name
    }
}

/// Riesz representation: convert a positive linear functional to a measure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RieszRepresentation {
    /// The representing measure.
    pub measure: Measure,
    /// The original functional.
    pub functional_name: String,
}

impl RieszRepresentation {
    /// Compute the Riesz representation of a positive linear functional.
    ///
    /// For finite spaces: μ({x}) = L(1_{x}).
    pub fn represent(
        functional: &LinearFunctional,
        sa: &crate::sigma_algebra::SigmaAlgebra,
    ) -> Result<Self, String> {
        if !functional.is_positive() {
            return Err("Functional must be positive for Riesz representation".to_string());
        }

        let elements: Vec<String> = sa.universal().0.iter().cloned().collect();
        let values: Vec<(Subset, f64)> = elements.iter().map(|elem| {
            let singleton = Subset::from_slice(&[elem]);
            let val = functional.basis_values.get(elem).copied().unwrap_or(0.0);
            (singleton, val)
        }).collect();

        let measure = Measure::new(
            &format!("μ_{}", functional.name),
            sa,
            &values,
        )?;

        Ok(RieszRepresentation {
            measure,
            functional_name: functional.name.clone(),
        })
    }

    /// Verify: L(f) should equal ∫f dμ for test functions.
    pub fn verify(
        &self,
        functional: &LinearFunctional,
        test_functions: &[std::collections::HashMap<String, f64>],
    ) -> bool {
        for f in test_functions {
            let l_val = functional.evaluate(f);
            let int_val: f64 = f.iter().map(|(elem, &f_val)| {
                let singleton = Subset::from_slice(&[elem]);
                f_val * self.measure.measure_of(&singleton)
            }).sum();
            if (l_val - int_val).abs() > 1e-10 {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sigma_algebra::SigmaAlgebra;

    fn setup_riesz() -> (SigmaAlgebra, LinearFunctional) {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let functional = LinearFunctional::new("L", vec![
            ("a", 0.5), ("b", 0.3), ("c", 0.2),
        ]);
        (sa, functional)
    }

    #[test]
    fn test_functional_is_positive() {
        let (_, f) = setup_riesz();
        assert!(f.is_positive());
    }

    #[test]
    fn test_functional_not_positive() {
        let f = LinearFunctional::new("L", vec![("a", 1.0), ("b", -0.5)]);
        assert!(!f.is_positive());
    }

    #[test]
    fn test_functional_evaluate() {
        let (_, f) = setup_riesz();
        let mut vals = std::collections::HashMap::new();
        vals.insert("a".into(), 2.0);
        vals.insert("b".into(), 3.0);
        vals.insert("c".into(), 1.0);
        // L(f) = 2*0.5 + 3*0.3 + 1*0.2 = 1.0 + 0.9 + 0.2 = 2.1
        assert!((f.evaluate(&vals) - 2.1).abs() < 1e-10);
    }

    #[test]
    fn test_functional_scale() {
        let (_, f) = setup_riesz();
        let g = f.scale(2.0);
        assert!((g.basis_values["a"] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_functional_add() {
        let f1 = LinearFunctional::new("L1", vec![("a", 1.0)]);
        let f2 = LinearFunctional::new("L2", vec![("a", 2.0)]);
        let sum = f1.add(&f2);
        assert!((sum.basis_values["a"] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_riesz_representation() {
        let (sa, functional) = setup_riesz();
        let rep = RieszRepresentation::represent(&functional, &sa).unwrap();
        assert!((rep.measure.measure_of(&Subset::from_slice(&["a"])) - 0.5).abs() < 1e-10);
        assert!((rep.measure.measure_of(&Subset::from_slice(&["b"])) - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_riesz_verify() {
        let (sa, functional) = setup_riesz();
        let rep = RieszRepresentation::represent(&functional, &sa).unwrap();
        
        let mut f1 = std::collections::HashMap::new();
        f1.insert("a".into(), 1.0); f1.insert("b".into(), 0.0); f1.insert("c".into(), 0.0);
        let mut f2 = std::collections::HashMap::new();
        f2.insert("a".into(), 2.0); f2.insert("b".into(), 3.0); f2.insert("c".into(), 1.0);
        
        assert!(rep.verify(&functional, &[f1, f2]));
    }

    #[test]
    fn test_riesz_negative_rejected() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a"]));
        let f = LinearFunctional::new("neg", vec![("a", -1.0)]);
        assert!(RieszRepresentation::represent(&f, &sa).is_err());
    }

    #[test]
    fn test_functional_norm() {
        let (_, f) = setup_riesz();
        assert!((f.norm() - 1.0).abs() < 1e-10); // 0.5 + 0.3 + 0.2
    }

    #[test]
    fn test_functional_name() {
        let f = LinearFunctional::new("my_L", vec![("a", 1.0)]);
        assert_eq!(f.name(), "my_L");
    }

    #[test]
    fn test_riesz_total_mass() {
        let (sa, functional) = setup_riesz();
        let rep = RieszRepresentation::represent(&functional, &sa).unwrap();
        // Total mass = L(1) = 0.5 + 0.3 + 0.2 = 1.0
        // Check by summing singleton measures
        let total: f64 = sa.universal().0.iter()
            .map(|e| rep.measure.measure_of(&Subset::from_slice(&[e])))
            .sum();
        assert!((total - 1.0).abs() < 1e-10);
    }
}
