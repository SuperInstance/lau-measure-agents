//! Radon-Nikodym theorem and derivative.
//!
//! If ν ≪ μ (ν absolutely continuous w.r.t. μ), then there exists a measurable function
//! f = dν/dμ such that ν(A) = ∫_A f dμ for all measurable A.
//! The function f is the Radon-Nikodym derivative.

use serde::{Serialize, Deserialize};
use crate::sigma_algebra::Subset;
use crate::measure::Measure;
use crate::integral::{RealValuedFunction, integrate};

/// Radon-Nikodym derivative dν/dμ.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RadonNikodymDerivative {
    /// The derivative function f where ν(A) = ∫_A f dμ.
    pub derivative: RealValuedFunction,
    /// Name of the reference measure μ.
    pub reference_name: String,
    /// Name of the target measure ν.
    pub target_name: String,
}

impl RadonNikodymDerivative {
    /// Compute the Radon-Nikodym derivative dν/dμ.
    ///
    /// For finite spaces, f(x) = ν({x}) / μ({x}) for each element x.
    /// Requires ν ≪ μ.
    pub fn compute(
        nu: &Measure,
        mu: &Measure,
        elements: &[String],
    ) -> Result<Self, String> {
        // Verify absolute continuity
        for elem in elements {
            let singleton = Subset::from_slice(&[elem]);
            let mu_val = mu.measure_of(&singleton);
            let nu_val = nu.measure_of(&singleton);
            if nu_val > 0.0 && mu_val == 0.0 {
                return Err(format!(
                    "ν is not absolutely continuous w.r.t. μ: ν({{{}}}) = {} but μ({{{}}}) = 0",
                    elem, nu_val, elem
                ));
            }
        }

        // Compute derivative pointwise
        let values: Vec<(&str, f64)> = elements.iter().map(|elem| {
            let singleton = Subset::from_slice(&[elem]);
            let mu_val = mu.measure_of(&singleton);
            let nu_val = nu.measure_of(&singleton);
            let f_val = if mu_val > 0.0 { nu_val / mu_val } else { 0.0 };
            (elem.as_str(), f_val)
        }).collect();

        Ok(RadonNikodymDerivative {
            derivative: RealValuedFunction::new(values),
            reference_name: mu.name().to_string(),
            target_name: nu.name().to_string(),
        })
    }

    /// Evaluate the derivative at a point.
    pub fn eval(&self, elem: &str) -> f64 {
        self.derivative.eval(elem)
    }

    /// Verify: ν(A) should equal ∫_A f dμ for all measurable A.
    pub fn verify(&self, nu: &Measure, mu: &Measure, test_sets: &[Subset]) -> bool {
        for set in test_sets {
            let nu_val = nu.measure_of(set);
            let integral_val = crate::integral::integrate_over(&self.derivative, mu, set);
            if (nu_val - integral_val).abs() > 1e-8 {
                return false;
            }
        }
        true
    }

    /// Chain rule: if dλ/dν = g and dν/dμ = f, then dλ/dμ = g·f.
    pub fn chain_rule(&self, other: &RadonNikodymDerivative) -> RadonNikodymDerivative {
        let mut values = std::collections::HashMap::new();
        for elem in self.derivative.elements() {
            values.insert(elem.to_string(), self.eval(elem) * other.eval(elem));
        }
        RadonNikodymDerivative {
            derivative: RealValuedFunction { values },
            reference_name: other.reference_name.clone(),
            target_name: self.target_name.clone(),
        }
    }

    /// Log-derivative (useful for belief updates).
    pub fn log_derivative(&self) -> RealValuedFunction {
        let mut values = std::collections::HashMap::new();
        for (k, &v) in &self.derivative.values {
            values.insert(k.clone(), if v > 0.0 { v.ln() } else { f64::NEG_INFINITY });
        }
        RealValuedFunction { values }
    }
}

/// Bayesian belief update via Radon-Nikodym derivative.
///
/// Given prior μ and observation likelihood, compute posterior ν.
/// The RN derivative dν/dμ represents the belief update ratio.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BeliefUpdate {
    /// The Radon-Nikodym derivative encoding the update.
    pub rn_derivative: RadonNikodymDerivative,
    /// Normalization constant.
    pub normalization: f64,
}

impl BeliefUpdate {
    /// Create a belief update from likelihood function.
    /// posterior(x) ∝ likelihood(x) · prior(x)
    pub fn from_likelihood(
        likelihood: &RealValuedFunction,
        prior: &Measure,
        elements: &[String],
    ) -> Result<Self, String> {
        // Compute unnormalized posterior weights
        let mut unnorm: Vec<(String, f64)> = elements.iter().map(|elem| {
            let singleton = Subset::from_slice(&[elem]);
            let prior_weight = prior.measure_of(&singleton);
            let like = likelihood.eval(elem);
            (elem.clone(), like * prior_weight)
        }).collect();

        // Normalize
        let total: f64 = unnorm.iter().map(|(_, w)| *w).sum();
        if total <= 0.0 {
            return Err("Total unnormalized weight is non-positive".to_string());
        }

        let sa = crate::sigma_algebra::SigmaAlgebra::power_set(
            Subset(elements.iter().cloned().collect())
        );
        let weights: Vec<(String, f64)> = unnorm.iter().map(|(e, w)| (e.clone(), w / total)).collect();
        let posterior = Measure::probability(&sa, &weights)?;

        let rn = RadonNikodymDerivative::compute(&posterior, prior, elements)?;

        Ok(BeliefUpdate {
            rn_derivative: rn,
            normalization: total,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sigma_algebra::SigmaAlgebra;

    fn setup_rn() -> (SigmaAlgebra, Measure, Measure) {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let mu = Measure::uniform(&sa).unwrap();
        let nu = Measure::probability(&sa, &[
            ("a".into(), 0.5), ("b".into(), 0.3), ("c".into(), 0.2),
        ]).unwrap();
        (sa, mu, nu)
    }

    #[test]
    fn test_rn_derivative_compute() {
        let (sa, mu, nu) = setup_rn();
        let elems = vec!["a".into(), "b".into(), "c".into()];
        let rn = RadonNikodymDerivative::compute(&nu, &mu, &elems).unwrap();
        // f(a) = 0.5 / (1/3) = 1.5
        assert!((rn.eval("a") - 1.5).abs() < 1e-10);
        // f(b) = 0.3 / (1/3) = 0.9
        assert!((rn.eval("b") - 0.9).abs() < 1e-10);
        // f(c) = 0.2 / (1/3) = 0.6
        assert!((rn.eval("c") - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_rn_verify() {
        let (sa, mu, nu) = setup_rn();
        let elems = vec!["a".into(), "b".into(), "c".into()];
        let rn = RadonNikodymDerivative::compute(&nu, &mu, &elems).unwrap();
        
        let test_sets = vec![
            Subset::from_slice(&["a"]),
            Subset::from_slice(&["b", "c"]),
            Subset::from_slice(&["a", "b", "c"]),
        ];
        assert!(rn.verify(&nu, &mu, &test_sets));
    }

    #[test]
    fn test_rn_not_absolutely_continuous() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let mu = Measure::dirac(&sa, "a").unwrap();
        let nu = Measure::dirac(&sa, "b").unwrap();
        let elems = vec!["a".into(), "b".into()];
        let result = RadonNikodymDerivative::compute(&nu, &mu, &elems);
        assert!(result.is_err());
    }

    #[test]
    fn test_rn_chain_rule() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let mu = Measure::uniform(&sa).unwrap();
        let nu = Measure::probability(&sa, &[
            ("a".into(), 0.6), ("b".into(), 0.4),
        ]).unwrap();
        let lambda = Measure::probability(&sa, &[
            ("a".into(), 0.8), ("b".into(), 0.2),
        ]).unwrap();

        let elems = vec!["a".into(), "b".into()];
        let rn1 = RadonNikodymDerivative::compute(&nu, &mu, &elems).unwrap();
        let rn2 = RadonNikodymDerivative::compute(&lambda, &nu, &elems).unwrap();
        let chain = rn1.chain_rule(&rn2);

        // dλ/dμ should give λ({a})/μ({a}) = 0.8/0.5 = 1.6
        assert!((chain.eval("a") - 1.6).abs() < 1e-8);
    }

    #[test]
    fn test_rn_log_derivative() {
        let (sa, mu, nu) = setup_rn();
        let elems = vec!["a".into(), "b".into(), "c".into()];
        let rn = RadonNikodymDerivative::compute(&nu, &mu, &elems).unwrap();
        let log_d = rn.log_derivative();
        assert!((log_d.eval("a") - 1.5_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_belief_update() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let prior = Measure::uniform(&sa).unwrap();
        let likelihood = RealValuedFunction::new(vec![("a", 0.9), ("b", 0.1), ("c", 0.0)]);
        let elems = vec!["a".into(), "b".into(), "c".into()];
        
        let update = BeliefUpdate::from_likelihood(&likelihood, &prior, &elems).unwrap();
        // Posterior should concentrate on "a"
        assert!(update.rn_derivative.eval("a") > update.rn_derivative.eval("b"));
    }

    #[test]
    fn test_belief_update_normalization() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let prior = Measure::uniform(&sa).unwrap();
        let likelihood = RealValuedFunction::new(vec![("a", 1.0), ("b", 1.0)]);
        let elems = vec!["a".into(), "b".into()];
        
        let update = BeliefUpdate::from_likelihood(&likelihood, &prior, &elems).unwrap();
        assert!((update.normalization - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rn_same_measure() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let mu = Measure::uniform(&sa).unwrap();
        let elems = vec!["a".into(), "b".into()];
        let rn = RadonNikodymDerivative::compute(&mu, &mu, &elems).unwrap();
        // dμ/dμ = 1 everywhere
        assert!((rn.eval("a") - 1.0).abs() < 1e-10);
        assert!((rn.eval("b") - 1.0).abs() < 1e-10);
    }
}
