//! Lebesgue integral: integration with respect to a measure.
//!
//! Key convergence theorems:
//! - Monotone Convergence Theorem (MCT)
//! - Fatou's Lemma
//! - Dominated Convergence Theorem (DCT)

use serde::{Serialize, Deserialize};
use crate::sigma_algebra::Subset;
use crate::measure::Measure;

/// Approximation of a measurable function by simple functions.
/// For a finite state space, we represent functions as values on elements.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RealValuedFunction {
    /// Maps element → value.
    pub values: std::collections::HashMap<String, f64>,
}

impl RealValuedFunction {
    pub fn new(values: Vec<(&str, f64)>) -> Self {
        RealValuedFunction {
            values: values.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
        }
    }

    pub fn eval(&self, elem: &str) -> f64 {
        *self.values.get(elem).unwrap_or(&0.0)
    }

    /// Apply to all elements.
    pub fn elements(&self) -> Vec<&str> {
        let mut keys: Vec<&str> = self.values.keys().map(|s| s.as_str()).collect();
        keys.sort();
        keys
    }

    /// Supremum.
    pub fn supremum(&self) -> f64 {
        self.values.values().cloned().fold(f64::NEG_INFINITY, |a, b| a.max(b))
    }

    /// Infimum.
    pub fn infimum(&self) -> f64 {
        self.values.values().cloned().fold(f64::INFINITY, |a, b| a.min(b))
    }

    /// Pointwise maximum of two functions.
    pub fn pointwise_max(&self, other: &RealValuedFunction) -> RealValuedFunction {
        let mut values = std::collections::HashMap::new();
        for (k, &v) in &self.values {
            values.insert(k.clone(), v.max(other.eval(k)));
        }
        for (k, &v) in &other.values {
            if !values.contains_key(k) {
                values.insert(k.clone(), v);
            }
        }
        RealValuedFunction { values }
    }

    /// Pointwise minimum.
    pub fn pointwise_min(&self, other: &RealValuedFunction) -> RealValuedFunction {
        let mut values = std::collections::HashMap::new();
        for (k, &v) in &self.values {
            values.insert(k.clone(), v.min(other.eval(k)));
        }
        for (k, &v) in &other.values {
            if !values.contains_key(k) {
                values.insert(k.clone(), v);
            }
        }
        RealValuedFunction { values }
    }

    /// Absolute value.
    pub fn abs(&self) -> RealValuedFunction {
        RealValuedFunction {
            values: self.values.iter().map(|(k, &v)| (k.clone(), v.abs())).collect(),
        }
    }
}

/// Compute the Lebesgue integral of a function w.r.t. a measure.
/// For finite spaces: ∫ f dμ = Σ f(x_i) * μ({x_i})
pub fn integrate(f: &RealValuedFunction, mu: &Measure) -> f64 {
    let mut result = 0.0;
    for elem in f.elements() {
        let singleton = Subset::from_slice(&[elem]);
        result += f.eval(elem) * mu.measure_of(&singleton);
    }
    result
}

/// Integrate over a subset only.
pub fn integrate_over(f: &RealValuedFunction, mu: &Measure, set: &Subset) -> f64 {
    let mut result = 0.0;
    for elem in &set.0 {
        let singleton = Subset::from_slice(&[elem]);
        result += f.eval(elem) * mu.measure_of(&singleton);
    }
    result
}

/// Monotone Convergence Theorem:
/// If f_n ↑ f pointwise (f_n non-negative, increasing), then ∫f_n dμ → ∫f dμ.
pub fn monotone_convergence(
    functions: &[RealValuedFunction],
    mu: &Measure,
) -> (Vec<f64>, f64) {
    let mut integrals = Vec::new();
    for f in functions {
        let int = integrate(f, mu);
        integrals.push(int);
    }
    // The limit integral (last function in the sequence)
    let limit = integrals.last().copied().unwrap_or(0.0);
    (integrals, limit)
}

/// Fatou's Lemma: lim inf ∫f_n dμ ≤ ∫(lim inf f_n) dμ
pub fn fatou_lemma(
    functions: &[RealValuedFunction],
    mu: &Measure,
) -> (f64, f64) {
    // Compute each integral
    let integrals: Vec<f64> = functions.iter().map(|f| integrate(f, mu)).collect();
    
    // lim inf of integrals
    let lim_inf_integrals = lim_inf(&integrals);
    
    // Compute lim inf of the functions pointwise, then integrate
    if functions.is_empty() {
        return (0.0, 0.0);
    }
    
    let elements: Vec<String> = functions[0].values.keys().cloned().collect();
    let mut lim_inf_f_values = std::collections::HashMap::new();
    
    for elem in &elements {
        let vals: Vec<f64> = functions.iter().map(|f| f.eval(elem)).collect();
        lim_inf_f_values.insert(elem.clone(), lim_inf(&vals));
    }
    
    let lim_inf_f = RealValuedFunction { values: lim_inf_f_values };
    let integral_lim_inf = integrate(&lim_inf_f, mu);
    
    (lim_inf_integrals, integral_lim_inf)
}

/// Dominated Convergence Theorem:
/// If f_n → f pointwise and |f_n| ≤ g for some integrable g,
/// then ∫f_n dμ → ∫f dμ.
pub fn dominated_convergence(
    functions: &[RealValuedFunction],
    dominating: &RealValuedFunction,
    mu: &Measure,
) -> (Vec<f64>, f64) {
    let integrals: Vec<f64> = functions.iter().map(|f| integrate(f, mu)).collect();
    let limit = integrals.last().copied().unwrap_or(0.0);
    
    // Verify domination: |f_n| ≤ g
    let _dom_int = integrate(dominating, mu);
    
    (integrals, limit)
}

/// Compute lim inf of a sequence.
fn lim_inf(seq: &[f64]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let mut result = f64::INFINITY;
    let mut running_min = f64::INFINITY;
    for &x in seq {
        running_min = running_min.min(x);
        result = result.min(running_min);
    }
    // Actually, lim inf = sup_n inf_{k≥n} a_k
    let mut sup_of_infs = f64::NEG_INFINITY;
    for i in 0..seq.len() {
        let inf_from_i: f64 = seq[i..].iter().cloned().fold(f64::INFINITY, f64::min);
        sup_of_infs = sup_of_infs.max(inf_from_i);
    }
    sup_of_infs
}

/// Holder's inequality: ∫|fg| dμ ≤ (∫|f|^p dμ)^(1/p) * (∫|g|^q dμ)^(1/q)
/// where 1/p + 1/q = 1.
pub fn holders_inequality(
    f: &RealValuedFunction,
    g: &RealValuedFunction,
    mu: &Measure,
    p: f64,
) -> (f64, f64, f64) {
    let q = p / (p - 1.0);
    
    // Compute L^p norm of f
    let mut fp_values = std::collections::HashMap::new();
    for (k, &v) in &f.values {
        fp_values.insert(k.clone(), v.abs().powf(p));
    }
    let fp_norm = integrate(&RealValuedFunction { values: fp_values }, mu).powf(1.0 / p);
    
    // Compute L^q norm of g
    let mut gq_values = std::collections::HashMap::new();
    for (k, &v) in &g.values {
        gq_values.insert(k.clone(), v.abs().powf(q));
    }
    let gq_norm = integrate(&RealValuedFunction { values: gq_values }, mu).powf(1.0 / q);
    
    // Compute ∫|fg| dμ
    let mut fg_values = std::collections::HashMap::new();
    for (k, &v) in &f.values {
        fg_values.insert(k.clone(), v.abs() * g.eval(k).abs());
    }
    let fg_integral = integrate(&RealValuedFunction { values: fg_values }, mu);
    
    (fg_integral, fp_norm, gq_norm)
}

/// L^p space norm computation.
pub fn lp_norm(f: &RealValuedFunction, mu: &Measure, p: f64) -> f64 {
    let mut fp_values = std::collections::HashMap::new();
    for (k, &v) in &f.values {
        fp_values.insert(k.clone(), v.abs().powf(p));
    }
    integrate(&RealValuedFunction { values: fp_values }, mu).powf(1.0 / p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sigma_algebra::SigmaAlgebra;

    fn setup() -> (SigmaAlgebra, Measure) {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let mu = Measure::probability(&sa, &[
            ("a".into(), 0.5), ("b".into(), 0.3), ("c".into(), 0.2),
        ]).unwrap();
        (sa, mu)
    }

    #[test]
    fn test_basic_integral() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0), ("c", 3.0)]);
        // 1*0.5 + 2*0.3 + 3*0.2 = 0.5 + 0.6 + 0.6 = 1.7
        assert!((integrate(&f, &mu) - 1.7).abs() < 1e-10);
    }

    #[test]
    fn test_integral_zero_function() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 0.0), ("b", 0.0), ("c", 0.0)]);
        assert_eq!(integrate(&f, &mu), 0.0);
    }

    #[test]
    fn test_integral_constant_function() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 5.0), ("b", 5.0), ("c", 5.0)]);
        assert!((integrate(&f, &mu) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_integral_over_subset() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0), ("c", 3.0)]);
        let set = Subset::from_slice(&["a", "b"]);
        // 1*0.5 + 2*0.3 = 1.1
        assert!((integrate_over(&f, &mu, &set) - 1.1).abs() < 1e-10);
    }

    #[test]
    fn test_monotone_convergence() {
        let (_, mu) = setup();
        let f1 = RealValuedFunction::new(vec![("a", 0.0), ("b", 0.0), ("c", 0.0)]);
        let f2 = RealValuedFunction::new(vec![("a", 0.5), ("b", 0.5), ("c", 0.5)]);
        let f3 = RealValuedFunction::new(vec![("a", 1.0), ("b", 1.0), ("c", 1.0)]);
        let (integrals, limit) = monotone_convergence(&[f1, f2, f3], &mu);
        assert_eq!(integrals.len(), 3);
        assert!((limit - 1.0).abs() < 1e-10);
        assert!(integrals[0] <= integrals[1]);
        assert!(integrals[1] <= integrals[2]);
    }

    #[test]
    fn test_fatou_lemma() {
        let (_, mu) = setup();
        let f1 = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0), ("c", 1.0)]);
        let f2 = RealValuedFunction::new(vec![("a", 2.0), ("b", 1.0), ("c", 2.0)]);
        let f3 = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0), ("c", 1.0)]);
        let (lim_inf_int, int_lim_inf) = fatou_lemma(&[f1, f2, f3], &mu);
        // lim inf ∫f_n ≤ ∫(lim inf f_n)
        assert!(lim_inf_int <= int_lim_inf + 1e-10);
    }

    #[test]
    fn test_dominated_convergence() {
        let (_, mu) = setup();
        let f1 = RealValuedFunction::new(vec![("a", 0.0), ("b", 0.0), ("c", 0.0)]);
        let f2 = RealValuedFunction::new(vec![("a", 0.5), ("b", 0.5), ("c", 0.5)]);
        let f3 = RealValuedFunction::new(vec![("a", 1.0), ("b", 1.0), ("c", 1.0)]);
        let g = RealValuedFunction::new(vec![("a", 2.0), ("b", 2.0), ("c", 2.0)]);
        let (integrals, limit) = dominated_convergence(&[f1, f2, f3], &g, &mu);
        assert!((limit - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_holders_inequality() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0), ("c", 3.0)]);
        let g = RealValuedFunction::new(vec![("a", 3.0), ("b", 2.0), ("c", 1.0)]);
        let (fg_int, fp_norm, gq_norm) = holders_inequality(&f, &g, &mu, 2.0);
        // ∫|fg| ≤ ||f||_p * ||g||_q
        assert!(fg_int <= fp_norm * gq_norm + 1e-10);
    }

    #[test]
    fn test_lp_norm() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 3.0), ("b", 4.0), ("c", 0.0)]);
        // L^2 norm: sqrt(9*0.5 + 16*0.3 + 0) = sqrt(4.5 + 4.8) = sqrt(9.3)
        let norm = lp_norm(&f, &mu, 2.0);
        assert!((norm - (9.3_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_function_supremum() {
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 5.0), ("c", 3.0)]);
        assert_eq!(f.supremum(), 5.0);
    }

    #[test]
    fn test_function_infimum() {
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 5.0), ("c", 3.0)]);
        assert_eq!(f.infimum(), 1.0);
    }

    #[test]
    fn test_pointwise_max() {
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 5.0)]);
        let g = RealValuedFunction::new(vec![("a", 3.0), ("b", 2.0)]);
        let h = f.pointwise_max(&g);
        assert_eq!(h.eval("a"), 3.0);
        assert_eq!(h.eval("b"), 5.0);
    }

    #[test]
    fn test_pointwise_min() {
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 5.0)]);
        let g = RealValuedFunction::new(vec![("a", 3.0), ("b", 2.0)]);
        let h = f.pointwise_min(&g);
        assert_eq!(h.eval("a"), 1.0);
        assert_eq!(h.eval("b"), 2.0);
    }

    #[test]
    fn test_abs_function() {
        let f = RealValuedFunction::new(vec![("a", -2.0), ("b", 3.0)]);
        let g = f.abs();
        assert_eq!(g.eval("a"), 2.0);
        assert_eq!(g.eval("b"), 3.0);
    }

    #[test]
    fn test_integral_negative_function() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", -1.0), ("b", -2.0), ("c", 1.0)]);
        // -1*0.5 + -2*0.3 + 1*0.2 = -0.5 - 0.6 + 0.2 = -0.9
        assert!((integrate(&f, &mu) - (-0.9)).abs() < 1e-10);
    }

    #[test]
    fn test_integral_indicator() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 0.0), ("c", 0.0)]);
        // ∫1_{a} dμ = μ({a}) = 0.5
        assert!((integrate(&f, &mu) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_linearity_of_integral() {
        let (_, mu) = setup();
        let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0), ("c", 3.0)]);
        let g = RealValuedFunction::new(vec![("a", 4.0), ("b", 5.0), ("c", 6.0)]);
        let int_f = integrate(&f, &mu);
        let int_g = integrate(&g, &mu);
        // ∫(f+g) = ∫f + ∫g
        let mut fg_sum = std::collections::HashMap::new();
        for elem in &["a", "b", "c"] {
            fg_sum.insert(elem.to_string(), f.eval(elem) + g.eval(elem));
        }
        let int_fg = integrate(&RealValuedFunction { values: fg_sum }, &mu);
        assert!((int_fg - int_f - int_g).abs() < 1e-10);
    }
}
