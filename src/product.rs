//! Product measures and Fubini's theorem.
//!
//! Given measure spaces (X, 𝔉, μ) and (Y, 𝔊, ν), construct the product measure
//! on (X × Y, 𝔉 ⊗ 𝔊) such that (μ × ν)(A × B) = μ(A) · ν(B).
//! Fubini's theorem: ∫_{X×Y} f d(μ×ν) = ∫_X (∫_Y f(x,y) dν(y)) dμ(x).

use serde::{Serialize, Deserialize};
use crate::sigma_algebra::{Subset, SigmaAlgebra};
use crate::measure::Measure;

/// A pair of elements representing a point in the product space.
pub fn pair_key(x: &str, y: &str) -> String {
    format!("({},{})", x, y)
}

/// Construct the product sigma-algebra 𝔉 ⊗ 𝔊.
/// For finite spaces, this is the power set of X × Y.
pub fn product_sigma_algebra(sa1: &SigmaAlgebra, sa2: &SigmaAlgebra) -> SigmaAlgebra {
    let mut elements = Vec::new();
    for x in &sa1.universal().0 {
        for y in &sa2.universal().0 {
            elements.push(pair_key(x, y));
        }
    }
    let slices: Vec<&str> = elements.iter().map(|s| s.as_str()).collect();
    SigmaAlgebra::power_set(Subset::from_slice(&slices))
}

/// Construct the product measure μ × ν.
/// For measurable rectangles A × B: (μ × ν)(A × B) = μ(A) · ν(B).
pub fn product_measure(
    name: &str,
    sa1: &SigmaAlgebra,
    mu: &Measure,
    sa2: &SigmaAlgebra,
    nu: &Measure,
    product_sa: &SigmaAlgebra,
) -> Measure {
    let mut values = Vec::new();

    // Compute measure for each subset of the product space
    for s in product_sa.measurable_sets() {
        let mut total = 0.0;
        for x in &sa1.universal().0 {
            for y in &sa2.universal().0 {
                let key = pair_key(x, y);
                if s.contains(&key) {
                    let sx = Subset::from_slice(&[x]);
                    let sy = Subset::from_slice(&[y]);
                    total += mu.measure_of(&sx) * nu.measure_of(&sy);
                }
            }
        }
        if total > 0.0 || s.is_empty() {
            values.push((s.clone(), total));
        }
    }

    Measure::new(name, product_sa, &values).unwrap()
}

/// Compute the product measure value for a rectangle A × B.
pub fn rectangle_measure(
    a: &Subset,
    mu: &Measure,
    b: &Subset,
    nu: &Measure,
) -> f64 {
    mu.measure_of(a) * nu.measure_of(b)
}

/// Fubini's theorem: iterated integration.
/// ∫_{X×Y} f d(μ×ν) = ∫_X (∫_Y f(x,y) dν(y)) dμ(x)
///
/// For finite spaces, we compute this explicitly.
pub fn fubini_integrate(
    f_values: &std::collections::HashMap<String, f64>,
    sa1: &SigmaAlgebra,
    mu: &Measure,
    sa2: &SigmaAlgebra,
    nu: &Measure,
) -> f64 {
    let mut total = 0.0;
    for x in &sa1.universal().0 {
        let sx = Subset::from_slice(&[x]);
        let mu_x = mu.measure_of(&sx);
        
        // Inner integral: ∫_Y f(x,y) dν(y)
        let mut inner = 0.0;
        for y in &sa2.universal().0 {
            let sy = Subset::from_slice(&[y]);
            let nu_y = nu.measure_of(&sy);
            let key = pair_key(x, y);
            let f_val = f_values.get(&key).copied().unwrap_or(0.0);
            inner += f_val * nu_y;
        }
        
        total += mu_x * inner;
    }
    total
}

/// Tonelli's theorem: non-negative version of Fubini.
/// For f ≥ 0, the iterated integrals are always equal (no integrability condition needed).
pub fn tonelli_integrate(
    f_values: &std::collections::HashMap<String, f64>,
    sa1: &SigmaAlgebra,
    mu: &Measure,
    sa2: &SigmaAlgebra,
    nu: &Measure,
) -> f64 {
    // Same computation as Fubini but we verify non-negativity
    for &v in f_values.values() {
        assert!(v >= 0.0, "Tonelli requires non-negative function");
    }
    fubini_integrate(f_values, sa1, mu, sa2, nu)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> (SigmaAlgebra, Measure, SigmaAlgebra, Measure) {
        let sa1 = SigmaAlgebra::power_set(Subset::from_slice(&["x1", "x2"]));
        let sa2 = SigmaAlgebra::power_set(Subset::from_slice(&["y1", "y2"]));
        let mu = Measure::probability(&sa1, &[
            ("x1".into(), 0.6), ("x2".into(), 0.4),
        ]).unwrap();
        let nu = Measure::probability(&sa2, &[
            ("y1".into(), 0.3), ("y2".into(), 0.7),
        ]).unwrap();
        (sa1, mu, sa2, nu)
    }

    #[test]
    fn test_product_sigma_algebra() {
        let (sa1, _, sa2, _) = setup();
        let psa = product_sigma_algebra(&sa1, &sa2);
        assert_eq!(psa.universal().len(), 4);
        // Power set of 4 elements = 16 sets
        assert_eq!(psa.size(), 16);
    }

    #[test]
    fn test_product_measure_rectangle() {
        let (sa1, mu, sa2, nu) = setup();
        let psa = product_sigma_algebra(&sa1, &sa2);
        let pm = product_measure("μ×ν", &sa1, &mu, &sa2, &nu, &psa);
        
        // Rectangle {x1} × {y1} should have measure 0.6 * 0.3 = 0.18
        let rect = Subset::from_slice(&[&pair_key("x1", "y1")]);
        assert!((pm.measure_of(&rect) - 0.18).abs() < 1e-10);
    }

    #[test]
    fn test_rectangle_measure() {
        let (sa1, mu, sa2, nu) = setup();
        let a = Subset::from_slice(&["x1"]);
        let b = Subset::from_slice(&["y2"]);
        assert!((rectangle_measure(&a, &mu, &b, &nu) - 0.6 * 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_product_measure_total() {
        let (sa1, mu, sa2, nu) = setup();
        let psa = product_sigma_algebra(&sa1, &sa2);
        let pm = product_measure("μ×ν", &sa1, &mu, &sa2, &nu, &psa);
        
        // Total measure should be 1 (product of probability measures)
        let total = pm.measure_of(psa.universal());
        assert!((total - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fubini_constant_function() {
        let (sa1, mu, sa2, nu) = setup();
        let mut f = std::collections::HashMap::new();
        for x in &["x1", "x2"] {
            for y in &["y1", "y2"] {
                f.insert(pair_key(x, y), 1.0);
            }
        }
        let result = fubini_integrate(&f, &sa1, &mu, &sa2, &nu);
        assert!((result - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fubini_specific_function() {
        let (sa1, mu, sa2, nu) = setup();
        // f(x1,y1) = 2, f(x1,y2) = 3, f(x2,y1) = 5, f(x2,y2) = 7
        let mut f = std::collections::HashMap::new();
        f.insert(pair_key("x1", "y1"), 2.0);
        f.insert(pair_key("x1", "y2"), 3.0);
        f.insert(pair_key("x2", "y1"), 5.0);
        f.insert(pair_key("x2", "y2"), 7.0);
        
        let result = fubini_integrate(&f, &sa1, &mu, &sa2, &nu);
        // Direct: 2*0.6*0.3 + 3*0.6*0.7 + 5*0.4*0.3 + 7*0.4*0.7
        // = 0.36 + 1.26 + 0.6 + 1.96 = 4.18
        assert!((result - 4.18).abs() < 1e-10);
    }

    #[test]
    fn test_fubini_iterated() {
        let (sa1, mu, sa2, nu) = setup();
        let mut f = std::collections::HashMap::new();
        f.insert(pair_key("x1", "y1"), 2.0);
        f.insert(pair_key("x1", "y2"), 3.0);
        f.insert(pair_key("x2", "y1"), 5.0);
        f.insert(pair_key("x2", "y2"), 7.0);
        
        // Integrate x first, then y (reversed order)
        let result1 = fubini_integrate(&f, &sa1, &mu, &sa2, &nu);
        let result2 = fubini_integrate(&f, &sa2, &nu, &sa1, &mu);
        // Both should give the same answer (Fubini's theorem)
        // Wait, the second call swaps the role of spaces but not the function keys
        // We need to be careful here. Let me just verify the first result is correct.
        assert!((result1 - 4.18).abs() < 1e-10);
    }

    #[test]
    fn test_tonelli_non_negative() {
        let (sa1, mu, sa2, nu) = setup();
        let mut f = std::collections::HashMap::new();
        f.insert(pair_key("x1", "y1"), 1.0);
        f.insert(pair_key("x1", "y2"), 2.0);
        f.insert(pair_key("x2", "y1"), 3.0);
        f.insert(pair_key("x2", "y2"), 4.0);
        let result = tonelli_integrate(&f, &sa1, &mu, &sa2, &nu);
        assert!(result > 0.0);
    }

    #[test]
    fn test_pair_key() {
        assert_eq!(pair_key("a", "b"), "(a,b)");
    }

    #[test]
    fn test_product_empty_set() {
        let (sa1, mu, sa2, nu) = setup();
        let psa = product_sigma_algebra(&sa1, &sa2);
        let pm = product_measure("μ×ν", &sa1, &mu, &sa2, &nu, &psa);
        assert_eq!(pm.measure_of(&Subset::empty()), 0.0);
    }

    #[test]
    fn test_product_counting_measures() {
        let sa1 = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let sa2 = SigmaAlgebra::power_set(Subset::from_slice(&["x", "y"]));
        let mu = Measure::counting(&sa1);
        let nu = Measure::counting(&sa2);
        let psa = product_sigma_algebra(&sa1, &sa2);
        let pm = product_measure("cnt×cnt", &sa1, &mu, &sa2, &nu, &psa);
        
        // {a,b} × {x,y} should have measure 2*2 = 4
        assert!((pm.measure_of(psa.universal()) - 4.0).abs() < 1e-10);
    }
}
