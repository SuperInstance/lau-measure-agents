//! Pushforward measure: measure induced by a measurable map.
//!
//! Given a measurable map T: (Ω, 𝔉, μ) → (Ω', 𝔊),
//! the pushforward measure T*μ on Ω' is defined by:
//! (T*μ)(B) = μ(T⁻¹(B)) for all B ∈ 𝔊.

use serde::{Serialize, Deserialize};
use crate::sigma_algebra::{Subset, SigmaAlgebra};
use crate::measure::Measure;
use crate::measurable_function::MeasurableFunction;

/// Compute the pushforward measure T*μ.
pub fn pushforward(
    name: &str,
    f: &MeasurableFunction,
    mu: &Measure,
    domain_algebra: &SigmaAlgebra,
    codomain_algebra: &SigmaAlgebra,
) -> Measure {
    let mut values = Vec::new();

    for b in codomain_algebra.measurable_sets() {
        let preimage = f.preimage(b);
        let val = mu.measure_of(&preimage);
        if val > 0.0 || b.is_empty() {
            values.push((b.clone(), val));
        }
    }

    Measure::new(name, codomain_algebra, &values).unwrap_or_else(|_| {
        // Fallback: create measure with just empty set
        Measure::new(name, codomain_algebra, &[]).unwrap()
    })
}

/// Compute the pushforward of a probability measure through a measurable function.
/// This gives the distribution of f(X) when X ~ μ.
pub fn pushforward_probability(
    name: &str,
    f: &MeasurableFunction,
    mu: &Measure,
    domain_algebra: &SigmaAlgebra,
    codomain_algebra: &SigmaAlgebra,
) -> Result<Measure, String> {
    let pf = pushforward(name, f, mu, domain_algebra, codomain_algebra);
    if !pf.is_finite() {
        return Err("Pushforward is not finite".to_string());
    }
    Ok(pf)
}

/// Pullback of a measure through a bijection (inverse of pushforward).
pub fn pullback(
    name: &str,
    f: &MeasurableFunction,
    nu: &Measure,
    domain_algebra: &SigmaAlgebra,
    codomain_algebra: &SigmaAlgebra,
) -> Measure {
    // If f is bijective, pullback is pushforward under f^{-1}
    // We construct inverse mapping
    let mut inv_map: std::collections::HashMap<String, Vec<String>> = std::collections::HashMap::new();
    for (src, dst) in f.mapping().iter() {
        inv_map.entry(dst.clone()).or_default().push(src.clone());
    }

    let mut values = Vec::new();
    for a in domain_algebra.measurable_sets() {
        let img = f.image(a);
        let val = nu.measure_of(&img);
        if val > 0.0 || a.is_empty() {
            values.push((a.clone(), val));
        }
    }

    Measure::new(name, domain_algebra, &values).unwrap_or_else(|_| {
        Measure::new(name, domain_algebra, &[]).unwrap()
    })
}

/// Change of variables formula:
/// ∫ g d(T*μ) = ∫ (g ∘ T) dμ
pub fn change_of_variables(
    g_values: &std::collections::HashMap<String, f64>,
    f: &MeasurableFunction,
    mu: &Measure,
    domain_elements: &[String],
) -> f64 {
    let mut total = 0.0;
    for elem in domain_elements {
        let singleton = Subset::from_slice(&[elem]);
        let mu_val = mu.measure_of(&singleton);
        if let Some(dst) = f.apply(elem) {
            let g_val = g_values.get(dst).copied().unwrap_or(0.0);
            total += g_val * mu_val;
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_pushforward() -> (SigmaAlgebra, SigmaAlgebra, Measure, MeasurableFunction) {
        let domain = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let codomain = SigmaAlgebra::power_set(Subset::from_slice(&["x", "y"]));
        let mu = Measure::probability(&domain, &[
            ("a".into(), 0.5), ("b".into(), 0.3), ("c".into(), 0.2),
        ]).unwrap();
        let f = MeasurableFunction::new(
            "T",
            &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into()), ("c".into(), "y".into())],
        ).unwrap();
        (domain, codomain, mu, f)
    }

    #[test]
    fn test_pushforward_basic() {
        let (domain, codomain, mu, f) = setup_pushforward();
        let pf = pushforward("T*μ", &f, &mu, &domain, &codomain);
        // T*μ({x}) = μ({a,b}) = 0.8
        assert!((pf.measure_of(&Subset::from_slice(&["x"])) - 0.8).abs() < 1e-10);
        // T*μ({y}) = μ({c}) = 0.2
        assert!((pf.measure_of(&Subset::from_slice(&["y"])) - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_pushforward_total() {
        let (domain, codomain, mu, f) = setup_pushforward();
        let pf = pushforward("T*μ", &f, &mu, &domain, &codomain);
        // Total should be 1 (probability measure)
        let total = pf.measure_of(codomain.universal());
        assert!((total - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pushforward_empty() {
        let (domain, codomain, mu, f) = setup_pushforward();
        let pf = pushforward("T*μ", &f, &mu, &domain, &codomain);
        assert_eq!(pf.measure_of(&Subset::empty()), 0.0);
    }

    #[test]
    fn test_change_of_variables() {
        let (domain, _, mu, f) = setup_pushforward();
        let g = std::collections::HashMap::from([
            ("x".into(), 3.0),
            ("y".into(), 5.0),
        ]);
        let elems = vec!["a".into(), "b".into(), "c".into()];
        let result = change_of_variables(&g, &f, &mu, &elems);
        // g(T(a))*μ(a) + g(T(b))*μ(b) + g(T(c))*μ(c)
        // = 3*0.5 + 3*0.3 + 5*0.2 = 1.5 + 0.9 + 1.0 = 3.4
        assert!((result - 3.4).abs() < 1e-10);
    }

    #[test]
    fn test_pushforward_identity() {
        let domain = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let mu = Measure::probability(&domain, &[
            ("a".into(), 0.6), ("b".into(), 0.4),
        ]).unwrap();
        let f = MeasurableFunction::new(
            "id", &domain, &domain,
            &[("a".into(), "a".into()), ("b".into(), "b".into())],
        ).unwrap();
        let pf = pushforward("id*μ", &f, &mu, &domain, &domain);
        assert!((pf.measure_of(&Subset::from_slice(&["a"])) - 0.6).abs() < 1e-10);
        assert!((pf.measure_of(&Subset::from_slice(&["b"])) - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_pushforward_constant() {
        let domain = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b"]));
        let codomain = SigmaAlgebra::power_set(Subset::from_slice(&["x"]));
        let mu = Measure::uniform(&domain).unwrap();
        let f = MeasurableFunction::new(
            "const", &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into())],
        ).unwrap();
        let pf = pushforward("const*μ", &f, &mu, &domain, &codomain);
        assert!((pf.measure_of(&Subset::from_slice(&["x"])) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pushforward_counting() {
        let domain = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let codomain = SigmaAlgebra::power_set(Subset::from_slice(&["x", "y"]));
        let mu = Measure::counting(&domain);
        let f = MeasurableFunction::new(
            "T", &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into()), ("c".into(), "y".into())],
        ).unwrap();
        let pf = pushforward("T*cnt", &f, &mu, &domain, &codomain);
        assert!((pf.measure_of(&Subset::from_slice(&["x"])) - 2.0).abs() < 1e-10);
        assert!((pf.measure_of(&Subset::from_slice(&["y"])) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pushforward_probability() {
        let (domain, codomain, mu, f) = setup_pushforward();
        let pf = pushforward_probability("T*μ", &f, &mu, &domain, &codomain).unwrap();
        assert!(pf.is_finite());
    }
}
