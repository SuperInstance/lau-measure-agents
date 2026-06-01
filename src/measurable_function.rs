//! Measurable functions: functions between measurable spaces.
//!
//! A function f: (Ω, 𝔉) → (Ω', 𝔉') is measurable if f⁻¹(A') ∈ 𝔉 for all A' ∈ 𝔉'.

use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use crate::sigma_algebra::{Subset, SigmaAlgebra};

/// A measurable function between finite measurable spaces.
/// Maps elements of the domain to elements of the codomain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeasurableFunction {
    name: String,
    /// Element-level mapping: domain element → codomain element.
    mapping: HashMap<String, String>,
}

impl MeasurableFunction {
    /// Create a new measurable function.
    /// Validates measurability: preimage of every measurable set in codomain must be measurable in domain.
    pub fn new(
        name: &str,
        domain: &SigmaAlgebra,
        codomain: &SigmaAlgebra,
        mapping: &[(String, String)],
    ) -> Result<Self, String> {
        let map: HashMap<String, String> = mapping.iter().cloned().collect();

        // Check that every domain element is mapped
        for elem in &domain.universal().0 {
            if !map.contains_key(elem) {
                return Err(format!("Domain element '{}' not mapped", elem));
            }
        }

        // Check measurability: for each measurable set B in codomain,
        // f⁻¹(B) must be measurable in domain
        for b in codomain.measurable_sets() {
            let preimage = Self::preimage_of_set(&map, b);
            if !domain.is_measurable(&preimage) {
                return Err(format!(
                    "Not measurable: preimage of {:?} is {:?}, which is not in the domain sigma-algebra",
                    b.0, preimage.0
                ));
            }
        }

        Ok(MeasurableFunction {
            name: name.to_string(),
            mapping: map,
        })
    }

    /// Create without validation (trusted).
    pub fn new_unchecked(name: &str, mapping: &[(String, String)]) -> Self {
        MeasurableFunction {
            name: name.to_string(),
            mapping: mapping.iter().cloned().collect(),
        }
    }

    /// Apply function to an element.
    pub fn apply(&self, elem: &str) -> Option<&str> {
        self.mapping.get(elem).map(|s| s.as_str())
    }

    /// Preimage of a subset of the codomain.
    pub fn preimage(&self, b: &Subset) -> Subset {
        Self::preimage_of_set(&self.mapping, b)
    }

    fn preimage_of_set(mapping: &HashMap<String, String>, b: &Subset) -> Subset {
        let mut result = std::collections::BTreeSet::new();
        for (src, dst) in mapping {
            if b.contains(dst) {
                result.insert(src.clone());
            }
        }
        Subset(result)
    }

    /// Image of a subset of the domain.
    pub fn image(&self, a: &Subset) -> Subset {
        let mut result = std::collections::BTreeSet::new();
        for elem in &a.0 {
            if let Some(dst) = self.mapping.get(elem) {
                result.insert(dst.clone());
            }
        }
        Subset(result)
    }

    /// Compose two measurable functions: g ∘ f.
    pub fn compose(&self, other: &MeasurableFunction) -> Result<MeasurableFunction, String> {
        let mut mapping = HashMap::new();
        for (elem, mid) in &self.mapping {
            if let Some(final_dst) = other.mapping.get(mid) {
                mapping.insert(elem.clone(), final_dst.clone());
            } else {
                return Err(format!("Intermediate value '{}' not in domain of second function", mid));
            }
        }
        Ok(MeasurableFunction {
            name: format!("{}_∘_{}", other.name, self.name),
            mapping,
        })
    }

    /// The name of this function.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get a reference to the internal mapping.
    pub fn mapping(&self) -> &HashMap<String, String> {
        &self.mapping
    }
}

/// An indicator function 1_A for a measurable set A.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndicatorFunction {
    set: Subset,
}

impl IndicatorFunction {
    pub fn new(set: Subset) -> Self {
        IndicatorFunction { set }
    }

    /// Evaluate: 1 if elem ∈ A, 0 otherwise.
    pub fn eval(&self, elem: &str) -> f64 {
        if self.set.contains(elem) { 1.0 } else { 0.0 }
    }

    /// Get the underlying set.
    pub fn set(&self) -> &Subset {
        &self.set
    }
}

/// A simple function: finite linear combination of indicator functions.
/// f = Σ c_i * 1_{A_i}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimpleFunction {
    /// (set, coefficient) pairs; sets should be disjoint for canonical form.
    pub pieces: Vec<(Subset, f64)>,
}

impl SimpleFunction {
    pub fn new(pieces: Vec<(Subset, f64)>) -> Self {
        SimpleFunction { pieces }
    }

    /// Evaluate at an element.
    pub fn eval(&self, elem: &str) -> f64 {
        self.pieces.iter()
            .filter(|(set, _)| set.contains(elem))
            .map(|(_, c)| *c)
            .sum()
    }

    /// Supremum of the function.
    pub fn supremum(&self) -> f64 {
        self.pieces.iter().map(|(_, c)| *c).fold(0.0_f64, f64::max)
    }

    /// Integral w.r.t. a measure.
    pub fn integrate(&self, mu: &crate::measure::Measure) -> f64 {
        self.pieces.iter()
            .map(|(set, c)| c * mu.measure_of(set))
            .sum()
    }

    /// Add two simple functions.
    pub fn add(&self, other: &SimpleFunction) -> SimpleFunction {
        let mut pieces = self.pieces.clone();
        pieces.extend(other.pieces.clone());
        SimpleFunction { pieces }
    }

    /// Scale a simple function.
    pub fn scale(&self, c: f64) -> SimpleFunction {
        SimpleFunction {
            pieces: self.pieces.iter().map(|(s, v)| (s.clone(), v * c)).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sigma_algebra::SigmaAlgebra;

    fn make_spaces() -> (SigmaAlgebra, SigmaAlgebra) {
        let domain = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let codomain = SigmaAlgebra::power_set(Subset::from_slice(&["x", "y"]));
        (domain, codomain)
    }

    #[test]
    fn test_measurable_function_constant() {
        let (domain, codomain) = make_spaces();
        let f = MeasurableFunction::new(
            "const",
            &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into()), ("c".into(), "x".into())],
        ).unwrap();
        assert_eq!(f.apply("a"), Some("x"));
        assert_eq!(f.apply("c"), Some("x"));
    }

    #[test]
    fn test_measurable_function_two_values() {
        let (domain, codomain) = make_spaces();
        let f = MeasurableFunction::new(
            "f",
            &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into()), ("c".into(), "y".into())],
        ).unwrap();
        assert_eq!(f.apply("a"), Some("x"));
        assert_eq!(f.apply("c"), Some("y"));
    }

    #[test]
    fn test_preimage() {
        let (domain, codomain) = make_spaces();
        let f = MeasurableFunction::new(
            "f",
            &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into()), ("c".into(), "y".into())],
        ).unwrap();
        let pre = f.preimage(&Subset::from_slice(&["x"]));
        assert_eq!(pre, Subset::from_slice(&["a", "b"]));
    }

    #[test]
    fn test_image() {
        let (domain, codomain) = make_spaces();
        let f = MeasurableFunction::new(
            "f",
            &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into()), ("c".into(), "y".into())],
        ).unwrap();
        let img = f.image(&Subset::from_slice(&["a", "c"]));
        assert_eq!(img, Subset::from_slice(&["x", "y"]));
    }

    #[test]
    fn test_compose() {
        let f = MeasurableFunction::new_unchecked("f", &[
            ("a".into(), "x".into()), ("b".into(), "y".into()),
        ]);
        let g = MeasurableFunction::new_unchecked("g", &[
            ("x".into(), "α".into()), ("y".into(), "β".into()),
        ]);
        let h = f.compose(&g).unwrap();
        assert_eq!(h.apply("a"), Some("α"));
        assert_eq!(h.apply("b"), Some("β"));
    }

    #[test]
    fn test_indicator_function() {
        let ind = IndicatorFunction::new(Subset::from_slice(&["a", "b"]));
        assert_eq!(ind.eval("a"), 1.0);
        assert_eq!(ind.eval("c"), 0.0);
    }

    #[test]
    fn test_simple_function() {
        let sf = SimpleFunction::new(vec![
            (Subset::from_slice(&["a"]), 1.0),
            (Subset::from_slice(&["b", "c"]), 2.0),
        ]);
        assert_eq!(sf.eval("a"), 1.0);
        assert_eq!(sf.eval("b"), 2.0);
        assert_eq!(sf.eval("c"), 2.0);
    }

    #[test]
    fn test_simple_function_integrate() {
        let sa = SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]));
        let mu = crate::measure::Measure::uniform(&sa).unwrap();
        let sf = SimpleFunction::new(vec![
            (Subset::from_slice(&["a"]), 3.0),
            (Subset::from_slice(&["b", "c"]), 6.0),
        ]);
        // 3 * (1/3) + 6 * (2/3) = 1 + 4 = 5
        let integral = sf.integrate(&mu);
        assert!((integral - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_simple_function_supremum() {
        let sf = SimpleFunction::new(vec![
            (Subset::from_slice(&["a"]), 1.0),
            (Subset::from_slice(&["b"]), 5.0),
        ]);
        assert_eq!(sf.supremum(), 5.0);
    }

    #[test]
    fn test_simple_function_add() {
        let f = SimpleFunction::new(vec![(Subset::from_slice(&["a"]), 1.0)]);
        let g = SimpleFunction::new(vec![(Subset::from_slice(&["a"]), 2.0)]);
        let h = f.add(&g);
        assert_eq!(h.eval("a"), 3.0);
    }

    #[test]
    fn test_simple_function_scale() {
        let f = SimpleFunction::new(vec![(Subset::from_slice(&["a"]), 2.0)]);
        let g = f.scale(3.0);
        assert_eq!(g.eval("a"), 6.0);
    }

    #[test]
    fn test_unmapped_element_rejected() {
        let (domain, codomain) = make_spaces();
        let result = MeasurableFunction::new(
            "bad",
            &domain, &codomain,
            &[("a".into(), "x".into()), ("b".into(), "x".into())], // missing "c"
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_function_name() {
        let f = MeasurableFunction::new_unchecked("my_func", &[]);
        assert_eq!(f.name(), "my_func");
    }
}
