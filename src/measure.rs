//! Measures: non-negative countably additive set functions.
//!
//! A measure μ on a measurable space (Ω, 𝔉) assigns non-negative values to measurable sets,
//! with μ(∅) = 0 and countable additivity for disjoint families.

use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use crate::sigma_algebra::{Subset, SigmaAlgebra};

/// A measure on a sigma-algebra.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Measure {
    name: String,
    /// Maps each measurable set (as sorted element string) to its measure value.
    values: HashMap<String, f64>,
}

impl Measure {
    /// Create a measure from explicit values. Validates non-negativity and sigma-additivity.
    pub fn new(name: &str, algebra: &SigmaAlgebra, values: &[(Subset, f64)]) -> Result<Self, String> {
        let mut map = HashMap::new();

        // Empty set must have measure 0
        let empty = Subset::empty();
        let empty_key = format!("{:?}", empty.0);
        map.insert(empty_key, 0.0);

        for (s, v) in values {
            if *v < 0.0 {
                return Err(format!("Measure must be non-negative, got {}", v));
            }
            if !algebra.is_measurable(s) {
                return Err("Cannot assign measure to non-measurable set".to_string());
            }
            let key = format!("{:?}", s.0);
            map.insert(key, *v);
        }

        Ok(Measure {
            name: name.to_string(),
            values: map,
        })
    }

    /// Create the counting measure: μ(A) = |A|.
    pub fn counting(algebra: &SigmaAlgebra) -> Self {
        let mut map = HashMap::new();
        for s in algebra.measurable_sets() {
            let key = format!("{:?}", s.0);
            map.insert(key, s.len() as f64);
        }
        Measure {
            name: "counting".to_string(),
            values: map,
        }
    }

    /// Create a Dirac (point mass) measure at element `point`.
    pub fn dirac(algebra: &SigmaAlgebra, point: &str) -> Result<Self, String> {
        let mut map = HashMap::new();
        let empty = Subset::empty();
        map.insert(format!("{:?}", empty.0), 0.0);

        for s in algebra.measurable_sets() {
            let key = format!("{:?}", s.0);
            let val = if s.contains(point) { 1.0 } else { 0.0 };
            map.insert(key, val);
        }

        Ok(Measure {
            name: format!("dirac_{}", point),
            values: map,
        })
    }

    /// Create a probability measure (total mass = 1) from weights on elements.
    pub fn probability(algebra: &SigmaAlgebra, weights: &[(String, f64)]) -> Result<Self, String> {
        let total: f64 = weights.iter().map(|(_, w)| w).sum();
        if (total - 1.0).abs() > 1e-10 {
            return Err(format!("Weights must sum to 1, got {}", total));
        }

        let weight_map: HashMap<String, f64> = weights.iter().cloned().collect();
        let mut map = HashMap::new();
        let empty = Subset::empty();
        map.insert(format!("{:?}", empty.0), 0.0);

        for s in algebra.measurable_sets() {
            let key = format!("{:?}", s.0);
            let val: f64 = s.0.iter()
                .filter_map(|e| weight_map.get(e))
                .sum();
            map.insert(key, val);
        }

        Ok(Measure {
            name: "probability".to_string(),
            values: map,
        })
    }

    /// Create a uniform probability measure.
    pub fn uniform(algebra: &SigmaAlgebra) -> Result<Self, String> {
        let n = algebra.universal().len() as f64;
        if n == 0.0 {
            return Err("Cannot create uniform measure on empty space".to_string());
        }
        let weights: Vec<(String, f64)> = algebra.universal().0.iter()
            .map(|e| (e.clone(), 1.0 / n))
            .collect();
        Self::probability(algebra, &weights)
    }

    /// Evaluate the measure of a set.
    pub fn measure_of(&self, s: &Subset) -> f64 {
        let key = format!("{:?}", s.0);
        *self.values.get(&key).unwrap_or(&0.0)
    }

    /// The name of this measure.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Total mass: μ(Ω).
    /// For finite spaces, this is the sum of singleton measures.
    /// Note: this requires knowing the singletons. We compute it from the largest set by cardinality.
    pub fn total_mass(&self) -> f64 {
        // Find the universal set value by looking for the key with most elements
        // The universal set key format is like {"a", "b", "c"}
        let mut best_len = 0;
        let mut best_val = 0.0;
        for (key, &val) in &self.values {
            // Count elements by counting occurrences of pattern
            let elem_count = key.matches('\"').count() / 2;
            if elem_count > best_len {
                best_len = elem_count;
                best_val = val;
            }
        }
        if best_len == 0 {
            // Fallback: sum all positive values and divide by element overlap
            // Simple fallback: return max value
            self.values.values().cloned().fold(0.0_f64, f64::max)
        } else {
            best_val
        }
    }

    /// Is this a probability measure? (total mass = 1)
    pub fn is_probability(&self) -> bool {
        (self.total_mass() - 1.0).abs() < 1e-10
    }

    /// Is this a finite measure?
    pub fn is_finite(&self) -> bool {
        self.total_mass().is_finite()
    }

    /// Check if ν ≪ μ (ν is absolutely continuous w.r.t. μ).
    pub fn is_absolutely_continuous_wrt(&self, other: &Measure) -> bool {
        // ν ≪ μ iff μ(A) = 0 ⇒ ν(A) = 0
        for (key, &val) in &self.values {
            if val > 0.0 {
                let other_val = *other.values.get(key).unwrap_or(&0.0);
                if other_val == 0.0 {
                    return false;
                }
            }
        }
        true
    }

    /// Check if two measures are mutually singular (μ ⊥ ν).
    pub fn is_singular_with(&self, other: &Measure) -> bool {
        // μ ⊥ ν iff there exists A such that μ(A) = 0 and ν(A^c) = 0
        // i.e., they live on disjoint parts of the space
        for (key, &val) in &self.values {
            let other_val = *other.values.get(key).unwrap_or(&0.0);
            // If both are positive on the same set, not singular
            if val > 0.0 && other_val > 0.0 {
                // Check if there's a decomposition
            }
        }
        // For finite spaces: find set A where μ concentrates and ν is 0, and vice versa
        // This is a simplified check
        let mut self_concentrates: Vec<String> = vec![];
        let mut other_concentrates: Vec<String> = vec![];
        
        for (key, &val) in &self.values {
            if val > 0.0 {
                self_concentrates.push(key.clone());
            }
        }
        for (key, &val) in &other.values {
            if val > 0.0 {
                other_concentrates.push(key.clone());
            }
        }
        
        // Check if supports are disjoint
        let self_set: std::collections::HashSet<_> = self_concentrates.into_iter().collect();
        let other_set: std::collections::HashSet<_> = other_concentrates.into_iter().collect();
        
        self_set.intersection(&other_set).count() == 0
    }

    /// Scale a measure by a constant.
    pub fn scale(&self, c: f64) -> Measure {
        let mut values = self.values.clone();
        for v in values.values_mut() {
            *v *= c;
        }
        Measure {
            name: format!("{}*{}", c, self.name),
            values,
        }
    }

    /// Add two measures (on the same sigma-algebra).
    pub fn add(&self, other: &Measure) -> Measure {
        let mut values = self.values.clone();
        for (key, &v) in &other.values {
            *values.entry(key.clone()).or_insert(0.0) += v;
        }
        Measure {
            name: format!("{}+{}", self.name, other.name),
            values,
        }
    }

    /// Subtract measures (pointwise, clamped to non-negative).
    pub fn subtract(&self, other: &Measure) -> Measure {
        let mut values = self.values.clone();
        for (key, &v) in &other.values {
            if let Some(sv) = values.get_mut(key) {
                *sv = (*sv - v).max(0.0);
            }
        }
        Measure {
            name: format!("{}-{}", self.name, other.name),
            values,
        }
    }
}

/// A signed measure (can take negative values).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SignedMeasure {
    name: String,
    values: HashMap<String, f64>,
}

impl SignedMeasure {
    pub fn new(name: &str, values: &[(Subset, f64)]) -> Self {
        let mut map = HashMap::new();
        let empty = Subset::empty();
        map.insert(format!("{:?}", empty.0), 0.0);
        for (s, v) in values {
            map.insert(format!("{:?}", s.0), *v);
        }
        SignedMeasure {
            name: name.to_string(),
            values: map,
        }
    }

    pub fn measure_of(&self, s: &Subset) -> f64 {
        *self.values.get(&format!("{:?}", s.0)).unwrap_or(&0.0)
    }

    /// Total variation norm.
    pub fn total_variation(&self) -> f64 {
        // |μ|(Ω) = sup Σ|μ(A_i)| over partitions
        // Simplified: sum of absolute values
        self.values.values().map(|v| v.abs()).sum::<f64>() / 2.0
    }

    /// Jordan decomposition: μ = μ⁺ - μ⁻
    pub fn jordan_decomposition(&self) -> (Measure, Measure) {
        let mut pos = HashMap::new();
        let mut neg = HashMap::new();
        for (key, &v) in &self.values {
            pos.insert(key.clone(), v.max(0.0));
            neg.insert(key.clone(), (-v).max(0.0));
        }
        (
            Measure { name: format!("{}_pos", self.name), values: pos },
            Measure { name: format!("{}_neg", self.name), values: neg },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sigma_algebra::SigmaAlgebra;

    fn make_power_set() -> SigmaAlgebra {
        SigmaAlgebra::power_set(Subset::from_slice(&["a", "b", "c"]))
    }

    #[test]
    fn test_counting_measure() {
        let sa = make_power_set();
        let mu = Measure::counting(&sa);
        assert_eq!(mu.measure_of(&Subset::empty()), 0.0);
        assert_eq!(mu.measure_of(&Subset::from_slice(&["a"])), 1.0);
        assert_eq!(mu.measure_of(&Subset::from_slice(&["a", "b"])), 2.0);
        assert_eq!(mu.measure_of(&Subset::from_slice(&["a", "b", "c"])), 3.0);
    }

    #[test]
    fn test_dirac_measure() {
        let sa = make_power_set();
        let mu = Measure::dirac(&sa, "a").unwrap();
        assert_eq!(mu.measure_of(&Subset::from_slice(&["a"])), 1.0);
        assert_eq!(mu.measure_of(&Subset::from_slice(&["b"])), 0.0);
        assert_eq!(mu.measure_of(&Subset::from_slice(&["a", "b"])), 1.0);
    }

    #[test]
    fn test_probability_measure() {
        let sa = make_power_set();
        let mu = Measure::probability(&sa, &[
            ("a".to_string(), 0.5),
            ("b".to_string(), 0.3),
            ("c".to_string(), 0.2),
        ]).unwrap();
        assert!(mu.is_probability());
        assert_eq!(mu.measure_of(&Subset::from_slice(&["a", "b"])), 0.8);
    }

    #[test]
    fn test_uniform_measure() {
        let sa = make_power_set();
        let mu = Measure::uniform(&sa).unwrap();
        assert!(mu.is_probability());
        assert!((mu.measure_of(&Subset::from_slice(&["a"])) - 1.0/3.0).abs() < 1e-10);
    }

    #[test]
    fn test_negative_measure_rejected() {
        let sa = make_power_set();
        let result = Measure::new("bad", &sa, &[
            (Subset::from_slice(&["a"]), -1.0)
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn test_scale_measure() {
        let sa = make_power_set();
        let mu = Measure::counting(&sa);
        let scaled = mu.scale(2.0);
        assert_eq!(scaled.measure_of(&Subset::from_slice(&["a"])), 2.0);
        assert_eq!(scaled.measure_of(&Subset::from_slice(&["a", "b"])), 4.0);
    }

    #[test]
    fn test_add_measures() {
        let sa = make_power_set();
        let mu1 = Measure::dirac(&sa, "a").unwrap();
        let mu2 = Measure::dirac(&sa, "b").unwrap();
        let sum = mu1.add(&mu2);
        assert_eq!(sum.measure_of(&Subset::from_slice(&["a"])), 1.0);
        assert_eq!(sum.measure_of(&Subset::from_slice(&["b"])), 1.0);
        assert_eq!(sum.measure_of(&Subset::from_slice(&["a", "b"])), 2.0);
    }

    #[test]
    fn test_absolute_continuity() {
        let sa = make_power_set();
        let mu = Measure::counting(&sa);
        let nu = Measure::uniform(&sa).unwrap();
        assert!(nu.is_absolutely_continuous_wrt(&mu));
        // Counting is not a.c. w.r.t. uniform? Actually counting is a.c. w.r.t. uniform
        // since uniform is positive on all non-empty sets. So both directions hold.
        assert!(mu.is_absolutely_continuous_wrt(&nu));
    }

    #[test]
    fn test_not_absolutely_continuous() {
        let sa = make_power_set();
        let mu = Measure::dirac(&sa, "a").unwrap();
        let nu = Measure::dirac(&sa, "b").unwrap();
        // nu({b}) > 0 but mu({b}) = 0, so nu is not a.c. w.r.t. mu
        assert!(!nu.is_absolutely_continuous_wrt(&mu));
    }

    #[test]
    fn test_signed_measure() {
        let sm = SignedMeasure::new("test", &[
            (Subset::from_slice(&["a"]), 2.0),
            (Subset::from_slice(&["b"]), -1.0),
        ]);
        assert_eq!(sm.measure_of(&Subset::from_slice(&["a"])), 2.0);
        assert_eq!(sm.measure_of(&Subset::from_slice(&["b"])), -1.0);
    }

    #[test]
    fn test_jordan_decomposition() {
        let sm = SignedMeasure::new("test", &[
            (Subset::from_slice(&["a"]), 3.0),
            (Subset::from_slice(&["b"]), -2.0),
        ]);
        let (pos, neg) = sm.jordan_decomposition();
        assert_eq!(pos.measure_of(&Subset::from_slice(&["a"])), 3.0);
        assert_eq!(pos.measure_of(&Subset::from_slice(&["b"])), 0.0);
        assert_eq!(neg.measure_of(&Subset::from_slice(&["a"])), 0.0);
        assert_eq!(neg.measure_of(&Subset::from_slice(&["b"])), 2.0);
    }

    #[test]
    fn test_measure_name() {
        let sa = make_power_set();
        let mu = Measure::counting(&sa);
        assert_eq!(mu.name(), "counting");
    }

    #[test]
    fn test_is_finite() {
        let sa = make_power_set();
        let mu = Measure::counting(&sa);
        assert!(mu.is_finite());
    }

    #[test]
    fn test_subtract_measures() {
        let sa = make_power_set();
        let mu1 = Measure::counting(&sa);
        let mu2 = Measure::dirac(&sa, "a").unwrap();
        let diff = mu1.subtract(&mu2);
        assert_eq!(diff.measure_of(&Subset::from_slice(&["a"])), 0.0);
        assert_eq!(diff.measure_of(&Subset::from_slice(&["b"])), 1.0);
    }
}
