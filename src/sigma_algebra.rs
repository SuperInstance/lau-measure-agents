//! Sigma-algebras: the foundational structure of measure theory.
//!
//! A σ-algebra on a set X is a collection of subsets containing ∅ and X,
//! closed under complements and countable unions (hence also countable intersections).

use serde::{Serialize, Deserialize};
use std::collections::{BTreeSet, HashSet};
use std::hash::{Hash, Hasher};

/// A set element represented as a string for generality.
pub type Element = String;

/// A subset of the underlying space, represented as a sorted set of elements.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord)]
pub struct Subset(pub BTreeSet<Element>);

impl Hash for Subset {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for elem in &self.0 {
            elem.hash(state);
        }
    }
}

impl Subset {
    /// Create an empty set.
    pub fn empty() -> Self {
        Subset(BTreeSet::new())
    }

    /// Create from a slice of elements.
    pub fn from_slice(elems: &[&str]) -> Self {
        Subset(elems.iter().map(|s| s.to_string()).collect())
    }

    /// Union of two subsets.
    pub fn union(&self, other: &Self) -> Self {
        Subset(self.0.union(&other.0).cloned().collect())
    }

    /// Intersection of two subsets.
    pub fn intersection(&self, other: &Self) -> Self {
        Subset(self.0.intersection(&other.0).cloned().collect())
    }

    /// Difference: elements in self but not in other.
    pub fn difference(&self, other: &Self) -> Self {
        Subset(self.0.difference(&other.0).cloned().collect())
    }

    /// Complement relative to the universal set.
    pub fn complement(&self, universal: &Self) -> Self {
        universal.difference(self)
    }

    /// Is this the empty set?
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Cardinality.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Does it contain an element?
    pub fn contains(&self, elem: &str) -> bool {
        self.0.contains(elem)
    }

    /// Is self a subset of other?
    pub fn is_subset_of(&self, other: &Self) -> bool {
        self.0.is_subset(&other.0)
    }
}

/// A sigma-algebra on a finite (or countably represented) set.
///
/// Stores the universal set and the collection of measurable subsets.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SigmaAlgebra {
    /// The universal set (the "space" Ω).
    universal: Subset,
    /// The measurable sets (must contain ∅ and Ω, closed under complement and countable union).
    measurable_sets: HashSet<Subset>,
}

impl SigmaAlgebra {
    /// Create the trivial sigma-algebra {∅, Ω}.
    pub fn trivial(universal: Subset) -> Self {
        let mut measurable_sets = HashSet::new();
        measurable_sets.insert(Subset::empty());
        measurable_sets.insert(universal.clone());
        SigmaAlgebra { universal, measurable_sets }
    }

    /// Create the power set sigma-algebra (all subsets are measurable).
    pub fn power_set(universal: Subset) -> Self {
        let elems: Vec<Element> = universal.0.iter().cloned().collect();
        let n = elems.len();
        let mut measurable_sets = HashSet::new();
        // Generate all 2^n subsets
        for mask in 0u64..(1u64 << n) {
            let mut s = BTreeSet::new();
            for i in 0..n {
                if mask & (1u64 << i) != 0 {
                    s.insert(elems[i].clone());
                }
            }
            measurable_sets.insert(Subset(s));
        }
        SigmaAlgebra { universal, measurable_sets }
    }

    /// Generate the sigma-algebra from a collection of generator sets.
    /// This computes the smallest sigma-algebra containing the generators.
    pub fn generate(universal: Subset, generators: &[Subset]) -> Self {
        let mut alg = HashSet::new();
        alg.insert(Subset::empty());
        alg.insert(universal.clone());

        // Add generators
        for g in generators {
            alg.insert(g.clone());
            alg.insert(g.complement(&universal));
        }

        // Close under finite unions, intersections, complements
        // For finite spaces, iterating to fixed point suffices
        let mut changed = true;
        while changed {
            changed = false;
            let current: Vec<Subset> = alg.iter().cloned().collect();
            for i in 0..current.len() {
                for j in i..current.len() {
                    let u = current[i].union(&current[j]);
                    if !alg.contains(&u) {
                        alg.insert(u);
                        changed = true;
                    }
                    let inter = current[i].intersection(&current[j]);
                    if !alg.contains(&inter) {
                        alg.insert(inter);
                        changed = true;
                    }
                }
                let comp = current[i].complement(&universal);
                if !alg.contains(&comp) {
                    alg.insert(comp);
                    changed = true;
                }
            }
        }

        SigmaAlgebra { universal, measurable_sets: alg }
    }

    /// The universal set.
    pub fn universal(&self) -> &Subset {
        &self.universal
    }

    /// Number of measurable sets.
    pub fn size(&self) -> usize {
        self.measurable_sets.len()
    }

    /// Is a given subset measurable?
    pub fn is_measurable(&self, s: &Subset) -> bool {
        self.measurable_sets.contains(s)
    }

    /// All measurable sets.
    pub fn measurable_sets(&self) -> &HashSet<Subset> {
        &self.measurable_sets
    }

    /// Validate sigma-algebra axioms.
    pub fn validate(&self) -> bool {
        // Must contain empty set
        if !self.measurable_sets.contains(&Subset::empty()) {
            return false;
        }
        // Must contain universal set
        if !self.measurable_sets.contains(&self.universal) {
            return false;
        }
        // Closed under complements
        for s in &self.measurable_sets {
            let comp = s.complement(&self.universal);
            if !self.measurable_sets.contains(&comp) {
                return false;
            }
        }
        // Closed under finite unions (implies countable for finite spaces)
        let sets: Vec<&Subset> = self.measurable_sets.iter().collect();
        for i in 0..sets.len() {
            for j in i..sets.len() {
                let u = sets[i].union(sets[j]);
                if !self.measurable_sets.contains(&u) {
                    return false;
                }
            }
        }
        true
    }

    /// Countable union of measurable sets.
    pub fn countable_union(&self, sets: &[Subset]) -> Result<Subset, String> {
        let mut result = Subset::empty();
        for s in sets {
            if !self.is_measurable(s) {
                return Err("Set is not measurable".to_string());
            }
            result = result.union(s);
        }
        Ok(result)
    }

    /// Countable intersection of measurable sets.
    pub fn countable_intersection(&self, sets: &[Subset]) -> Result<Subset, String> {
        if sets.is_empty() {
            return Ok(self.universal.clone());
        }
        let mut result = self.universal.clone();
        for s in sets {
            if !self.is_measurable(s) {
                return Err("Set is not measurable".to_string());
            }
            result = result.intersection(s);
        }
        Ok(result)
    }

    /// Complement of a measurable set.
    pub fn complement(&self, s: &Subset) -> Result<Subset, String> {
        if !self.is_measurable(s) {
            return Err("Set is not measurable".to_string());
        }
        Ok(s.complement(&self.universal))
    }
}

/// A measurable space: a set equipped with a sigma-algebra.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeasurableSpace {
    algebra: SigmaAlgebra,
}

impl MeasurableSpace {
    pub fn new(algebra: SigmaAlgebra) -> Self {
        MeasurableSpace { algebra }
    }

    pub fn sigma_algebra(&self) -> &SigmaAlgebra {
        &self.algebra
    }

    pub fn universal(&self) -> &Subset {
        self.algebra.universal()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_and_universal() {
        let u = Subset::from_slice(&["a", "b", "c"]);
        let sa = SigmaAlgebra::trivial(u.clone());
        assert!(sa.validate());
        assert_eq!(sa.size(), 2);
        assert!(sa.is_measurable(&Subset::empty()));
        assert!(sa.is_measurable(&u));
    }

    #[test]
    fn test_power_set() {
        let u = Subset::from_slice(&["a", "b"]);
        let sa = SigmaAlgebra::power_set(u);
        assert!(sa.validate());
        assert_eq!(sa.size(), 4); // 2^2
    }

    #[test]
    fn test_power_set_three() {
        let u = Subset::from_slice(&["a", "b", "c"]);
        let sa = SigmaAlgebra::power_set(u);
        assert!(sa.validate());
        assert_eq!(sa.size(), 8);
    }

    #[test]
    fn test_generated_algebra() {
        let u = Subset::from_slice(&["a", "b", "c"]);
        let gen = Subset::from_slice(&["a"]);
        let sa = SigmaAlgebra::generate(u.clone(), &[gen]);
        assert!(sa.validate());
        // Should contain: ∅, {a}, {b,c}, {a,b,c}
        assert_eq!(sa.size(), 4);
    }

    #[test]
    fn test_countable_union() {
        let u = Subset::from_slice(&["a", "b", "c"]);
        let sa = SigmaAlgebra::power_set(u);
        let s1 = Subset::from_slice(&["a"]);
        let s2 = Subset::from_slice(&["b"]);
        let result = sa.countable_union(&[s1, s2]).unwrap();
        assert_eq!(result, Subset::from_slice(&["a", "b"]));
    }

    #[test]
    fn test_countable_intersection() {
        let u = Subset::from_slice(&["a", "b", "c"]);
        let sa = SigmaAlgebra::power_set(u);
        let s1 = Subset::from_slice(&["a", "b"]);
        let s2 = Subset::from_slice(&["b", "c"]);
        let result = sa.countable_intersection(&[s1, s2]).unwrap();
        assert_eq!(result, Subset::from_slice(&["b"]));
    }

    #[test]
    fn test_complement() {
        let u = Subset::from_slice(&["a", "b", "c"]);
        let sa = SigmaAlgebra::power_set(u);
        let s = Subset::from_slice(&["a", "b"]);
        let comp = sa.complement(&s).unwrap();
        assert_eq!(comp, Subset::from_slice(&["c"]));
    }

    #[test]
    fn test_subset_operations() {
        let s1 = Subset::from_slice(&["a", "b"]);
        let s2 = Subset::from_slice(&["b", "c"]);
        assert_eq!(s1.union(&s2), Subset::from_slice(&["a", "b", "c"]));
        assert_eq!(s1.intersection(&s2), Subset::from_slice(&["b"]));
        assert_eq!(s1.difference(&s2), Subset::from_slice(&["a"]));
    }

    #[test]
    fn test_subset_contains() {
        let s = Subset::from_slice(&["a", "b"]);
        assert!(s.contains("a"));
        assert!(!s.contains("c"));
    }

    #[test]
    fn test_subset_is_subset_of() {
        let s1 = Subset::from_slice(&["a"]);
        let s2 = Subset::from_slice(&["a", "b"]);
        assert!(s1.is_subset_of(&s2));
        assert!(!s2.is_subset_of(&s1));
    }

    #[test]
    fn test_trivial_algebra_two_elements() {
        let u = Subset::from_slice(&["x", "y"]);
        let sa = SigmaAlgebra::trivial(u);
        assert!(sa.validate());
        assert!(!sa.is_measurable(&Subset::from_slice(&["x"])));
    }

    #[test]
    fn test_measurable_space() {
        let u = Subset::from_slice(&["a", "b"]);
        let sa = SigmaAlgebra::power_set(u);
        let ms = MeasurableSpace::new(sa);
        assert_eq!(ms.universal().len(), 2);
    }
}
