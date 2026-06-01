//! Lebesgue measure on R^n (discrete approximation).
//!
//! The Lebesgue measure is the standard "volume" measure on R^n.
//! For computational purposes, we work with intervals/boxes and their measures.

use serde::{Serialize, Deserialize};
use nalgebra::DVector;

/// An interval [a, b] on the real line.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Interval {
    pub left: f64,
    pub right: f64,
}

impl Interval {
    pub fn new(left: f64, right: f64) -> Result<Self, String> {
        if left > right {
            return Err(format!("left ({}) > right ({})", left, right));
        }
        Ok(Interval { left, right })
    }

    /// Lebesgue measure of this interval: b - a.
    pub fn length(&self) -> f64 {
        self.right - self.left
    }

    /// Is point in interval?
    pub fn contains(&self, x: f64) -> bool {
        x >= self.left && x <= self.right
    }

    /// Intersection of two intervals.
    pub fn intersection(&self, other: &Interval) -> Option<Interval> {
        let left = self.left.max(other.left);
        let right = self.right.min(other.right);
        if left <= right {
            Some(Interval { left, right })
        } else {
            None
        }
    }

    /// Union (if overlapping or adjacent).
    pub fn union(&self, other: &Interval) -> Option<Interval> {
        if self.right < other.left || other.right < self.left {
            return None; // disjoint
        }
        Some(Interval {
            left: self.left.min(other.left),
            right: self.right.max(other.right),
        })
    }

    /// Complement in [min_val, max_val].
    pub fn complement(&self, min_val: f64, max_val: f64) -> Vec<Interval> {
        let mut result = Vec::new();
        if self.left > min_val {
            result.push(Interval { left: min_val, right: self.left });
        }
        if self.right < max_val {
            result.push(Interval { left: self.right, right: max_val });
        }
        result
    }
}

/// A box in R^n (product of intervals).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Box {
    pub intervals: Vec<Interval>,
}

impl Box {
    pub fn new(intervals: Vec<Interval>) -> Self {
        Box { intervals }
    }

    /// Lebesgue measure (volume) of this box.
    pub fn volume(&self) -> f64 {
        self.intervals.iter().map(|i| i.length()).product()
    }

    /// Dimension.
    pub fn dimension(&self) -> usize {
        self.intervals.len()
    }

    /// Does this box contain a point?
    pub fn contains(&self, point: &DVector<f64>) -> bool {
        if point.len() != self.intervals.len() {
            return false;
        }
        self.intervals.iter().enumerate().all(|(i, interval)| {
            interval.contains(point[i])
        })
    }

    /// Intersection of two boxes.
    pub fn intersection(&self, other: &Box) -> Option<Box> {
        if self.intervals.len() != other.intervals.len() {
            return None;
        }
        let intervals: Vec<Interval> = self.intervals.iter()
            .zip(&other.intervals)
            .filter_map(|(a, b)| a.intersection(b))
            .collect();
        if intervals.len() == self.intervals.len() {
            Some(Box::new(intervals))
        } else {
            None
        }
    }
}

/// Lebesgue measure on R (using interval representation).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LebesgueMeasure {
    dimension: usize,
}

impl LebesgueMeasure {
    pub fn new(dimension: usize) -> Self {
        LebesgueMeasure { dimension }
    }

    /// Measure of a box.
    pub fn measure_box(&self, b: &Box) -> f64 {
        if b.dimension() != self.dimension {
            return 0.0;
        }
        b.volume()
    }

    /// Measure of an interval (1D case).
    pub fn measure_interval(&self, i: &Interval) -> f64 {
        if self.dimension != 1 {
            panic!("Use measure_box for higher dimensions");
        }
        i.length()
    }

    /// Measure of a countable union of disjoint intervals.
    pub fn measure_disjoint_intervals(&self, intervals: &[Interval]) -> f64 {
        intervals.iter().map(|i| i.length()).sum()
    }

    /// Outer measure of a set of points (approximated by covering with intervals).
    pub fn outer_measure_approx(&self, points: &[f64], epsilon: f64) -> f64 {
        if points.is_empty() {
            return 0.0;
        }
        let mut sorted = points.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let mut total = 0.0;
        let mut i = 0;
        while i < sorted.len() {
            let mut j = i + 1;
            while j < sorted.len() && sorted[j] - sorted[j - 1] < epsilon {
                j += 1;
            }
            total += (sorted[j - 1] - sorted[i]) + epsilon;
            i = j;
        }
        total
    }

    /// Dimension.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
}

/// Simple function: a finite linear combination of indicator functions.
/// f = Σ a_i * 1_{A_i} where A_i are disjoint measurable sets (intervals).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimpleIntervalFunction {
    /// (interval, value) pairs
    pub pieces: Vec<(Interval, f64)>,
}

impl SimpleIntervalFunction {
    pub fn new(pieces: Vec<(Interval, f64)>) -> Self {
        SimpleIntervalFunction { pieces }
    }

    /// Evaluate at a point.
    pub fn eval(&self, x: f64) -> f64 {
        for (interval, val) in &self.pieces {
            if interval.contains(x) {
                return *val;
            }
        }
        0.0
    }

    /// Lebesgue integral of this simple function w.r.t. Lebesgue measure.
    pub fn integrate(&self) -> f64 {
        self.pieces.iter()
            .map(|(interval, val)| interval.length() * val)
            .sum()
    }

    /// Supremum.
    pub fn supremum(&self) -> f64 {
        self.pieces.iter().map(|(_, v)| *v).fold(0.0_f64, f64::max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_length() {
        let i = Interval::new(0.0, 1.0).unwrap();
        assert_eq!(i.length(), 1.0);
    }

    #[test]
    fn test_interval_contains() {
        let i = Interval::new(0.0, 1.0).unwrap();
        assert!(i.contains(0.5));
        assert!(!i.contains(1.5));
    }

    #[test]
    fn test_interval_intersection() {
        let a = Interval::new(0.0, 2.0).unwrap();
        let b = Interval::new(1.0, 3.0).unwrap();
        let c = a.intersection(&b).unwrap();
        assert_eq!(c.left, 1.0);
        assert_eq!(c.right, 2.0);
    }

    #[test]
    fn test_disjoint_intersection() {
        let a = Interval::new(0.0, 1.0).unwrap();
        let b = Interval::new(2.0, 3.0).unwrap();
        assert!(a.intersection(&b).is_none());
    }

    #[test]
    fn test_interval_union() {
        let a = Interval::new(0.0, 2.0).unwrap();
        let b = Interval::new(1.0, 3.0).unwrap();
        let c = a.union(&b).unwrap();
        assert_eq!(c.left, 0.0);
        assert_eq!(c.right, 3.0);
    }

    #[test]
    fn test_interval_complement() {
        let i = Interval::new(1.0, 2.0).unwrap();
        let comp = i.complement(0.0, 3.0);
        assert_eq!(comp.len(), 2);
        assert_eq!(comp[0], Interval::new(0.0, 1.0).unwrap());
        assert_eq!(comp[1], Interval::new(2.0, 3.0).unwrap());
    }

    #[test]
    fn test_box_volume() {
        let b = Box::new(vec![
            Interval::new(0.0, 2.0).unwrap(),
            Interval::new(0.0, 3.0).unwrap(),
        ]);
        assert_eq!(b.volume(), 6.0);
    }

    #[test]
    fn test_box_contains() {
        let b = Box::new(vec![
            Interval::new(0.0, 1.0).unwrap(),
            Interval::new(0.0, 1.0).unwrap(),
        ]);
        assert!(b.contains(&DVector::from_vec(vec![0.5, 0.5])));
        assert!(!b.contains(&DVector::from_vec(vec![1.5, 0.5])));
    }

    #[test]
    fn test_lebesgue_1d() {
        let m = LebesgueMeasure::new(1);
        let i = Interval::new(0.0, 5.0).unwrap();
        assert_eq!(m.measure_interval(&i), 5.0);
    }

    #[test]
    fn test_lebesgue_2d() {
        let m = LebesgueMeasure::new(2);
        let b = Box::new(vec![
            Interval::new(0.0, 3.0).unwrap(),
            Interval::new(0.0, 4.0).unwrap(),
        ]);
        assert_eq!(m.measure_box(&b), 12.0);
    }

    #[test]
    fn test_disjoint_intervals_measure() {
        let m = LebesgueMeasure::new(1);
        let intervals = vec![
            Interval::new(0.0, 1.0).unwrap(),
            Interval::new(2.0, 3.0).unwrap(),
            Interval::new(5.0, 7.0).unwrap(),
        ];
        assert_eq!(m.measure_disjoint_intervals(&intervals), 4.0);
    }

    #[test]
    fn test_simple_function_integrate() {
        let f = SimpleIntervalFunction::new(vec![
            (Interval::new(0.0, 1.0).unwrap(), 2.0),
            (Interval::new(1.0, 3.0).unwrap(), 1.0),
        ]);
        // 2*1 + 1*2 = 4
        assert_eq!(f.integrate(), 4.0);
    }

    #[test]
    fn test_simple_function_eval() {
        let f = SimpleIntervalFunction::new(vec![
            (Interval::new(0.0, 1.0).unwrap(), 5.0),
            (Interval::new(1.0, 2.0).unwrap(), 10.0),
        ]);
        assert_eq!(f.eval(0.5), 5.0);
        assert_eq!(f.eval(1.5), 10.0);
        assert_eq!(f.eval(3.0), 0.0);
    }

    #[test]
    fn test_simple_function_supremum() {
        let f = SimpleIntervalFunction::new(vec![
            (Interval::new(0.0, 1.0).unwrap(), 3.0),
            (Interval::new(1.0, 2.0).unwrap(), 7.0),
        ]);
        assert_eq!(f.supremum(), 7.0);
    }

    #[test]
    fn test_box_intersection() {
        let b1 = Box::new(vec![
            Interval::new(0.0, 2.0).unwrap(),
            Interval::new(0.0, 2.0).unwrap(),
        ]);
        let b2 = Box::new(vec![
            Interval::new(1.0, 3.0).unwrap(),
            Interval::new(1.0, 3.0).unwrap(),
        ]);
        let bi = b1.intersection(&b2).unwrap();
        assert_eq!(bi.volume(), 1.0);
    }

    #[test]
    fn test_invalid_interval() {
        assert!(Interval::new(2.0, 1.0).is_err());
    }

    #[test]
    fn test_lebesgue_dimension() {
        let m = LebesgueMeasure::new(3);
        assert_eq!(m.dimension(), 3);
    }
}
