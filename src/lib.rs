//! # lau-measure-agents
//!
//! Measure theory for agents — the rigorous foundation of integration and probability.
//!
//! This crate provides:
//! - **Sigma-algebras**: countable unions/intersections/complements
//! - **Measures**: non-negative countably additive set functions
//! - **Lebesgue measure**: the standard measure on R^n
//! - **Measurable functions**: functions between measurable spaces
//! - **Lebesgue integral**: with monotone/dominated convergence and Fatou's lemma
//! - **Product measure**: Fubini's theorem
//! - **Radon-Nikodym theorem**: densities when ν ≪ μ
//! - **Lebesgue decomposition**: ν = ν_ac + ν_singular
//! - **Riesz representation**: positive linear functionals as integrals
//! - **Pushforward measure**: measure induced by a measurable map
//! - **Agent applications**: observations as measurable functions, belief updates via Radon-Nikodym

pub mod sigma_algebra;
pub mod measure;
pub mod lebesgue;
pub mod measurable_function;
pub mod integral;
pub mod product;
pub mod radon_nikodym;
pub mod decomposition;
pub mod riesz;
pub mod pushforward;
pub mod agent;

pub use sigma_algebra::*;
pub use measure::*;
pub use lebesgue::*;
pub use measurable_function::*;
pub use integral::*;
pub use product::*;
pub use radon_nikodym::*;
pub use decomposition::*;
pub use riesz::*;
pub use pushforward::*;
pub use agent::*;
