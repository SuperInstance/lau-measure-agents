//! Agent applications: measure-theoretic foundations for agent systems.
//!
//! - Agent observations as measurable functions
//! - Agent state space as a measurable space
//! - Belief states as probability measures
//! - Belief updates as Radon-Nikodym derivatives
//! - Expected utility as Lebesgue integral

use serde::{Serialize, Deserialize};
use crate::sigma_algebra::{Subset, SigmaAlgebra, MeasurableSpace};
use crate::measure::Measure;
use crate::measurable_function::MeasurableFunction;
use crate::integral::{RealValuedFunction, integrate};
use crate::radon_nikodym::{RadonNikodymDerivative, BeliefUpdate};
use crate::pushforward::pushforward;
use crate::lebesgue::Interval;

/// An agent's state space (a measurable space).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentStateSpace {
    name: String,
    space: MeasurableSpace,
    elements: Vec<String>,
}

impl AgentStateSpace {
    /// Create a discrete agent state space.
    pub fn discrete(name: &str, elements: &[&str]) -> Self {
        let universal = Subset::from_slice(elements);
        let sa = SigmaAlgebra::power_set(universal);
        AgentStateSpace {
            name: name.to_string(),
            space: MeasurableSpace::new(sa),
            elements: elements.iter().map(|s| s.to_string()).collect(),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn measurable_space(&self) -> &MeasurableSpace {
        &self.space
    }

    pub fn sigma_algebra(&self) -> &SigmaAlgebra {
        self.space.sigma_algebra()
    }

    pub fn elements(&self) -> &[String] {
        &self.elements
    }

    pub fn universal(&self) -> &Subset {
        self.space.universal()
    }
}

/// An agent's belief state (a probability measure on the state space).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BeliefState {
    state_space: String,
    measure: Measure,
}

impl BeliefState {
    /// Create a uniform belief over all states.
    pub fn uniform(space: &AgentStateSpace) -> Self {
        let mu = Measure::uniform(space.sigma_algebra()).unwrap();
        BeliefState {
            state_space: space.name().to_string(),
            measure: mu,
        }
    }

    /// Create a belief from explicit weights.
    pub fn from_weights(space: &AgentStateSpace, weights: &[(String, f64)]) -> Result<Self, String> {
        let mu = Measure::probability(space.sigma_algebra(), weights)?;
        Ok(BeliefState {
            state_space: space.name().to_string(),
            measure: mu,
        })
    }

    /// Create a point belief (Dirac at a single state).
    pub fn point_mass(space: &AgentStateSpace, point: &str) -> Result<Self, String> {
        let mu = Measure::dirac(space.sigma_algebra(), point)?;
        Ok(BeliefState {
            state_space: space.name().to_string(),
            measure: mu,
        })
    }

    /// Probability of a particular state.
    pub fn probability_of(&self, state: &str) -> f64 {
        let singleton = Subset::from_slice(&[state]);
        self.measure.measure_of(&singleton)
    }

    /// Probability of a set of states.
    pub fn probability_of_set(&self, set: &Subset) -> f64 {
        self.measure.measure_of(set)
    }

    /// Update belief given an observation (likelihood function).
    pub fn update(
        &self,
        likelihood: &RealValuedFunction,
        space: &AgentStateSpace,
    ) -> Result<BeliefState, String> {
        let update = BeliefUpdate::from_likelihood(
            likelihood,
            &self.measure,
            &space.elements,
        )?;
        
        // Extract the posterior measure from the RN derivative
        let mut weights = Vec::new();
        for elem in &space.elements {
            let mu_singleton = self.measure.measure_of(&Subset::from_slice(&[elem]));
            let rn_val = update.rn_derivative.eval(elem);
            weights.push((elem.clone(), mu_singleton * rn_val));
        }

        let posterior = Measure::probability(space.sigma_algebra(), &weights)?;
        Ok(BeliefState {
            state_space: self.state_space.clone(),
            measure: posterior,
        })
    }

    /// Get the underlying measure.
    pub fn measure(&self) -> &Measure {
        &self.measure
    }
}

/// An observation model: a measurable function from states to observations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObservationModel {
    name: String,
    /// Maps state → observation.
    mapping: std::collections::HashMap<String, String>,
}

impl ObservationModel {
    pub fn new(name: &str, mappings: &[(String, String)]) -> Self {
        ObservationModel {
            name: name.to_string(),
            mapping: mappings.iter().cloned().collect(),
        }
    }

    /// Observe: what observation does a state produce?
    pub fn observe(&self, state: &str) -> Option<&str> {
        self.mapping.get(state).map(|s| s.as_str())
    }

    /// Get the likelihood function given an observation.
    /// P(observation | state) = 1 if observation matches, 0 otherwise (deterministic).
    pub fn likelihood(&self, observation: &str, space: &AgentStateSpace) -> RealValuedFunction {
        let values: Vec<(&str, f64)> = space.elements.iter()
            .map(|s| {
                let obs = self.mapping.get(s).map(|s| s.as_str()).unwrap_or("");
                (s.as_str(), if obs == observation { 1.0 } else { 0.0 })
            })
            .collect();
        RealValuedFunction::new(values)
    }

    /// Convert to a measurable function.
    pub fn as_measurable_function(
        &self,
        state_space: &AgentStateSpace,
        observation_space: &AgentStateSpace,
    ) -> Result<MeasurableFunction, String> {
        let mapping_refs: Vec<(String, String)> = self.mapping.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        MeasurableFunction::new(
            &self.name,
            state_space.sigma_algebra(),
            observation_space.sigma_algebra(),
            &mapping_refs,
        )
    }
}

/// A reward/utility function on states.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UtilityFunction {
    values: std::collections::HashMap<String, f64>,
}

impl UtilityFunction {
    pub fn new(values: Vec<(&str, f64)>) -> Self {
        UtilityFunction {
            values: values.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
        }
    }

    pub fn utility(&self, state: &str) -> f64 {
        *self.values.get(state).unwrap_or(&0.0)
    }

    /// Expected utility under a belief state.
    pub fn expected_utility(&self, belief: &BeliefState, space: &AgentStateSpace) -> f64 {
        let f = RealValuedFunction {
            values: self.values.clone(),
        };
        integrate(&f, belief.measure())
    }
}

/// A sequential agent that maintains beliefs and makes decisions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Agent {
    name: String,
    state_space: String,
    observation_model: ObservationModel,
    utility: UtilityFunction,
}

impl Agent {
    pub fn new(
        name: &str,
        observation_model: ObservationModel,
        utility: UtilityFunction,
    ) -> Self {
        Agent {
            name: name.to_string(),
            state_space: String::new(),
            observation_model,
            utility,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// Compute expected utility given a belief.
    pub fn expected_utility(&self, belief: &BeliefState, space: &AgentStateSpace) -> f64 {
        self.utility.expected_utility(belief, space)
    }

    /// Get observation model.
    pub fn observation_model(&self) -> &ObservationModel {
        &self.observation_model
    }

    /// Get utility function.
    pub fn utility_function(&self) -> &UtilityFunction {
        &self.utility
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn world_space() -> AgentStateSpace {
        AgentStateSpace::discrete("world", &["sunny", "cloudy", "rainy"])
    }

    fn obs_space() -> AgentStateSpace {
        AgentStateSpace::discrete("obs", &["bright", "dim"])
    }

    #[test]
    fn test_state_space_creation() {
        let space = world_space();
        assert_eq!(space.elements().len(), 3);
        assert_eq!(space.name(), "world");
    }

    #[test]
    fn test_uniform_belief() {
        let space = world_space();
        let belief = BeliefState::uniform(&space);
        assert!((belief.probability_of("sunny") - 1.0/3.0).abs() < 1e-10);
        assert!((belief.probability_of("rainy") - 1.0/3.0).abs() < 1e-10);
    }

    #[test]
    fn test_weighted_belief() {
        let space = world_space();
        let belief = BeliefState::from_weights(&space, &[
            ("sunny".into(), 0.6), ("cloudy".into(), 0.3), ("rainy".into(), 0.1),
        ]).unwrap();
        assert!((belief.probability_of("sunny") - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_point_mass_belief() {
        let space = world_space();
        let belief = BeliefState::point_mass(&space, "sunny").unwrap();
        assert_eq!(belief.probability_of("sunny"), 1.0);
        assert_eq!(belief.probability_of("rainy"), 0.0);
    }

    #[test]
    fn test_belief_update() {
        let space = world_space();
        let belief = BeliefState::uniform(&space);
        let likelihood = RealValuedFunction::new(vec![
            ("sunny", 0.9), ("cloudy", 0.1), ("rainy", 0.0),
        ]);
        let updated = belief.update(&likelihood, &space).unwrap();
        assert!(updated.probability_of("sunny") > belief.probability_of("sunny"));
        assert_eq!(updated.probability_of("rainy"), 0.0);
    }

    #[test]
    fn test_observation_model() {
        let model = ObservationModel::new("weather_obs", &[
            ("sunny".into(), "bright".into()),
            ("cloudy".into(), "dim".into()),
            ("rainy".into(), "dim".into()),
        ]);
        assert_eq!(model.observe("sunny"), Some("bright"));
        assert_eq!(model.observe("rainy"), Some("dim"));
    }

    #[test]
    fn test_observation_likelihood() {
        let space = world_space();
        let model = ObservationModel::new("weather_obs", &[
            ("sunny".into(), "bright".into()),
            ("cloudy".into(), "dim".into()),
            ("rainy".into(), "dim".into()),
        ]);
        let like = model.likelihood("bright", &space);
        assert_eq!(like.eval("sunny"), 1.0);
        assert_eq!(like.eval("rainy"), 0.0);
    }

    #[test]
    fn test_utility_function() {
        let util = UtilityFunction::new(vec![
            ("sunny", 10.0), ("cloudy", 5.0), ("rainy", -5.0),
        ]);
        assert_eq!(util.utility("sunny"), 10.0);
    }

    #[test]
    fn test_expected_utility() {
        let space = world_space();
        let belief = BeliefState::from_weights(&space, &[
            ("sunny".into(), 0.5), ("cloudy".into(), 0.3), ("rainy".into(), 0.2),
        ]).unwrap();
        let util = UtilityFunction::new(vec![
            ("sunny", 10.0), ("cloudy", 5.0), ("rainy", -5.0),
        ]);
        let eu = util.expected_utility(&belief, &space);
        // 10*0.5 + 5*0.3 + (-5)*0.2 = 5 + 1.5 - 1 = 5.5
        assert!((eu - 5.5).abs() < 1e-10);
    }

    #[test]
    fn test_agent() {
        let obs_model = ObservationModel::new("obs", &[
            ("sunny".into(), "bright".into()),
            ("cloudy".into(), "dim".into()),
            ("rainy".into(), "dim".into()),
        ]);
        let util = UtilityFunction::new(vec![
            ("sunny", 10.0), ("cloudy", 5.0), ("rainy", -5.0),
        ]);
        let agent = Agent::new("weather_agent", obs_model, util);
        assert_eq!(agent.name(), "weather_agent");
    }

    #[test]
    fn test_agent_expected_utility() {
        let space = world_space();
        let belief = BeliefState::uniform(&space);
        let obs_model = ObservationModel::new("obs", &[
            ("sunny".into(), "bright".into()),
            ("cloudy".into(), "dim".into()),
            ("rainy".into(), "dim".into()),
        ]);
        let util = UtilityFunction::new(vec![
            ("sunny", 10.0), ("cloudy", 5.0), ("rainy", -5.0),
        ]);
        let agent = Agent::new("agent", obs_model, util);
        let eu = agent.expected_utility(&belief, &space);
        // (10 + 5 + (-5))/3 = 10/3 ≈ 3.333
        assert!((eu - 10.0/3.0).abs() < 1e-10);
    }

    #[test]
    fn test_belief_probability_of_set() {
        let space = world_space();
        let belief = BeliefState::from_weights(&space, &[
            ("sunny".into(), 0.6), ("cloudy".into(), 0.3), ("rainy".into(), 0.1),
        ]).unwrap();
        let set = Subset::from_slice(&["sunny", "cloudy"]);
        assert!((belief.probability_of_set(&set) - 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_agent_sequential_update() {
        let space = world_space();
        let mut belief = BeliefState::uniform(&space);
        
        // First observation: bright → likelihood favors sunny
        let like1 = RealValuedFunction::new(vec![
            ("sunny", 0.9), ("cloudy", 0.1), ("rainy", 0.0),
        ]);
        belief = belief.update(&like1, &space).unwrap();
        assert!(belief.probability_of("sunny") > 0.5);
        
        // Second observation: also bright → even more confident
        belief = belief.update(&like1, &space).unwrap();
        assert!(belief.probability_of("sunny") > 0.8);
    }
}
