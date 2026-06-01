# lau-measure-agents

> Measure theory for agents: sigma-algebras, Lebesgue integration, Radon-Nikodym, Fubini's theorem, and Riesz representation

## What This Does

This crate implements the foundations of measure theory — the mathematical framework underlying probability, integration, and functional analysis — applied to agent systems. It provides sigma-algebras, measures, Lebesgue measure on R^n, measurable functions, Lebesgue integration with convergence theorems (MCT, Fatou, DCT), product measures and Fubini's theorem, Radon-Nikodym derivatives, Lebesgue decomposition, Riesz representation, pushforward measures, and agent-specific constructs (belief states as measures, observations as measurable functions).

## The Key Idea

Measure theory is the operating system of probability and analysis. Everything built on top — distributions, expectations, Bayesian updates — reduces to measures on sigma-algebras. This crate builds that foundation explicitly: a sigma-algebra is a collection of sets closed under complement and countable union, a measure is a countably additive set function, and everything else follows. Agent beliefs are measures, observations are measurable functions, and Bayesian updates are Radon-Nikodym derivatives.

## Install

```toml
[dependencies]
lau-measure-agents = { git = "https://github.com/SuperInstance/lau-measure-agents" }
```

## Quick Start

```rust
use lau_measure_agents::*;
use lau_measure_agents::sigma_algebra::{Subset, SigmaAlgebra};

// Create a measurable space
let universal = Subset::from_slice(&["sunny", "cloudy", "rainy"]);
let sa = SigmaAlgebra::power_set(universal);

// Define a probability measure
let measure = Measure::new("weather", &sa, &[
    (Subset::from_slice(&["sunny"]), 0.4),
    (Subset::from_slice(&["cloudy"]), 0.3),
    (Subset::from_slice(&["rainy"]), 0.3),
])?;

// Measure of a set
let rainy = Subset::from_slice(&["rainy"]);
println!("P(rainy) = {}", measure.measure_of(&rainy));

// Measurable function: weather → mood
let f = MeasurableFunction::new(
    "weather_to_mood",
    &sa,
    &target_sa,
    &[
        ("sunny".into(), "happy".into()),
        ("cloudy".into(), "neutral".into()),
        ("rainy".into(), "sad".into()),
    ],
)?;

// Lebesgue integral: E[f(X)]
let g = RealValuedFunction::new(vec![("happy", 1.0), ("neutral", 0.0), ("sad", -1.0)]);
let integral = integrate(&g, &measure, &sa);
println!("Expected mood: {}", integral);

// Radon-Nikodym derivative (Bayesian update)
let prior = /* ... */;
let posterior = RadonNikodymDerivative::compute(&likelihood_measure, &prior, &elements, &sa);

// Agent state space
let state_space = AgentStateSpace::discrete("robot", &["idle", "working", "charging"]);
let belief = AgentBelief::uniform(&state_space);
```

## API Reference

### `sigma_algebra`

| Type | Description |
|------|-------------|
| `Subset` | Finite subset of elements. Supports union, intersection, complement. |
| `SigmaAlgebra::power_set(universal)` | Full power set (all subsets measurable). |
| `is_measurable(set)` | Check if a set is in the algebra. |
| `measurable_sets()` | Enumerate all measurable sets. |

### `measure`

| Type | Description |
|------|-------------|
| `Measure::new(name, algebra, values)` | Define a measure with explicit values. |
| `Measure::counting(algebra)` | Counting measure μ(A) = \|A\|. |
| `Measure::dirac(algebra, point)` | Point mass δₓ. |
| `measure_of(set)` | Measure of a measurable set. |
| `is_probability()` | Check μ(Ω) = 1. |

### `lebesgue`

| Type | Description |
|------|-------------|
| `Interval::new(left, right)` | Closed interval on R. |
| `length()` | Lebesgue measure of interval. |
| `intersection()`, `union()`, `complement()` | Set operations. |

### `measurable_function`

| Type | Description |
|------|-------------|
| `MeasurableFunction::new(name, domain, codomain, mapping)` | Verified measurable map. |
| `apply(element)` | Evaluate at a point. |
| `preimage(set)` | Inverse image of a measurable set. |

### `integral`

| Type | Description |
|------|-------------|
| `RealValuedFunction::new(values)` | Function on the measurable space. |
| `integrate(f, measure, algebra)` | Lebesgue integral ∫ f dμ. |
| `monotone_convergence(f_n)` | Verify MCT. |
| `fatous_lemma(f_n)` | Verify Fatou's lemma. |
| `dominated_convergence(f_n, g)` | Verify DCT. |

### `product`

| Function | Description |
|----------|-------------|
| `product_sigma_algebra(sa1, sa2)` | Product σ-algebra 𝔉⊗𝔊. |
| `product_measure(name, sa1, μ, sa2, ν, product_sa)` | Product measure μ×ν. |
| `fubini(f, μ, ν, sa1, sa2)` | Fubini: ∫∫ = ∫(∫). |

### `radon_nikodym`

| Type | Description |
|------|-------------|
| `RadonNikodymDerivative::compute(ν, μ, elements, sa)` | dν/dμ when ν ≪ μ. |
| `BeliefUpdate` | Bayesian update as Radon-Nikodym derivative. |

### `decomposition`

| Type | Description |
|------|-------------|
| `LebesgueDecomposition::compute(ν, μ, elements, sa)` | ν = ν_ac + ν_singular. |

### `riesz`

| Type | Description |
|------|-------------|
| `RieszRepresentation` | Positive linear functional ↔ integral against a measure. |

### `pushforward`

| Function | Description |
|----------|-------------|
| `pushforward(name, f, μ, domain_sa, codomain_sa)` | T*μ(B) = μ(T⁻¹(B)). |

### `agent`

| Type | Description |
|------|-------------|
| `AgentStateSpace::discrete(name, elements)` | Agent state as measurable space. |
| `AgentBelief` | Belief state as probability measure. |
| `update_belief(observation)` | Bayesian update via Radon-Nikodym. |

## How It Works

1. **Sigma-Algebra**: For finite spaces, the power set gives all subsets as measurable.
2. **Measure**: Assigns non-negative values to measurable sets. Validated for σ-additivity.
3. **Integration**: Lebesgue integral via simple function approximation (exact for finite spaces).
4. **Convergence Theorems**: MCT (monotone ↑ → integral ↑), Fatou (liminf integral ≤ integral of liminf), DCT (dominated convergence).
5. **Radon-Nikodym**: When ν ≪ μ, the derivative dν/dμ exists and ν(A) = ∫_A (dν/dμ) dμ.
6. **Agent Belief Updates**: Observation updates the belief measure via Radon-Nikodym derivative.

## The Math

- **Sigma-Algebra**: Σ ⊆ 2^Ω closed under complement and countable union.
- **Measure**: μ: Σ → [0,∞] with μ(∅) = 0, μ(∪ Aₙ) = Σ μ(Aₙ) for disjoint Aₙ.
- **Lebesgue Integral**: ∫ f dμ = sup{∫ s dμ : s ≤ f, s simple}.
- **Radon-Nikodym**: If ν ≪ μ, ∃ dν/dμ such that ν(A) = ∫_A (dν/dμ) dμ.
- **Fubini**: ∫_{X×Y} f d(μ×ν) = ∫_X(∫_Y f dν) dμ.
- **Riesz**: Every positive linear functional on C(X) is integration against a measure.

## Testing

134 tests covering:
- Sigma-algebra closure properties
- Measure construction and validation
- Lebesgue integral computation
- Monotone/Fatou/Dominated convergence theorems
- Product measure and Fubini's theorem
- Radon-Nikodym derivative computation
- Lebesgue decomposition
- Riesz representation
- Pushforward measure
- Agent belief updates

## License

MIT
