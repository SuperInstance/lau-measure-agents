# lau-measure-agents

> Measure theory for agents — the rigorous foundation of integration and probability.

A Rust crate implementing **measure theory** from the ground up: σ-algebras, measures, Lebesgue measure, measurable functions, the Lebesgue integral with all three convergence theorems, product measures and Fubini's theorem, the Radon-Nikodym theorem, Lebesgue decomposition, the Riesz representation theorem, pushforward measures, and Bayesian belief updates for agents.

Every result is verified by **134 unit tests**.

---

## What This Does

This crate gives you concrete, computable implementations of the core structures of measure theory on finite (or countably represented) spaces:

- **σ-algebras** — trivial, power-set, and generated from seed sets; validated against axioms
- **Measures** — counting, Dirac, probability, uniform, with signed measures and Jordan decomposition
- **Lebesgue measure** — intervals, boxes, simple interval functions on ℝⁿ
- **Measurable functions** — with preimage/image, composition, indicator and simple functions
- **Lebesgue integral** — with the Monotone Convergence Theorem, Fatou's Lemma, and the Dominated Convergence Theorem; also Lᵖ norms and Hölder's inequality
- **Product measures** — product σ-algebras, product measures on rectangles, and Fubini/Tonelli iterated integration
- **Radon-Nikodym theorem** — derivatives dν/dμ, chain rule, verification of the representation
- **Lebesgue decomposition** — split ν = ν_ac + ν_singular; Hahn decomposition for signed measures
- **Riesz representation** — positive linear functionals ↔ measures
- **Pushforward measures** — measure transport through measurable maps, change of variables
- **Agent applications** — Bayesian belief states, observation models, utility functions, expected utility

All structures serialize via `serde`.

---

## Key Idea

**Measure theory is the foundation of probability, integration, and information theory.** This crate treats an agent's *beliefs* as probability measures, observations as measurable functions, and belief updates as Radon-Nikodym derivatives.

When an agent sees evidence, it updates via Bayes' rule:

```
P(state | evidence) ∝ P(evidence | state) · P(state)
```

In measure-theoretic terms, the posterior dν/dμ is the Radon-Nikodym derivative of the joint measure with respect to the prior. The crate makes this explicit: `BeliefState::update` computes the posterior as a new probability measure, verifying that the update is well-defined (no division by zero = no absolute continuity violation).

The rest of the crate provides the mathematical machinery (σ-algebras, integrals, product measures) that makes this rigorous.

---

## Install

```toml
[dependencies]
lau-measure-agents = "0.1"
```

Or clone directly:

```bash
git clone https://github.com/SuperInstance/lau-measure-agents.git
cargo build
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `nalgebra` 0.33 | Linear algebra for Lebesgue measure on ℝⁿ |
| `serde` 1 | Serialization of all mathematical structures |
| `approx` 0.5 (dev) | Floating-point assertions in tests |

---

## Quick Start

### σ-algebra and measures

```rust
use lau_measure_agents::{Subset, SigmaAlgebra, Measure};

let omega = Subset::from_slice(&["sunny", "cloudy", "rainy"]);
let sa = SigmaAlgebra::power_set(omega);
assert!(sa.validate());
assert_eq!(sa.size(), 8); // 2³

let prior = Measure::probability(&sa, &[
    ("sunny".into(), 0.5),
    ("cloudy".into(), 0.3),
    ("rainy".into(), 0.2),
]).unwrap();
assert!(prior.is_probability());
```

### Lebesgue integral and convergence theorems

```rust
use lau_measure_agents::{RealValuedFunction, integrate, monotone_convergence};

let f1 = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0)]);
let f2 = RealValuedFunction::new(vec![("a", 1.5), ("b", 2.5)]);
let f3 = RealValuedFunction::new(vec![("a", 2.0), ("b", 3.0)]);

// Monotone convergence: lim ∫fₙ = ∫(lim fₙ)
let (lim_int, int_lim) = monotone_convergence(&[f1, f2, f3], &mu);
assert!((lim_int - int_lim).abs() < 1e-10);
```

### Radon-Nikodym derivative

```rust
use lau_measure_agents::{RadonNikodymDerivative, LebesgueDecomposition};

let rn = RadonNikodymDerivative::compute(&posterior, &prior, &elements).unwrap();
assert!(rn.verify(&posterior, &prior, &test_sets));

// Lebesgue decomposition: ν = ν_ac + ν_singular
let decomp = LebesgueDecomposition::compute(&nu, &mu, &elements, &sa);
```

### Agent belief update

```rust
use lau_measure_agents::{BeliefState, ObservationModel, UtilityFunction, Agent};

let belief = BeliefState::from_weights(&space, &[
    ("sunny".into(), 0.5), ("cloudy".into(), 0.3), ("rainy".into(), 0.2),
]).unwrap();

// Bayes update: see "bright" → likelihood favors sunny
let likelihood = RealValuedFunction::new(vec![
    ("sunny", 0.9), ("cloudy", 0.1), ("rainy", 0.0),
]);
let posterior = belief.update(&likelihood, &space).unwrap();
assert!(posterior.probability_of("sunny") > 0.5);
```

---

## API Reference

### `sigma_algebra` — The Foundation of Measure Theory

| Type / Function | Description |
|----------------|-------------|
| `Element` | Alias for `String`; an element of the underlying set |
| `Subset` | Sorted set of elements. Methods: `empty`, `from_slice`, `union`, `intersection`, `difference`, `complement`, `is_empty`, `len`, `contains`, `is_subset_of` |
| `SigmaAlgebra` | Collection of measurable subsets closed under complement and countable union. Methods: `trivial(ω)`, `power_set(ω)`, `generate(ω, seeds)`, `universal`, `size`, `is_measurable`, `measurable_sets`, `validate`, `countable_union`, `countable_intersection`, `complement` |
| `MeasurableSpace` | Pair (Ω, 𝔉). Methods: `new`, `sigma_algebra`, `universal` |

`generate` computes the smallest σ-algebra containing the seed sets by closing under unions, intersections, and complements to a fixed point.

### `measure` — Set Functions

| Type / Function | Description |
|----------------|-------------|
| `Measure` | Non-negative countably additive set function. Constructors: `new`, `counting`, `dirac(point)`, `probability(weights)`, `uniform`. Methods: `measure_of`, `name`, `total_mass`, `is_probability`, `is_finite`, `is_absolutely_continuous_wrt`, `is_singular_with`, `scale`, `add`, `subtract` |
| `SignedMeasure` | Can take negative values. Methods: `new`, `measure_of`, `total_variation`, `jordan_decomposition` |

**Absolute continuity**: ν ≪ μ iff μ(A) = 0 ⟹ ν(A) = 0. **Mutual singularity**: μ ⊥ ν iff they concentrate on disjoint sets.

### `lebesgue` — Lebesgue Measure on ℝⁿ

| Type / Function | Description |
|----------------|-------------|
| `Interval` | [a, b] ⊂ ℝ. Methods: `new`, `length`, `contains`, `intersection`, `union`, `complement` |
| `Box` | Product of intervals ⊂ ℝⁿ. Methods: `new`, `volume`, `contains`, `intersection` |
| `LebesgueMeasure` | Standard measure on ℝⁿ. Methods: `new(dim)`, `measure_interval`, `measure_box`, `measure_disjoint_intervals`, `dimension` |
| `SimpleIntervalFunction` | Piecewise constant on intervals. Methods: `new`, `eval`, `integrate`, `supremum` |

### `measurable_function` — Maps Between Measurable Spaces

| Type / Function | Description |
|----------------|-------------|
| `MeasurableFunction` | f: (Ω, 𝔉) → (Ω', 𝔊') with preimages of measurable sets being measurable. Methods: `new` (validates), `new_unchecked`, `apply`, `preimage`, `image`, `compose`, `name`, `mapping` |
| `IndicatorFunction` | 𝟙_A(x) = 1 if x ∈ A, 0 otherwise. Methods: `new`, `eval`, `set` |
| `SimpleFunction` | f = Σ cᵢ · 𝟙_{Aᵢ}. Methods: `new`, `eval`, `supremum`, `integrate(μ)`, `add`, `scale` |

Measurability is validated on construction: for every measurable B in the codomain, f⁻¹(B) must be measurable in the domain.

### `integral` — Lebesgue Integration and Convergence Theorems

| Type / Function | Description |
|----------------|-------------|
| `RealValuedFunction` | Element → ℝ. Methods: `new`, `eval`, `elements`, `supremum`, `infimum`, `pointwise_max`, `pointwise_min`, `abs` |
| `integrate(f, μ)` | ∫f dμ = Σ f(xᵢ) · μ({xᵢ}) |
| `integrate_over(f, μ, A)` | ∫_A f dμ (integral restricted to a measurable set) |
| `lp_norm(f, μ, p)` | ‖f‖_p = (∫|f|ᵖ dμ)^{1/p} |
| `monotone_convergence(fₙ, μ)` | MCT: lim ∫fₙ = ∫(lim fₙ) for 0 ≤ f₁ ≤ f₂ ≤ ⋯ |
| `fatou_lemma(fₙ, μ)` | ∫liminf fₙ ≤ liminf ∫fₙ |
| `dominated_convergence(fₙ, g, μ)` | DCT: if |fₙ| ≤ g and fₙ → f, then ∫fₙ → ∫f |
| `holders_inequality(f, g, μ, p)` | ∫|fg| ≤ ‖f‖_p · ‖g‖_q |

### `product` — Product Measures and Fubini's Theorem

| Type / Function | Description |
|----------------|-------------|
| `pair_key(x, y)` | String key "(x,y)" for product space elements |
| `product_sigma_algebra(𝔉, 𝔊)` | 𝔉 ⊗ 𝔊 (power set of X × Y for finite spaces) |
| `product_measure(name, sa₁, μ, sa₂, ν, product_sa)` | μ × ν with (μ × ν)(A × B) = μ(A) · ν(B) |
| `rectangle_measure(A, μ, B, ν)` | Shortcut for μ(A) · ν(B) |
| `fubini_integrate(f, sa₁, μ, sa₂, ν)` | ∫_{X×Y} f d(μ×ν) = ∫_X (∫_Y f dν) dμ |
| `tonelli_integrate(f, sa₁, μ, sa₂, ν)` | Same for non-negative f (no integrability condition) |

### `radon_nikodym` — Densities and Derivatives

| Type / Function | Description |
|----------------|-------------|
| `RadonNikodymDerivative` | dν/dμ: the density f where ν(A) = ∫_A f dμ. Methods: `compute(ν, μ, elements)`, `eval`, `verify(ν, μ, test_sets)`, `chain_rule` |

`compute` verifies ν ≪ μ, then sets f(x) = ν({x}) / μ({x}) for each element. `chain_rule` composes derivatives: dλ/dμ = (dλ/dν) · (dν/dμ).

### `decomposition` — Lebesgue and Hahn Decompositions

| Type / Function | Description |
|----------------|-------------|
| `LebesgueDecomposition` | ν = ν_ac + ν_singular. Methods: `compute(ν, μ, elements, sa)`, `verify`, `is_pure` → `PureType` (AC / Singular / Mixed) |
| `PureType` | Enum: `AbsolutelyContinuous`, `Singular`, `Mixed` |
| `HahnDecomposition` | Ω = P ∪ N for signed measure (P = positive, N = negative). Method: `compute(values)` |

### `riesz` — Riesz Representation Theorem

| Type / Function | Description |
|----------------|-------------|
| `LinearFunctional` | L: functions → ℝ via L(f) = Σ f(xᵢ) · L(𝟙_{xᵢ}). Methods: `new`, `evaluate`, `is_positive`, `norm`, `scale`, `add` |
| `RieszRepresentation` | Converts positive linear functional to measure: μ({x}) = L(𝟙_x). Methods: `represent(L, sa)`, `verify` |

### `pushforward` — Measure Transport

| Type / Function | Description |
|----------------|-------------|
| `pushforward(name, f, μ, dom, cod)` | T\*μ where (T\*μ)(B) = μ(f⁻¹(B)) |
| `pushforward_probability(name, f, μ, dom, cod)` | Same, with finiteness check |
| `pullback(name, f, ν, dom, cod)` | Inverse of pushforward for bijections |
| `change_of_variables(g, f, μ, elements)` | ∫g d(T\*μ) = ∫(g ∘ T) dμ |

### `agent` — Bayesian Agents

| Type / Function | Description |
|----------------|-------------|
| `BeliefState` | Probability measure over world states. Constructors: `from_weights`, `uniform`. Methods: `probability_of`, `probability_of_set`, `entropy`, `update(likelihood, space)`, `measure` |
| `ObservationModel` | Maps world states to observations. Methods: `new`, `observe`, `likelihood(obs, space)` |
| `UtilityFunction` | Maps states to utility values. Methods: `new`, `utility`, `expected_utility(belief, space)` |
| `Agent` | Bundles observation model + utility. Methods: `new`, `name`, `expected_utility` |

`BeliefState::update` implements Bayes' rule: P'(s) ∝ likelihood(s) · P(s), then normalizes.

---

## How It Works

### Architecture

```
sigma_algebra ──→ measure ──→ lebesgue
     │               │
     ├── measurable_function ──→ integral (MCT, Fatou, DCT, Hölder)
     │                                      │
     │                                      ├── product (Fubini, Tonelli)
     │                                      ├── radon_nikodym
     │                                      ├── decomposition (Lebesgue, Hahn)
     │                                      ├── riesz
     │                                      └── pushforward
     │
     └────────────────────────────→ agent (Bayesian belief update)
```

### Finite Spaces

This crate works on **finite** (or countably represented) spaces. This is deliberate:

- Every computation is exact (no infinite sums or limits to approximate)
- σ-algebras can be explicitly enumerated (power sets)
- Measures are stored as HashMaps
- All convergence theorems can be verified directly

The mathematical structure extends to infinite-dimensional measure theory conceptually; here we compute the finite-dimensional cases exactly.

### String-Based Elements

Elements are represented as `String`s rather than numeric indices. This makes the API more intuitive (`"sunny"`, `"rainy"`) and lets the crate handle arbitrary discrete spaces without requiring a fixed dimension at compile time.

### Integration Strategy

For finite spaces, the Lebesgue integral simplifies to a weighted sum:

```
∫f dμ = Σᵢ f(xᵢ) · μ({xᵢ})
```

The convergence theorems (MCT, Fatou, DCT) are verified by comparing sequences of such sums — still rigorous, but computationally trivial.

---

## The Math

### σ-Algebras

A **σ-algebra** 𝔉 on a set Ω is a collection of subsets satisfying:

1. ∅ ∈ 𝔉 and Ω ∈ 𝔉
2. Closed under complement: A ∈ 𝔉 ⟹ Aᶜ ∈ 𝔉
3. Closed under countable union: A₁, A₂, … ∈ 𝔉 ⟹ ∪Aᵢ ∈ 𝔉

The pair (Ω, 𝔉) is a **measurable space**. The smallest σ-algebra containing a collection of sets is the **generated** σ-algebra. The largest is the **power set** 2^Ω.

### Measures

A **measure** μ on (Ω, 𝔉) assigns [0, ∞] to each measurable set with:

1. μ(∅) = 0
2. **Countable additivity**: for disjoint {Aᵢ}, μ(∪Aᵢ) = Σμ(Aᵢ)

A **probability measure** has total mass 1. The **counting measure** gives μ(A) = |A|. The **Dirac measure** δₓ concentrates all mass at one point.

### Lebesgue Integral

The integral is built in stages:

1. **Indicator functions**: ∫𝟙_A dμ = μ(A)
2. **Simple functions**: ∫(Σ cᵢ · 𝟙_{Aᵢ}) dμ = Σ cᵢ · μ(Aᵢ)
3. **Non-negative functions**: ∫f dμ = sup{∫s dμ : s ≤ f, s simple}
4. **General functions**: ∫f dμ = ∫f⁺ dμ − ∫f⁻ dμ

### Convergence Theorems

Three pillars of integration theory:

- **Monotone Convergence (MCT)**: If 0 ≤ f₁ ≤ f₂ ≤ ⋯ and fₙ → f pointwise, then ∫fₙ → ∫f.
- **Fatou's Lemma**: ∫liminf fₙ ≤ liminf ∫fₙ (always holds; gives a lower bound).
- **Dominated Convergence (DCT)**: If fₙ → f pointwise and |fₙ| ≤ g for some integrable g, then ∫fₙ → ∫f.

### Product Measures and Fubini

Given (X, μ) and (Y, ν), the **product measure** μ × ν satisfies (μ × ν)(A × B) = μ(A) · ν(B).

**Fubini's theorem**: ∫_{X×Y} f d(μ×ν) = ∫_X (∫_Y f(x,y) dν(y)) dμ(x)

**Tonelli's theorem**: Same conclusion for non-negative f, without the integrability hypothesis.

### Radon-Nikodym Theorem

If ν ≪ μ (ν is absolutely continuous w.r.t. μ), there exists a unique (up to μ-null sets) measurable function f = dν/dμ such that:

```
ν(A) = ∫_A f dμ   for all A ∈ 𝔉
```

This is the rigorous foundation of **probability densities**. The chain rule dλ/dμ = (dλ/dν) · (dν/dμ) follows from the theorem.

### Lebesgue Decomposition

Any σ-finite measure ν can be uniquely decomposed as ν = ν_ac + ν_singular where:

- ν_ac ≪ μ (absolutely continuous part)
- ν_singular ⊥ μ (singular part, concentrates on μ-null sets)

The **Hahn decomposition** splits Ω into positive and negative sets for signed measures, enabling the Jordan decomposition μ = μ⁺ − μ⁻.

### Riesz Representation

Every **positive linear functional** L on C_c(X) (continuous functions with compact support) is integration against a unique regular Borel measure:

```
L(f) = ∫f dμ
```

In finite dimensions, this is trivially μ({x}) = L(𝟙_x), but the structural correspondence between functionals and measures is profound.

### Pushforward Measures

For a measurable map T: (Ω, 𝔉, μ) → (Ω', 𝔊), the **pushforward** T\*μ is defined by:

```
(T*μ)(B) = μ(T⁻¹(B))   for all B ∈ 𝔊
```

The **change of variables formula**: ∫g d(T\*μ) = ∫(g ∘ T) dμ. This is the measure-theoretic version of substitution in integrals.

---

## License

MIT
