# lau-measure-agents

**Measure theory for agents — the rigorous foundation of integration and probability**

A Rust library implementing measure theory from the ground up: σ-algebras, measures, Lebesgue measure, measurable functions, the Lebesgue integral with convergence theorems, product measures with Fubini's theorem, Radon-Nikodym derivatives, Lebesgue decomposition, Riesz representation, pushforward measures, and agent belief-state systems built on these foundations.

134 tests · 12 modules · ~3,600 LOC

---

## What This Does

Measure theory is the **rigorous mathematical framework** underlying probability, integration, and analysis. This library implements it computationally for finite/discrete spaces:

- **σ-algebras** — collections of measurable sets, closed under complements and countable unions; generated from arbitrary subsets; validated for closure properties
- **Measures** — non-negative countably additive set functions: counting, Dirac (point mass), uniform, probability measures; signed measures with Jordan decomposition
- **Lebesgue measure** — intervals, boxes in ℝⁿ, outer measure approximation, simple function integration
- **Measurable functions** — maps between measurable spaces with validated preimages; indicator functions, simple functions
- **Lebesgue integral** — integration of real-valued functions over measure spaces; convergence theorems (MCT, Fatou, DCT); Hölder's inequality; Lᵖ norms
- **Product measures** — Fubini's theorem for iterated integration; Tonelli's theorem for non-negative functions
- **Radon-Nikodym theorem** — density functions dν/dμ when ν ≪ μ; chain rule; log-derivative; Bayesian belief updates
- **Lebesgue decomposition** — ν = ν_ac + ν_singular relative to a reference measure; Hahn decomposition for signed measures
- **Riesz representation** — every positive linear functional ↔ integration against a measure
- **Pushforward measure** — measure induced by a measurable map; change of variables
- **Agent systems** — state spaces as measurable spaces, beliefs as probability measures, observations as measurable functions, belief updates via Radon-Nikodym, expected utility as Lebesgue integral

---

## Key Idea

Probability is just measure theory where μ(Ω) = 1. Every probabilistic concept — random variables, expectations, Bayes' theorem, conditional probability — has a precise measure-theoretic definition. This library makes those definitions concrete and computational:

1. **Belief states = probability measures** on a σ-algebra of states
2. **Observations = measurable functions** mapping states to observations
3. **Belief updates = Radon-Nikodym derivatives** (Bayes' rule is dπ(·|obs)/dπ)
4. **Expected utility = Lebesgue integral** of utility against belief

This gives agents a **mathematically rigorous** foundation for reasoning under uncertainty.

---

## Install

Add to your `Cargo.toml`:

```toml
[dependencies]
lau-measure-agents = { git = "https://github.com/SuperInstance/lau-measure-agents" }
```

### Dependencies

- `nalgebra` 0.33 — linear algebra (Lebesgue measure in ℝⁿ)
- `serde` 1.x (with `derive`) — serialization
- `approx` 0.5 (dev) — floating-point comparison

---

## Quick Start

```rust
use lau_measure_agents::*;

// --- Sigma-algebra ---
let universal = Subset::from_slice(&["sunny", "cloudy", "rainy"]);
let algebra = SigmaAlgebra::power_set(universal);

// --- Measure (probability) ---
let weather = Measure::probability(&algebra, &[
    ("sunny".into(), 0.5),
    ("cloudy".into(), 0.3),
    ("rainy".into(), 0.2),
]).unwrap();

// --- Measurable function (observation) ---
let obs = MeasurableFunction::new(
    "sky_obs",
    &algebra, &algebra,
    &[
        ("sunny".into(), "bright".into()),
        ("cloudy".into(), "dim".into()),
        ("rainy".into(), "dim".into()),
    ],
).unwrap();

// Preimage: obs⁻¹({"dim"}) = {"cloudy", "rainy"}
let dim_set = Subset::from_slice(&["dim"]);
let preimage = obs.preimage(&dim_set);

// --- Lebesgue integral ---
let utility = RealValuedFunction::new(vec![
    ("sunny", 10.0), ("cloudy", 3.0), ("rainy", -2.0),
]);
let expected = integrate(&utility, &weather);
// = 10*0.5 + 3*0.3 + (-2)*0.2 = 5.5

// --- Agent with beliefs ---
let space = AgentStateSpace::discrete("weather", &["sunny", "cloudy", "rainy"]);
let belief = BeliefState::uniform(&space);
let utility_fn = UtilityFunction::new(vec![
    ("sunny", 10.0), ("cloudy", 3.0), ("rainy", -2.0),
]);
let agent = Agent::new("planner", &space, &utility_fn, &[
    ("sunny", "bright"), ("cloudy", "dim"), ("rainy", "dim"),
]);
let eu = agent.expected_utility(&belief, &space);

// --- Pushforward measure ---
let pushed = pushforward("sky_dist", &obs, &weather, &algebra, &algebra);

// --- Radon-Nikodym (Bayes' rule) ---
let prior = Measure::uniform(&algebra).unwrap();
let posterior = Measure::probability(&algebra, &[
    ("sunny".into(), 0.7), ("cloudy".into(), 0.2), ("rainy".into(), 0.1),
]).unwrap();
let rn = RadonNikodymDerivative::compute(&posterior, &prior, &algebra);
// dν/dμ is the Bayes factor
```

---

## API Reference

### `sigma_algebra` — σ-Algebras

| Type / Function | Description |
|---|---|
| `Subset(BTreeSet<String>)` | A measurable subset of the universal set |
| `Subset::empty()` | The empty set |
| `Subset::from_slice(elems)` | Construct from string slice |
| `a.union(b)` | Set union |
| `a.intersection(b)` | Set intersection |
| `a.complement(universal)` | Set complement |
| `a.difference(b)` | Set difference |
| `SigmaAlgebra` | A σ-algebra: collection of measurable sets |
| `SigmaAlgebra::trivial(universal)` | {∅, Ω} — the trivial σ-algebra |
| `SigmaAlgebra::power_set(universal)` | 2^Ω — the finest σ-algebra |
| `SigmaAlgebra::generate(universal, generators)` | Generated σ-algebra (smallest containing generators) |
| `sa.is_measurable(s)` | Check if a set is in the σ-algebra |
| `sa.validate()` | Verify closure under complement and countable union |
| `MeasurableSpace` | A set equipped with a σ-algebra |

### `measure` — Measures

| Type / Function | Description |
|---|---|
| `Measure` | A non-negative countably additive set function |
| `Measure::new(name, algebra, values)` | Create from explicit set-measure pairs |
| `Measure::counting(algebra)` | Counting measure: μ(A) = \|A\| |
| `Measure::dirac(algebra, point)` | Dirac measure: δ_x(A) = 1 if x ∈ A |
| `Measure::uniform(algebra)` | Uniform probability measure |
| `Measure::probability(algebra, weights)` | Weighted probability measure |
| `m.measure_of(set)` | μ(A) — measure of a set |
| `m.total_mass()` | μ(Ω) |
| `m.is_probability()` | Check μ(Ω) = 1 |
| `m.is_finite()` | Check μ(Ω) < ∞ |
| `m.is_absolutely_continuous_wrt(other)` | Check ν ≪ μ |
| `m.is_singular_with(other)` | Check ν ⊥ μ |
| `m.scale(c)` | c · μ |
| `m.add(other)` | μ + ν |
| `SignedMeasure` | A signed measure (allows negative values) |
| `sm.jordan_decomposition()` | μ = μ⁺ - μ⁻ (Jordan decomposition) |

### `lebesgue` — Lebesgue Measure on ℝⁿ

| Type / Function | Description |
|---|---|
| `Interval` | A closed interval [a, b] |
| `Interval::length()` | b - a (Lebesgue measure of interval) |
| `Interval::intersection(other)` | Interval intersection |
| `Interval::complement(lo, hi)` | Complement within bounds |
| `Box` | A box in ℝⁿ (Cartesian product of intervals) |
| `Box::volume()` | ∏(b_i - a_i) (Lebesgue measure of box) |
| `LebesgueMeasure` | Approximate Lebesgue measure on ℝⁿ |
| `lm.measure_box(b)` | Volume of a box |
| `lm.outer_measure_approx(points, ε)` | Outer measure approximation |
| `SimpleIntervalFunction` | Piecewise constant on intervals |
| `sif.integrate()` | ∫ f dλ over its support |

### `measurable_function` — Measurable Functions

| Type / Function | Description |
|---|---|
| `MeasurableFunction` | A map f: (Ω, 𝔉) → (Ω', 𝔊') validated for measurability |
| `MeasurableFunction::new(name, domain, codomain, mapping)` | Create with measurability check |
| `mf.apply(elem)` | f(x) |
| `mf.preimage(set)` | f⁻¹(B) — preimage of a measurable set |
| `mf.image(set)` | f(A) — direct image |
| `mf.compose(other)` | g ∘ f — composition |
| `IndicatorFunction` | 𝟙_A(x) = 1 if x ∈ A, 0 otherwise |
| `SimpleFunction` | Σ cᵢ · 𝟙_{Aᵢ} — simple function on measurable sets |
| `sf.integrate(mu)` | ∫ f dμ (exact for simple functions) |

### `integral` — Lebesgue Integral

| Type / Function | Description |
|---|---|
| `RealValuedFunction` | A function f: Ω → ℝ on a finite space |
| `integrate(f, mu)` | ∫ f dμ — Lebesgue integral |
| `integrate_over(f, mu, set)` | ∫_A f dμ — integral over a subset |
| `monotone_convergence(f_n, mu)` | Verify MCT: lim ∫ fₙ dμ = ∫ lim fₙ dμ |
| `fatou_lemma(f_n, mu)` | Fatou: ∫ lim inf fₙ ≤ lim inf ∫ fₙ |
| `dominated_convergence(f_n, mu, g)` | DCT: if \|fₙ\| ≤ g, swap limit and integral |
| `holders_inequality(f, g, mu, p, q)` | Hölder: \|fg\|₁ ≤ \|f\|_p · \|g\|_q |
| `lp_norm(f, mu, p)` | Lᵖ norm (\|f\|_p = (∫ \|f\|ᵖ dμ)^{1/p}) |

### `product` — Product Measures & Fubini's Theorem

| Function | Description |
|---|---|
| `product_sigma_algebra(sa1, sa2)` | Construct 𝔉 ⊗ 𝔊 |
| `product_measure(mu1, mu2)` | Construct μ × ν |
| `rectangle_measure(mu1, mu2, A, B)` | (μ × ν)(A × B) = μ(A) · ν(B) |
| `fubini_integrate(f, mu1, mu2)` | ∫∫ f d(μ × ν) = ∫ (∫ f(x,y) dν) dμ |
| `tonelli_integrate(f, mu1, mu2)` | Tonelli's theorem (non-negative case) |

### `radon_nikodym` — Radon-Nikodym Theorem

| Type / Function | Description |
|---|---|
| `RadonNikodymDerivative` | dν/dμ — the density function |
| `RadonNikodymDerivative::compute(nu, mu, algebra)` | Compute dν/dμ element-by-element |
| `rn.eval(elem)` | Evaluate dν/dμ at a point |
| `rn.verify(nu, mu, sets)` | Verify ν(A) = ∫_A (dν/dμ) dμ for test sets |
| `rn.chain_rule(other)` | Chain rule: dλ/dμ = (dλ/dν)(dν/dμ) |
| `rn.log_derivative()` | log(dν/dμ) — log-likelihood ratio |
| `BeliefUpdate` | Bayesian update: prior → posterior via likelihood |
| `BeliefUpdate::from_likelihood(prior, likelihood, algebra)` | Bayes' rule as Radon-Nikodym |

### `decomposition` — Lebesgue & Hahn Decomposition

| Type / Function | Description |
|---|---|
| `LebesgueDecomposition` | ν = ν_ac + ν_singular |
| `LebesgueDecomposition::compute(nu, mu)` | Decompose ν relative to μ |
| `ld.verify(nu, elements)` | Verify decomposition correctness |
| `ld.is_pure()` | Check if purely AC, purely singular, or mixed |
| `HahnDecomposition` | Ω = P ∪ N for a signed measure |
| `HahnDecomposition::compute(signed_measure)` | Find positive and negative sets |

### `riesz` — Riesz Representation Theorem

| Type / Function | Description |
|---|---|
| `LinearFunctional` | L: functions → ℝ (linear map) |
| `LinearFunctional::new(name, values)` | Create from values on basis functions |
| `lf.evaluate(f_values)` | L(f) |
| `lf.is_positive()` | Check L(f) ≥ 0 for f ≥ 0 |
| `lf.norm()` | Operator norm of L |
| `RieszRepresentation` | The representing measure μ where L(f) = ∫ f dμ |
| `RieszRepresentation::represent(functional)` | Recover μ from L |
| `RieszRepresentation::verify(functional, measure)` | Verify L(f) = ∫ f dμ |

### `pushforward` — Pushforward Measure

| Function | Description |
|---|---|
| `pushforward(name, f, mu, domain, codomain)` | T\*μ where (T\*μ)(B) = μ(f⁻¹(B)) |
| `pushforward_probability(name, f, mu, domain, codomain)` | Pushforward of a probability measure |
| `pullback(set, f)` | f⁻¹(A) — pullback of a measurable set |
| `change_of_variables(f, mu, domain, codomain, g)` | ∫ g d(T\*μ) = ∫ (g ∘ T) dμ |

### `agent` — Agent Systems

| Type / Function | Description |
|---|---|
| `AgentStateSpace` | A measurable space for agent states |
| `AgentStateSpace::discrete(name, elements)` | Create finite state space |
| `BeliefState` | A probability measure on the state space |
| `BeliefState::uniform(space)` | Uniform belief |
| `BeliefState::from_weights(space, weights)` | Weighted belief |
| `BeliefState::point_mass(space, point)` | Certainty at a single state |
| `bs.probability_of(state)` | P(state) |
| `bs.probability_of_set(set)` | P(A) |
| `bs.update(likelihood, space)` | Bayesian update via Radon-Nikodym |
| `ObservationModel` | Maps states → observations (measurable function) |
| `om.likelihood(observation, space)` | Likelihood function L(obs \| state) |
| `UtilityFunction` | f: states → ℝ |
| `uf.expected_utility(belief, space)` | E[u] = ∫ u dμ (Lebesgue integral) |
| `Agent` | Complete agent: state space + utility + observation model |

---

## How It Works

### Architecture

```
sigma_algebra (measurable sets, σ-algebra generation)
    ├── measure (non-negative measures, signed measures)
    │   ├── lebesgue (intervals, boxes, Lebesgue measure on ℝⁿ)
    │   ├── product (product measures, Fubini/Tonelli)
    │   ├── radon_nikodym (densities, Bayes updates)
    │   ├── decomposition (Lebesgue/Hahn decomposition)
    │   └── pushforward (pushforward/pullback measures)
    ├── measurable_function (measurable maps, simple functions)
    │   └── integral (Lebesgue integral, convergence theorems)
    │       └── riesz (Riesz representation)
    └── agent (state spaces, beliefs, utilities, observations)
```

### Finite Measure Theory

All constructions work over **finite sets** with explicit enumeration. A σ-algebra on a finite set is just a collection of subsets closed under complement and union. This makes every operation computable:

- Measure: stored as a HashMap from sets to values
- Integral: sum over elements weighted by measure
- Radon-Nikodym: ratio of measure values element-by-element
- Fubini: double sum (iterating over product space)

### Sigma-Algebra Generation

Given generators {A₁, A₂, ...}, the generated σ-algebra is the smallest collection containing them that's closed under complement and countable union. For finite sets, we compute this by iterative closure: start with generators + ∅ + Ω, then repeatedly add complements and unions until stable.

### Bayesian Updates as Radon-Nikodym

Bayes' rule says: P(A|obs) ∝ P(obs|A) · P(A). In measure-theoretic terms, the posterior dπ(·|obs) is the Radon-Nikodym derivative of the pushforward times the prior, normalized. This is exactly `BeliefUpdate::from_likelihood`.

---

## The Math

### σ-Algebra

A **σ-algebra** 𝔉 on Ω satisfies: ∅ ∈ 𝔉, A ∈ 𝔉 ⟹ Aᶜ ∈ 𝔉, A₁, A₂, ... ∈ 𝔉 ⟹ ∪Aᵢ ∈ 𝔉.

### Measure

A **measure** μ on (Ω, 𝔉) satisfies: μ(∅) = 0, μ(A) ≥ 0, and for disjoint {Aᵢ}: μ(∪Aᵢ) = Σμ(Aᵢ).

### Lebesgue Integral

For f ≥ 0 measurable:

$$\int f \, d\mu = \sup\left\{\int s \, d\mu : s \text{ simple}, s \le f\right\}$$

For general f: ∫ f dμ = ∫ f⁺ dμ - ∫ f⁻ dμ.

### Convergence Theorems

- **Monotone Convergence**: If 0 ≤ f₁ ≤ f₂ ≤ ... ↗ f, then ∫ fₙ dμ ↗ ∫ f dμ
- **Fatou's Lemma**: ∫ lim inf fₙ dμ ≤ lim inf ∫ fₙ dμ
- **Dominated Convergence**: If fₙ → f and |fₙ| ≤ g ∈ L¹, then ∫ fₙ → ∫ f

### Radon-Nikodym Theorem

If ν ≪ μ (ν absolutely continuous w.r.t. μ), then ∃ f ≥ 0 measurable such that:

$$\nu(A) = \int_A f \, d\mu \quad \text{for all measurable } A$$

The function f = dν/dμ is the **Radon-Nikodym derivative**. Bayes' rule is a special case.

### Fubini's Theorem

For f ∈ L¹(μ × ν):

$$\int_{X \times Y} f \, d(\mu \times \nu) = \int_X \left(\int_Y f(x,y) \, d\nu(y)\right) d\mu(x)$$

### Riesz Representation

Every positive linear functional L on C_c(X) has a unique representing measure:

$$L(f) = \int f \, d\mu$$

This connects functional analysis to measure theory.

---

## Running Tests

```bash
cargo test
```

134 tests: σ-algebra closure validation, measure additivity, Lebesgue measure on intervals/boxes, measurability verification, integral computation, convergence theorem verification, product measure correctness, Radon-Nikodym computation, decomposition verification, Riesz representation, pushforward correctness, and agent belief updates.

---

## License

MIT
