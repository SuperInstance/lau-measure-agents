# lau-measure-agents

> Rigorous measure theory: sigma-algebras, Lebesgue integration, Radon-Nikodym derivatives, and convergence theorems.

## What This Does

This crate implements the foundations of measure theory — the mathematical framework underlying integration, probability, and functional analysis. It covers sigma-algebras and measurable spaces, measures (counting, Dirac, probability, Lebesgue), measurable functions, the Lebesgue integral with all three major convergence theorems (Monotone Convergence, Fatou's Lemma, Dominated Convergence), product measures and Fubini's theorem, the Radon-Nikodym theorem, Lebesgue decomposition, Riesz representation, pushforward measures, and an agent application layer that models observations as measurable functions and belief updates via Radon-Nikodym derivatives.

Use this when you need mathematically rigorous integration beyond Riemann, or when building probabilistic agent systems that require proper measure-theoretic foundations.

## The Key Idea

Measure theory generalizes "length," "area," and "volume" to arbitrary sets. A **measure** μ assigns non-negative numbers to sets in a sigma-algebra, with countable additivity for disjoint unions. The **Lebesgue integral** extends integration to functions that Riemann can't handle — it works by partitioning the *range* of the function rather than the domain. This gives you powerful convergence theorems: under mild conditions, the limit of integrals equals the integral of the limit.

## Install

```bash
cargo add lau-measure-agents
```

## Quick Start

```rust
use lau_measure_agents::*;
use nalgebra::DVector;

fn main() {
    // Build a sigma-algebra (power set) on {a, b, c}
    let universe = Subset::from_slice(&["a", "b", "c"]);
    let algebra = SigmaAlgebra::power_set(universe);

    // Create a probability measure
    let mu = Measure::probability(&algebra, &[
        ("a".into(), 0.5), ("b".into(), 0.3), ("c".into(), 0.2),
    ]).unwrap();
    println!("μ({{a}}) = {}", mu.measure_of(&Subset::from_slice(&["a"])));
    println!("Is probability: {}", mu.is_probability());

    // Integrate a function: f(a)=2, f(b)=4, f(c)=6
    let f = RealValuedFunction::new(vec![("a", 2.0), ("b", 4.0), ("c", 6.0)]);
    let integral = integrate(&f, &mu);
    println!("∫ f dμ = {}", integral); // 0.5*2 + 0.3*4 + 0.2*6 = 3.4

    // Lebesgue measure on R
    let leb = LebesgueMeasure::new(1);
    let interval = Interval::new(0.0, 5.0).unwrap();
    println!("λ([0,5]) = {}", leb.measure_interval(&interval)); // 5.0

    // Simple function integration
    let sf = SimpleIntervalFunction::new(vec![
        (Interval::new(0.0, 1.0).unwrap(), 3.0),
        (Interval::new(1.0, 4.0).unwrap(), 2.0),
    ]);
    println!("∫ f dλ = {}", sf.integrate()); // 3*1 + 2*3 = 9
}
```

## API Reference

### Sigma-Algebras

#### `Subset`
A finite subset of a universal set.

```rust
let s = Subset::from_slice(&["a", "b", "c"]);
let empty = Subset::empty();

s.contains("a");
s.len();
s.union(&other);
s.intersection(&other);
s.difference(&other);
s.complement(&universe);
s.is_subset_of(&other);
```

#### `SigmaAlgebra`
A collection of measurable sets closed under complement and countable union.

```rust
// Construction
SigmaAlgebra::trivial(universe)           // {∅, Ω}
SigmaAlgebra::power_set(universe)          // all subsets
SigmaAlgebra::generate(universe, &[seed]) // smallest σ-algebra containing seeds

// Operations
sa.is_measurable(&set);
sa.countable_union(&[s1, s2]);
sa.countable_intersection(&[s1, s2]);
sa.complement(&set);
sa.measurable_sets();
sa.universal();
sa.size();
sa.validate();  // check σ-algebra axioms
```

#### `MeasurableSpace`
A set equipped with a sigma-algebra.

```rust
let ms = MeasurableSpace::new(algebra);
ms.sigma_algebra();
ms.universal();
```

### Measures

#### `Measure`
A non-negative countably additive set function.

```rust
// Construction
Measure::new("my_mu", &algebra, &[(set, value)])?;
Measure::counting(&algebra);                        // μ(A) = |A|
Measure::dirac(&algebra, "a")?;                     // point mass at "a"
Measure::probability(&algebra, &[("a", 0.5)])?;     // total mass = 1
Measure::uniform(&algebra)?;                        // equal weights

// Evaluation
mu.measure_of(&set);
mu.total_mass();
mu.is_probability();
mu.is_finite();

// Relationships
mu.is_absolutely_continuous_wrt(&nu);  // ν ≪ μ?
mu.is_singular_with(&nu);              // μ ⊥ ν?

// Arithmetic
mu.scale(2.0);
mu.add(&other);
mu.subtract(&other);
```

#### `SignedMeasure`
A measure that can take negative values.

```rust
let sm = SignedMeasure::new("signed", &[(set, -1.0)]);
sm.measure_of(&set);
sm.total_variation();              // |μ|(Ω)
let (pos, neg) = sm.jordan_decomposition();  // μ = μ⁺ - μ⁻
```

### Lebesgue Measure

#### `Interval`
A closed interval [a, b] on the real line.

```rust
let i = Interval::new(0.0, 1.0)?;
i.length();
i.contains(0.5);
i.intersection(&other);
i.union(&other);
i.complement(0.0, 10.0);
```

#### `Box`
A product of intervals in R^n.

```rust
let b = Box::new(vec![Interval::new(0.0, 1.0)?, Interval::new(0.0, 2.0)?]);
b.volume();          // 2.0
b.dimension();       // 2
b.contains(&point);
b.intersection(&other);
```

#### `LebesgueMeasure`
Standard volume measure on R^n.

```rust
let m = LebesgueMeasure::new(2);  // R^2
m.measure_box(&b);
m.measure_interval(&i);                    // 1D only
m.measure_disjoint_intervals(&intervals);
m.outer_measure_approx(&points, epsilon);
```

#### `SimpleIntervalFunction`
A finite linear combination of indicator functions.

```rust
let f = SimpleIntervalFunction::new(vec![
    (Interval::new(0.0, 1.0)?, 3.0),
    (Interval::new(1.0, 4.0)?, 1.0),
]);
f.eval(0.5);       // 3.0
f.integrate();      // 3*1 + 1*3 = 6.0
f.supremum();       // 3.0
```

### Measurable Functions

#### `MeasurableFunction`
A function between measurable spaces where preimages of measurable sets are measurable.

```rust
let mf = MeasurableFunction::new("f", &domain, &codomain, mapping);
mf.is_measurable();
mf.apply(&element);
mf.preimage(&measurable_set);
mf.compose(&other);
```

### Lebesgue Integral

#### `RealValuedFunction`
A function on a finite space.

```rust
let f = RealValuedFunction::new(vec![("a", 1.0), ("b", 2.0)]);
f.eval("a");
f.supremum();
f.infimum();
f.pointwise_max(&g);
f.pointwise_min(&g);
f.abs();
```

#### Integration

```rust
integrate(&f, &mu);                       // ∫ f dμ
integrate_over(&f, &mu, &set);            // ∫_A f dμ
```

#### Convergence Theorems

```rust
// Monotone Convergence: f_n ↑ f ⟹ ∫f_n → ∫f
let (integrals, limit) = monotone_convergence(&functions, &mu);

// Fatou's Lemma: lim inf ∫f_n ≤ ∫(lim inf f_n)
let (lim_inf_integrals, integral_of_lim_inf) = fatou_lemma(&functions, &mu);

// Dominated Convergence: |f_n| ≤ g, f_n → f ⟹ ∫f_n → ∫f
let (integrals, limit, converged) = dominated_convergence(&functions, &mu, &dominator);
```

### Product Measures

```rust
let (product, product_sa) = product_measure(&mu1, &sa1, &mu2, &sa2);

// Fubini's theorem
let result = fubini(&f, &mu1, &mu2);
// ∫∫ f d(μ₁×μ₂) = ∫ [∫ f dμ₂] dμ₁
```

### Radon-Nikodym Theorem

```rust
// If ν ≪ μ, find dν/dμ
let density = radon_nikodym_derivative(&nu, &mu);
// Verify: ν(A) = ∫_A (dν/dμ) dμ for all measurable A
verify_radon_nikodym(&nu, &mu, &density);
```

### Lebesgue Decomposition

```rust
// Decompose ν = ν_ac + ν_singular
let (abs_continuous, singular) = lebesgue_decomposition(&nu, &mu);
// ν_ac ≪ μ, ν_singular ⊥ μ
```

### Riesz Representation

```rust
// Every positive linear functional on C(X) is an integral against a measure
let (representing_measure, integral_value) = riesz_representation(&functional, &space);
```

### Pushforward Measure

```rust
// μ∘f⁻¹: the measure induced on the codomain by a measurable map
let pushfwd = pushforward(&mu, &f);
// pushfwd(B) = μ(f⁻¹(B))
```

### Agent Application

```rust
// Model observations as measurable functions
let obs = Observation::new(measurable_function);

// Belief update via Radon-Nikodym
let updated_belief = bayesian_update(&prior, &observation, &evidence);
// posterior density ∝ likelihood × prior density
```

## How It Works

**Sigma-algebras** are built by enumerating all subsets of a finite universal set (power set) or generated from seed sets by closing under complement and union. The `validate()` method checks the three axioms: contains ∅, closed under complement, closed under countable union.

**Measures** are stored as hash maps from set keys to values. Counting, Dirac, probability, and uniform measures are pre-built constructors. Absolute continuity (ν ≪ μ) is checked by verifying μ(A) = 0 ⟹ ν(A) = 0. Singularity checks for disjoint supports.

**Lebesgue integration** on finite spaces reduces to Σ f(xᵢ)μ({xᵢ}). For real-valued functions on R, simple interval functions approximate integrands, with integral computed as Σ aᵢ·λ(Aᵢ).

**Convergence theorems** are demonstrated by evaluating sequences of functions:
- **MCT**: monotone increasing non-negative functions, integrals converge upward.
- **Fatou**: lim inf of integrals ≤ integral of lim inf.
- **DCT**: dominated sequences converge in integral.

**Product measures** are constructed as Cartesian products of sigma-algebras with μ₁×μ₂(A×B) = μ₁(A)·μ₂(B). Fubini's theorem iterates integrals.

**Radon-Nikodym** computes dν/dμ element-wise for finite spaces: density(x) = ν({x})/μ({x}).

**Lebesgue decomposition** separates ν into absolutely continuous and singular parts with respect to μ.

## The Math

### Sigma-Algebra

A collection 𝔉 of subsets of Ω satisfying:
1. ∅ ∈ 𝔉
2. A ∈ 𝔉 ⟹ Aᶜ ∈ 𝔉
3. A₁, A₂, ... ∈ 𝔉 ⟹ ∪Aᵢ ∈ 𝔉

### Measure

A function μ: 𝔉 → [0, ∞] with:
1. μ(∅) = 0
2. **Countable additivity**: Aᵢ disjoint ⟹ μ(∪Aᵢ) = Σμ(Aᵢ)

### Lebesgue Integral

For a non-negative simple function φ = Σ aᵢ·1_{Aᵢ}:

$$\int \varphi \, d\mu = \sum a_i \, \mu(A_i)$$

For general f ≥ 0:

$$\int f \, d\mu = \sup\left\{\int \varphi \, d\mu : \varphi \leq f, \varphi \text{ simple}\right\}$$

### Monotone Convergence Theorem

If 0 ≤ f₁ ≤ f₂ ≤ ... and fₙ → f pointwise, then:

$$\lim_{n \to \infty} \int f_n \, d\mu = \int f \, d\mu$$

### Fatou's Lemma

$$\int \liminf_{n} f_n \, d\mu \leq \liminf_{n} \int f_n \, d\mu$$

### Dominated Convergence Theorem

If fₙ → f pointwise and |fₙ| ≤ g with ∫g dμ < ∞, then:

$$\lim_{n \to \infty} \int f_n \, d\mu = \int f \, d\mu$$

### Radon-Nikodym Theorem

If ν ≪ μ, there exists a unique (up to μ-a.e.) density f = dν/dμ such that:

$$\nu(A) = \int_A f \, d\mu \quad \forall A \in \mathfrak{F}$$

### Lebesgue Decomposition

Any σ-finite ν can be uniquely written as ν = ν_ac + ν_s where ν_ac ≪ μ and ν_s ⊥ μ.

### Fubini's Theorem

If f is integrable on (X × Y, μ × ν):

$$\int_{X \times Y} f \, d(\mu \times \nu) = \int_X \left[\int_Y f(x,y) \, d\nu(y)\right] d\mu(x)$$

### Riesz Representation Theorem

Every positive linear functional Λ on C_c(X) corresponds to a unique regular Borel measure μ:

$$\Lambda(f) = \int f \, d\mu$$

## License

MIT
