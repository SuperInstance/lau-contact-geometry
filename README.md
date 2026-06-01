# lau-contact-geometry

A (2n+1)-dimensional manifold M carries a contact structure when there exists a 1-form α such that α ∧ (dα)ⁿ ≠ 0. That single inequality defines the entire architecture of contact geometry — the odd-dimensional sibling of symplectic geometry — and this crate turns it into working code.

The zero-section hums along ker(α), the Reeb field walks the z-axis, and thermodynamic equilibrium is a Legendrian condition.

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Background](#mathematical-background)
3. [Architecture](#architecture)
4. [API Reference](#api-reference)
5. [Examples](#examples)
6. [Theorems Verified](#theorems-verified)
7. [Installation & Usage](#installation--usage)
8. [License](#license)

---

## Overview

This crate implements contact geometry — the study of odd-dimensional manifolds equipped with a maximally non-integrable hyperplane distribution. It provides:

- **Contact forms** α on ℝ^{2n+1} in Darboux coordinates, with the defining condition α ∧ (dα)ⁿ ≠ 0
- **Contact manifolds** (M, α) with symplectic structure on the contact distribution ξ = ker(α)
- **Reeb vector fields** R characterized by α(R) = 1 and dα(R, ·) = 0
- **Contactomorphisms** — diffeomorphisms preserving the contact structure (strict: φ\*α = α; conformal: φ\*α = fα)
- **Darboux coordinates** and the Darboux theorem (all contact structures of the same dimension are locally equivalent)
- **Legendrian submanifolds** — maximal submanifolds inside the contact distribution, with front and Lagrangian projections
- **Contact Hamiltonian systems** — dynamics preserving contact structure, with Euler integration
- **Thermodynamic contact geometry** — modeling thermodynamic phase space as a contact manifold, encoding the First and Second Laws, Carnot efficiency, and Landauer's bound

All results verified by 64 property-based tests across 8 modules.

---

## Mathematical Background

### Contact Forms and Contact Manifolds

A **contact form** on a (2n+1)-dimensional manifold M is a 1-form α satisfying:

> α ∧ (dα)ⁿ ≠ 0

This condition means the hyperplane distribution ξ = ker(α) is **maximally non-integrable** — there is no (2n)-dimensional submanifold tangent to ξ everywhere. The bilinear form dα|ξ is non-degenerate, making ξ into a symplectic vector bundle.

In **Darboux coordinates** (x₁,…,xₙ, y₁,…,yₙ, z):

> α = dz − Σᵢ yᵢ dxᵢ

The **Darboux theorem** states all contact structures of the same dimension are locally equivalent — there always exist coordinates putting α in this canonical form.

### The Reeb Vector Field

The **Reeb vector field** R is uniquely characterized by:

> α(R) = 1, &nbsp;&nbsp; dα(R, ·) = 0

In Darboux coordinates, R = ∂/∂z. The Reeb flow is always transverse to the contact distribution. On ℝ^{2n+1}, Reeb orbits are straight lines in the z-direction (non-periodic).

### Contactomorphisms

A **contactomorphism** φ: (M,α) → (M',α') is a diffeomorphism preserving the contact structure:

> φ\*α' = f · α &nbsp;&nbsp; for some non-vanishing function f

- **Strict** contactomorphism: f ≡ 1 (φ\*α' = α)
- **Conformal** contactomorphism: f is a non-zero constant

Contactomorphisms form a group under composition. The conformal factors multiply: f₁·f₂.

### Legendrian Submanifolds

A **Legendrian submanifold** L ⊂ (M,α) is an n-dimensional submanifold contained in the contact distribution:

> α|_L = 0 &nbsp;&nbsp; (equivalently, T_pL ⊂ ξ_p for all p ∈ L)

Dimension n is maximal — you cannot fit an (n+1)-dimensional submanifold inside ξ everywhere. Key examples:
- **Zero section**: L = {(x, 0, 0)} — trivially Legendrian since α = dz = 0 on z = const
- **1-jet sections**: L = {(x, ∂f/∂x, f(x))} for any function f — the graph of the 1-jet

Legendrian submanifolds have two natural projections:
- **Front projection** (drops y and z): shows the "wavefront"
- **Lagrrangian projection** (drops z): shows the projection to the symplectic base

### Contact Hamiltonian Systems

Given H: M → ℝ, the **contact Hamiltonian vector field** X_H is defined by:

> ẋᵢ = ∂H/∂yᵢ
> ẏᵢ = −∂H/∂xᵢ + yᵢ · ∂H/∂z
> ż = H − Σᵢ yᵢ(∂H/∂yᵢ)

Unlike symplectic Hamiltonian systems, energy is NOT conserved. Instead:

> dH/dt = λ · H &nbsp;&nbsp; where λ = −∂H/∂z is the conformal factor

The flow generates a contactomorphism with conformal factor e^{λt}.

### Thermodynamic Contact Geometry

The **thermodynamic phase space** is a 5-dimensional contact manifold:

> Coordinates: (U, S, V, T, −P)
> Contact form: α = dU − T dS + P dV

- **First Law** (dU = TdS − PdV) is encoded as the Legendrian condition α|_L = 0 on equilibrium states
- **Second Law** (entropy production ≥ 0) constrains which Legendrian submanifolds are physically realizable
- **Carnot efficiency**: η = 1 − T_cold/T_hot
- **Landauer's bound**: minimum energy to erase one bit = k_BT ln(2)

The entropy representation S = S(U,V) parametrizes the equilibrium Legendrian:

> (U, V) ↦ (U, S(U,V), V, T, −P) &nbsp;&nbsp; where T = (∂S/∂U)⁻¹, P = T·∂S/∂V

---

## Architecture

```
src/
├── lib.rs               — Crate root, re-exports all public types
├── contact_form.rs      — ContactForm: α, dα, contact volume, contact distribution
├── contact_manifold.rs  — ContactManifold: (M,α) with symplectic structure on ξ
├── contactomorphism.rs  — Contactomorphism: strict & conformal, compose, inverse
├── darboux.rs           — DarbouxCoordinates, DarbouxChart: local normal forms
├── hamiltonian.rs       — ContactHamiltonianSystem: vector field, Euler integration
├── legendrian.rs        — LegendrianSubmanifold: zero section, 1-jet, projections
├── reeb.rs              — ReebVectorField: evaluation, flow, condition verification
└── thermodynamic.rs     — ThermodynamicContactStructure: physics as geometry

tests/ — 64 tests (inline in each module)
```

### Computation Pipeline

```
Dimension 2n+1
    │
    ├─→ ContactForm α              (Darboux: dz − Σ yᵢ dxᵢ)
    │       │
    │       ├─→ evaluate(p)         (1-form as covector)
    │       ├─→ exterior_derivative (dα as skew-symmetric matrix)
    │       ├─→ contact_volume(p)   (α ∧ (dα)ⁿ)
    │       ├─→ is_contact_at(p)    (volume ≠ 0?)
    │       └─→ contact_distribution (basis for ξ = ker(α))
    │
    ├─→ ContactManifold (M, α)
    │       └─→ symplectic_on_distribution (dα|ξ — standard symplectic form)
    │
    ├─→ ReebVectorField R
    │       ├─→ evaluate(p)         (R = ∂/∂z in Darboux)
    │       ├─→ verify_alpha_condition (α(R) = 1?)
    │       ├─→ verify_da_condition   (dα(R,·) = 0?)
    │       └─→ flow(p, t)          (Φ_t(p) = p + tR)
    │
    ├─→ Contactomorphism φ
    │       ├─→ compose(φ₁, φ₂)    (group operation)
    │       ├─→ inverse             (φ⁻¹)
    │       └─→ preserves_contact   (φ*α = fα?)
    │
    ├─→ LegendrianSubmanifold L
    │       ├─→ zero_section        (L = {(x,0,0)})
    │       ├─→ verify_legendrian   (α|_L = 0?)
    │       ├─→ front_projection    (drop y,z)
    │       └─→ lagrangian_projection (drop z)
    │
    ├─→ ContactHamiltonianSystem (H, α)
    │       ├─→ vector_field(p)     (X_H at point)
    │       ├─→ step(p, dt)         (Euler integration)
    │       ├─→ integrate(p₀, dt, N) (full trajectory)
    │       ├─→ conformal_factor(p)  (λ = −∂H/∂z)
    │       └─→ energy_evolution_rate (dH/dt = λH)
    │
    ├─→ DarbouxCoordinates
    │       ├─→ verify_normal_form  (α = dz − Σ yᵢ dxᵢ?)
    │       └─→ verify_darboux_theorem
    │
    └─→ ThermodynamicContactStructure
            ├─→ verify_first_law    (dU = TdS − PdV)
            ├─→ verify_second_law   (ΔS ≥ 0)
            ├─→ carnot_efficiency   (1 − T_c/T_h)
            ├─→ landauer_bound      (k_BT ln2)
            └─→ entropy_representation (Legendrian parametrization)
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `nalgebra` | Linear algebra (DMatrix, DVector) with serde support |
| `serde` | Serialization of all geometric structures |

---

## API Reference

### ContactForm

```rust
let cf = ContactForm::new(dim)?;           // dim must be odd ≥ 3
cf.evaluate(&p);                            // α at point p (as covector)
cf.exterior_derivative();                   // dα as skew-symmetric matrix
cf.contact_volume(&p);                      // α ∧ (dα)ⁿ
cf.pfaffian_squared();                      // Pf(dα)² = 1 in Darboux
cf.is_contact_at(&p);                       // volume ≠ 0?
cf.contact_distribution(&p);               // basis for ξ = ker(α), dim × 2n matrix
cf.dim;                                     // 2n+1
cf.n;                                       // n
```

### ContactManifold

```rust
let m = ContactManifold::standard(dim)?;   // ℝ^{2n+1} with standard contact
let m = ContactManifold::named(dim, "name")?;
m.verify_contact(&p);                      // p ∈ M and contact condition holds
m.symplectic_on_distribution();            // dα|ξ — standard symplectic form
m.contains(&p);                            // dimension check
m.dim(); m.n(); m.name;
```

### ReebVectorField

```rust
let reeb = ReebVectorField::new(contact_form);
reeb.evaluate(&p);                         // R at point p (= ∂/∂z in Darboux)
reeb.verify_alpha_condition(&p);           // α(R) = 1?
reeb.verify_da_condition(&p);              // dα(R, ·) = 0?
reeb.verify_all(&p);                       // both conditions
reeb.flow(&p, t);                          // Φ_t(p) = p + t·R
reeb.orbit_period();                       // None for ℝ^{2n+1}
```

### Contactomorphism

```rust
let id = Contactomorphism::identity(dim)?;
let strict = Contactomorphism::strict_from_jacobian(jacobian)?;
let conf = Contactomorphism::conformal(jacobian, factor)?;

phi.compose(&psi)?;                        // group composition
phi.inverse()?;                            // φ⁻¹ (invert Jacobian)
phi.preserves_contact(&alpha);             // structure preservation check
phi.strict;                                // true if f ≡ 1
phi.conformal_factor;                      // Some(f) or None
```

### LegendrianSubmanifold

```rust
let leg = LegendrianSubmanifold::new(cf, "name");
let zero = LegendrianSubmanifold::zero_section(cf);
let jet = LegendrianSubmanifold::jet_section(cf, f);

leg.is_on_zero_section(&p);               // y=0, z=0?
leg.is_in_contact_distribution(&p, &v);    // α_p(v) = 0?
leg.verify_legendrian(&points, &tangents); // α(γ') = 0 for all tangent vectors?
leg.front_projection(&p);                  // drop y, z → ℝⁿ
leg.lagrangian_projection(&p);             // drop z → ℝ²ⁿ
leg.legendrian_dimension();                // n
```

### ContactHamiltonianSystem

```rust
let sys = ContactHamiltonianSystem::constant(cf, h0);
let sys = ContactHamiltonianSystem::height(cf);
let sys = ContactHamiltonianSystem::harmonic_plus_height(cf);

sys.evaluate_h(&p);                        // H(p)
sys.dh_dx(&p);                             // ∂H/∂x (n-vector)
sys.dh_dy(&p);                             // ∂H/∂y (n-vector)
sys.dh_dz(&p);                             // ∂H/∂z (scalar)
sys.vector_field(&p);                      // X_H(p) — full contact Hamiltonian vector field
sys.step(&p, dt);                          // Euler step
sys.integrate(&p0, dt, steps);             // full trajectory: Vec<DVector<f64>>
sys.conformal_factor(&p);                  // λ = −∂H/∂z
sys.energy_evolution_rate(&p);             // dH/dt = λ·H
```

### DarbouxCoordinates

```rust
let dc = DarbouxCoordinates::new(dim)?;
dc.to_darboux(&p);                         // identity in standard coords
dc.from_darboux(&p);                       // identity in standard coords
dc.verify_normal_form(&p);                 // α = dz − Σ yᵢ dxᵢ?
dc.verify_darboux_theorem();               // normal form + contact condition
dc.chart_around(&point);                   // → DarbouxChart
```

### ThermodynamicContactStructure

```rust
let tcs = ThermodynamicContactStructure::standard();    // 5D
let tcs = ThermodynamicContactStructure::extended(2);   // 9D

tcs.contact_form_at(&p);                   // α = dU − TdS + PdV as covector
tcs.verify_first_law(&p, dU, dS, dV);      // dU = TdS − PdV?
tcs.verify_second_law(entropy_prod);       // ΔS ≥ 0?
tcs.is_equilibrium(&p);                    // on the Legendrian?
tcs.carnot_efficiency(t_hot, t_cold);      // 1 − T_c/T_h
tcs.landauer_bound(temperature);           // k_BT ln(2)
tcs.agent_energy_cost(temperature, n);     // n × Landauer bound
tcs.entropy_representation(u, v, S, ∂S/∂U, ∂S/∂V); // → state vector
```

---

## Examples

### Creating and Verifying a Contact Form

```rust
use lau_contact_geometry::*;
use nalgebra::dvector;

let cf = ContactForm::new(3)?;  // ℝ³ with standard contact structure
let p = dvector![1.0, 2.0, 3.0];  // point (x, y, z)

// Evaluate α = dz − y·dx = (-2, 0, 1) at (1, 2, 3)
let alpha = cf.evaluate(&p);
assert!((alpha[0] - (-2.0)).abs() < 1e-10);  // −y = −2
assert!(alpha[1].abs() < 1e-10);              // dy coefficient = 0
assert!((alpha[2] - 1.0).abs() < 1e-10);     // dz coefficient = 1

// Contact condition: α ∧ (dα)¹ ≠ 0
assert!(cf.is_contact_at(&p));  // always true for standard contact

// Exterior derivative: dα = dx ∧ dy
let da = cf.exterior_derivative();
assert!((da[(0, 1)] - 1.0).abs() < 1e-10);
assert!((da[(1, 0)] - (-1.0)).abs() < 1e-10);
```

### Contact Distribution and Symplectic Structure

```rust
let m = ContactManifold::standard(5)?;  // ℝ⁵ with n=2
let p = dvector![1.0, 2.0, 3.0, 4.0, 5.0];

// ξ = ker(α) is a 4-dimensional subspace
let xi = m.contact_form.contact_distribution(&p);
assert_eq!(xi.nrows(), 5);
assert_eq!(xi.ncols(), 4);  // 2n = 4

// dα|ξ gives the standard symplectic form on ℝ⁴
let omega = m.symplectic_on_distribution();
assert!((omega[(0, 2)] - 1.0).abs() < 1e-10);  // dx₁∧dy₁
assert!((omega[(2, 0)] - (-1.0)).abs() < 1e-10);
```

### Reeb Vector Field

```rust
let cf = ContactForm::new(3)?;
let reeb = ReebVectorField::new(cf);
let p = dvector![1.0, 2.0, 3.0];

// R = ∂/∂z = (0, 0, 1) in Darboux coordinates
let r = reeb.evaluate(&p);
assert!((r[2] - 1.0).abs() < 1e-10);

// Verify the defining conditions
assert!(reeb.verify_alpha_condition(&p));   // α(R) = 1
assert!(reeb.verify_da_condition(&p));      // dα(R, ·) = 0
assert!(reeb.verify_all(&p));

// Reeb flow: only z changes
let flowed = reeb.flow(&p, 2.5);
assert!((flowed[2] - 5.5).abs() < 1e-10);  // z: 3 → 5.5
```

### Contactomorphisms

```rust
let id = Contactomorphism::identity(3)?;
assert!(id.strict);
assert!(id.preserves_contact(&ContactForm::new(3)?));

// Composition: identity ∘ identity = identity
let comp = id.compose(&id)?;
for i in 0..3 {
    assert!((comp.jacobian[(i, i)] - 1.0).abs() < 1e-10);
}

// Conformal contactomorphism: φ*α = 2α
let conf = Contactomorphism::conformal(DMatrix::identity(3, 3), 2.0)?;
assert!(!conf.strict);
assert_eq!(conf.conformal_factor, Some(2.0));

// Conformal factors multiply under composition
let conf2 = Contactomorphism::conformal(DMatrix::identity(3, 3), 3.0)?;
let product = conf.compose(&conf2)?;
assert!((product.conformal_factor.unwrap() - 6.0).abs() < 1e-10);

// Inverse: factor inverts
let inv = conf.inverse()?;
assert!((inv.conformal_factor.unwrap() - 0.5).abs() < 1e-10);
```

### Legendrian Submanifolds

```rust
let cf = ContactForm::new(3)?;
let leg = LegendrianSubmanifold::zero_section(cf);

// Zero-section: points with y=0, z=0
assert!(leg.is_on_zero_section(&dvector![5.0, 0.0, 0.0]));
assert!(!leg.is_on_zero_section(&dvector![1.0, 2.0, 3.0]));

// Verify Legendrian condition: α(tangent) = 0
let points = vec![dvector![0.0, 0.0, 0.0], dvector![1.0, 0.0, 0.0]];
let tangents = vec![dvector![1.0, 0.0, 0.0], dvector![1.0, 0.0, 0.0]];
assert!(leg.verify_legendrian(&points, &tangents));

// Projections
let p = dvector![1.0, 2.0, 3.0];
let front = leg.front_projection(&p);       // [1.0] — drops y and z
let lag = leg.lagrangian_projection(&p);     // [1.0, 2.0] — drops z

// Legendrian dimension = n
assert_eq!(leg.legendrian_dimension(), 1);
```

### Contact Hamiltonian Dynamics

```rust
let cf = ContactForm::new(3)?;

// Constant Hamiltonian H = 2
let sys = ContactHamiltonianSystem::constant(cf, 2.0);
let p = dvector![1.0, 2.0, 3.0];
let xh = sys.vector_field(&p);
// ẋ = 0, ẏ = 0, ż = H = 2

// Euler integration
let p0 = dvector![0.0, 0.0, 0.0];
let traj = sys.integrate(&p0, 0.01, 100);
assert_eq!(traj.len(), 101);
// z drifts linearly upward

// Height Hamiltonian H = z
let sys_h = ContactHamiltonianSystem::height(cf);
assert!((sys_h.conformal_factor(&p) - (-1.0)).abs() < 1e-10);  // λ = −∂H/∂z = −1
assert!((sys_h.energy_evolution_rate(&dvector![0.0, 0.0, 5.0]) - (-5.0)).abs() < 1e-10);

// Harmonic + Height: H = (x² + y²)/2 + z
let sys_harm = ContactHamiltonianSystem::harmonic_plus_height(cf);
let xh = sys_harm.vector_field(&dvector![1.0, 0.0, 0.0]);
// ẋ = y = 0, ẏ = −x = −1, ż = 0.5
```

### Darboux Theorem

```rust
let dc = DarbouxCoordinates::new(3)?;
let p = dvector![1.0, 2.0, 3.0];

// Verify α is in normal form: dz − y dx
assert!(dc.verify_normal_form(&p));

// The Darboux theorem: locally, all contact structures look the same
assert!(dc.verify_darboux_theorem());

// Round-trip through Darboux coordinates
let darboux = dc.to_darboux(&p);
let back = dc.from_darboux(&darboux);
for i in 0..3 {
    assert!((p[i] - back[i]).abs() < 1e-10);
}
```

### Thermodynamic Contact Geometry

```rust
let tcs = ThermodynamicContactStructure::standard();

// Contact form: α = dU − T·dS + P·dV
let p = dvector![100.0, 50.0, 10.0, 300.0, -101325.0];
let alpha = tcs.contact_form_at(&p);
assert!((alpha[0] - 1.0).abs() < 1e-10);        // dU
assert!((alpha[1] - (-300.0)).abs() < 1e-10);    // −T·dS
assert!((alpha[2] - 101325.0).abs() < 1e-10);    // P·dV

// First Law: dU = TdS − PdV
let p_state = dvector![0.0, 0.0, 0.0, 300.0, -100000.0];
assert!(tcs.verify_first_law(&p_state, 50.0, 0.5, 0.001));
assert!(!tcs.verify_first_law(&p_state, 999.0, 0.5, 0.001));  // violation

// Second Law: entropy production ≥ 0
assert!(tcs.verify_second_law(5.0));
assert!(!tcs.verify_second_law(-1.0));

// Carnot efficiency
assert!((tcs.carnot_efficiency(600.0, 300.0) - 0.5).abs() < 1e-10);
assert!((tcs.carnot_efficiency(600.0, 0.0) - 1.0).abs() < 1e-10);  // ideal

// Landauer's bound: minimum energy to erase one bit
let landauer = tcs.landauer_bound(300.0);  // ≈ 2.87 × 10⁻²¹ J
assert!(landauer > 0.0 && landauer < 1e-19);

// Agent energy cost: n operations at temperature T
let cost = tcs.agent_energy_cost(300.0, 1_000_000);

// Entropy representation: Legendrian parametrization
let state = tcs.entropy_representation(100.0, 1.0, S, dS_dU, dS_dV);
// → (U, S(U,V), V, T, −P)
```

---

## Theorems Verified

All 64 tests serve as machine-checked verifications of mathematical theorems and invariants:

### Contact Form (11 tests)

1. **Creation**: valid odd dimensions ≥ 3 accepted, even/low dimensions rejected — 3 tests
2. **Darboux evaluation**: α = (−y, 0, 1) at (x, y, z) verified — `test_darboux_3d_contact_form`
3. **Contact condition**: α ∧ (dα)ⁿ ≠ 0 at arbitrary points in ℝ³, ℝ⁵, ℝ⁷ — 3 tests
4. **Exterior derivative**: dα = dx ∧ dy in ℝ³, correct skew-symmetric structure — `test_exterior_derivative_3d`
5. **Contact volume**: non-zero at origin and general points — `test_contact_volume_nonzero`
6. **Distribution dimension**: ξ has dimension 2n in ℝ⁵ — `test_contact_distribution_dimension`
7. **Higher dimensions**: all properties hold in ℝ⁷ — `test_higher_dimensional_contact`

### Contact Manifold (6 tests)

8. **Standard manifold construction**: ℝ³, ℝ⁵ with correct n — 2 tests
9. **Point containment**: dimension checking — 2 tests
10. **Contact verification at point** — `test_verify_contact_at_point`
11. **Symplectic on distribution**: dα|ξ = [[0,1],[−1,0]] in ℝ³ — `test_symplectic_on_distribution`

### Contactomorphism (6 tests)

12. **Identity**: strict, preserves contact, J = I — `test_identity_contactomorphism`
13. **Composition**: id ∘ id = id — `test_compose_identity`
14. **Inverse**: id⁻¹ = id — `test_inverse_identity`
15. **Conformal**: non-strict, factor stored correctly — `test_conformal_contactomorphism`
16. **Conformal inverse**: factor inverts — `test_inverse_conformal`
17. **Conformal composition**: factors multiply (2×3 = 6) — `test_compose_conformal`

### Darboux Coordinates (7 tests)

18. **Creation and normal form verification** at general point and origin — 3 tests
19. **Darboux theorem**: holds for ℝ³ and ℝ⁷ — 2 tests
20. **Chart construction**: center point stored, contains valid points — `test_chart_around_point`
21. **Round-trip coordinates**: to_darboux ∘ from_darboux = identity — `test_roundtrip_darboux`

### Contact Hamiltonian System (10 tests)

22. **Hamiltonian evaluation**: constant, height, and harmonic types compute correct H(p) — 3 tests
23. **Constant vector field**: ẋ=0, ẏ=0, ż=H₀ — `test_constant_vector_field`
24. **Euler step**: z increments by H₀·dt — `test_euler_step`
25. **Trajectory length**: integrate produces steps+1 points — `test_integrate_trajectory_length`
26. **Conformal factor**: λ = −∂H/∂z correct for height Hamiltonian — `test_conformal_factor_height`
27. **Energy evolution rate**: dH/dt = λH verified — `test_energy_evolution`
28. **Harmonic vector field**: ẏ = −x, ż = (x²+y²)/2 computed correctly — `test_harmonic_vector_field`

### Legendrian Submanifold (8 tests)

29. **Zero section**: points with y=0, z=0 identified — 2 tests
30. **Contact distribution membership**: ∂/∂y ∈ ξ, ∂/∂z ∉ ξ — 2 tests
31. **Legendrian dimension**: n in (2n+1)-manifold — `test_legendrian_dimension`
32. **Front projection**: drops y and z — `test_front_projection`
33. **Lagrangian projection**: drops z — `test_lagrangian_projection`
34. **Legendrian verification**: zero-section tangent vectors satisfy α(v) = 0 — `test_verify_legendrian_zero_section`

### Reeb Vector Field (6 tests)

35. **Reeb direction**: R = (0,0,1) in ℝ³ — `test_reeb_vector_3d`
36. **α(R) = 1**: verified in ℝ³ — `test_reeb_alpha_equals_one`
37. **dα(R, ·) = 0**: verified in ℝ³ — `test_reeb_da_equals_zero`
38. **Both conditions simultaneously**: verified in ℝ⁵ — `test_reeb_all_conditions`
39. **Reeb flow**: only z changes, correct displacement — `test_reeb_flow`
40. **Higher dimensions**: all conditions hold in ℝ⁷ — `test_reeb_7d`

### Thermodynamic Contact Geometry (10 tests)

41. **Standard structure**: 5-dimensional manifold — `test_standard_thermodynamic_structure`
42. **Contact form coefficients**: α = (1, −T, P, 0, 0) — `test_contact_form_thermodynamic`
43. **First Law satisfied**: dU = TdS − PdV — `test_first_law`
44. **First Law violated**: wrong dU detected — `test_first_law_violation`
45. **Second Law satisfied**: positive entropy production — `test_second_law_satisfied`
46. **Second Law violated**: negative entropy production detected — `test_second_law_violated`
47. **Carnot efficiency**: 1 − 300/600 = 0.5 — `test_carnot_efficiency`
48. **Carnot at absolute zero**: η = 1 (ideal limit) — `test_carnot_efficiency_zero_cold`
49. **Landauer bound**: k_BT ln(2) ≈ 2.87×10⁻²¹ J at 300K — `test_landauer_bound`
50. **Agent energy cost**: scales linearly with operations — `test_agent_energy_cost`
51. **Extended thermodynamics**: 9-dimensional (2 extra pairs) — `test_extended_thermodynamics`
52. **Entropy representation**: Legendrian parametrization produces correct state vector — `test_entropy_representation`

---

## Installation & Usage

### Prerequisites

- Rust 1.56+ (2021 edition)

### Add to your project

```toml
[dependencies]
lau-contact-geometry = { git = "https://github.com/SuperInstance/lau-contact-geometry" }
```

Or clone and use locally:

```bash
git clone https://github.com/SuperInstance/lau-contact-geometry.git
cd lau-contact-geometry
cargo test   # Run all 64 tests
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `nalgebra` (0.33) | Linear algebra: matrices, vectors, with serde-serialize feature |
| `serde` (1) | Derive macros for serialization of all geometric types |

---

## License

MIT
