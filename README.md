# lau-contact-geometry

Contact geometry is the odd-dimensional cousin of symplectic geometry. A contact form α on a (2n+1)-dimensional manifold satisfies α∧(dα)ⁿ ≠ 0 — it's maximally non-integrable. The kernel of α defines a hyperplane distribution that twists so much it can never be tangent to a hypersurface.

Contact geometry is where thermodynamics becomes geometry: entropy is the contact form, the first and second laws are its properties.

## The math in 60 seconds

A **contact manifold** (M, α) has dimension 2n+1. The **Reeb vector field** R is defined by α(R)=1 and dα(R,·)=0 — it's the unique direction the contact form doesn't "see." Key results:

- **Darboux theorem:** all contact structures look locally the same — α = dz + Σxᵢdyᵢ
- **Legendrian submanifolds:** n-dimensional submanifolds inside ker(α) — the "constrained" subspaces
- **Contact Hamiltonians:** equations of motion that preserve the contact structure
- **Contactomorphisms:** diffeomorphisms φ with φ*α = f·α — structure-preserving maps
- **Thermodynamic phase space:** (U, S, V, T, P) as a contact manifold with α = dU - TdS + PdV

References: Geiges, *An Introduction to Contact Topology* (2008)

## Quick start

```rust
use lau_contact_geometry::{ContactForm, ContactManifold, ReebField, ThermodynamicPhase};

// Standard contact structure on R³
let alpha = ContactForm::standard_r3();
assert!(alpha.verify_contact_condition()); // α∧(dα)¹ ≠ 0

// Reeb vector field
let reeb = ReebField::from_form(&alpha);

// Contact Hamilton's equations
let hamiltonian = 1.0;
let flow = reeb.contact_flow(hamiltonian, 100);

// Thermodynamic phase space
let phase = ThermodynamicPhase::ideal_gas();
let carnot = phase.carnot_efficiency(400.0, 300.0); // hot=400K, cold=300K
assert!((carnot - 0.25).abs() < 1e-10);

// Landauer bound: minimum energy to erase one bit
let landauer = phase.landauer_bound(300.0); // kT·ln(2)
```

## Key types

| Type | What it is |
|------|-----------|
| `ContactForm` | A 1-form α satisfying α∧(dα)ⁿ ≠ 0 |
| `ContactManifold` | A (2n+1)-manifold with contact structure |
| `ReebField` | The unique vector field R with α(R)=1, dα(R,·)=0 |
| `Contactomorphism` | Structure-preserving map φ*α = f·α |
| `Legendrian` | Maximal submanifold inside ker(α) |
| `ContactHamiltonian` | Equations of motion preserving the contact structure |
| `ThermodynamicPhase` | Thermodynamic quantities as contact geometry |

## Contributing

[Open an issue](https://github.com/SuperInstance/lau-contact-geometry/issues) or PR.
