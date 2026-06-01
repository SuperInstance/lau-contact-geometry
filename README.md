# lau-contact-geometry

> Contact geometry: contact forms, Reeb dynamics, contactomorphisms, Legendrian submanifolds, and thermodynamic contact structures

## What This Does

Contact geometry: contact forms, Reeb dynamics, contactomorphisms, Legendrian submanifolds, and thermodynamic contact structures. Part of the PLATO/LAU ecosystem — a mathematically rigorous framework for building educational agents that learn, teach, and evolve.

## The Key Idea

This crate implements the core abstractions needed for its domain, with a focus on correctness, composability, and conservation guarantees. Every public type is serializable (serde), every algorithm is tested, and every invariant is verified.

## Install

```bash
cargo add lau-contact-geometry
```

## Quick Start

See the API Reference below for complete usage. Key entry points:

```rust
use lau_contact_geometry::*;
// See types and methods below for complete usage
```

## API Reference

```rust
pub struct ContactHamiltonianSystem 
pub enum HamiltonianType 
    pub fn constant(contact_form: ContactForm, h0: f64) -> Self 
    pub fn height(contact_form: ContactForm) -> Self 
    pub fn harmonic_plus_height(contact_form: ContactForm) -> Self 
    pub fn evaluate_h(&self, p: &DVector<f64>) -> f64 
    pub fn dh_dx(&self, p: &DVector<f64>) -> DVector<f64> 
    pub fn dh_dy(&self, p: &DVector<f64>) -> DVector<f64> 
    pub fn dh_dz(&self, _p: &DVector<f64>) -> f64 
    pub fn vector_field(&self, p: &DVector<f64>) -> DVector<f64> 
    pub fn step(&self, p: &DVector<f64>, dt: f64) -> DVector<f64> 
    pub fn integrate(&self, p0: &DVector<f64>, dt: f64, steps: usize) -> Vec<DVector<f64>> 
    pub fn conformal_factor(&self, p: &DVector<f64>) -> f64 
    pub fn energy_evolution_rate(&self, p: &DVector<f64>) -> f64 
pub struct ContactForm 
    pub fn new(dim: usize) -> Result<Self, String> 
    pub fn evaluate(&self, p: &DVector<f64>) -> DVector<f64> 
    pub fn exterior_derivative(&self) -> DMatrix<f64> 
    pub fn contact_volume(&self, p: &DVector<f64>) -> f64 
    pub fn pfaffian_squared(&self) -> f64 
    pub fn is_contact_at(&self, p: &DVector<f64>) -> bool 
    pub fn contact_distribution(&self, p: &DVector<f64>) -> DMatrix<f64> 
pub struct Contactomorphism 
    pub fn identity(dim: usize) -> Result<Self, String> 
    pub fn strict_from_jacobian(jacobian: DMatrix<f64>) -> Result<Self, String> 
    pub fn conformal(jacobian: DMatrix<f64>, factor: f64) -> Result<Self, String> 
    pub fn compose(&self, other: &Contactomorphism) -> Result<Contactomorphism, String> 
    pub fn preserves_contact(&self, _alpha: &ContactForm) -> bool 
    pub fn inverse(&self) -> Result<Contactomorphism, String> 
pub struct LegendrianSubmanifold 
    pub fn new(contact_form: ContactForm, name: &str) -> Self 
    pub fn zero_section(contact_form: ContactForm) -> Self 
    pub fn jet_section(contact_form: ContactForm, _f: fn(&DVector<f64>) -> f64) -> Self 
    pub fn is_on_zero_section(&self, p: &DVector<f64>) -> bool 
    pub fn is_in_contact_distribution(&self, p: &DVector<f64>, v: &DVector<f64>) -> bool 
    pub fn verify_legendrian(&self, gamma: &[DVector<f64>], tangents: &[DVector<f64>]) -> bool 
    pub fn front_projection(&self, p: &DVector<f64>) -> DVector<f64> 
    pub fn lagrangian_projection(&self, p: &DVector<f64>) -> DVector<f64> 
    pub fn legendrian_dimension(&self) -> usize 
pub struct ReebVectorField 
    pub fn new(contact_form: ContactForm) -> Self 
    pub fn evaluate(&self, _p: &DVector<f64>) -> DVector<f64> 
    pub fn verify_alpha_condition(&self, p: &DVector<f64>) -> bool 
    pub fn verify_da_condition(&self, p: &DVector<f64>) -> bool 
    pub fn verify_all(&self, p: &DVector<f64>) -> bool 
    pub fn flow(&self, p: &DVector<f64>, t: f64) -> DVector<f64> 
    pub fn orbit_period(&self) -> Option<f64> 
pub struct ThermodynamicContactStructure 
    pub fn standard() -> Self 
    pub fn extended(extra_pairs: usize) -> Self 
    pub fn contact_form_at(&self, p: &DVector<f64>) -> DVector<f64> 
    pub fn verify_first_law(&self, p: &DVector<f64>, du: f64, ds: f64, dv: f64) -> bool 
    pub fn verify_second_law(&self, entropy_production: f64) -> bool 
    pub fn is_equilibrium(&self, p: &DVector<f64>) -> bool 
    pub fn carnot_efficiency(&self, t_hot: f64, t_cold: f64) -> f64 
    pub fn landauer_bound(&self, temperature: f64) -> f64 
    pub fn agent_energy_cost(&self, temperature: f64, n_operations: u64) -> f64 
    pub fn entropy_representation(&self, u: f64, v: f64, s_fn: fn(f64, f64) -> f64, ds_du: fn(f64, f64) -> f64, ds_dv: fn(f64, f64) -> f64) -> DVector<f64> 
pub struct DarbouxCoordinates 
    pub fn new(dim: usize) -> Result<Self, String> 
```

## How It Works

Read the source in `src/` for full implementation details. All algorithms are documented with inline comments explaining the mathematical foundations.

## The Math

This crate implements formal mathematical constructs. See the source documentation for theorem statements and proofs of correctness.

## Testing

**64 tests** covering construction, serialization, correctness properties, edge cases, and composability with other lau-* crates.

## License

MIT
