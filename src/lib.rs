#![deny(unsafe_code)]
//! # lau-contact-geometry
//!
//! Contact geometry — the odd-dimensional cousin of symplectic geometry.
//!
//! A contact manifold is a (2n+1)-dimensional manifold equipped with a contact form α
//! satisfying the maximally non-integrable condition α ∧ (dα)ⁿ ≠ 0.
//!
//! This crate provides:
//! - Contact forms and contact manifolds
//! - Reeb vector field dynamics
//! - Contactomorphisms (structure-preserving maps)
//! - Darboux theorem (local normal forms)
//! - Legendrian submanifolds
//! - Contact Hamiltonian systems
//! - Thermodynamic contact geometry

pub mod contact_form;
pub mod contact_manifold;
pub mod contactomorphism;
pub mod darboux;
pub mod hamiltonian;
pub mod legendrian;
pub mod reeb;
pub mod thermodynamic;

pub use contact_form::ContactForm;
pub use contact_manifold::ContactManifold;
pub use contactomorphism::Contactomorphism;
pub use darboux::DarbouxCoordinates;
pub use hamiltonian::ContactHamiltonianSystem;
pub use legendrian::LegendrianSubmanifold;
pub use reeb::ReebVectorField;
pub use thermodynamic::ThermodynamicContactStructure;
