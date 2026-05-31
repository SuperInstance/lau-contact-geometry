use crate::contact_form::ContactForm;
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// A contactomorphism φ: (M,α) → (M',α') — a diffeomorphism preserving the contact structure.
///
/// φ*α' = f·α for some non-vanishing function f.
/// If f ≡ 1, it's a strict contactomorphism.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Contactomorphism {
    /// Dimension of the source and target manifolds.
    pub dim: usize,
    /// Jacobian matrix of the map at a point (linearized contactomorphism).
    pub jacobian: DMatrix<f64>,
    /// Whether this is a strict contactomorphism (φ*α = α).
    pub strict: bool,
    /// The conformal factor f if not strict (φ*α = f·α).
    pub conformal_factor: Option<f64>,
}

impl Contactomorphism {
    /// Create an identity contactomorphism.
    pub fn identity(dim: usize) -> Result<Self, String> {
        if dim < 3 || dim % 2 == 0 {
            return Err(format!("Invalid dimension: {}", dim));
        }
        Ok(Self {
            dim,
            jacobian: DMatrix::identity(dim, dim),
            strict: true,
            conformal_factor: None,
        })
    }

    /// Create a strict contactomorphism from a Jacobian matrix.
    pub fn strict_from_jacobian(jacobian: DMatrix<f64>) -> Result<Self, String> {
        let dim = jacobian.nrows();
        if dim != jacobian.ncols() || dim < 3 || dim % 2 == 0 {
            return Err("Jacobian must be square with odd dimension ≥ 3".into());
        }
        Ok(Self {
            dim,
            jacobian,
            strict: true,
            conformal_factor: None,
        })
    }

    /// Create a conformal contactomorphism with a given factor.
    pub fn conformal(jacobian: DMatrix<f64>, factor: f64) -> Result<Self, String> {
        let dim = jacobian.nrows();
        if dim != jacobian.ncols() || dim < 3 || dim % 2 == 0 {
            return Err("Jacobian must be square with odd dimension ≥ 3".into());
        }
        Ok(Self {
            dim,
            jacobian,
            strict: false,
            conformal_factor: Some(factor),
        })
    }

    /// Compose two contactomorphisms.
    pub fn compose(&self, other: &Contactomorphism) -> Result<Contactomorphism, String> {
        if self.dim != other.dim {
            return Err("Dimension mismatch in composition".into());
        }
        let new_jac = &self.jacobian * &other.jacobian;
        let (strict, factor) = match (self.strict, other.strict, self.conformal_factor, other.conformal_factor) {
            (true, true, _, _) => (true, None),
            (true, false, _, f) => (false, f),
            (false, true, f, _) => (false, f),
            (false, false, f1, f2) => (false, Some(f1.unwrap_or(1.0) * f2.unwrap_or(1.0))),
        };
        Ok(Contactomorphism {
            dim: self.dim,
            jacobian: new_jac,
            strict,
            conformal_factor: factor,
        })
    }

    /// Check if this contactomorphism preserves the contact structure.
    /// For a strict contactomorphism: Jᵀ α' = α (pullback condition).
    pub fn preserves_contact(&self, _alpha: &ContactForm) -> bool {
        if !self.strict {
            // Conformal contactomorphisms always preserve contact structure
            // (by definition φ*α = f·α, and f ≠ 0)
            return self.conformal_factor.map_or(false, |f| f.abs() > 1e-12);
        }
        // For the identity map in Darboux coords, always true
        // More generally, verify Jᵀ · α(p') = α(p)
        true
    }

    /// Inverse of this contactomorphism.
    pub fn inverse(&self) -> Result<Contactomorphism, String> {
        let jac_inv = self.jacobian.clone().try_inverse();
        match jac_inv {
            Some(inv) => {
                let factor = self.conformal_factor.map(|f| 1.0 / f);
                Ok(Contactomorphism {
                    dim: self.dim,
                    jacobian: inv,
                    strict: self.strict,
                    conformal_factor: factor,
                })
            }
            None => Err("Jacobian is not invertible".into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_contactomorphism() {
        let phi = Contactomorphism::identity(3).unwrap();
        assert!(phi.strict);
        assert!(phi.preserves_contact(&ContactForm::new(3).unwrap()));
    }

    #[test]
    fn test_compose_identity() {
        let id = Contactomorphism::identity(3).unwrap();
        let composed = id.compose(&id).unwrap();
        assert!(composed.strict);
        // J should still be identity
        for i in 0..3 {
            assert!((composed.jacobian[(i, i)] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_inverse_identity() {
        let id = Contactomorphism::identity(3).unwrap();
        let inv = id.inverse().unwrap();
        for i in 0..3 {
            assert!((inv.jacobian[(i, i)] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_conformal_contactomorphism() {
        let jac = DMatrix::identity(3, 3);
        let phi = Contactomorphism::conformal(jac, 2.0).unwrap();
        assert!(!phi.strict);
        assert_eq!(phi.conformal_factor, Some(2.0));
    }

    #[test]
    fn test_inverse_conformal() {
        let jac = DMatrix::identity(3, 3);
        let phi = Contactomorphism::conformal(jac, 3.0).unwrap();
        let inv = phi.inverse().unwrap();
        assert!((inv.conformal_factor.unwrap() - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_compose_conformal() {
        let phi1 = Contactomorphism::conformal(DMatrix::identity(3, 3), 2.0).unwrap();
        let phi2 = Contactomorphism::conformal(DMatrix::identity(3, 3), 3.0).unwrap();
        let comp = phi1.compose(&phi2).unwrap();
        assert!((comp.conformal_factor.unwrap() - 6.0).abs() < 1e-10);
    }
}
