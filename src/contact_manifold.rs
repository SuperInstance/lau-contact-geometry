use crate::contact_form::ContactForm;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// A contact manifold (M, α) — a (2n+1)-dimensional manifold with a contact form.
///
/// The fundamental structure is the hyperplane distribution ξ = ker(α)
/// together with the symplectic form dα|ᵏᵉʳ(α).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactManifold {
    /// The contact form defining the structure.
    pub contact_form: ContactForm,
    /// Name/label for this manifold.
    pub name: String,
}

impl ContactManifold {
    /// Create a standard Darboux contact manifold R^{2n+1}.
    pub fn standard(dim: usize) -> Result<Self, String> {
        Ok(Self {
            contact_form: ContactForm::new(dim)?,
            name: format!("R^{} (standard contact)", dim),
        })
    }

    /// Create a named contact manifold.
    pub fn named(dim: usize, name: &str) -> Result<Self, String> {
        Ok(Self {
            contact_form: ContactForm::new(dim)?,
            name: name.to_string(),
        })
    }

    /// Dimension of the manifold.
    pub fn dim(&self) -> usize {
        self.contact_form.dim
    }

    /// Half-dimension n.
    pub fn n(&self) -> usize {
        self.contact_form.n
    }

    /// Check if a point is in the manifold (trivially true for R^{2n+1}).
    pub fn contains(&self, p: &DVector<f64>) -> bool {
        p.len() == self.dim()
    }

    /// Verify contact condition at a point.
    pub fn verify_contact(&self, p: &DVector<f64>) -> bool {
        self.contains(p) && self.contact_form.is_contact_at(p)
    }

    /// The symplectic structure on the contact distribution.
    /// Returns dα restricted to ξ = ker(α), which is non-degenerate.
    pub fn symplectic_on_distribution(&self) -> nalgebra::DMatrix<f64> {
        let n = self.n();
        // In Darboux coordinates, dα|ᵏᵉʳ(α) = Σ dxᵢ ∧ dyᵢ
        // This is just the standard symplectic form on R^{2n}
        let mut omega = nalgebra::DMatrix::zeros(2 * n, 2 * n);
        for i in 0..n {
            omega[(i, n + i)] = 1.0;
            omega[(n + i, i)] = -1.0;
        }
        omega
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn test_standard_contact_manifold_3d() {
        let m = ContactManifold::standard(3).unwrap();
        assert_eq!(m.dim(), 3);
        assert_eq!(m.n(), 1);
    }

    #[test]
    fn test_manifold_contains_point() {
        let m = ContactManifold::standard(3).unwrap();
        let p = dvector![1.0, 2.0, 3.0];
        assert!(m.contains(&p));
    }

    #[test]
    fn test_manifold_rejects_wrong_dimension_point() {
        let m = ContactManifold::standard(3).unwrap();
        let p = dvector![1.0, 2.0];
        assert!(!m.contains(&p));
    }

    #[test]
    fn test_verify_contact_at_point() {
        let m = ContactManifold::standard(3).unwrap();
        let p = dvector![1.0, 2.0, 3.0];
        assert!(m.verify_contact(&p));
    }

    #[test]
    fn test_symplectic_on_distribution() {
        let m = ContactManifold::standard(3).unwrap();
        let omega = m.symplectic_on_distribution();
        // 2x2 matrix [[0,1],[-1,0]]
        assert!((omega[(0, 1)] - 1.0).abs() < 1e-10);
        assert!((omega[(1, 0)] - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_named_manifold() {
        let m = ContactManifold::named(5, "Thermodynamic Phase Space").unwrap();
        assert_eq!(m.name, "Thermodynamic Phase Space");
        assert_eq!(m.dim(), 5);
    }
}
