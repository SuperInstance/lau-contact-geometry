use crate::contact_form::ContactForm;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// Reeb vector field R on a contact manifold (M, α).
///
/// Characterized by:
///   α(R) = 1
///   dα(R, ·) = 0
///
/// In Darboux coordinates (x₁,...,xₙ, y₁,...,yₙ, z):
///   R = ∂/∂z
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReebVectorField {
    /// The contact form.
    pub contact_form: ContactForm,
}

impl ReebVectorField {
    /// Create the Reeb vector field for a given contact form.
    pub fn new(contact_form: ContactForm) -> Self {
        Self { contact_form }
    }

    /// Evaluate the Reeb vector field at a point.
    ///
    /// In Darboux coordinates: R = (0, ..., 0, 1) — only the z-component.
    pub fn evaluate(&self, _p: &DVector<f64>) -> DVector<f64> {
        let mut reeb = DVector::zeros(self.contact_form.dim);
        reeb[2 * self.contact_form.n] = 1.0;
        reeb
    }

    /// Verify α(R) = 1 at a point.
    pub fn verify_alpha_condition(&self, p: &DVector<f64>) -> bool {
        let alpha = self.contact_form.evaluate(p);
        let reeb = self.evaluate(p);
        let ar = alpha.dot(&reeb);
        (ar - 1.0).abs() < 1e-10
    }

    /// Verify dα(R, ·) = 0 at a point (Reeb is in kernel of dα).
    pub fn verify_da_condition(&self, p: &DVector<f64>) -> bool {
        let da = self.contact_form.exterior_derivative();
        let reeb = self.evaluate(p);
        // dα(R, ·) = da * R (as matrix-vector product)
        let da_r = &da * &reeb;
        da_r.iter().all(|&x| x.abs() < 1e-10)
    }

    /// Verify both Reeb conditions simultaneously.
    pub fn verify_all(&self, p: &DVector<f64>) -> bool {
        self.verify_alpha_condition(p) && self.verify_da_condition(p)
    }

    /// Compute the Reeb flow: Φ_t(p) = p + t·R(p) (in Darboux coords, just z → z + t).
    pub fn flow(&self, p: &DVector<f64>, t: f64) -> DVector<f64> {
        let reeb = self.evaluate(p);
        p + t * reeb
    }

    /// The Reeb orbits are curves along the z-direction.
    /// Period of a Reeb orbit (None = non-periodic / infinite).
    /// For the standard structure on R^{2n+1}, Reeb orbits are lines (no period).
    pub fn orbit_period(&self) -> Option<f64> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn test_reeb_vector_3d() {
        let cf = ContactForm::new(3).unwrap();
        let reeb = ReebVectorField::new(cf);
        let p = dvector![1.0, 2.0, 3.0];
        let r = reeb.evaluate(&p);
        // R = (0, 0, 1) in Darboux
        assert!(r[0].abs() < 1e-10);
        assert!(r[1].abs() < 1e-10);
        assert!((r[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reeb_alpha_equals_one() {
        let cf = ContactForm::new(3).unwrap();
        let reeb = ReebVectorField::new(cf);
        let p = dvector![1.0, 2.0, 3.0];
        assert!(reeb.verify_alpha_condition(&p));
    }

    #[test]
    fn test_reeb_da_equals_zero() {
        let cf = ContactForm::new(3).unwrap();
        let reeb = ReebVectorField::new(cf);
        let p = dvector![1.0, 2.0, 3.0];
        assert!(reeb.verify_da_condition(&p));
    }

    #[test]
    fn test_reeb_all_conditions() {
        let cf = ContactForm::new(5).unwrap();
        let reeb = ReebVectorField::new(cf);
        let p = dvector![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(reeb.verify_all(&p));
    }

    #[test]
    fn test_reeb_flow() {
        let cf = ContactForm::new(3).unwrap();
        let reeb = ReebVectorField::new(cf);
        let p = dvector![1.0, 2.0, 3.0];
        let flowed = reeb.flow(&p, 2.5);
        // Only z changes
        assert!((flowed[0] - 1.0).abs() < 1e-10);
        assert!((flowed[1] - 2.0).abs() < 1e-10);
        assert!((flowed[2] - 5.5).abs() < 1e-10);
    }

    #[test]
    fn test_reeb_7d() {
        let cf = ContactForm::new(7).unwrap();
        let reeb = ReebVectorField::new(cf);
        let p = dvector![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let r = reeb.evaluate(&p);
        assert!((r[6] - 1.0).abs() < 1e-10);
        assert!(reeb.verify_all(&p));
    }
}
