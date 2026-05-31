use crate::contact_form::ContactForm;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// A Legendrian submanifold L ⊂ (M, α).
///
/// Legendrian submanifolds are the maximal submanifolds contained in
/// the contact distribution ξ = ker(α). They have dimension n in a
/// (2n+1)-dimensional contact manifold.
///
/// Equivalently: L is Legendrian if α|_L = 0.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LegendrianSubmanifold {
    /// The ambient contact form.
    pub contact_form: ContactForm,
    /// Name/label.
    pub name: String,
}

impl LegendrianSubmanifold {
    /// Create a new Legendrian submanifold in the given contact manifold.
    pub fn new(contact_form: ContactForm, name: &str) -> Self {
        Self {
            contact_form,
            name: name.to_string(),
        }
    }

    /// The zero-section: L = {(x, 0, 0)} in Darboux coordinates.
    /// This is Legendrian since α|_L = dz - 0·dx = dz = 0 on L (z=const).
    /// More precisely: for any tangent vector to L, α(v) = 0.
    pub fn zero_section(contact_form: ContactForm) -> Self {
        Self {
            contact_form,
            name: "Zero section".into(),
        }
    }

    /// A fiber: L = {(x₀, y, z₀)} for fixed x₀, z₀ — NOT Legendrian in general.
    /// The standard Legendrian examples are graphs: z = f(x), y = ∂f/∂x (1-jet).
    pub fn jet_section(contact_form: ContactForm, _f: fn(&DVector<f64>) -> f64) -> Self {
        Self {
            contact_form,
            name: "1-jet section".into(),
        }
    }

    /// Check if a point lies on the Legendrian zero-section.
    /// Zero-section: all yᵢ = 0 and z = 0.
    pub fn is_on_zero_section(&self, p: &DVector<f64>) -> bool {
        let n = self.contact_form.n;
        for i in 0..n {
            if p[n + i].abs() > 1e-10 {
                return false;
            }
        }
        p[2 * n].abs() < 1e-10
    }

    /// Check if a tangent vector v at point p lies in the contact distribution.
    /// α_p(v) = 0.
    pub fn is_in_contact_distribution(&self, p: &DVector<f64>, v: &DVector<f64>) -> bool {
        let alpha = self.contact_form.evaluate(p);
        alpha.dot(v).abs() < 1e-10
    }

    /// Verify that a given parametric submanifold is Legendrian.
    /// The submanifold is given by a map γ: Rⁿ → R^{2n+1}.
    /// Check: α(∂γ/∂uᵢ) = 0 for all i.
    pub fn verify_legendrian(&self, gamma: &[DVector<f64>], tangents: &[DVector<f64>]) -> bool {
        for (p, v) in gamma.iter().zip(tangents.iter()) {
            if !self.is_in_contact_distribution(p, v) {
                return false;
            }
        }
        true
    }

    /// Compute the Legendrian projection π: R^{2n+1} → Rⁿ (drop y and z).
    pub fn front_projection(&self, p: &DVector<f64>) -> DVector<f64> {
        let n = self.contact_form.n;
        p.rows(0, n).into()
    }

    /// Compute the Lagrangian projection π: R^{2n+1} → R^{2n} (drop z).
    pub fn lagrangian_projection(&self, p: &DVector<f64>) -> DVector<f64> {
        let n = self.contact_form.n;
        let mut proj = DVector::zeros(2 * n);
        for i in 0..(2 * n) {
            proj[i] = p[i];
        }
        proj
    }

    /// The expected dimension of a Legendrian submanifold.
    pub fn legendrian_dimension(&self) -> usize {
        self.contact_form.n
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn test_zero_section_is_legendrian() {
        let cf = ContactForm::new(3).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        // Points on zero-section: (x, 0, 0)
        let p = dvector![1.0, 0.0, 0.0];
        assert!(leg.is_on_zero_section(&p));
    }

    #[test]
    fn test_nonzero_not_on_zero_section() {
        let cf = ContactForm::new(3).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        let p = dvector![1.0, 2.0, 3.0];
        assert!(!leg.is_on_zero_section(&p));
    }

    #[test]
    fn test_tangent_in_contact_distribution() {
        let cf = ContactForm::new(3).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        let p = dvector![1.0, 2.0, 3.0];
        // ∂/∂y is in ξ since α(∂/∂y) = 0 (α = -y dx + dz, coefficient of dy is 0)
        let v = dvector![0.0, 1.0, 0.0];
        assert!(leg.is_in_contact_distribution(&p, &v));
    }

    #[test]
    fn test_reeb_not_in_contact_distribution() {
        let cf = ContactForm::new(3).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        let p = dvector![1.0, 2.0, 3.0];
        // ∂/∂z is NOT in ξ since α(∂/∂z) = 1
        let v = dvector![0.0, 0.0, 1.0];
        assert!(!leg.is_in_contact_distribution(&p, &v));
    }

    #[test]
    fn test_legendrian_dimension() {
        let cf = ContactForm::new(5).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        assert_eq!(leg.legendrian_dimension(), 2);
    }

    #[test]
    fn test_front_projection() {
        let cf = ContactForm::new(3).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        let p = dvector![1.0, 2.0, 3.0];
        let front = leg.front_projection(&p);
        assert_eq!(front.len(), 1);
        assert!((front[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lagrangian_projection() {
        let cf = ContactForm::new(3).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        let p = dvector![1.0, 2.0, 3.0];
        let lag = leg.lagrangian_projection(&p);
        assert_eq!(lag.len(), 2);
        assert!((lag[0] - 1.0).abs() < 1e-10);
        assert!((lag[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_verify_legendrian_zero_section() {
        let cf = ContactForm::new(3).unwrap();
        let leg = LegendrianSubmanifold::zero_section(cf);
        // Zero section: γ(t) = (t, 0, 0), tangent = (1, 0, 0)
        // α(1, 0, 0) at (t, 0, 0): α = -0·1 + 0 = 0 ✓
        let points = vec![dvector![0.0, 0.0, 0.0], dvector![1.0, 0.0, 0.0]];
        let tangents = vec![dvector![1.0, 0.0, 0.0], dvector![1.0, 0.0, 0.0]];
        assert!(leg.verify_legendrian(&points, &tangents));
    }
}
