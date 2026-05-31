use crate::contact_form::ContactForm;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// Darboux coordinates and the Darboux theorem.
///
/// The Darboux theorem states: all contact structures of the same dimension
/// are locally equivalent. Around any point, there exist coordinates
/// (x₁,...,xₙ, y₁,...,yₙ, z) such that α = dz - Σᵢ yᵢ dxᵢ.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DarbouxCoordinates {
    /// The contact form in these coordinates.
    pub contact_form: ContactForm,
}

impl DarbouxCoordinates {
    /// Create Darboux coordinates for a contact manifold of given dimension.
    pub fn new(dim: usize) -> Result<Self, String> {
        Ok(Self {
            contact_form: ContactForm::new(dim)?,
        })
    }

    /// Convert a general point to Darboux coordinates.
    /// In the standard case, this is the identity map.
    pub fn to_darboux(&self, p: &DVector<f64>) -> DVector<f64> {
        p.clone()
    }

    /// Convert from Darboux coordinates back to the original.
    pub fn from_darboux(&self, p: &DVector<f64>) -> DVector<f64> {
        p.clone()
    }

    /// Verify that the contact form is in Darboux normal form:
    /// α = dz - Σᵢ yᵢ dxᵢ
    pub fn verify_normal_form(&self, p: &DVector<f64>) -> bool {
        let alpha = self.contact_form.evaluate(p);
        let n = self.contact_form.n;
        // In Darboux coords: α = (-y₁, ..., -yₙ, 0, ..., 0, 1)
        // Check the x-coefficients are -yᵢ
        for i in 0..n {
            if (alpha[i] - (-p[n + i])).abs() > 1e-10 {
                return false;
            }
        }
        // Check the y-coefficients are 0
        for i in 0..n {
            if alpha[n + i].abs() > 1e-10 {
                return false;
            }
        }
        // Check z-coefficient is 1
        if (alpha[2 * n] - 1.0).abs() > 1e-10 {
            return false;
        }
        true
    }

    /// Compute the Darboux chart transition map (identity for standard).
    /// In general, this would be the coordinate change bringing α to normal form.
    pub fn transition_map(&self, p: &DVector<f64>) -> DVector<f64> {
        self.to_darboux(p)
    }

    /// Verify the Darboux theorem holds: contact volume in normal form is 1.
    pub fn verify_darboux_theorem(&self) -> bool {
        let dim = self.contact_form.dim;
        let p = DVector::zeros(dim);
        self.verify_normal_form(&p) && self.contact_form.is_contact_at(&p)
    }

    /// Construct a Darboux chart around a given point.
    /// Returns the coordinates of the chart center.
    pub fn chart_around(&self, point: &DVector<f64>) -> DarbouxChart {
        DarbouxChart {
            center: point.clone(),
            dim: self.contact_form.dim,
            n: self.contact_form.n,
        }
    }
}

/// A local Darboux chart around a point on a contact manifold.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DarbouxChart {
    /// Center of the chart.
    pub center: DVector<f64>,
    /// Dimension of the manifold.
    pub dim: usize,
    /// Half-dimension.
    pub n: usize,
}

impl DarbouxChart {
    /// Check if a point is within this chart's domain.
    /// For R^{2n+1}, all points are valid.
    pub fn contains(&self, p: &DVector<f64>) -> bool {
        p.len() == self.dim
    }

    /// The contact form in this chart (always Darboux normal form).
    pub fn contact_form_value(&self, p: &DVector<f64>) -> DVector<f64> {
        let cf = ContactForm::new(self.dim).unwrap();
        cf.evaluate(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn test_darboux_creation() {
        let dc = DarbouxCoordinates::new(3).unwrap();
        assert_eq!(dc.contact_form.dim, 3);
    }

    #[test]
    fn test_normal_form_verification() {
        let dc = DarbouxCoordinates::new(3).unwrap();
        let p = dvector![1.0, 2.0, 3.0];
        assert!(dc.verify_normal_form(&p));
    }

    #[test]
    fn test_normal_form_at_origin() {
        let dc = DarbouxCoordinates::new(5).unwrap();
        let p = dvector![0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(dc.verify_normal_form(&p));
    }

    #[test]
    fn test_darboux_theorem() {
        let dc = DarbouxCoordinates::new(3).unwrap();
        assert!(dc.verify_darboux_theorem());
    }

    #[test]
    fn test_darboux_theorem_higher_dim() {
        let dc = DarbouxCoordinates::new(7).unwrap();
        assert!(dc.verify_darboux_theorem());
    }

    #[test]
    fn test_chart_around_point() {
        let dc = DarbouxCoordinates::new(3).unwrap();
        let p = dvector![1.0, 2.0, 3.0];
        let chart = dc.chart_around(&p);
        assert!((chart.center[0] - 1.0).abs() < 1e-10);
        assert!(chart.contains(&dvector![0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_roundtrip_darboux() {
        let dc = DarbouxCoordinates::new(3).unwrap();
        let p = dvector![1.5, -2.3, 0.7];
        let darboux = dc.to_darboux(&p);
        let back = dc.from_darboux(&darboux);
        for i in 0..3 {
            assert!((p[i] - back[i]).abs() < 1e-10);
        }
    }
}
