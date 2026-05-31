use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};

/// A contact form α on a (2n+1)-dimensional manifold.
///
/// The defining property: α ∧ (dα)ⁿ ≠ 0 (maximal non-integrability).
///
/// In Darboux coordinates (x₁,...,xₙ, y₁,...,yₙ, z):
///   α = dz - Σᵢ yᵢ dxᵢ
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactForm {
    /// Dimension of the contact manifold (must be odd: 2n+1).
    pub dim: usize,
    /// Half-dimension n where dim = 2n+1.
    pub n: usize,
}

impl ContactForm {
    /// Create a new contact form on a (2n+1)-dimensional manifold.
    pub fn new(dim: usize) -> Result<Self, String> {
        if dim < 3 || dim.is_multiple_of(2) {
            return Err(format!(
                "Contact manifold dimension must be odd and ≥ 3, got {}",
                dim
            ));
        }
        let n = (dim - 1) / 2;
        Ok(Self { dim, n })
    }

    /// Evaluate the contact form α at a point p in Darboux coordinates.
    ///
    /// α = dz - Σᵢ yᵢ dxᵢ
    ///
    /// Point layout: p = (x₁,...,xₙ, y₁,...,yₙ, z)
    pub fn evaluate(&self, p: &DVector<f64>) -> DVector<f64> {
        assert_eq!(p.len(), self.dim, "Point dimension mismatch");
        // α is a 1-form; its action on a tangent vector v gives a scalar.
        // In coordinates, α = -y₁ dx₁ - ... - yₙ dxₙ + dz
        // Represented as a covector (row vector):
        let mut alpha = DVector::zeros(self.dim);
        // Coefficients of dxᵢ: -yᵢ
        for i in 0..self.n {
            alpha[i] = -p[self.n + i];
        }
        // Coefficient of dz: 1
        alpha[2 * self.n] = 1.0;
        alpha
    }

    /// Compute dα (the exterior derivative of α).
    ///
    /// In Darboux coordinates: dα = -Σᵢ dyᵢ ∧ dxᵢ = Σᵢ dxᵢ ∧ dyᵢ
    ///
    /// Returns the matrix of dα as a 2-form (skew-symmetric).
    pub fn exterior_derivative(&self) -> DMatrix<f64> {
        let d = self.dim;
        let mut da = DMatrix::zeros(d, d);
        // dα = Σᵢ dxᵢ ∧ dyᵢ
        for i in 0..self.n {
            // dxᵢ ∧ dyᵢ means da[xᵢ, yᵢ] = 1, da[yᵢ, xᵢ] = -1
            da[(i, self.n + i)] = 1.0;
            da[(self.n + i, i)] = -1.0;
        }
        da
    }

    /// Compute α ∧ (dα)ⁿ — the contact volume form.
    ///
    /// For a valid contact form, this must be non-zero everywhere.
    /// Returns the Pfaffian-like scalar measuring the volume.
    pub fn contact_volume(&self, p: &DVector<f64>) -> f64 {
        let alpha = self.evaluate(p);
        // α ∧ (dα)ⁿ evaluated as: α_z × Pfaffian²(dα|ξ)
        // In Darboux coordinates, α ∧ (dα)ⁿ = dz ∧ (dx₁∧dy₁) ∧ ... ∧ (dxₙ∧dyₙ) = 1
        let alpha_z = alpha[2 * self.n];
        alpha_z * self.pfaffian_squared()
    }

    /// Pfaffian squared of dα restricted to the contact distribution.
    /// In Darboux coordinates this is 1.
    pub fn pfaffian_squared(&self) -> f64 {
        // dα on the 2n-subspace is block diagonal with [[0,1],[-1,0]] blocks
        // Pfaffian of each block is 1, so Pf(dα)² = 1
        1.0
    }

    /// Verify the contact condition at a point: α ∧ (dα)ⁿ ≠ 0.
    pub fn is_contact_at(&self, p: &DVector<f64>) -> bool {
        self.contact_volume(p).abs() > 1e-12
    }

    /// The contact distribution ξ = ker(α) at point p.
    /// Returns a basis for the (2n)-dimensional subspace.
    pub fn contact_distribution(&self, p: &DVector<f64>) -> DMatrix<f64> {
        // ξ = {v : α(v) = 0} = {v : αᵀv = 0}
        // In Darboux coordinates:
        // α = -y₁ dx₁ - ... - yₙ dxₙ + dz
        // So ξ is spanned by ∂/∂yᵢ and ∂/∂xᵢ + yᵢ ∂/∂z
        let mut basis = DMatrix::zeros(self.dim, 2 * self.n);
        for i in 0..self.n {
            // ∂/∂yᵢ
            basis[(self.n + i, i)] = 1.0;
            // ∂/∂xᵢ + yᵢ ∂/∂z
            basis[(i, self.n + i)] = 1.0;
            basis[(2 * self.n, self.n + i)] = p[self.n + i];
        }
        basis
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn test_contact_form_creation() {
        let cf = ContactForm::new(3).unwrap();
        assert_eq!(cf.dim, 3);
        assert_eq!(cf.n, 1);
    }

    #[test]
    fn test_reject_even_dimension() {
        assert!(ContactForm::new(2).is_err());
        assert!(ContactForm::new(4).is_err());
    }

    #[test]
    fn test_reject_dimension_less_than_3() {
        assert!(ContactForm::new(1).is_err());
        assert!(ContactForm::new(0).is_err());
    }

    #[test]
    fn test_darboux_3d_contact_form() {
        let cf = ContactForm::new(3).unwrap();
        // Point (x, y, z) = (1, 2, 3)
        let p = dvector![1.0, 2.0, 3.0];
        let alpha = cf.evaluate(&p);
        // α = -y dx + dz → covector (-2, 0, 1)
        assert!((alpha[0] - (-2.0)).abs() < 1e-10);
        assert!(alpha[1].abs() < 1e-10);
        assert!((alpha[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_contact_condition_3d() {
        let cf = ContactForm::new(3).unwrap();
        let p = dvector![1.0, 2.0, 3.0];
        assert!(cf.is_contact_at(&p));
    }

    #[test]
    fn test_contact_condition_5d() {
        let cf = ContactForm::new(5).unwrap();
        let p = dvector![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(cf.is_contact_at(&p));
    }

    #[test]
    fn test_exterior_derivative_3d() {
        let cf = ContactForm::new(3).unwrap();
        let da = cf.exterior_derivative();
        // dα = dx ∧ dy → matrix [[0,1,0],[-1,0,0],[0,0,0]]
        assert!((da[(0, 1)] - 1.0).abs() < 1e-10);
        assert!((da[(1, 0)] - (-1.0)).abs() < 1e-10);
        assert!(da[(0, 2)].abs() < 1e-10);
        assert!(da[(2, 2)].abs() < 1e-10);
    }

    #[test]
    fn test_contact_volume_nonzero() {
        let cf = ContactForm::new(3).unwrap();
        let p = dvector![0.0, 0.0, 0.0];
        // At origin: α = dz, so contact volume = 1.0
        assert!((cf.contact_volume(&p) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_contact_distribution_dimension() {
        let cf = ContactForm::new(5).unwrap();
        let p = dvector![1.0, 2.0, 3.0, 4.0, 5.0];
        let xi = cf.contact_distribution(&p);
        assert_eq!(xi.nrows(), 5);
        assert_eq!(xi.ncols(), 4);
    }

    #[test]
    fn test_higher_dimensional_contact() {
        let cf = ContactForm::new(7).unwrap();
        assert_eq!(cf.n, 3);
        let p = dvector![1.0, 0.0, -1.0, 2.0, 0.5, -0.5, 3.0];
        assert!(cf.is_contact_at(&p));
    }
}
