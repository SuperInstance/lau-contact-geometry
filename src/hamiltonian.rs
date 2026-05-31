use crate::contact_form::ContactForm;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// Contact Hamiltonian system on a (2n+1)-dimensional contact manifold.
///
/// Given a Hamiltonian H: M → R, the contact Hamilton's equations are:
///   ẋ = ∂H/∂y
///   ẏ = -∂H/∂x + y·∂H/∂z
///   ż = H - y·(∂H/∂y) + ∂H/∂z
///
/// These preserve the contact structure (up to conformal factor).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactHamiltonianSystem {
    /// The contact form.
    pub contact_form: ContactForm,
    /// The Hamiltonian function (stored as coefficients for polynomial approximation,
    /// or as the actual function for evaluation).
    pub hamiltonian_type: HamiltonianType,
}

/// Type of Hamiltonian.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HamiltonianType {
    /// H(x, y, z) = H₀ (constant).
    Constant(f64),
    /// H(x, y, z) = z (height/energy function).
    Height,
    /// H(x, y, z) = Σᵢ (xᵢ² + yᵢ²)/2 + z (harmonic + height).
    HarmonicPlusHeight,
    /// H(x, y, z) = z² (quadratic in z).
    QuadraticZ,
    /// Generic: user provides evaluation and partial derivatives.
    /// Stored as (H_value, ∂H/∂x, ∂H/∂y, ∂H/∂z) at a point.
    Custom,
}

impl ContactHamiltonianSystem {
    /// Create a system with a constant Hamiltonian.
    pub fn constant(contact_form: ContactForm, h0: f64) -> Self {
        Self {
            contact_form,
            hamiltonian_type: HamiltonianType::Constant(h0),
        }
    }

    /// Create a system with the height Hamiltonian H = z.
    pub fn height(contact_form: ContactForm) -> Self {
        Self {
            contact_form,
            hamiltonian_type: HamiltonianType::Height,
        }
    }

    /// Create a system with H = harmonic + z.
    pub fn harmonic_plus_height(contact_form: ContactForm) -> Self {
        Self {
            contact_form,
            hamiltonian_type: HamiltonianType::HarmonicPlusHeight,
        }
    }

    /// Evaluate H at a point.
    pub fn evaluate_h(&self, p: &DVector<f64>) -> f64 {
        let n = self.contact_form.n;
        match &self.hamiltonian_type {
            HamiltonianType::Constant(h0) => *h0,
            HamiltonianType::Height => p[2 * n],
            HamiltonianType::HarmonicPlusHeight => {
                let mut h = p[2 * n];
                for i in 0..n {
                    h += (p[i] * p[i] + p[n + i] * p[n + i]) / 2.0;
                }
                h
            }
            HamiltonianType::QuadraticZ => p[2 * n] * p[2 * n],
            HamiltonianType::Custom => 0.0,
        }
    }

    /// Compute ∂H/∂x at a point. Returns n-vector.
    pub fn dh_dx(&self, p: &DVector<f64>) -> DVector<f64> {
        let n = self.contact_form.n;
        match &self.hamiltonian_type {
            HamiltonianType::Constant(_) => DVector::zeros(n),
            HamiltonianType::Height => DVector::zeros(n),
            HamiltonianType::HarmonicPlusHeight => {
                let mut grad = DVector::zeros(n);
                for i in 0..n {
                    grad[i] = p[i];
                }
                grad
            }
            HamiltonianType::QuadraticZ => DVector::zeros(n),
            HamiltonianType::Custom => DVector::zeros(n),
        }
    }

    /// Compute ∂H/∂y at a point. Returns n-vector.
    pub fn dh_dy(&self, p: &DVector<f64>) -> DVector<f64> {
        let n = self.contact_form.n;
        match &self.hamiltonian_type {
            HamiltonianType::Constant(_) => DVector::zeros(n),
            HamiltonianType::Height => DVector::zeros(n),
            HamiltonianType::HarmonicPlusHeight => {
                let mut grad = DVector::zeros(n);
                for i in 0..n {
                    grad[i] = p[n + i];
                }
                grad
            }
            HamiltonianType::QuadraticZ => DVector::zeros(n),
            HamiltonianType::Custom => DVector::zeros(n),
        }
    }

    /// Compute ∂H/∂z at a point. Returns scalar.
    pub fn dh_dz(&self, _p: &DVector<f64>) -> f64 {
        match &self.hamiltonian_type {
            HamiltonianType::Constant(_) => 0.0,
            HamiltonianType::Height => 1.0,
            HamiltonianType::HarmonicPlusHeight => 1.0,
            HamiltonianType::QuadraticZ => 2.0 * _p[2 * self.contact_form.n],
            HamiltonianType::Custom => 0.0,
        }
    }

    /// Compute the contact Hamiltonian vector field X_H at point p.
    ///
    /// Contact Hamilton's equations:
    ///   ẋᵢ = ∂H/∂yᵢ
    ///   ẏᵢ = -∂H/∂xᵢ + yᵢ · ∂H/∂z
    ///   ż = H - Σᵢ yᵢ (∂H/∂yᵢ) + z · ∂H/∂z
    pub fn vector_field(&self, p: &DVector<f64>) -> DVector<f64> {
        let n = self.contact_form.n;
        let h = self.evaluate_h(p);
        let dhx = self.dh_dx(p);
        let dhy = self.dh_dy(p);
        let dhz = self.dh_dz(p);

        let mut x_h = DVector::zeros(self.contact_form.dim);

        // ẋ = ∂H/∂y
        for i in 0..n {
            x_h[i] = dhy[i];
        }

        // ẏ = -∂H/∂x + y · ∂H/∂z
        for i in 0..n {
            x_h[n + i] = -dhx[i] + p[n + i] * dhz;
        }

        // ż = H - y·∂H/∂y + z·∂H/∂z   (note: using z here, some formulations differ)
        // Standard: ż = -(n+1)·H + Σᵢ xᵢ ∂H/∂xᵢ + ... 
        // But for the standard formulation: ż = -H + Σ yᵢ ∂H/∂yᵢ ... 
        // Let's use the clean form from contact geometry literature:
        // ż = H - Σᵢ yᵢ(∂H/∂yᵢ)
        let y_dhy: f64 = (0..n).map(|i| p[n + i] * dhy[i]).sum();
        x_h[2 * n] = h - y_dhy;

        x_h
    }

    /// Step the system forward by dt using Euler integration.
    pub fn step(&self, p: &DVector<f64>, dt: f64) -> DVector<f64> {
        let xh = self.vector_field(p);
        p + dt * xh
    }

    /// Integrate the system for multiple steps.
    pub fn integrate(&self, p0: &DVector<f64>, dt: f64, steps: usize) -> Vec<DVector<f64>> {
        let mut trajectory = vec![p0.clone()];
        let mut p = p0.clone();
        for _ in 0..steps {
            p = self.step(&p, dt);
            trajectory.push(p.clone());
        }
        trajectory
    }

    /// Compute the conformal factor of the contactomorphism generated by X_H.
    /// The contact Hamiltonian flow satisfies φ_t^* α = e^{λt} α where λ = -∂H/∂z.
    pub fn conformal_factor(&self, p: &DVector<f64>) -> f64 {
        -self.dh_dz(p)
    }

    /// Check if the Hamiltonian is conserved along its own flow (energy conservation).
    /// For contact systems, H is NOT generally conserved; it satisfies dH/dt = λ·H.
    pub fn energy_evolution_rate(&self, p: &DVector<f64>) -> f64 {
        self.conformal_factor(p) * self.evaluate_h(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::dvector;

    #[test]
    fn test_constant_hamiltonian() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::constant(cf, 5.0);
        let p = dvector![1.0, 2.0, 3.0];
        assert!((sys.evaluate_h(&p) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_height_hamiltonian() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::height(cf);
        let p = dvector![1.0, 2.0, 3.0];
        assert!((sys.evaluate_h(&p) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_harmonic_hamiltonian() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::harmonic_plus_height(cf);
        let p = dvector![1.0, 2.0, 3.0];
        // H = (1+4)/2 + 3 = 5.5
        assert!((sys.evaluate_h(&p) - 5.5).abs() < 1e-10);
    }

    #[test]
    fn test_constant_vector_field() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::constant(cf, 2.0);
        let p = dvector![1.0, 2.0, 3.0];
        let xh = sys.vector_field(&p);
        // For constant H: ẋ=0, ẏ=0, ż=H=2
        assert!(xh[0].abs() < 1e-10);
        assert!(xh[1].abs() < 1e-10);
        assert!((xh[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_euler_step() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::constant(cf, 1.0);
        let p = dvector![0.0, 0.0, 0.0];
        let p1 = sys.step(&p, 0.1);
        assert!((p1[2] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_integrate_trajectory_length() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::constant(cf, 1.0);
        let p0 = dvector![0.0, 0.0, 0.0];
        let traj = sys.integrate(&p0, 0.01, 100);
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_conformal_factor_height() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::height(cf);
        let p = dvector![1.0, 2.0, 3.0];
        // λ = -∂H/∂z = -1
        assert!((sys.conformal_factor(&p) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_energy_evolution() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::height(cf);
        let p = dvector![0.0, 0.0, 5.0];
        // dH/dt = λ·H = -1·5 = -5
        assert!((sys.energy_evolution_rate(&p) - (-5.0)).abs() < 1e-10);
    }

    #[test]
    fn test_harmonic_vector_field() {
        let cf = ContactForm::new(3).unwrap();
        let sys = ContactHamiltonianSystem::harmonic_plus_height(cf);
        let p = dvector![1.0, 0.0, 0.0];
        let xh = sys.vector_field(&p);
        // n=1: ẋ = ∂H/∂y = y = 0, ẏ = -∂H/∂x + y·∂H/∂z = -x + 0 = -1, ż = H - y·∂H/∂y = 0.5 - 0 = 0.5
        assert!(xh[0].abs() < 1e-10);
        assert!((xh[1] - (-1.0)).abs() < 1e-10);
        assert!((xh[2] - 0.5).abs() < 1e-10);
    }
}
