// test function for importing
pub fn test(){
    println!("gamer")
}

use num_complex::Complex64;
use scilib::math::bessel::h1_nu; // H^{(1)}_\nu(z) 

// Functions for the BEM resonance scan

// Kernel function
#[inline(always)]
pub fn kernel(x: [f64; 2], n: [f64; 2], y: [f64; 2], k: f64) -> Complex64 {

    // dx = x - y
    let dx0 = x[0] - y[0];
    let dx1 = x[1] - y[1];

    // r = norm(dx)
    let r = (dx0.mul_add(dx0, dx1 * dx1)).sqrt();

    // jump term at singular point (dx -> 0)
    const THRSH: f64 = 1e-16; // (1e-8)^2
    if r < THRSH {
        return Complex64::new(-0.5, 0.0);
    }

    // dot(x - y, n)
    let dot = dx0.mul_add(n[0], dx1 * n[1]);
    let scale = dot / r;

    // H1^{(1)}(k r)
    let kr = k * r; 
    let h1 = h1_nu(1.0, Complex64::new(kr, 0.0)); // H^{(1)}_1(kr)

    // imaginary coefficient
    let coeff = Complex64::new(0.0, -(k * 0.25));

    // imaginary coefficient * Hankel function * scaling factor
    coeff * h1 * scale
}



pub fn construct_matrix() {
    

}

pub fn min_sin_value() {


}

pub fn compute_spectrum() {


}