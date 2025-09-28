// test function for importing
pub fn test(){
    println!("gamer")
}

use num_complex::Complex64; // not getting far without this
use scilib::math::bessel::h1_nu; // H^{(1)}_\nu(z): Hankel function of order 1 and type nu (here 1)
use rayon::prelude::*; // concurrency fun
use crate::boundarydata::{self, Nodes}; // need nodes from boundarydata.rs
use plotters::prelude::*;
use colorous::VIRIDIS;
use std::error::Error;
use bessel;

// Functions for the BEM resonance scan

// Kernel function
#[inline(always)]
pub fn kernel(nodes: &Nodes, i: usize, j: usize, k: f64) -> Complex64 {
    let xi  = nodes.x[i];
    let yi  = nodes.y[i];
    let nxi = nodes.nx[i];
    let nyi = nodes.ny[i];
    let xj = nodes.x[j];
    let yj = nodes.y[j];

    let dx0 = xi - xj;
    let dx1 = yi - yj;
    let r2  = dx0.mul_add(dx0, dx1 * dx1);

    const EPS2: f64 = 1.0e-16; // (1e-8)^2
    if r2 < EPS2 {
        return Complex64::new(-0.5, 0.0);
    }

    let dot   = dx0.mul_add(nxi, dx1 * nyi);
    let r     = r2.sqrt();
    let scale = dot / r;

    let kr = k * r;
    // Fast real Bessels:
    let a = bessel::j1(kr); // Re H1
    let b = bessel::y1(kr); // Im H1

    // -(i k / 4)*(a + i b) = (k/4)*(b - i a)
    let c = k * 0.25 * scale;
    Complex64::new(c * b, -c * a)
}
// Constructs the BEM matrix operator
pub fn construct_matrix(
    k: f64,
    nodes: &Nodes,
    kernel_fn: fn(&Nodes, usize, usize, f64) -> Complex64,
) -> Vec<Complex64> {
    let n = nodes.x.len();
    assert_eq!(nodes.y.len(), n);
    assert_eq!(nodes.nx.len(), n);
    assert_eq!(nodes.ny.len(), n);
    assert_eq!(nodes.w.len(), n);

    let mut A = vec![Complex64::new(0.0, 0.0); n * n];

    A.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
        for j in 0..n {
            // jump term + kernel * weight
            let mut aij = if i == j {
                Complex64::new(-0.5, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            aij += kernel_fn(nodes, i, j, k) * nodes.w[j];
            row[j] = aij;
        }
    });

    A
}
// these next functions will be visualisation stuff 

// the matrix Helmholtz equation is A * phi = 0
// for a given A, extract an approximate nontrivial solution phi on the boundary. 
// choose the singular vector corresponding to the smallest singular value.
