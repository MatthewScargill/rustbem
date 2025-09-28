mod bem;
mod boundarydata;
use num_complex::Complex64;
use std::time::Instant;


fn main() {
    let n: usize = 24;
    let nodes = boundarydata::square_billiard_nodes(1f64, n);

    let start = Instant::now();
    let A = bem::construct_matrix(2.5, &nodes,bem::kernel);
    print_matrix(&A, n);

    println!("V1 done in {:?}", start.elapsed());
}

pub fn print_matrix(a: &[Complex64], n: usize) {
    assert_eq!(a.len(), n * n, "length does not match n*n");

    for i in 0..n {
        for j in 0..n {
            let v = a[i * n + j];
            // Print as "a+bi" with scientific formatting
            print!("{:>12.4e}{:+12.4e}i  ", v.re, v.im);
        }
        println!();
    }
}