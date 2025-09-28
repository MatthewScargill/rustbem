use plotters::prelude::*;
use std::error::Error;

// Welcome to the Nodes structure 
// This (zero cost wow) abstraction holds all of the boundary data for our billiard
#[derive(Debug)]
pub struct Nodes {
    pub x: Vec<f64>,   // x-coordinates of boundary nodes
    pub y: Vec<f64>,   // y-coordinates
    pub nx: Vec<f64>,  // outward normal x
    pub ny: Vec<f64>,  // outward normal y
    pub w: Vec<f64>,   // quadrature weights (panel lengths)
    pub s: Vec<f64>,   // arclength parameters used (optional but handy)
    pub l_total: f64,  // total perimeter (will become important for Weyl estimation stuff later)
}

// Computes the Nodes for the square boundary of length "a" and node count "n"
#[inline(always)]
pub fn square_billiard_nodes(a: f64, n: usize) -> Nodes {

    // check length and number of nodes for safety
    assert!(a > 0.0, "side length must be > 0");
    assert!(n >= 4 && n % 4 == 0, "N must be a multiple of 4 (got {n})");

    // compute total length, weights, and midpoints (to avoid corners)
    let l_total = 4.0 * a;
    let h = l_total / n as f64;
    let half_h = 0.5 * h;

    // initialise the nodes entries
    let mut x  = Vec::with_capacity(n);
    let mut y  = Vec::with_capacity(n);
    let mut nx = Vec::with_capacity(n);
    let mut ny = Vec::with_capacity(n);
    let mut w  = Vec::with_capacity(n);
    let mut s  = Vec::with_capacity(n);

    // midpoint nodes: s_i = (i + 0.5) * h, i=0..N-1
    for i in 0..n {
        let si = (i as f64) * h + half_h;
        let mut s_wrap = si;
        // wrap just in case (for safety with float ops)
        if s_wrap >= l_total {
            s_wrap -= l_total * (s_wrap / l_total).floor();
        }

        // organise the normals 
        // side = 0: bottom (0,0) -> (a,0),   n = (0,-1)
        // side = 1: right  (a,0) -> (a,a),   n = (1, 0)
        // side = 2: top    (a,a) -> (0,a),   n = (0, 1)
        // side = 3: left   (0,a) -> (0,0),   n = (-1,0)
        let side = (s_wrap / a).floor() as i32; // 0..=3
        let t = s_wrap - (side as f64) * a;     // position along the current side in [0,a)

        let (xi, yi, nxi, nyi) = match side {
            0 => (t,       0.0,   0.0, -1.0),
            1 => (a,       t,     1.0,  0.0),
            2 => (a - t,   a,     0.0,  1.0),
            _ => (0.0,     a - t, -1.0, 0.0), 
        };

        // pop them onto their respective piles
        x.push(xi);
        y.push(yi);
        nx.push(nxi);
        ny.push(nyi);
        w.push(h);
        s.push(s_wrap);
    }

    Nodes { x, y, nx, ny, w, s, l_total }
}




// some visualisation stuff to sanity check the boundaries (need to fix and generalise)
pub fn plot_nodes_svg(
    path: &str,
    nodes: &Nodes,
    normal_scale: f64,
) -> Result<(), Box<dyn Error>> {
    let n = nodes.x.len();
    assert!(n > 1);

    // ---- bounds with padding ----
    let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..n {
        xmin = xmin.min(nodes.x[i]); xmax = xmax.max(nodes.x[i]);
        ymin = ymin.min(nodes.y[i]); ymax = ymax.max(nodes.y[i]);
    }
    let dx = (xmax - xmin).max(1e-12);
    let dy = (ymax - ymin).max(1e-12);
    let pad_x = 0.05 * dx;
    let pad_y = 0.05 * dy;

    // ---- canvas ----
    let root = SVGBackend::new(path, (900, 900)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption("Boundary Nodes & Normals", ("Arial", 25))
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d((xmin - pad_x)..(xmax + pad_x), (ymin - pad_y)..(ymax + pad_y))?;

    chart.configure_mesh().x_labels(10).y_labels(10).draw()?;

    // ---- boundary polyline (closed visually) ----
    {
        let mut poly: Vec<(f64, f64)> = (0..n).map(|i| (nodes.x[i], nodes.y[i])).collect();
        poly.push((nodes.x[0], nodes.y[0]));
        chart.draw_series(LineSeries::new(poly, &BLACK.mix(0.4)))?;
    }

    // ---- nodes as dots ----
    chart.draw_series((0..n).map(|i| Circle::new((nodes.x[i], nodes.y[i]), 3, BLUE.filled())))?;

    // ---- normals as three simple series (shaft + 2 heads) ----
    let arrow_color = &RED.mix(0.9);
    let head_scale = normal_scale * 0.25;

    // 1) shafts
    chart.draw_series((0..n).map(|i| {
        let (x0, y0) = (nodes.x[i], nodes.y[i]);
        let (x1, y1) = (x0 + normal_scale * nodes.nx[i], y0 + normal_scale * nodes.ny[i]);
        PathElement::new(vec![(x0, y0), (x1, y1)], arrow_color)
    }))?;

    // 2) heads
    chart.draw_series((0..n).map(|i| {
        let (x0, y0) = (nodes.x[i], nodes.y[i]);
        let (nx, ny) = (nodes.nx[i], nodes.ny[i]);
        let (x1, y1) = (x0 + normal_scale * nx, y0 + normal_scale * ny);

        // unit direction
        let dx = x1 - x0;
        let dy = y1 - y0;
        let len = (dx * dx + dy * dy).sqrt().max(1e-16);
        let (ux, uy) = (dx / len, dy / len);

        // rotate (-ux,-uy) by ±30°
        let cos_th = (30f64).to_radians().cos();
        let sin_th = (30f64).to_radians().sin();

        let rx1 = -ux * cos_th - (-uy) * sin_th;
        let ry1 = -ux * sin_th + (-uy) * cos_th;
        let rx2 = -ux * cos_th - ( uy) * sin_th;
        let ry2 = -ux * sin_th + ( uy) * cos_th;

        let h1 = (x1 + head_scale * rx1, y1 + head_scale * ry1);
        let h2 = (x1 + head_scale * rx2, y1 + head_scale * ry2);

        // We'll emit TWO elements from this iterator position by returning a dummy;
        // instead, just build a short polyline for one head here,
        // and the second head in another series right below.
        PathElement::new(vec![(x1, y1), h1], arrow_color)
    }))?;

    chart.draw_series((0..n).map(|i| {
        let (x0, y0) = (nodes.x[i], nodes.y[i]);
        let (nx, ny) = (nodes.nx[i], nodes.ny[i]);
        let (x1, y1) = (x0 + normal_scale * nx, y0 + normal_scale * ny);

        let dx = x1 - x0;
        let dy = y1 - y0;
        let len = (dx * dx + dy * dy).sqrt().max(1e-16);
        let (ux, uy) = (dx / len, dy / len);

        let cos_th = (30f64).to_radians().cos();
        let sin_th = (30f64).to_radians().sin();

        let rx2 = -ux * cos_th - ( uy) * sin_th;
        let ry2 = -ux * sin_th + ( uy) * cos_th;
        let h2 = (x1 + head_scale * rx2, y1 + head_scale * ry2);

        PathElement::new(vec![(x1, y1), h2], arrow_color)
    }))?;

    root.present()?;
    Ok(())
}