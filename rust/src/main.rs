mod algo;

use algo::p::p;
use algo::hat::Hat;
use algo::point::Point;

use num::complex::Complex;
use std::io::Write;

fn main() {
    let z = 0.16;
    let hat = Hat::new(256.0, z, false);

    let cps = vec![
        Point::new(0.08, 0.08, 0.08 + 0.00428 - 0.00428 / 2.0),
        Point::new(0.08, 0.08, 0.08 - 0.00428 - 0.00428 / 2.0)
    ];

    let phases = hat.run_hat(&cps);

    let transducers: Vec<Point> = hat.transducers;
    let reflected_transducers: Vec<Point> = transducers
        .iter()
        .map(|p| Point {
            x: p.x,
            y: p.y,
            z: -p.z,
        })
        .collect();

    let x0 = 0.02;
    let y0 = 0.08;
    let z0 = 0.02;

    let xsize = 0.14 - 0.02;
    let ysize = 0.0;
    let zsize = 0.14 - 0.02;

    let nx = 160;
    let ny = 1;
    let nz = 160;

    let nx_sep = xsize / nx as f32;
    let ny_sep = ysize / ny as f32;
    let nz_sep = zsize / nz as f32;
    let mut field = vec![vec![vec![Complex::new(0.0, 0.0); nz]; ny]; nx];

    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                let point = Point::new(
                    x as f32 * nx_sep + nx_sep / 2.0 + x0,
                    y as f32 * ny_sep + ny_sep / 2.0 + y0,
                    z as f32 * nz_sep + nz_sep / 2.0 + z0,
                );

                for i in 0..transducers.len() {
                    let vec_r = point - transducers[i];
                    let r = vec_r.norm();
                    let theta = (vec_r.z / r).acos();
                    field[x][y][z] += phases[i] * p(r, theta, 0.0);
                }

                for i in 0..reflected_transducers.len() {
                    let vec_r = point - reflected_transducers[i];
                    let r = vec_r.norm();
                    let theta = (vec_r.z / r).acos();
                    field[x][y][z] += phases[i] * p(r, theta, 0.0);
                }
            }
        }
    }

    println!(
        "\nx: {}, {}",
        nx_sep / 2.0 + x0,
        (nx as f32 + nx_sep / 2.0) * nx_sep + x0
    );
    println!(
        "y: {}, {}",
        ny_sep / 2.0 + y0,
        (ny as f32 + ny_sep / 2.0) * ny_sep + y0
    );
    println!(
        "z: {}, {}",
        nz_sep / 2.0 + z0,
        (nz as f32 + nz_sep / 2.0) * nz_sep + z0
    );

    let s = serde_pickle::to_vec(&field, Default::default()).unwrap();
    let mut file = std::fs::File::create("field.pickle").unwrap();
    file.write_all(&s).unwrap();
}