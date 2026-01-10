mod algo;

use algo::p::p;
use algo::hat::Hat;
use algo::point::Point;

use num::complex::Complex;
use std::io::Write;
use csv::ReaderBuilder;

fn main() {
    let distance = 0.16; // Set height between transducers and reflective plate
    let hat = Hat::new(256.0, distance, false); // Initialize phase solver

    // Load in .csv file
    let base = env!("CARGO_MANIFEST_DIR");
    let filename = format!("{}/pos/{}.csv", base, "test");

    // Read trajectory CSV file
    let mut reader = ReaderBuilder::new().has_headers(true).from_path(&filename).expect("Could not open trajectory data file!");

    println!("Loading trajectory data...");

    // Define a "frame", containing a time value and transducer phase data
    #[derive(Debug)]
    struct Frame {
        t: f32,
        x: f32,
        y: f32,
        z: f32,
        phases: Vec<f32>,
    }

    let mut frames: Vec<Frame> = Vec::new();

    for result in reader.records() {
        let rec = result.expect("Could not parse trajectory data!");

        let t: f32 = rec[0].parse().unwrap();
        let x: f32 = rec[1].parse().unwrap();
        let y: f32 = rec[2].parse().unwrap();
        let z_raw: f32 = rec[3].parse().unwrap();

        let z: f32 = z_raw * distance;

        let cps = vec![
            Point::new(x, y, z + 0.00428 - 0.00428 / 2.0),
            Point::new(x, y, z - 0.00428 - 0.00428 / 2.0)
        ];

        let phases = hat.run_hat(&cps);

        frames.push(Frame { t, x, y, z, phases });
    }

    println!("Loaded {} frames, with solved phases", frames.len());

    // Set simulation range and spacing
    let x0 = 0.02;
    let y0 = 0.08;
    let z0 = 0.02;

    let xsize = 0.14 - 0.02;
    let ysize = 0.0;
    let zsize = 0.14 - 0.02;

    let nx = 80;
    let ny = 1;
    let nz = 80;

    let nx_sep = xsize / nx as f32;
    let ny_sep = ysize / ny as f32;
    let nz_sep = zsize / nz as f32;

    // Print simulation range (optional, for verification)
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

    // Initialize transducers
    let transducers: Vec<Point> = hat.transducers;
    let reflected_transducers: Vec<Point> = transducers
        .iter()
        .map(|p| Point {
            x: p.x,
            y: p.y,
            z: -p.z,
        })
        .collect();

    let mut all_field = vec![vec![vec![vec![Complex::new(0.0, 0.0); nz]; ny]; nx]; frames.len()];

    // Iterate over frames
    for (index, frame) in frames.iter().enumerate() {
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
                        // field[x][y][z] += phases[i] * p(r, theta, 0.0);
                        field[x][y][z] += p(r, theta, 0.0) * Complex::from_polar(1.0, frame.phases[i]);
                    }

                    for i in 0..reflected_transducers.len() {
                        let vec_r = point - reflected_transducers[i];
                        let r = vec_r.norm();
                        let theta = (vec_r.z / r).acos();
                        // field[x][y][z] += phases[i] * p(r, theta, 0.0);
                        field[x][y][z] += p(r, theta, 0.0) * Complex::from_polar(1.0, frame.phases[i]);
                    }
                }
            }
        }

        all_field[index] = field; // Note: times not saved (need to adjust to incorporate this)
    }

    let s = serde_pickle::to_vec(&all_field, Default::default()).unwrap();
    let mut file = std::fs::File::create("field.pickle").unwrap();
    file.write_all(&s).unwrap();
}