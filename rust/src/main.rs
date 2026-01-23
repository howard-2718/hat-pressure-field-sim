mod algo;

use algo::p::p_modified;
use algo::hat::Hat;
use algo::point::Point;

use num::complex::Complex;
use std::io::Write;
use csv::ReaderBuilder;

use std::time::Instant;

fn main() {
    let timer = Instant::now();

    let distance = 0.16; // Set height between transducers and reflective plate
    let hat = Hat::new(256.0, distance, false); // Initialize phase solver

    // Load in .csv file
    let base = env!("CARGO_MANIFEST_DIR");
    let filename = format!("{}/pos/{}.csv", base, "test2");

    // Read trajectory CSV file
    let mut reader = ReaderBuilder::new().has_headers(true).from_path(&filename).expect("Could not open trajectory data file!");

    println!("\nLoading trajectory data...");

    // Define a "frame", containing a time value and transducer phase data
    #[derive(Debug)]
    struct Frame {
        t: f32,
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

        frames.push(Frame { t, phases });
    }

    println!("Loaded {} frames, with solved phases", frames.len());
    println!("Current elapsed time: {:?}", timer.elapsed());

    // Set simulation range and spacing
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

    // Print simulation range (optional, for verification)
    println!("\nSimulation spatial range:");

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
        "z: {}, {}\n",
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

    // Initialize points at which to construct pressure field
    let xs: Vec<f32> = (0..nx)
        .map(|x| x as f32 * nx_sep + nx_sep / 2.0 + x0)
        .collect();

    let ys: Vec<f32> = (0..ny)
        .map(|y| y as f32 * ny_sep + ny_sep / 2.0 + y0)
        .collect();

    let zs: Vec<f32> = (0..nz)
        .map(|z| z as f32 * nz_sep + nz_sep / 2.0 + z0)
        .collect();

    let mut all_field = vec![vec![vec![vec![Complex::new(0.0, 0.0); nz]; ny]; nx]; frames.len()];
    let mut all_times = vec![0.0; frames.len()];

    // Iterate over frames
    for (index, frame) in frames.iter().enumerate() {
        println!("Currently solving for frame {}", index + 1);

        let mut field = vec![vec![vec![Complex::new(0.0, 0.0); nz]; ny]; nx];

        // Precompute complex phase factors
        let phase_factors: Vec<Complex<f32>> =
            frame.phases.iter()
            .map(|&phi| Complex::from_polar(1.0, phi))
            .collect();

        for x in 0..nx {
            let px = xs[x];

            for y in 0..ny {
                let py = ys[y];

                for z in 0..nz {
                    let pz = zs[z];

                    let point = Point::new(px, py, pz);

                    let mut acc = Complex::new(0.0, 0.0);

                    for i in 0..transducers.len() {
                        let vec_r = point - transducers[i];
                        let r = vec_r.norm();
                        let cos_theta = vec_r.z / r;
                        acc += p_modified(r, cos_theta, 0.0) * phase_factors[i];
                    }

                    for i in 0..reflected_transducers.len() {
                        let vec_r = point - reflected_transducers[i];
                        let r = vec_r.norm();
                        let cos_theta = vec_r.z / r;
                        acc += p_modified(r, cos_theta, 0.0) * phase_factors[i];
                    }

                    field[x][y][z] = acc;
                }
            }
        }

        all_field[index] = field;
        all_times[index] = frame.t;
    }

    println!("Current elapsed time: {:?}", timer.elapsed());
    println!("\nSaving fields to pickles...");

    // Save fields to field.pickle
    let s = serde_pickle::to_vec(&all_field, Default::default()).unwrap();
    let mut file = std::fs::File::create("field.pickle").unwrap();
    file.write_all(&s).unwrap();

    // Save times to time.pickle
    let s_2 = serde_pickle::to_vec(&all_times, Default::default()).unwrap();
    let mut file_2 = std::fs::File::create("time.pickle").unwrap();
    file_2.write_all(&s_2).unwrap();

    println!("\nAll done!\nCurrent elapsed time: {:?}", timer.elapsed());
}