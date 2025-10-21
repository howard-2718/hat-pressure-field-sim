use num::complex::Complex;
use std::f32::consts::PI;

use scilib::math::bessel;

const WAVE_LENGTH: f32 = 343.0 / 40000.0; // Speed of sound: 343 m/s
const OMEGA: f32 = 2.0 * PI * WAVE_LENGTH;
const K: f32 = 2.0 * PI / WAVE_LENGTH;

const P_0: f32 = 1.293; // Density of air
const EMITTER_RADIUS: f32 = 0.005; // Radius of the emitter: 5 mm

// Far field piston-source model (https://jontallen.ece.illinois.edu/uploads/473.F18/Lectures/Chapter_7b.pdf)
pub fn p(r: f32, theta: f32, t: f32) -> Complex<f32> {
    let sin_theta = f32::sin(theta);

    let amp: Complex<f32> = Complex::new(0.0, 1.0) * OMEGA * P_0 * EMITTER_RADIUS.powi(2) / (2.0 * r) * Complex::new(0.0, OMEGA * t + K * r).exp();
    
    if sin_theta == 0.0 {
        amp
    } 
    else {
        amp * 2.0 * bessel::j_n(1, (K * EMITTER_RADIUS * sin_theta) as f64) as f32 / (K * EMITTER_RADIUS * sin_theta)
    }
}