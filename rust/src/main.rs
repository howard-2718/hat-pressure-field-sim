mod algo;

use algo::p::p;
use algo::hat::Hat;
use algo::point::Point;
fn main() {
    let h = Hat::new(5.0);
    let pt = Point::new(1.0, 2.0);
    let val = p(pt.x as f32, 0.0, 0.0);

    println!("Hat size: {}, p(x): {}", h.size, val);
}