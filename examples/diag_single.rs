use mfs::prelude::*;
fn main() -> Result<()> {
    for ord in [2usize, 3, 4] {
        for rl in [20.0, 22.0] {
            let d = FilterDesign::prototype(ord, rl).synthesize()?;
            let m = BandPassMapping::new(15.0e9, 0.2e9)?;
            let grid = FrequencyGrid::linspace(14.8e9, 15.2e9, 401)?;
            let resp = ResponseSolver.evaluate(d.matrix(), &grid, &m)?;
            let peak = resp.samples.iter().map(|s| s.s21_db()).fold(f64::NEG_INFINITY, f64::max);
            println!("order {ord} rl {rl}: peak S21 = {peak:.4} dB");
        }
    }
    Ok(())
}
