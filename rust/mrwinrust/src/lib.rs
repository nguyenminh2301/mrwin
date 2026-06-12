//! extendr bindings exposing the sub-quadratic win/loss kernel to R.
//!
//! Build requires an R toolchain (libR) and the extendr crates; it is not
//! compiled in R-less environments. The numerical core lives in the
//! dependency-free `mrwinkernel` crate, which is tested independently.
//!
//! R sees one function, `rust_pair_win_loss(th, dh, wh, tl, dl, wl, k)`,
//! returning `c(wins, losses)`. Times/weights arrive as doubles; status
//! vectors arrive as doubles in {0,1} and are cast to bytes.

use extendr_api::prelude::*;
use mrwinkernel::pair_win_loss;

/// Weighted win/loss totals for one high-vs-low stratum pair.
///
/// `th`/`dh` are length `n_high * k` (row-major), `tl`/`dl` length
/// `n_low * k`, `wh`/`wl` per-individual weights, `k` the number of priorities.
#[extendr]
fn rust_pair_win_loss(
    th: &[f64],
    dh: &[f64],
    wh: &[f64],
    tl: &[f64],
    dl: &[f64],
    wl: &[f64],
    k: i32,
) -> Doubles {
    let k = k as usize;
    let dh_u: Vec<u8> = dh.iter().map(|&v| v as u8).collect();
    let dl_u: Vec<u8> = dl.iter().map(|&v| v as u8).collect();
    let (wins, losses) = pair_win_loss(th, &dh_u, wh, tl, &dl_u, wl, k);
    Doubles::from_values([wins, losses])
}

extendr_module! {
    mod mrwinrust;
    fn rust_pair_win_loss;
}
