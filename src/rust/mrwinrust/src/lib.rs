//! extendr bindings exposing the sub-quadratic win/loss kernel to R.
//!
//! Build requires an R toolchain (libR) and the extendr crates; it is not
//! compiled in R-less environments. The numerical core lives in the
//! dependency-free `mrwinkernel` crate, which is tested independently.
//!
//! R sees one function, `rust_pair_win_loss(th, dh, wh, tl, dl, wl, k)`,
//! returning `c(wins, losses)`. All numeric inputs arrive as R doubles; the
//! status vectors hold 0/1 and are cast to bytes.

use extendr_api::prelude::*;
use mrwinkernel::pair_win_loss;

/// Weighted win/loss totals for one high-vs-low stratum pair.
///
/// `th`/`dh` are length `n_high * k` (row-major), `tl`/`dl` length
/// `n_low * k`, `wh`/`wl` per-individual weights, `k` the number of priorities.
/// Returns `c(wins, losses)`.
#[extendr]
fn rust_pair_win_loss(
    th: Vec<f64>,
    dh: Vec<f64>,
    wh: Vec<f64>,
    tl: Vec<f64>,
    dl: Vec<f64>,
    wl: Vec<f64>,
    k: i32,
) -> Vec<f64> {
    let k = k as usize;
    let dh_u: Vec<u8> = dh.iter().map(|&v| v as u8).collect();
    let dl_u: Vec<u8> = dl.iter().map(|&v| v as u8).collect();
    let (wins, losses) = pair_win_loss(&th, &dh_u, &wh, &tl, &dl_u, &wl, k);
    vec![wins, losses]
}

extendr_module! {
    mod mrwinrust;
    fn rust_pair_win_loss;
}
