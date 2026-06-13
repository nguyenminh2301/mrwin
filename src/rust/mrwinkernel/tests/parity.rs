//! Cross-language parity: the Rust core must reproduce the Python reference
//! (`p1_engine_v5.fast_kernel`) win/loss totals exactly (bit-exact unweighted,
//! <= 1e-9 weighted). Regenerate the fixture with:
//!   PYTHONPATH=python python3 -m p1_engine_v5.gen_parity_fixtures \
//!       src/rust/mrwinkernel/tests/parity_fixtures.txt

use mrwinkernel::pair_win_loss;
use std::fs;
use std::path::Path;

fn parse_floats(line: &str) -> Vec<f64> {
    line.split_whitespace()
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>().unwrap())
        .collect()
}

fn parse_ints(line: &str) -> Vec<u8> {
    line.split_whitespace()
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<u8>().unwrap())
        .collect()
}

#[test]
fn matches_python_reference() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/parity_fixtures.txt");
    let text = fs::read_to_string(&path).expect("fixture file present");
    let mut lines = text.lines();
    let n_cases: usize = lines.next().unwrap().trim().parse().unwrap();

    let mut checked = 0;
    for case in 0..n_cases {
        let header = parse_floats(lines.next().unwrap());
        let (k, n_h, n_l, weighted) = (
            header[0] as usize,
            header[1] as usize,
            header[2] as usize,
            header[3] as usize == 1,
        );
        let th = parse_floats(lines.next().unwrap());
        let dh = parse_ints(lines.next().unwrap());
        let wh = parse_floats(lines.next().unwrap());
        let tl = parse_floats(lines.next().unwrap());
        let dl = parse_ints(lines.next().unwrap());
        let wl = parse_floats(lines.next().unwrap());
        let expect = parse_floats(lines.next().unwrap());
        let (exp_w, exp_l) = (expect[0], expect[1]);

        assert_eq!(th.len(), n_h * k);
        assert_eq!(tl.len(), n_l * k);

        let (w, l) = pair_win_loss(&th, &dh, &wh, &tl, &dl, &wl, k);
        if weighted {
            assert!((w - exp_w).abs() < 1e-9, "case {case}: wins {w} vs {exp_w}");
            assert!((l - exp_l).abs() < 1e-9, "case {case}: loss {l} vs {exp_l}");
        } else {
            assert_eq!(w, exp_w, "case {case}: wins (bit-exact)");
            assert_eq!(l, exp_l, "case {case}: loss (bit-exact)");
        }
        checked += 1;
    }
    assert_eq!(checked, n_cases);
    assert!(n_cases >= 20, "expected a reasonable number of parity cases");
}
