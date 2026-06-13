"""
Generate cross-language parity fixtures for the Rust core crate.

Writes a plain-text fixture (std-parseable, no serde needed) of random
two-group inputs and the win/loss totals from the Python reference
(`fast_kernel.fast_pair_win_loss`), which is itself pinned to the dense kernel.
The Rust core reads this file and must reproduce the same totals
(bit-exact unweighted; <=1e-9 weighted).

Format:
    <n_cases>
    repeated per case:
        K nH nL weighted(0/1)
        TH: nH*K floats (row-major)
        DH: nH*K ints
        wH: nH floats
        TL: nL*K floats
        DL: nL*K ints
        wL: nL floats
        wins losses
Run: PYTHONPATH=python python3 -m p1_engine_v5.gen_parity_fixtures <out>
"""
from __future__ import annotations

import sys
import numpy as np

from .fast_kernel import fast_pair_win_loss


def _fmt_floats(a):
    return " ".join(repr(float(x)) for x in np.asarray(a).ravel())


def _fmt_ints(a):
    return " ".join(str(int(x)) for x in np.asarray(a).ravel())


def main(out_path):
    rng = np.random.default_rng(20260612)
    cases = []
    cid = 0
    for K in (1, 2, 3, 4):
        for _ in range(6):
            weighted = cid % 2
            nH = int(rng.integers(5, 60))
            nL = int(rng.integers(5, 60))
            ties = cid % 3 == 0
            if ties:
                TH = rng.integers(1, 5, size=(nH, K)).astype(float)
                TL = rng.integers(1, 5, size=(nL, K)).astype(float)
            else:
                TH = rng.exponential(1.0, size=(nH, K))
                TL = rng.exponential(1.0, size=(nL, K))
            DH = rng.integers(0, 2, size=(nH, K))
            DL = rng.integers(0, 2, size=(nL, K))
            if weighted:
                wH = rng.exponential(1.0, size=nH)
                wL = rng.exponential(1.0, size=nL)
                wins, losses, _ = fast_pair_win_loss(TH, DH, TL, DL, wH, wL)
            else:
                wH = np.ones(nH)
                wL = np.ones(nL)
                wins, losses, _ = fast_pair_win_loss(TH, DH, TL, DL)
            cases.append((K, nH, nL, weighted, TH, DH, wH, TL, DL, wL,
                          float(wins), float(losses)))
            cid += 1

    lines = [str(len(cases))]
    for (K, nH, nL, w, TH, DH, wH, TL, DL, wL, wins, losses) in cases:
        lines.append(f"{K} {nH} {nL} {w}")
        lines.append(_fmt_floats(TH))
        lines.append(_fmt_ints(DH))
        lines.append(_fmt_floats(wH))
        lines.append(_fmt_floats(TL))
        lines.append(_fmt_ints(DL))
        lines.append(_fmt_floats(wL))
        lines.append(f"{repr(wins)} {repr(losses)}")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Wrote {len(cases)} parity cases to {out_path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "parity_fixtures.txt")
