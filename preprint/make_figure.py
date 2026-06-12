"""Generate the scaling figure for the preprint from scaling_k3.csv."""
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).parent
rows = list(csv.DictReader((here / "scaling_k3.csv").open()))

n = [int(r["n_per_group"]) for r in rows]
fast = [float(r["fast_seconds"]) for r in rows]
brute_n = [int(r["n_per_group"]) for r in rows if r["brute_seconds"]]
brute = [float(r["brute_seconds"]) for r in rows if r["brute_seconds"]]

fig, ax = plt.subplots(figsize=(6.4, 4.6))

ax.loglog(brute_n, brute, "s-", color="#c0392b", label="dense kernel  $O(N^2)$", zorder=3)
ax.loglog(n, fast, "o-", color="#2c6fbb", label="fast algorithm  $O(N\\log^2 N)$", zorder=3)

# Reference guide lines anchored at the first fast point.
n0, t0 = n[0], fast[0]
xs = [n[0], n[-1]]
lin = [t0 * (x / n0) for x in xs]                          # slope 1
quad = [t0 * (x / n0) ** 2 for x in xs]                    # slope 2
ax.loglog(xs, lin, "--", color="gray", lw=1, label="slope 1 (linear)")
ax.loglog(xs, quad, ":", color="gray", lw=1, label="slope 2 (quadratic)")

# Fitted exponent annotation.
logn = [math.log(x) for x in n]
logt = [math.log(x) for x in fast]
nbar = sum(logn) / len(logn)
tbar = sum(logt) / len(logt)
slope = sum((a - nbar) * (b - tbar) for a, b in zip(logn, logt)) / \
        sum((a - nbar) ** 2 for a in logn)

ax.set_xlabel("sample size per group, $N$")
ax.set_ylabel("wall-clock time (s)")
ax.set_title("Hierarchical win/loss: dense vs. fast ($K=3$, weighted, censored)")
ax.text(0.04, 0.93, f"fitted fast-path exponent: {slope:.2f}",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.8))
ax.grid(True, which="both", ls=":", alpha=0.4)
ax.legend(frameon=False, fontsize=9, loc="lower right")
fig.tight_layout()

out = here / "fig_scaling.png"
fig.savefig(out, dpi=200)
print(f"Wrote {out} (fitted exponent {slope:.3f})")
