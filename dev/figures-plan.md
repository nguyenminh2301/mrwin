# Figure plan for the `mrwin` papers (01–04) + a new figure type

Status: plan (2026-06-30). Companion to `dev/research-frontier-roadmap.md`. Goal: design
figures that make the **novel** content of each paper legible, and invent **one wholly new
figure type** for the programme's unifying object — *causal inference on a degree-2 (pairwise)
U-statistic functional of potential-outcome distributions, on a hierarchical composite
endpoint, identified by an instrument.*

## Why the standard figures are not enough

The off-the-shelf figures each paper would otherwise reach for — the MR scatter
(β_GX vs β_GY), the forest plot, the funnel/Egger plot, the Kaplan–Meier curve — were all
built for **degree-1 (mean) effects on a single endpoint**. None of them can show the four
things that make this programme new:

1. **The pairwise (win) structure** — the estimand is built from `h(O_i,O_j)`, a *comparison*,
   not a mean. A KM curve or a forest plot has no place to put the pair.
2. **The hierarchical resolution** — a win is decided at the first priority tier (death →
   hospitalization → recurrence) that separates the pair. Which tier drives the effect is a
   first-order clinical question that no MR/survival figure expresses.
3. **IV identification of a pairwise functional** — the reduced form is `E[h | instrument
   contrast]`, a curve no MR scatter draws.
4. **The strong → weak → degenerate identification continuum** (Paper 03) — the AR set's width
   and *shape* (bounded / unbounded / disjoint) is the whole point and has no standard glyph.

The per-paper figures below cover 1–3 case by case; the new **Concordance–Contrast (C²)
portrait** (§5) unifies all four in a single figure, and the **Identification Phase Diagram**
(§6) is a second new figure dedicated to point 4.

## 1. Paper 01 — scalable, calibrated DS-CWR (one-sample)

The novel claims are *calibration* (the R²/Fieller lesson: SE-ratio ≈ 1, coverage at nominal)
and *scalability* (subquadratic kernel). Figures:

- **F1.1 Calibration band.** Across the validation grid (effect size × N × censoring) plot, per
  cell, the **estimated-SE / empirical-SD ratio** and the **CI coverage** as a paired 2-row
  heatmap (ratio on top, coverage below) with the nominal lines (1.0, 0.95) marked. This is the
  figure that would have caught the mis-calibrated bivariate-Δ CI; it *shows* the fix.
  Source: `inst/spec/validation-findings.md` + `tools/validation-scripts/paper01_mc_grid.R`.
- **F1.2 Wall-clock scaling.** Runtime vs N on log–log, fast kernel vs the O(n²) reference, with
  the fitted slopes annotated (≈1 vs 2). Source: the kernel benchmark.
- **F1.3 Operating characteristics.** Type-I (at the null) and power vs effect size at fixed N,
  one panel per regime (weak-IV, pleiotropy). Standard but needed for a methods paper.

## 2. Paper 02 — two-sample win-odds MR (IVW)

Novel claim: a per-SNP **win-odds** coefficient pooled by IVW recovers the causal win gradient.

- **F2.1 Win-odds MR scatter (augmented).** Each point a SNP: x = β_GX, y = δ (win-odds slope),
  whisker = ±SE on both axes, IVW slope line = pooled gradient, oracle gradient as a dashed
  reference. **New twist standard MR scatter can't do:** render each SNP as a small *win/loss/tie
  stacked bar* (its outcome-side composition) instead of a dot — so the reader sees that the
  pooled slope is built from pairwise wins, not a mean shift.
- **F2.2 Heterogeneity / pleiotropy funnel.** Per-SNP ratio δ/β_GX vs instrument strength
  (β_GX/SE), with the IVW estimate and the Cochran-Q implied funnel; asymmetry = directional
  pleiotropy. Source: `R/twosample.R` + `tools/twosample-*.R`.

## 3. Paper 03 — weak-IV-robust AR (degenerate-U-statistic)

Novel claims: AR set **uniformly valid** across strong→weak→degenerate; the degeneracy boundary.

- **F3.1 AR curve gallery.** AR(b) vs b with the χ²_L threshold line; the accepted set shaded.
  Three stacked panels — strong (tight bounded set), weak (wide/half-line), degenerate
  (set = whole line) — so the reader sees the CI *honestly* widen rather than the Wald CI
  falsely staying tight. Source: `R/twosample.R::mrwin_winmr_ar`, `tools/.../p3-*-probe.R`.
- **F3.2 Identification Phase Diagram** — see §6 (new figure).

## 4. Paper 04 — within-family win-MR

Novel claims: the **statistic-equals-design identity**; restricted-U (sibship) vs complete-U
(trio); robustness to the **PC-irreducible** confounders (dynastic, assortative mating); the
intrinsic weak-instrument regime that mandates AR.

- **F4.1 Confounder-robustness bars.** Bias vs oracle for {naive, PC-adjusted, within-family}
  grouped by {none, stratification, dynastic, AM}. Reproduces the corrected-scope story: PCs fix
  stratification but not dynastic/AM; within-family fixes all. Source:
  `tools/validation-scripts/within-family-winmr-probe.R` + `within-family-ar-inference.R`.
- **F4.2 Concordance–Contrast (C²) portrait** — see §5 (the flagship new figure). This is the
  natural Paper-04 figure: it draws the statistic-equals-design identity directly.
- **F4.3 Weak-IV ladder.** Per-SNP F-statistic vs IVW bias vs AR coverage across N, showing the
  attenuation vanish and the AR hold. Source: `within-family-twosample-winmr.R` (the F-ladder).

---

## 5. NEW FIGURE TYPE — the Concordance–Contrast (C²) portrait

**One glyph for "an instrument identifying a pairwise win on a hierarchical outcome."** It does
not exist in the MR literature (which plots means) or the win-statistics literature (which plots
win *ratios*, never against an instrument). It is the reduced form of a U-statistic-valued IV.

### Construction
Take the pairs that define the estimator (within-family sib pairs for Paper 04; any pairs for
02/03). For each pair (i,j):
- **x = instrument contrast** `ΔZ_ij` (sib genetic contrast `Z_i−Z_j`, or the trio mid-parent
  residual contrast `g_i−g_j`, or the per-SNP dosage contrast).
- **y = win outcome** `h(O_i,O_j) ∈ {−1,0,+1}`.
Then plot the **binned reduced form** `Ê[h | ΔZ]` (running mean ± MC band) as the spine.

### The three overlays (what makes it a *portrait*, not a scatter)
1. **Identification spine + AR fan.** The slope of `Ê[h|ΔZ]` through the origin is the outcome
   reduced form; divided by the exposure reduced form `Ê[ΔX|ΔZ]` it is the causal win gradient.
   Draw the point-estimate line, then a **fan of admissible slopes** = the AR confidence set
   mapped into this plane. Fan **angular width is the identification strength**, read directly:
   a pencil = strong instrument, a wide wedge = weak, a fan spanning all positive slopes =
   degenerate/unidentified. This is the strong→weak→degenerate continuum *in one picture*.
2. **De-confounding overlay.** Plot the **naive / between-family** pairs' `Ê[h|ΔZ]` in a second
   colour. Under dynastic/AM confounding it is tilted/shifted off the oracle slope (dashed); the
   **within-family** spine sits on it. The gap between the two clouds *is* the confounding the
   design removes — the statistic-equals-design identity made visual.
3. **Hierarchical resolution ribbon.** Below the main panel, a stacked area across the same ΔZ
   bins showing the fraction of *decided* wins resolved at tier 1 (death) / tier 2
   (hospitalization) / tier 3 (recurrence). Reveals **which endpoint drives the genetic win
   gradient** — a clinical reading no existing MR or win figure offers.

### Why it is useful and new
- It is the **only** figure that simultaneously shows the pairwise estimand, the IV identifying
  it, the identification *strength* (AR fan), the de-confounding, and the hierarchical mechanism.
- It degrades gracefully: drop overlay 2 → a generic two-sample win-MR reduced form (Paper 02);
  keep only the fan → a Paper-03 identification view; keep only the ribbon → a clinical
  decomposition. So one new figure type serves Papers 02–04.
- Prototype: `tools/figures/c2-portrait.R` (public package + simulator only).

## 6. SECOND new figure — the Identification Phase Diagram (Paper 03)

A 2-D map: x = instrument strength (concentration parameter `μ²` or mean per-SNP F-stat),
y = degeneracy (`ζ₂ / (n·ζ₁)`, distance to the U-statistic boundary). Colour each cell by the
**state of the AR set** (bounded / unbounded half-line / disjoint / empty) and overlay
**coverage contours** (the 0.95 isocline). It shows, at a glance, the region where the Wald CI
under-covers and the AR remains valid — the "uniform validity across the continuum" theorem as a
phase portrait. Source: the P3 degeneracy probes (`tools/.../p3-degeneracy-*-probe.R`).

## Implementation notes
- Base R graphics + cairo PNG/PDF (ggplot2 is not installed in the dev container); reproducible,
  fixed seeds, public package + simulator only — no confidential data, consistent with §5 of
  `CLAUDE.md`. Each figure ships a self-contained generator under `tools/figures/`.
