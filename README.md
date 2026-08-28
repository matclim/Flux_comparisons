# muonflux_compare
 
**Reference vs test.** The first sample is the reference; everything is
reported as `(test - reference) / reference`.
 
 
Are the muons delivered by two FairShip productions compatible — in their
relative amounts and in their characteristics — when the overall rate is a
generation setting and must not be tested?
 
Three tests, one multi-page PDF, one ROOT file.
 
| # | test | what it answers |
|---|------|-----------------|
| 1 | μ⁺/μ⁻ ratio | a relative amount with no rate scale in it |
| 2 | flux fractions in phase-space regions | what share of the flux sits where |
| 3 | weighted two-sample AD / CvM / KS with a permutation null | do the 1-D shapes agree, tails included |
 
## Build
 
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
ctest --test-dir build --output-on-failure
```
 
Needs ROOT ≥ 6.22 with RDataFrame. If ROOT was built with C++20, pass
`-DCMAKE_CXX_STANDARD=20`; mismatched standards show up as missing symbols.
 
## Run
 
```bash
build/muonflux_compare A.root B.root
build/muonflux_compare @listA.txt @listB.txt --label-a nominal --label-b variant
build/muonflux_compare A.root --closure          # run this first
./scripts/run_compare.sh 'prodA/*.root' 'prodB/*.root'
```
 
`--help` lists every option.
 
### Outputs
 
Everything lands in `--outdir` (default `outputs/`):
 
| file | contents |
|------|----------|
| `PREFIX_NN_<name>.pdf` | **one file per plot**, numbered and named |
| `PREFIX_all.pdf` | the same pages as one book -- only with `--combined-pdf` |
| `PREFIX.root` | `A_<var>_<tag>`, `B_<var>_<tag>`, `R_<var>_<tag>` for restyling |
| `PREFIX.txt` | the console output verbatim, so runs are reproducible on paper |
 
`--png` additionally writes a `.png` beside every `.pdf`.
 
### Calibration: `--closure-scan N`
 
`--closure` tells you whether the pipeline invents differences. `--closure-scan
N` tells you whether the p-values are *correctly calibrated*, which is a
stronger and more useful statement.
 
```bash
build/muonflux_compare A.root --closure-scan 20 --nperm 200 --max-events-per-file 50000
```
 
It splits sample A into halves N times with different seeds. Each split is an
independent draw, so it yields an independent p-value per variable. Under the
null those p-values must be uniform on (0,1), and the P-P plot shows whether
they are: sorted p-values on the y-axis against `k/(n+1)` on the x-axis, with
the diagonal marked.
 
* **on the diagonal** — calibrated, trust the p-values
* **above the diagonal** — conservative; real differences will be missed
* **below the diagonal** — anti-conservative; the p-values are not trustworthy
 
The table alongside gives mean p (should be ~0.50), the fraction below 0.05
(should be ~0.05), and the KS distance from uniform.
 
Note the p-values are independent **across splits**, not across variables:
`log10_p`, `log10_pz` and `asinh_pT` all derive from the same momenta, so they
move together. That is why the scan replicates over seeds rather than counting
variables.
 
### Species and input branch
 
The tool is not muon-specific. Two options select what is compared:
 
| flag | meaning |
|------|---------|
| `--pdg LIST` | `\|PDG\|` codes to accept. Default `13` (muons). `--pdg 14` for muon neutrinos, `--pdg 12,14,16` for all three flavours. |
| `--branch NAME` | the scoring-plane branch to read. Default `PlaneHAPoint`. Any branch exposing `fX, fY, fPx, fPy, fPz, fPdgCode, fTrackID` works, since those member names are assumed common to all such branches. |
 
Neutrinos are recorded at a scoring plane like any other track, so the
plane-hit path needs no change. Two things do adapt automatically:
 
* **`q·x` is dropped.** For a neutral species it is identically zero, the
  equal-weight binning degenerates and `A²` is undefined. Requesting it
  explicitly via `--vars` is an error rather than a silent no-op.
* **The sign split becomes `nu` / `nubar`.** The two classes are defined by the
  sign of the PDG code, not by charge, so `--split-sign` (formerly
  `--split-charge`, still accepted) works for both families. Index 0 is the
  antilepton: `mu+` for muons, `nubar` for neutrinos. The charge ratio
  generalises to `nubar / nu`.
 
Mixing charged and neutral codes in one `--pdg` list is rejected: the sign
classes and `q·x` would mean different things for different members.
 
Two caveats for neutrinos. Flavour composition (`νe` vs `νμ` vs `ντ`) is not a
single-run comparison — the accumulator carries a sign index, not a flavour
index — so compare flavours with separate runs. And the weights of a
muon-focused production may not be meaningful for neutrinos if the biasing was
applied to muons only; run `muonflux_weight_audit --pdg 14` first, since the
per-entry `sum(w)` invariant applies to any species selection.
 
### Scale
 
All files of a sample are handed to one RDataFrame and read as a single dataset
with implicit MT -- not opened one after another. With a few hundred files that
is the difference between minutes and hours.
 
Nothing is stored per muon. Muons stream into fixed fine-binned histograms, one
per (unit, charge, variable), so **memory does not depend on how much data goes
in**: 256 units x 2 charges x 8 variables x 512 bins x 8 B = 8 MB per thread.
The driver prints the estimate before reading anything. A single split100 file
holds 36M muons; keeping them would have cost 2.3 GB, and a 300-file sample
would have been impossible.
 
Two approximations come with streaming, both reported at runtime:
 
* bin edges snap to the fine grid, so a cut at 100 GeV lands within one fine
  bin of it (raise `--nfine` to tighten);
* values outside the recorded range pile into the end bins. The range is fixed
  by a probe pass over the first entries of a few files and padded by 35%; the
  clamp count is printed and should be ~0.
 
Selftest case 7 checks the streamed path against the exact per-muon path on the
same data and units: `A^2` agrees to 0.2% with identical permutation p-values.
 
### Limiting the input
 
For quick iteration, or when one production is far larger than the other:
 
| flag | effect |
|------|--------|
| `--sample-frac F` | keep a deterministic random fraction of events |
| `--max-files N` | use only the first N files of each list |
 
`--sample-frac` hashes (file index, entry) rather than truncating, so the subset
is random, unbiased, reproducible across runs and thread counts, and does not
depend on what `rdfentry_` counts from across a chain boundary.
 
`0` means no limit. Entries are taken in file order, and files are read in list
order; once `--max-events` is reached the remaining files are not opened at
all, so the flag saves wall time as well as memory.
 
```bash
build/muonflux_compare @listA.txt @listB.txt --max-events-per-file 2000
build/muonflux_compare @listA.txt B.root --max-events-a 500000   # balance the two
```
 
Truncation is unbiased — nothing in a FairShip file orders events by physics —
but it cuts `N_eff`, and `N_eff` is the ceiling on every test here. A run with
limits prints `input is TRUNCATED` and records the settings on the summary
page. Treat those p-values as a preview and re-run unrestricted before quoting
anything.
 
One interaction to watch: `--max-files 3` drops you below `--min-files` (5), so
the bootstrap silently switches from files to random event blocks, and the
error bars stop seeing per-file correlations. The run header prints which unit
was used for each sample.
 
**Run `--closure` before you trust anything.** It splits sample A in half by
event and pushes the two halves through the identical pipeline. Every p-value
should be flat and every fraction pull should sit inside ±3σ. Anything that
fires there is a pipeline artefact — duplicated muons crossing event
boundaries, a per-file correlation, a binning pathology — not physics.
 
## The invariant everything rests on
 
The two productions are generated with deliberately different muon rates, so
`sum(w)` is a design choice and carries no physics. Every statistic here is
invariant under `w_A → a·w_A` and `w_B → b·w_B`:
 
* fractions and ratios are quotients of weighted sums **within one sample**;
* the shape tests normalise each sample to unit total weight internally, and
  the permutation null does the same for **every** relabelling, so no
  permutation can see the rate difference.
 
The total weight is printed once, for the record, and then used only as a
denominator. `test/core_selftest.cxx` asserts the invariant by rescaling one
sample by 1000 and requiring bit-identical statistics.
 
This matters more than it sounds. A classifier trained on the two samples with
raw weights will find the normalisation offset, because it is the easiest thing
in the data, and report a large separation that means nothing.
 
## Why Anderson–Darling and not the Fisher separation
 
The Fisher discriminant reports `d = (μ_B − μ_A)·c / √(cᵀ S_p c)` — a
difference-of-means statistic. Two distributions with identical means *and*
identical variances but different tails give `d ≈ 0`. Selftest case 3 builds
exactly that pair:
 
```
Fisher-style mean separation d = +0.0002   (blind)
AD  =   111.344   p < 0.005   z = +176.5
CvM =    15.720   p < 0.005   z = +124.4
KS  =     5.880   p < 0.005   z = +19.8
```
 
The AD weighting `1/(H(1−H))` is what buys the tail sensitivity; KS, a
supremum statistic, is nearly blind there. Both are computed, so you can see
which one is carrying the signal.
 
The usual objection — "AD assumes a Gaussian" — is about the *one-sample* AD
test against a named distribution. The **two-sample** AD statistic used here is
distribution-free, and in any case the null is built by permutation, so the
p-values are correct whatever the parent distribution is.
 
The real caveat is different: AD's tail weighting diverges at the edges, so a
handful of extreme-weight muons can dominate it. Check that the p-value is
stable across `--seed`, and if it is not, raise `--ad-clip` (try `1e-4`, then
`1e-3`) to drop the singular edge bins, or read the CvM column instead.
 
## Resampling units
 
Two different units, on purpose:
 
* **Permutation unit = event** (`--perm-unit`). Many units, so the null is well
  resolved. Muons from one event move together, which keeps splitmult copies
  from narrowing the null — provided the copies stay inside one tree entry.
  Selftest case 5 checks this. If your production writes copies across entries,
  switch to `--perm-unit file`.
* **Bootstrap unit = file**, or a random block of whole events when there are
  fewer than `--min-files` (default 5). Files are the natural correlated block
  in a production; bootstrapping over them propagates per-file correlations
  into the error bars automatically. That is what replaces the hand-tuned error
  inflation factor.
 
## Binning
 
Plots use the ROOT `kCool` palette; the two sample colours are drawn from 12%
and 88% of it rather than hard-coded, so overlays match anything else on the
page.
 
Bins hold equal shares of the **pooled weight**, not equal width. The tails get
wide bins with real statistics in them instead of a row of near-empty fine bins
that a ratio panel has to blank out. Plots divide by bin width, so what you see
is a density and the shape is undistorted; the ratio panel is unaffected,
because both samples are divided by the same widths.
 
Every variable also gets an **N_eff-per-bin page**. Read it first. If N_eff in
the high-momentum tail is a few tens, no statistic will resolve anything there
and the answer is more MC, not a better test.
 
## Reading the profile pages
 
`*_profile_<var>_<case>.pdf` answers *where* two samples differ, which no
p-value can.
 
* **Solid line, left axis** -- the signed CDF difference `F_test - F_ref`.
  Positive means the test sample has more of its flux below that value, i.e.
  it is softer there. Grey bands are the pointwise 68% / 95% envelopes from the
  same permutations that produce the p-value. Where the solid line leaves the
  band, the samples differ.
* **Dashed line, right axis** -- the cumulative share of `A^2`, running 0 to 1.
  The steep part is where the statistic is being generated. This localises
  better than the CDF difference, which is cumulative and therefore smears a
  local difference into a broad hump.
 
A single broad positive hump with a steadily rising cumulative means a uniform
shift; a step in the cumulative at one value means a localised feature. Those
two look identical in `A^2` alone.
 
The bands are **pointwise**, not simultaneous: with 40 bins about two will
stray outside the 95% band under the null, so read coherent runs of bins rather
than isolated excursions.
 
## Variables
 
All measured on the muon at the reference plane. No truth, no ancestry, no
production vertex.
 
`log10_p`, `log10_pz`, `asinh_pT`, `log10_theta`, `x`, `y`, `qx`, `r`
 
`qx = q·x` is the bending variable — in a magnetised shield the sign of the
horizontal displacement relative to the charge is where the information is, and
unsigned `x` throws it away. `--split-sign` additionally runs every test for each sign class
separately.
 
## Layout
 
```
include/muonflux/Core.h     ROOT-free numerics: EDF statistics, permutation
                            null, block bootstrap, quantile binning
include/muonflux/Config.h   variables, regions, reported quantities
include/muonflux/Load.h     RDataFrame loader (plane hits or MCTrack)
include/muonflux/Plot.h     overlay + ratio panels, nulls, pull summary
test/core_selftest.cxx      six assertions, no ROOT required
```
 
`Core.h` includes no ROOT, so `g++ -O2 -std=c++17 -Iinclude
test/core_selftest.cxx -o selftest -pthread` builds and runs the tests on any
machine. Same arrangement as `fisher_core.h`.
 
## Extending it
 
Add a region to `defaultRegions()` in `Config.h` and a matching entry to
`defaultQuantities()` (indices there are 1-based into the region list; column 0
is the sample total). A plain fraction is `{name, k, 0}`; a ratio of two
regions is `{name, k, m}`.
 
Adding a statistic: check it against the scale invariance above, and add the
check to the selftest.
