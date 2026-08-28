// Core.h -- ROOT-free numerics for the muon-flux two-sample comparison.
//
// Nothing in this header includes ROOT.  That is deliberate: it lets the whole
// statistical machinery be compiled and unit-tested with a plain compiler
// (see test/core_selftest.cxx), the same way fisher_core.h is arranged.
//
// SCALE INVARIANCE -- THE CENTRAL INVARIANT OF THIS CODE
// -----------------------------------------------------
// The two productions are generated with deliberately different muon rates, so
// sum(w) is a design choice and carries no physics.  Every statistic here is
// therefore invariant under  w_A -> a*w_A  and  w_B -> b*w_B  for any a,b > 0:
//
//   * fractions and ratios are quotients of weighted sums within one sample;
//   * the shape tests normalise each sample to unit total weight internally,
//     and the permutation null does the same for every relabelling, so no
//     permutation can ever "see" the rate difference.
//
// If you add a statistic, check it against that invariant.  The selftest does.
//
// RESAMPLING UNITS
// ----------------
// Two different units are used, on purpose:
//
//   * PERMUTATION unit = event.  Many units, so the null is well resolved.
//     Muons from one event (including splitmult copies, if they stay inside
//     one entry) move together, which is what keeps the null honest.
//   * BOOTSTRAP unit = file, or a random block of whole events when there are
//     too few files.  Files are the natural correlated block in a FairShip
//     production; using them propagates per-file correlations into the error
//     bars automatically, instead of the ad-hoc error inflation factor.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace fluxval {

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

inline double sumW(const std::vector<double>& w) {
  double s = 0.0;
  for (double v : w) s += v;
  return s;
}

/// Kish effective sample size (sum w)^2 / sum w^2 -- the real statistics.
inline double nEff(const std::vector<double>& w) {
  double s = 0.0, s2 = 0.0;
  for (double v : w) { s += v; s2 += v * v; }
  return (s2 > 0.0) ? s * s / s2 : 0.0;
}

inline double nEffFrom(double s, double s2) {
  return (s2 > 0.0) ? s * s / s2 : 0.0;
}

/// Run f(i) for i in [0,n) across `nthreads` threads (1 => serial, 0 => auto).
template <typename F>
inline void parallelFor(std::size_t n, int nthreads, F&& f) {
  unsigned nt = (nthreads > 0) ? static_cast<unsigned>(nthreads)
                               : std::max(1u, std::thread::hardware_concurrency());
  if (nt <= 1 || n < 2) {
    for (std::size_t i = 0; i < n; ++i) f(i);
    return;
  }
  nt = static_cast<unsigned>(std::min<std::size_t>(nt, n));
  std::vector<std::thread> pool;
  pool.reserve(nt);
  for (unsigned t = 0; t < nt; ++t) {
    pool.emplace_back([&, t] {
      for (std::size_t i = t; i < n; i += nt) f(i);
    });
  }
  for (auto& th : pool) th.join();
}

// ---------------------------------------------------------------------------
// One sample, flattened to muon rows
// ---------------------------------------------------------------------------

/// Column-major muon table.  `var[j]` holds variable j for every row.
/// `event` and `file` are the two candidate resampling units; both are
/// compacted to contiguous 0-based ids by compactUnits().
struct Sample {
  std::string label;
  int nvar = 0;
  std::vector<std::vector<double>> var;  // [nvar][rows]
  std::vector<double> w;                 // per-row weight, as read from the file
  std::vector<int> q;                    // +1 / -1
  std::vector<int> event;                // contiguous event id
  std::vector<int> file;                 // contiguous file id
  int nEvents = 0, nFiles = 0;

  std::size_t rows() const { return w.size(); }

  void compactUnits() {
    auto compact = [](std::vector<int>& v) {
      if (v.empty()) return 0;
      std::vector<int> u(v);
      std::sort(u.begin(), u.end());
      u.erase(std::unique(u.begin(), u.end()), u.end());
      for (int& x : v)
        x = static_cast<int>(std::lower_bound(u.begin(), u.end(), x) - u.begin());
      return static_cast<int>(u.size());
    };
    nEvents = compact(event);
    nFiles = compact(file);
  }

  /// Rescale every weight so that sum(w) == 1.  Called before any pooling:
  /// after this, an A row and a B row are on the same scale and may legally be
  /// exchanged by a permutation.  Returns the original total.
  double normaliseWeights() {
    const double s = sumW(w);
    if (!(s > 0.0)) throw std::runtime_error("sample " + label + " has sum(w) <= 0");
    for (double& x : w) x /= s;
    return s;
  }
};

/// Copy the rows of one charge.  Units are recompacted; the caller must
/// re-normalise the weights afterwards.
inline Sample filterCharge(const Sample& S, int charge) {
  if (charge == 0) return S;
  Sample O;
  O.label = S.label;
  O.nvar = S.nvar;
  O.var.assign(S.nvar, {});
  for (std::size_t i = 0; i < S.rows(); ++i) {
    if (S.q[i] != charge) continue;
    for (int j = 0; j < S.nvar; ++j) O.var[j].push_back(S.var[j][i]);
    O.w.push_back(S.w[i]);
    O.q.push_back(S.q[i]);
    O.event.push_back(S.event[i]);
    O.file.push_back(S.file[i]);
  }
  O.compactUnits();
  return O;
}

/// Split a sample into two halves along whole events.  Used for the closure
/// test: running A1 vs A2 through the identical pipeline must give a flat
/// p-value.  Any structure found there is a pipeline artefact, not physics.
inline std::pair<Sample, Sample> splitEvents(const Sample& S, uint64_t seed) {
  std::vector<uint8_t> side(S.nEvents);
  std::mt19937_64 rng(seed);
  for (int e = 0; e < S.nEvents; ++e) side[e] = static_cast<uint8_t>(rng() & 1u);

  Sample A, B;
  for (Sample* p : {&A, &B}) {
    p->nvar = S.nvar;
    p->var.assign(S.nvar, {});
  }
  A.label = S.label + " (half 1)";
  B.label = S.label + " (half 2)";
  for (std::size_t i = 0; i < S.rows(); ++i) {
    Sample& T = side[S.event[i]] ? B : A;
    for (int j = 0; j < S.nvar; ++j) T.var[j].push_back(S.var[j][i]);
    T.w.push_back(S.w[i]);
    T.q.push_back(S.q[i]);
    T.event.push_back(S.event[i]);
    T.file.push_back(S.file[i]);
  }
  A.compactUnits();
  B.compactUnits();
  return {A, B};
}

// ---------------------------------------------------------------------------
// Weighted quantiles and binning
// ---------------------------------------------------------------------------

/// Weighted quantiles of v.  `probs` must be sorted ascending in (0,1).
inline std::vector<double> weightedQuantiles(const std::vector<double>& v,
                                             const std::vector<double>& w,
                                             const std::vector<double>& probs) {
  const std::size_t n = v.size();
  std::vector<std::size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return v[a] < v[b]; });

  double tot = 0.0;
  for (double x : w) tot += x;
  std::vector<double> out(probs.size(), v.empty() ? 0.0 : v[idx.front()]);
  if (!(tot > 0.0) || n == 0) return out;

  double cum = 0.0;
  std::size_t p = 0;
  for (std::size_t k = 0; k < n && p < probs.size(); ++k) {
    cum += w[idx[k]];
    while (p < probs.size() && cum >= probs[p] * tot) out[p++] = v[idx[k]];
  }
  for (; p < probs.size(); ++p) out[p] = v[idx.back()];
  return out;
}

/// Bin edges placed so that each bin holds an equal share of the POOLED
/// weight.  Equal-N_eff-ish bins: the tails get wide bins with real statistics
/// in them instead of a row of empty fine bins that a ratio plot cannot use.
inline std::vector<double> equalWeightEdges(const std::vector<double>& v,
                                            const std::vector<double>& w, int nbins) {
  if (nbins < 2) nbins = 2;
  std::vector<double> probs;
  probs.reserve(nbins - 1);
  for (int k = 1; k < nbins; ++k) probs.push_back(double(k) / nbins);
  std::vector<double> inner = weightedQuantiles(v, w, probs);

  const auto mm = std::minmax_element(v.begin(), v.end());
  std::vector<double> e;
  e.reserve(nbins + 1);
  e.push_back(*mm.first);
  for (double x : inner) e.push_back(x);
  e.push_back(std::nextafter(*mm.second, *mm.second + 1.0));

  // Heavy ties (x exactly 0, say) collapse edges; drop the duplicates rather
  // than carry zero-width bins that make dH == 0.
  e.erase(std::unique(e.begin(), e.end()), e.end());
  if (e.size() < 3) {  // degenerate variable, fall back to a 2-bin split
    e = {*mm.first, 0.5 * (*mm.first + *mm.second),
         std::nextafter(*mm.second, *mm.second + 1.0)};
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }
  return e;
}

inline int binOf(const std::vector<double>& edges, double x) {
  const int nb = static_cast<int>(edges.size()) - 1;
  if (x < edges.front()) return 0;
  if (x >= edges.back()) return nb - 1;
  int k = static_cast<int>(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
  return std::max(0, std::min(nb - 1, k));
}

// ---------------------------------------------------------------------------
// Pooled, binned, unit-compressed representation for the shape tests
// ---------------------------------------------------------------------------

/// Compressed sparse rows: for each permutation unit, the list of (bin, weight)
/// it contributes.  Entries sharing a (unit,bin) are merged once, up front, so
/// a permutation costs one pass over the compressed entries rather than over
/// all muons.
struct PooledVar {
  int nbins = 0;
  int nUnits = 0;
  std::vector<double> edges;
  std::vector<int> off;       // nUnits+1
  std::vector<uint16_t> bin;  // nnz
  std::vector<double> wgt;    // nnz
  std::vector<double> uSw;    // per unit: sum w
  std::vector<double> uSw2;   // per unit: sum w^2
  std::vector<uint8_t> label; // per unit: 0 = A, 1 = B  (observed labelling)
};

/// Build the pooled representation of variable `j`.
/// `unitsAreFiles` picks the permutation unit; events are the default.
/// PRECONDITION: A.w and B.w have each been normalised to unit sum.
inline PooledVar buildPooled(const Sample& A, const Sample& B, int j, int nbins,
                             bool unitsAreFiles = false) {
  const auto& vA = A.var.at(j);
  const auto& vB = B.var.at(j);

  std::vector<double> vAll;
  std::vector<double> wAll;
  vAll.reserve(vA.size() + vB.size());
  wAll.reserve(vA.size() + vB.size());
  vAll.insert(vAll.end(), vA.begin(), vA.end());
  vAll.insert(vAll.end(), vB.begin(), vB.end());
  wAll.insert(wAll.end(), A.w.begin(), A.w.end());
  wAll.insert(wAll.end(), B.w.begin(), B.w.end());

  PooledVar P;
  P.edges = equalWeightEdges(vAll, wAll, nbins);
  P.nbins = static_cast<int>(P.edges.size()) - 1;

  const int nUA = unitsAreFiles ? A.nFiles : A.nEvents;
  const int nUB = unitsAreFiles ? B.nFiles : B.nEvents;
  P.nUnits = nUA + nUB;
  P.label.assign(P.nUnits, 0);
  for (int u = nUA; u < P.nUnits; ++u) P.label[u] = 1;

  // (unit, bin, w) triples, then sort+merge.
  struct Ent { int u; uint16_t b; double w; };
  std::vector<Ent> ent;
  ent.reserve(vA.size() + vB.size());
  const auto& uA = unitsAreFiles ? A.file : A.event;
  const auto& uB = unitsAreFiles ? B.file : B.event;
  for (std::size_t i = 0; i < vA.size(); ++i)
    ent.push_back({uA[i], static_cast<uint16_t>(binOf(P.edges, vA[i])), A.w[i]});
  for (std::size_t i = 0; i < vB.size(); ++i)
    ent.push_back({nUA + uB[i], static_cast<uint16_t>(binOf(P.edges, vB[i])), B.w[i]});

  std::sort(ent.begin(), ent.end(), [](const Ent& a, const Ent& b) {
    return a.u != b.u ? a.u < b.u : a.b < b.b;
  });

  P.off.assign(P.nUnits + 1, 0);
  P.uSw.assign(P.nUnits, 0.0);
  P.uSw2.assign(P.nUnits, 0.0);
  P.bin.reserve(ent.size());
  P.wgt.reserve(ent.size());

  // sum w^2 must be accumulated over MUONS, not over merged (unit,bin) cells:
  // N_eff is a property of the weights, not of the binning.
  for (const Ent& e : ent) { P.uSw[e.u] += e.w; P.uSw2[e.u] += e.w * e.w; }

  int cur = -1;
  for (std::size_t i = 0; i < ent.size(); ++i) {
    if (ent[i].u != cur) {
      for (int u = cur + 1; u <= ent[i].u; ++u) P.off[u] = static_cast<int>(P.bin.size());
      cur = ent[i].u;
    }
    if (!P.bin.empty() && P.off[cur] < static_cast<int>(P.bin.size()) &&
        P.bin.back() == ent[i].b) {
      P.wgt.back() += ent[i].w;
    } else {
      P.bin.push_back(ent[i].b);
      P.wgt.push_back(ent[i].w);
    }
  }
  for (int u = cur + 1; u <= P.nUnits; ++u) P.off[u] = static_cast<int>(P.bin.size());
  return P;
}

// ---------------------------------------------------------------------------
// Weighted two-sample EDF statistics
// ---------------------------------------------------------------------------

struct Stat {
  double AD = 0.0;   // Anderson-Darling: 1/(H(1-H)) tail weighting
  double CvM = 0.0;  // Cramer-von Mises: flat weighting, tame in the tails
  double KS = 0.0;   // Kolmogorov-Smirnov: supremum, blind to the tails
};

/// Two-sample EDF statistics from binned weighted histograms.
///
///   A^2   = (nA nB / N) * sum (F_A - F_B)^2 / (H (1-H)) dH
///   omega^2 = (nA nB / N) * sum (F_A - F_B)^2 dH
///   D     = sqrt(nA nB / N) * max |F_A - F_B|
///
/// with H the pooled CDF weighted by the EFFECTIVE sizes nA, nB.  Both
/// histograms are renormalised internally, which is what makes the whole thing
/// invariant to the arbitrary production normalisation.
///
/// `clip` drops bins with H(1-H) <= clip.  clip = 0 keeps everything except the
/// mathematically singular last bin (H == 1).  Raise it (1e-4, 1e-3) if the
/// bootstrap shows AD is being driven by a handful of extreme-weight muons.
/// `dcdf` and `contrib`, when non-null, receive the per-bin signed CDF
/// difference F_test - F_ref and the per-bin contribution to A^2.  Those two
/// curves are what turn "the samples differ" into "they differ here, by this
/// much, in this direction".
inline Stat statFromHists(const double* hA, const double* hB, int nb, double nA,
                          double nB, double clip = 0.0, double* dcdf = nullptr,
                          double* contrib = nullptr) {
  Stat s;
  double sA = 0.0, sB = 0.0;
  for (int k = 0; k < nb; ++k) { sA += hA[k]; sB += hB[k]; }
  if (!(sA > 0.0) || !(sB > 0.0) || !(nA > 0.0) || !(nB > 0.0)) return s;

  const double pref = nA * nB / (nA + nB);
  const double rA = nA / (nA + nB), rB = nB / (nA + nB);
  double FA = 0.0, FB = 0.0, Hprev = 0.0, maxd = 0.0;
  for (int k = 0; k < nb; ++k) {
    FA += hA[k] / sA;
    FB += hB[k] / sB;
    const double H = rA * FA + rB * FB;
    const double dH = H - Hprev;
    Hprev = H;
    const double d = FA - FB;
    maxd = std::max(maxd, std::fabs(d));
    if (dcdf) dcdf[k] = FB - FA;  // signed: positive = test is softer here
    if (contrib) contrib[k] = 0.0;
    if (!(dH > 0.0)) continue;
    s.CvM += d * d * dH;
    const double den = H * (1.0 - H);
    if (den > clip && den > 0.0) {
      const double term = d * d / den * dH;
      s.AD += term;
      if (contrib) contrib[k] = pref * term;
    }
  }
  s.AD *= pref;
  s.CvM *= pref;
  s.KS = std::sqrt(pref) * maxd;
  return s;
}

/// Fill hA/hB (size nbins each) for a given per-unit labelling.
inline void fillByLabel(const PooledVar& P, const uint8_t* lab, double* hA, double* hB,
                        double& nAout, double& nBout) {
  std::fill(hA, hA + P.nbins, 0.0);
  std::fill(hB, hB + P.nbins, 0.0);
  double sA = 0, s2A = 0, sB = 0, s2B = 0;
  for (int u = 0; u < P.nUnits; ++u) {
    double* h = lab[u] ? hB : hA;
    for (int e = P.off[u]; e < P.off[u + 1]; ++e) h[P.bin[e]] += P.wgt[e];
    if (lab[u]) { sB += P.uSw[u]; s2B += P.uSw2[u]; }
    else { sA += P.uSw[u]; s2A += P.uSw2[u]; }
  }
  nAout = nEffFrom(sA, s2A);
  nBout = nEffFrom(sB, s2B);
}

struct ShapeResult {
  std::string var;
  Stat obs;
  double pAD = 1.0, pCvM = 1.0, pKS = 1.0;
  double zAD = 0.0, zCvM = 0.0, zKS = 0.0;
  double neffRef = 0.0, neffTest = 0.0;
  int nbins = 0, nUnits = 0, nperm = 0;
  std::vector<double> nullAD, nullCvM, nullKS;

  // Discrepancy profile: where along the variable the difference sits.
  std::vector<double> edges;      // coarse bin edges, nbins+1
  std::vector<double> dCdf;       // F_test - F_ref at each bin's upper edge
  std::vector<double> a2Contrib;  // per-bin contribution to A^2
  std::vector<double> cumA2;      // running sum of the above, normalised to 1
  std::vector<double> lo68, hi68, lo95, hi95;  // pointwise permutation bands
};

/// Permutation test.  Labels are shuffled at the unit level, preserving the
/// number of A-units and B-units, and each pseudo-sample is renormalised
/// exactly like the observed one -- so the null is a pure shape null and knows
/// nothing about the rate.
inline ShapeResult shapeTest(const PooledVar& P, const std::string& varName, int nperm,
                             uint64_t seed, int nthreads = 0, double clip = 0.0) {
  ShapeResult R;
  R.var = varName;
  R.nbins = P.nbins;
  R.nUnits = P.nUnits;
  R.nperm = nperm;

  std::vector<double> hA(P.nbins), hB(P.nbins);
  fillByLabel(P, P.label.data(), hA.data(), hB.data(), R.neffRef, R.neffTest);
  R.edges = P.edges;
  R.dCdf.assign(P.nbins, 0.0);
  R.a2Contrib.assign(P.nbins, 0.0);
  R.obs = statFromHists(hA.data(), hB.data(), P.nbins, R.neffRef, R.neffTest, clip,
                        R.dCdf.data(), R.a2Contrib.data());
  R.cumA2.assign(P.nbins, 0.0);
  {
    double run = 0.0, tot = 0.0;
    for (double x : R.a2Contrib) tot += x;
    for (int k = 0; k < P.nbins; ++k) {
      run += R.a2Contrib[k];
      R.cumA2[k] = (tot > 0.0) ? run / tot : 0.0;
    }
  }

  if (nperm <= 0) return R;
  R.nullAD.resize(nperm);
  R.nullCvM.resize(nperm);
  R.nullKS.resize(nperm);
  // Per-bin CDF differences under relabelling: the envelope of these is the
  // band the observed curve must escape to mean anything.
  std::vector<std::vector<double>> nullD(nperm);

  parallelFor(static_cast<std::size_t>(nperm), nthreads, [&](std::size_t k) {
    std::mt19937_64 rng(seed + 0x9E3779B97F4A7C15ULL * (k + 1));
    std::vector<uint8_t> lab(P.label);
    std::shuffle(lab.begin(), lab.end(), rng);
    std::vector<double> a(P.nbins), b(P.nbins);
    double nA = 0, nB = 0;
    fillByLabel(P, lab.data(), a.data(), b.data(), nA, nB);
    nullD[k].assign(P.nbins, 0.0);
    const Stat s = statFromHists(a.data(), b.data(), P.nbins, nA, nB, clip,
                                 nullD[k].data(), nullptr);
    R.nullAD[k] = s.AD;
    R.nullCvM[k] = s.CvM;
    R.nullKS[k] = s.KS;
  });

  // Pointwise quantiles.  Pointwise, not simultaneous: with n bins, about 5%
  // of them stray outside the 95% band under the null, so read coherent runs
  // of bins rather than isolated excursions.
  R.lo68.assign(P.nbins, 0.0);
  R.hi68.assign(P.nbins, 0.0);
  R.lo95.assign(P.nbins, 0.0);
  R.hi95.assign(P.nbins, 0.0);
  {
    std::vector<double> col(nperm);
    for (int k = 0; k < P.nbins; ++k) {
      for (int t = 0; t < nperm; ++t) col[t] = nullD[t][k];
      std::sort(col.begin(), col.end());
      auto q = [&](double f) {
        const int i = std::max(0, std::min(nperm - 1,
                                           static_cast<int>(f * (nperm - 1) + 0.5)));
        return col[i];
      };
      R.lo95[k] = q(0.025);
      R.lo68[k] = q(0.16);
      R.hi68[k] = q(0.84);
      R.hi95[k] = q(0.975);
    }
  }

  auto finish = [](const std::vector<double>& null, double obs, double& p, double& z) {
    int ge = 0;
    double m = 0.0, m2 = 0.0;
    for (double x : null) {
      if (x >= obs) ++ge;
      m += x;
      m2 += x * x;
    }
    const double n = static_cast<double>(null.size());
    m /= n;
    const double sd = std::sqrt(std::max(0.0, m2 / n - m * m));
    p = (ge + 1.0) / (n + 1.0);           // one-sided: large statistic = different
    z = (sd > 0.0) ? (obs - m) / sd : 0.0;
  };
  finish(R.nullAD, R.obs.AD, R.pAD, R.zAD);
  finish(R.nullCvM, R.obs.CvM, R.pCvM, R.zCvM);
  finish(R.nullKS, R.obs.KS, R.pKS, R.zKS);
  return R;
}

// ---------------------------------------------------------------------------
// Scale-free fractions and ratios, with a block bootstrap
// ---------------------------------------------------------------------------

/// A phase-space region: var in [lo,hi), optionally restricted to one charge.
struct Region {
  std::string name;
  int var = 0;
  double lo = -1e300, hi = 1e300;
  int charge = 0;  // 0 = both, +1 = mu+, -1 = mu-
};

/// A reported number: the ratio of two accumulator columns.  Column 0 is
/// always the sample total, so a plain fraction is {name, k, 0}.
struct Quantity {
  std::string name;
  int num = 0;
  int den = 0;
};

/// Per-unit weighted sums, one row per bootstrap unit, one column per region
/// (+ column 0 = total).  Small: nUnits is ~ n_files or ~200 blocks.
struct UnitSums {
  int nUnits = 0, ncol = 0;
  std::vector<double> m;  // [nUnits * ncol]
  double at(int u, int c) const { return m[std::size_t(u) * ncol + c]; }
  double& at(int u, int c) { return m[std::size_t(u) * ncol + c]; }
};

/// Assign each row to a bootstrap unit: files when there are enough of them,
/// otherwise random blocks of whole events.
inline std::vector<int> bootstrapUnits(const Sample& S, int minFiles, int nBlocks,
                                       uint64_t seed, int& nUnitsOut, bool& usedFiles) {
  if (S.nFiles >= minFiles) {
    usedFiles = true;
    nUnitsOut = S.nFiles;
    return S.file;
  }
  usedFiles = false;
  const int nb = std::max(2, std::min(nBlocks, std::max(2, S.nEvents)));
  std::vector<int> ev2blk(S.nEvents);
  std::mt19937_64 rng(seed ^ 0xD1B54A32D192ED03ULL);
  for (int e = 0; e < S.nEvents; ++e) ev2blk[e] = static_cast<int>(rng() % nb);
  std::vector<int> u(S.rows());
  for (std::size_t i = 0; i < S.rows(); ++i) u[i] = ev2blk[S.event[i]];
  nUnitsOut = nb;
  return u;
}

inline UnitSums accumulate(const Sample& S, const std::vector<Region>& regs,
                           const std::vector<int>& unit, int nUnits) {
  UnitSums U;
  U.nUnits = nUnits;
  U.ncol = static_cast<int>(regs.size()) + 1;
  U.m.assign(std::size_t(nUnits) * U.ncol, 0.0);
  for (std::size_t i = 0; i < S.rows(); ++i) {
    const int u = unit[i];
    const double w = S.w[i];
    U.at(u, 0) += w;
    for (std::size_t r = 0; r < regs.size(); ++r) {
      const Region& R = regs[r];
      if (R.charge != 0 && S.q[i] != R.charge) continue;
      const double x = S.var[R.var][i];
      if (x >= R.lo && x < R.hi) U.at(u, int(r) + 1) += w;
    }
  }
  return U;
}

inline std::vector<double> evalQuantities(const std::vector<double>& col,
                                          const std::vector<Quantity>& qs) {
  std::vector<double> out(qs.size(), std::nan(""));
  for (std::size_t i = 0; i < qs.size(); ++i) {
    const double d = col[qs[i].den];
    if (d > 0.0) out[i] = col[qs[i].num] / d;
  }
  return out;
}

struct QuantityResult {
  std::string name;
  double a = 0, ea = 0, b = 0, eb = 0;
  double diff = 0, ediff = 0, pull = 0;
  double relDiff = 0;  // (b-a)/a, the scale-free headline number
  // A region can be empty, or a ratio can have an empty denominator.  Those
  // cases are not differences, they are absences of measurement, and must not
  // be printed as "-100%".
  bool okA = false, okB = false, comparable = false;
};

/// Bootstrap over units, independently for A and B, then compare.
inline std::vector<QuantityResult> compareQuantities(
    const Sample& A, const Sample& B, const std::vector<Region>& regs,
    const std::vector<Quantity>& qs, int nboot, uint64_t seed, int minFiles = 5,
    int nBlocks = 200, bool* usedFilesA = nullptr, bool* usedFilesB = nullptr) {
  int nUA = 0, nUB = 0;
  bool fA = false, fB = false;
  const std::vector<int> uA = bootstrapUnits(A, minFiles, nBlocks, seed, nUA, fA);
  const std::vector<int> uB = bootstrapUnits(B, minFiles, nBlocks, seed + 7, nUB, fB);
  if (usedFilesA) *usedFilesA = fA;
  if (usedFilesB) *usedFilesB = fB;

  const UnitSums UA = accumulate(A, regs, uA, nUA);
  const UnitSums UB = accumulate(B, regs, uB, nUB);

  auto total = [](const UnitSums& U) {
    std::vector<double> c(U.ncol, 0.0);
    for (int u = 0; u < U.nUnits; ++u)
      for (int j = 0; j < U.ncol; ++j) c[j] += U.at(u, j);
    return c;
  };
  const std::vector<double> obsA = evalQuantities(total(UA), qs);
  const std::vector<double> obsB = evalQuantities(total(UB), qs);

  auto boot = [&](const UnitSums& U, uint64_t sd) {
    std::vector<std::vector<double>> rep(qs.size());
    std::mt19937_64 rng(sd);
    std::uniform_int_distribution<int> pick(0, U.nUnits - 1);
    std::vector<double> c(U.ncol);
    for (int b = 0; b < nboot; ++b) {
      std::fill(c.begin(), c.end(), 0.0);
      for (int k = 0; k < U.nUnits; ++k) {
        const int u = pick(rng);
        for (int j = 0; j < U.ncol; ++j) c[j] += U.at(u, j);
      }
      const std::vector<double> v = evalQuantities(c, qs);
      for (std::size_t i = 0; i < qs.size(); ++i)
        if (std::isfinite(v[i])) rep[i].push_back(v[i]);
    }
    std::vector<double> sd_out(qs.size(), 0.0);
    for (std::size_t i = 0; i < qs.size(); ++i) {
      if (rep[i].size() < 2) continue;
      double m = 0.0;
      for (double x : rep[i]) m += x;
      m /= rep[i].size();
      double s2 = 0.0;
      for (double x : rep[i]) s2 += (x - m) * (x - m);
      sd_out[i] = std::sqrt(s2 / (rep[i].size() - 1));
    }
    return sd_out;
  };
  const std::vector<double> eA = boot(UA, seed + 11);
  const std::vector<double> eB = boot(UB, seed + 13);

  std::vector<QuantityResult> out(qs.size());
  for (std::size_t i = 0; i < qs.size(); ++i) {
    QuantityResult& r = out[i];
    r.name = qs[i].name;
    r.a = obsA[i];
    r.b = obsB[i];
    r.ea = eA[i];
    r.eb = eB[i];
    r.diff = r.b - r.a;
    r.ediff = std::hypot(r.ea, r.eb);
    r.pull = (r.ediff > 0.0) ? r.diff / r.ediff : 0.0;
    r.okA = std::isfinite(obsA[i]);
    r.okB = std::isfinite(obsB[i]);
    r.comparable = r.okA && r.okB && r.a > 0.0;
    r.relDiff = r.comparable ? r.diff / r.a : std::nan("");
  }
  return out;
}

// ---------------------------------------------------------------------------
// Shape-ratio band for the ratio panels
// ---------------------------------------------------------------------------

/// Per-bin SHAPE ratio (B/A, each normalised to unit area) with a bootstrap
/// error band.  Bootstrapping over files or event blocks -- rather than using
/// sqrt(sum w^2) -- automatically accounts for correlations between muons in
/// the same event or file, including splitmult copies.  That removes the need
/// for a hand-tuned error inflation factor.
struct RatioBand {
  std::vector<double> ratio, err;   // per bin
  std::vector<double> shapeA, shapeB, errA, errB;
};

inline RatioBand shapeRatioBootstrap(const Sample& A, const Sample& B, int j,
                                     const std::vector<double>& edges, int nboot,
                                     uint64_t seed, int minFiles = 5, int nBlocks = 200) {
  const int nb = static_cast<int>(edges.size()) - 1;
  int nUA = 0, nUB = 0;
  bool fA = false, fB = false;
  const std::vector<int> uA = bootstrapUnits(A, minFiles, nBlocks, seed, nUA, fA);
  const std::vector<int> uB = bootstrapUnits(B, minFiles, nBlocks, seed + 7, nUB, fB);

  auto build = [&](const Sample& S, const std::vector<int>& u, int nU) {
    std::vector<double> m(std::size_t(nU) * nb, 0.0);
    for (std::size_t i = 0; i < S.rows(); ++i)
      m[std::size_t(u[i]) * nb + binOf(edges, S.var[j][i])] += S.w[i];
    return m;
  };
  const std::vector<double> MA = build(A, uA, nUA);
  const std::vector<double> MB = build(B, uB, nUB);

  auto shapeOf = [&](const std::vector<double>& M, int nU, const int* pickIdx) {
    std::vector<double> h(nb, 0.0);
    for (int k = 0; k < nU; ++k) {
      const int u = pickIdx ? pickIdx[k] : k;
      for (int b = 0; b < nb; ++b) h[b] += M[std::size_t(u) * nb + b];
    }
    double s = 0.0;
    for (double x : h) s += x;
    if (s > 0.0) for (double& x : h) x /= s;
    return h;
  };

  RatioBand R;
  R.shapeA = shapeOf(MA, nUA, nullptr);
  R.shapeB = shapeOf(MB, nUB, nullptr);
  R.ratio.assign(nb, std::nan(""));
  for (int b = 0; b < nb; ++b)
    if (R.shapeA[b] > 0.0) R.ratio[b] = R.shapeB[b] / R.shapeA[b];

  std::vector<std::vector<double>> rep(nb);
  std::vector<std::vector<double>> repA(nb), repB(nb);
  std::mt19937_64 rng(seed + 17);
  std::vector<int> ia(nUA), ib(nUB);
  for (int t = 0; t < nboot; ++t) {
    for (int k = 0; k < nUA; ++k) ia[k] = static_cast<int>(rng() % nUA);
    for (int k = 0; k < nUB; ++k) ib[k] = static_cast<int>(rng() % nUB);
    const std::vector<double> ha = shapeOf(MA, nUA, ia.data());
    const std::vector<double> hb = shapeOf(MB, nUB, ib.data());
    for (int b = 0; b < nb; ++b) {
      repA[b].push_back(ha[b]);
      repB[b].push_back(hb[b]);
      if (ha[b] > 0.0) rep[b].push_back(hb[b] / ha[b]);
    }
  }
  auto sdev = [](const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    double m = 0.0;
    for (double x : v) m += x;
    m /= v.size();
    double s2 = 0.0;
    for (double x : v) s2 += (x - m) * (x - m);
    return std::sqrt(s2 / (v.size() - 1));
  };
  R.err.resize(nb);
  R.errA.resize(nb);
  R.errB.resize(nb);
  for (int b = 0; b < nb; ++b) {
    R.err[b] = sdev(rep[b]);
    R.errA[b] = sdev(repA[b]);
    R.errB[b] = sdev(repB[b]);
  }
  return R;
}

// ---------------------------------------------------------------------------
// Per-bin effective statistics -- the ceiling on everything above
// ---------------------------------------------------------------------------

inline std::vector<double> nEffPerBin(const Sample& S, int j,
                                      const std::vector<double>& edges) {
  const int nb = static_cast<int>(edges.size()) - 1;
  std::vector<double> s(nb, 0.0), s2(nb, 0.0), out(nb, 0.0);
  for (std::size_t i = 0; i < S.rows(); ++i) {
    const int b = binOf(edges, S.var[j][i]);
    s[b] += S.w[i];
    s2[b] += S.w[i] * S.w[i];
  }
  for (int b = 0; b < nb; ++b) out[b] = nEffFrom(s[b], s2[b]);
  return out;
}


// ===========================================================================
// STREAMING ACCUMULATOR
// ===========================================================================
//
// The in-memory Sample above keeps every muon.  One split100 file holds 36M
// muons; a 300-file reference sample would need hundreds of GB.  So the real
// analysis never stores muons at all.  It streams them straight into a fixed
// set of fine-binned histograms, one per (unit, charge, variable), and every
// statistic downstream is computed from those.
//
// Memory is then O(nUnits * 2 * nvar * nfine) and completely independent of
// how much data goes in: 200 units x 2 x 8 x 512 x 8 B = 13 MB.
//
// Two approximations come with it, both controllable and both reported:
//   * bin edges snap to the fine grid, so a cut at 100 GeV lands within one
//     fine bin of 100 GeV.  The effective cut is printed.
//   * values outside the recorded range pile into the end bins.  The clamp
//     count is printed; if it is not ~0, widen the range.
// The selftest checks the streamed path against the exact per-muon path.

struct HistSpec {
  int nvar = 0;
  int nfine = 512;
  std::vector<double> lo, hi;  // per variable, in transformed units

  double width(int v) const { return (hi[v] - lo[v]) / nfine; }
  int bin(int v, double x, bool& clamped) const {
    const double f = (x - lo[v]) / (hi[v] - lo[v]);
    clamped = (f < 0.0 || f >= 1.0);
    int b = static_cast<int>(f * nfine);
    return std::max(0, std::min(nfine - 1, b));
  }
  double edge(int v, int b) const { return lo[v] + b * width(v); }
  /// Fine-grid index whose lower edge is nearest to x.
  int snap(int v, double x) const {
    if (x <= lo[v]) return 0;
    if (x >= hi[v]) return nfine;
    return static_cast<int>(std::lround((x - lo[v]) / width(v)));
  }
};

/// Charge index: 0 = mu+, 1 = mu-.  Pass chg = -1 to any accessor for "both".
inline int chgIndex(int q) { return q > 0 ? 0 : 1; }

struct UnitHists {
  std::string label;
  HistSpec spec;
  int nUnits = 0;
  bool unitsAreFiles = false;
  std::vector<double> h;     // [unit][chg][var][bin], size nUnits*2*nvar*nfine
  std::vector<double> hw2;   // [chg][var][bin]  -- sum w^2, for N_eff per bin
  std::vector<double> uSw, uSw2, uSwP, uSwM;  // per unit
  long long rows = 0, clamped = 0;

  std::size_t idx(int u, int c, int v, int b) const {
    return ((static_cast<std::size_t>(u) * 2 + c) * spec.nvar + v) * spec.nfine + b;
  }
  std::size_t idx2(int c, int v, int b) const {
    return (static_cast<std::size_t>(c) * spec.nvar + v) * spec.nfine + b;
  }
  void alloc(const HistSpec& sp, int nU) {
    spec = sp;
    nUnits = nU;
    h.assign(static_cast<std::size_t>(nU) * 2 * sp.nvar * sp.nfine, 0.0);
    hw2.assign(static_cast<std::size_t>(2) * sp.nvar * sp.nfine, 0.0);
    uSw.assign(nU, 0.0);
    uSw2.assign(nU, 0.0);
    uSwP.assign(nU, 0.0);
    uSwM.assign(nU, 0.0);
  }
  /// One particle, addressed by sign class index (0 = PDG < 0).
  void fillSign(int unit, int c, const double* vals, double w) {
    bool cl = false, any = false;
    for (int v = 0; v < spec.nvar; ++v) {
      const int b = spec.bin(v, vals[v], cl);
      any = any || cl;
      h[idx(unit, c, v, b)] += w;
      hw2[idx2(c, v, b)] += w * w;
    }
    uSw[unit] += w;
    uSw2[unit] += w * w;
    (c == 0 ? uSwP : uSwM)[unit] += w;
    ++rows;
    if (any) ++clamped;
  }
  /// Convenience overload used by the tests, which speak in charges.
  void fill(int unit, int q, const double* vals, double w) {
    fillSign(unit, q > 0 ? 0 : 1, vals, w);
  }
  void merge(const UnitHists& o) {
    if (h.empty()) { *this = o; return; }
    for (std::size_t i = 0; i < h.size(); ++i) h[i] += o.h[i];
    for (std::size_t i = 0; i < hw2.size(); ++i) hw2[i] += o.hw2[i];
    for (int u = 0; u < nUnits; ++u) {
      uSw[u] += o.uSw[u];
      uSw2[u] += o.uSw2[u];
      uSwP[u] += o.uSwP[u];
      uSwM[u] += o.uSwM[u];
    }
    rows += o.rows;
    clamped += o.clamped;
  }
  double totalW() const { return sumW(uSw); }
  double nEffTotal() const {
    double s = 0, s2 = 0;
    for (int u = 0; u < nUnits; ++u) { s += uSw[u]; s2 += uSw2[u]; }
    return nEffFrom(s, s2);
  }
  /// Fine-bin spectrum summed over units, for one charge selection.
  std::vector<double> fine(int var, int chg) const {
    std::vector<double> out(spec.nfine, 0.0);
    for (int u = 0; u < nUnits; ++u)
      for (int c = 0; c < 2; ++c) {
        if (chg >= 0 && c != chg) continue;
        for (int b = 0; b < spec.nfine; ++b) out[b] += h[idx(u, c, var, b)];
      }
    return out;
  }
  /// Per-unit fine spectrum, for one charge selection.
  std::vector<double> fineByUnit(int var, int chg) const {
    std::vector<double> m(static_cast<std::size_t>(nUnits) * spec.nfine, 0.0);
    for (int u = 0; u < nUnits; ++u)
      for (int c = 0; c < 2; ++c) {
        if (chg >= 0 && c != chg) continue;
        for (int b = 0; b < spec.nfine; ++b)
          m[static_cast<std::size_t>(u) * spec.nfine + b] += h[idx(u, c, var, b)];
      }
    return m;
  }
};

/// Merge every `factor` adjacent fine bins.  This is EXACT: summing adjacent
/// bins gives precisely what accumulating on the coarser grid would have
/// given, because the coarse edges are a subset of the fine ones.  So a whole
/// binning-stability scan costs one read, not several.
inline UnitHists coarsenFine(const UnitHists& U, int factor) {
  if (factor < 2 || U.spec.nfine % factor != 0) return U;
  HistSpec sp = U.spec;
  sp.nfine = U.spec.nfine / factor;
  UnitHists O;
  O.alloc(sp, U.nUnits);
  O.label = U.label;
  O.unitsAreFiles = U.unitsAreFiles;
  O.rows = U.rows;
  O.clamped = U.clamped;
  O.uSw = U.uSw;
  O.uSw2 = U.uSw2;
  O.uSwP = U.uSwP;
  O.uSwM = U.uSwM;
  for (int u = 0; u < U.nUnits; ++u)
    for (int c = 0; c < 2; ++c)
      for (int v = 0; v < U.spec.nvar; ++v)
        for (int b = 0; b < U.spec.nfine; ++b)
          O.h[O.idx(u, c, v, b / factor)] += U.h[U.idx(u, c, v, b)];
  for (int c = 0; c < 2; ++c)
    for (int v = 0; v < U.spec.nvar; ++v)
      for (int b = 0; b < U.spec.nfine; ++b)
        O.hw2[O.idx2(c, v, b / factor)] += U.hw2[U.idx2(c, v, b)];
  return O;
}

/// Split the units into two halves.  This is the closure test in the streamed
/// world: no re-reading, and the halves are exactly disjoint sets of units.
inline std::pair<UnitHists, UnitHists> splitUnits(const UnitHists& U, uint64_t seed) {
  std::vector<int> order(U.nUnits);
  std::iota(order.begin(), order.end(), 0);
  std::mt19937_64 rng(seed);
  std::shuffle(order.begin(), order.end(), rng);
  const int half = U.nUnits / 2;

  UnitHists A, B;
  A.label = U.label + " (half 1)";
  B.label = U.label + " (half 2)";
  A.alloc(U.spec, half);
  B.alloc(U.spec, U.nUnits - half);
  A.unitsAreFiles = B.unitsAreFiles = U.unitsAreFiles;
  for (int k = 0; k < U.nUnits; ++k) {
    UnitHists& T = (k < half) ? A : B;
    const int t = (k < half) ? k : k - half;
    const int u = order[k];
    for (int c = 0; c < 2; ++c)
      for (int v = 0; v < U.spec.nvar; ++v)
        for (int b = 0; b < U.spec.nfine; ++b) {
          const double w = U.h[U.idx(u, c, v, b)];
          T.h[T.idx(t, c, v, b)] += w;
          T.hw2[T.idx2(c, v, b)] += 0.0;  // w^2 is not separable per unit
        }
    T.uSw[t] = U.uSw[u];
    T.uSw2[t] = U.uSw2[u];
    T.uSwP[t] = U.uSwP[u];
    T.uSwM[t] = U.uSwM[u];
    T.rows += 0;
  }
  // hw2 is only used for the N_eff diagnostic; approximate it by sharing the
  // parent's sum in proportion to the weight in each bin.
  for (int c = 0; c < 2; ++c)
    for (int v = 0; v < U.spec.nvar; ++v)
      for (int b = 0; b < U.spec.nfine; ++b) {
        double tot = 0.0;
        for (int u = 0; u < U.nUnits; ++u) tot += U.h[U.idx(u, c, v, b)];
        if (!(tot > 0.0)) continue;
        const double parent = U.hw2[U.idx2(c, v, b)];
        double sa = 0.0;
        for (int t = 0; t < A.nUnits; ++t) sa += A.h[A.idx(t, c, v, b)];
        A.hw2[A.idx2(c, v, b)] = parent * (sa / tot);
        B.hw2[B.idx2(c, v, b)] = parent * (1.0 - sa / tot);
      }
  A.rows = U.rows / 2;
  B.rows = U.rows - A.rows;
  return {A, B};
}

/// Coarse bin boundaries, as indices into the fine grid, placed so that each
/// coarse bin holds an equal share of the pooled weight.
inline std::vector<int> equalWeightFineEdges(const UnitHists& R, const UnitHists& T,
                                             int var, int chg, int nbins) {
  const int nf = R.spec.nfine;
  const std::vector<double> fr = R.fine(var, chg), ft = T.fine(var, chg);
  double sr = 0, st = 0;
  for (int b = 0; b < nf; ++b) { sr += fr[b]; st += ft[b]; }
  std::vector<double> pooled(nf, 0.0);
  for (int b = 0; b < nf; ++b)
    pooled[b] = (sr > 0 ? fr[b] / sr : 0.0) + (st > 0 ? ft[b] / st : 0.0);

  double tot = 0.0;
  for (double x : pooled) tot += x;
  std::vector<int> edges{0};
  if (!(tot > 0.0)) { edges.push_back(nf); return edges; }
  double cum = 0.0;
  int k = 1;
  for (int b = 0; b < nf && k < nbins; ++b) {
    cum += pooled[b];
    while (k < nbins && cum >= k * tot / nbins) {
      if (b + 1 > edges.back()) edges.push_back(b + 1);
      ++k;
    }
  }
  if (edges.back() != nf) edges.push_back(nf);
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  return edges;
}

/// Physical coordinates of coarse edges.
inline std::vector<double> edgeValues(const HistSpec& sp, int var,
                                      const std::vector<int>& fineEdges) {
  std::vector<double> e;
  e.reserve(fineEdges.size());
  for (int b : fineEdges) e.push_back(sp.edge(var, b));
  return e;
}

/// CSR over units, coarse-binned -- the input the permutation test wants.
inline PooledVar buildPooledFromHists(const UnitHists& R, const UnitHists& T, int var,
                                      int chg, const std::vector<int>& fineEdges) {
  PooledVar P;
  P.nbins = static_cast<int>(fineEdges.size()) - 1;
  P.nUnits = R.nUnits + T.nUnits;
  P.edges = edgeValues(R.spec, var, fineEdges);
  P.label.assign(P.nUnits, 0);
  for (int u = R.nUnits; u < P.nUnits; ++u) P.label[u] = 1;
  P.off.assign(P.nUnits + 1, 0);
  P.uSw.assign(P.nUnits, 0.0);
  P.uSw2.assign(P.nUnits, 0.0);

  // Each sample is normalised to unit total weight before pooling, which is
  // what makes every statistic blind to the production rate.
  const double sR = R.totalW(), sT = T.totalW();
  std::vector<int> f2c(R.spec.nfine, 0);
  for (int cb = 0; cb < P.nbins; ++cb)
    for (int b = fineEdges[cb]; b < fineEdges[cb + 1]; ++b) f2c[b] = cb;

  std::vector<double> tmp(P.nbins);
  for (int side = 0; side < 2; ++side) {
    const UnitHists& S = side ? T : R;
    const double norm = (side ? sT : sR) > 0 ? 1.0 / (side ? sT : sR) : 0.0;
    for (int u = 0; u < S.nUnits; ++u) {
      const int gu = side ? R.nUnits + u : u;
      std::fill(tmp.begin(), tmp.end(), 0.0);
      double sw = 0.0;
      for (int c = 0; c < 2; ++c) {
        if (chg >= 0 && c != chg) continue;
        for (int b = 0; b < S.spec.nfine; ++b) {
          const double w = S.h[S.idx(u, c, var, b)] * norm;
          if (w != 0.0) { tmp[f2c[b]] += w; sw += w; }
        }
      }
      P.off[gu] = static_cast<int>(P.bin.size());
      for (int cb = 0; cb < P.nbins; ++cb)
        if (tmp[cb] != 0.0) {
          P.bin.push_back(static_cast<uint16_t>(cb));
          P.wgt.push_back(tmp[cb]);
        }
      P.uSw[gu] = sw;
      // N_eff of a unit is a property of its weights, and rescaling by `norm`
      // scales sum w^2 by norm^2.  Charge selection is ignored here: it would
      // need per-charge sum w^2, and the effect on N_eff is second order.
      P.uSw2[gu] = S.uSw2[u] * norm * norm;
    }
  }
  P.off[P.nUnits] = static_cast<int>(P.bin.size());
  for (int u = P.nUnits - 1; u >= 0; --u)
    if (P.off[u] == 0 && u > 0 && P.off[u - 1] > 0) P.off[u] = P.off[u + 1];
  return P;
}

/// Weighted sums per unit for the region/fraction machinery.
inline UnitSums unitSumsFromHists(const UnitHists& S, const std::vector<Region>& regs) {
  UnitSums U;
  U.nUnits = S.nUnits;
  U.ncol = static_cast<int>(regs.size()) + 1;
  U.m.assign(static_cast<std::size_t>(U.nUnits) * U.ncol, 0.0);
  for (int u = 0; u < S.nUnits; ++u) {
    U.at(u, 0) = S.uSw[u];
    for (std::size_t r = 0; r < regs.size(); ++r) {
      const Region& R = regs[r];
      const int b0 = S.spec.snap(R.var, R.lo);
      const int b1 = S.spec.snap(R.var, R.hi);
      double acc = 0.0;
      for (int c = 0; c < 2; ++c) {
        if (R.charge != 0 && c != chgIndex(R.charge)) continue;
        for (int b = b0; b < b1 && b < S.spec.nfine; ++b)
          acc += S.h[S.idx(u, c, R.var, b)];
      }
      U.at(u, static_cast<int>(r) + 1) = acc;
    }
  }
  return U;
}

/// The cut actually applied, after snapping to the fine grid.
inline std::pair<double, double> effectiveCut(const HistSpec& sp, const Region& R) {
  return {sp.edge(R.var, sp.snap(R.var, R.lo)), sp.edge(R.var, sp.snap(R.var, R.hi))};
}

/// Bootstrap comparison driven by two pre-built UnitSums.
inline std::vector<QuantityResult> compareQuantitySums(const UnitSums& UA,
                                                       const UnitSums& UB,
                                                       const std::vector<Quantity>& qs,
                                                       int nboot, uint64_t seed) {
  auto total = [](const UnitSums& U) {
    std::vector<double> c(U.ncol, 0.0);
    for (int u = 0; u < U.nUnits; ++u)
      for (int j = 0; j < U.ncol; ++j) c[j] += U.at(u, j);
    return c;
  };
  const std::vector<double> obsA = evalQuantities(total(UA), qs);
  const std::vector<double> obsB = evalQuantities(total(UB), qs);

  auto boot = [&](const UnitSums& U, uint64_t sd) {
    std::vector<std::vector<double>> rep(qs.size());
    std::mt19937_64 rng(sd);
    std::vector<double> c(U.ncol);
    for (int b = 0; b < nboot; ++b) {
      std::fill(c.begin(), c.end(), 0.0);
      for (int k = 0; k < U.nUnits; ++k) {
        const int u = static_cast<int>(rng() % U.nUnits);
        for (int j = 0; j < U.ncol; ++j) c[j] += U.at(u, j);
      }
      const std::vector<double> v = evalQuantities(c, qs);
      for (std::size_t i = 0; i < qs.size(); ++i)
        if (std::isfinite(v[i])) rep[i].push_back(v[i]);
    }
    std::vector<double> sd_out(qs.size(), 0.0);
    for (std::size_t i = 0; i < qs.size(); ++i) {
      if (rep[i].size() < 2) continue;
      double m = 0.0;
      for (double x : rep[i]) m += x;
      m /= rep[i].size();
      double s2 = 0.0;
      for (double x : rep[i]) s2 += (x - m) * (x - m);
      sd_out[i] = std::sqrt(s2 / (rep[i].size() - 1));
    }
    return sd_out;
  };
  const std::vector<double> eA = boot(UA, seed + 11), eB = boot(UB, seed + 13);

  std::vector<QuantityResult> out(qs.size());
  for (std::size_t i = 0; i < qs.size(); ++i) {
    QuantityResult& r = out[i];
    r.name = qs[i].name;
    r.a = obsA[i];
    r.b = obsB[i];
    r.ea = eA[i];
    r.eb = eB[i];
    r.diff = r.b - r.a;
    r.ediff = std::hypot(r.ea, r.eb);
    r.pull = (r.ediff > 0.0) ? r.diff / r.ediff : 0.0;
    r.okA = std::isfinite(obsA[i]) && (r.a == 0.0 || r.ea < 0.5 * std::fabs(r.a));
    r.okB = std::isfinite(obsB[i]) && (r.b == 0.0 || r.eb < 0.5 * std::fabs(r.b));
    // A fraction known only to +-100% cannot support a comparison; say so
    // rather than printing a spurious "-100%".
    r.comparable = r.okA && r.okB && r.a > 0.0;
    r.relDiff = r.comparable ? r.diff / r.a : std::nan("");
  }
  return out;
}

/// Shape ratio (test/reference, each normalised to unit area) with a bootstrap
/// band, from per-unit coarse histograms.
inline RatioBand shapeRatioFromHists(const UnitHists& R, const UnitHists& T, int var,
                                     int chg, const std::vector<int>& fineEdges,
                                     int nboot, uint64_t seed) {
  const int nb = static_cast<int>(fineEdges.size()) - 1;
  auto coarse = [&](const UnitHists& S) {
    std::vector<double> m(static_cast<std::size_t>(S.nUnits) * nb, 0.0);
    for (int u = 0; u < S.nUnits; ++u)
      for (int cb = 0; cb < nb; ++cb) {
        double acc = 0.0;
        for (int c = 0; c < 2; ++c) {
          if (chg >= 0 && c != chg) continue;
          for (int b = fineEdges[cb]; b < fineEdges[cb + 1]; ++b)
            acc += S.h[S.idx(u, c, var, b)];
        }
        m[static_cast<std::size_t>(u) * nb + cb] = acc;
      }
    return m;
  };
  const std::vector<double> MR = coarse(R), MT = coarse(T);

  auto shapeOf = [&](const std::vector<double>& M, int nU, const int* pick) {
    std::vector<double> hh(nb, 0.0);
    for (int k = 0; k < nU; ++k) {
      const int u = pick ? pick[k] : k;
      for (int b = 0; b < nb; ++b) hh[b] += M[static_cast<std::size_t>(u) * nb + b];
    }
    double s = 0.0;
    for (double x : hh) s += x;
    if (s > 0.0) for (double& x : hh) x /= s;
    return hh;
  };

  RatioBand B;
  B.shapeA = shapeOf(MR, R.nUnits, nullptr);
  B.shapeB = shapeOf(MT, T.nUnits, nullptr);
  B.ratio.assign(nb, std::nan(""));
  for (int b = 0; b < nb; ++b)
    if (B.shapeA[b] > 0.0) B.ratio[b] = B.shapeB[b] / B.shapeA[b];

  std::vector<std::vector<double>> rep(nb), repA(nb), repB(nb);
  std::mt19937_64 rng(seed + 17);
  std::vector<int> ia(R.nUnits), ib(T.nUnits);
  for (int t = 0; t < nboot; ++t) {
    for (int k = 0; k < R.nUnits; ++k) ia[k] = static_cast<int>(rng() % R.nUnits);
    for (int k = 0; k < T.nUnits; ++k) ib[k] = static_cast<int>(rng() % T.nUnits);
    const std::vector<double> ha = shapeOf(MR, R.nUnits, ia.data());
    const std::vector<double> hb = shapeOf(MT, T.nUnits, ib.data());
    for (int b = 0; b < nb; ++b) {
      repA[b].push_back(ha[b]);
      repB[b].push_back(hb[b]);
      if (ha[b] > 0.0) rep[b].push_back(hb[b] / ha[b]);
    }
  }
  auto sdev = [](const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    double m = 0.0;
    for (double x : v) m += x;
    m /= v.size();
    double s2 = 0.0;
    for (double x : v) s2 += (x - m) * (x - m);
    return std::sqrt(s2 / (v.size() - 1));
  };
  B.err.resize(nb);
  B.errA.resize(nb);
  B.errB.resize(nb);
  for (int b = 0; b < nb; ++b) {
    B.err[b] = sdev(rep[b]);
    B.errA[b] = sdev(repA[b]);
    B.errB[b] = sdev(repB[b]);
  }
  return B;
}

inline std::vector<double> nEffPerBinFromHists(const UnitHists& S, int var, int chg,
                                               const std::vector<int>& fineEdges) {
  const int nb = static_cast<int>(fineEdges.size()) - 1;
  std::vector<double> out(nb, 0.0);
  for (int cb = 0; cb < nb; ++cb) {
    double s = 0.0, s2 = 0.0;
    for (int c = 0; c < 2; ++c) {
      if (chg >= 0 && c != chg) continue;
      for (int b = fineEdges[cb]; b < fineEdges[cb + 1]; ++b) {
        for (int u = 0; u < S.nUnits; ++u) s += S.h[S.idx(u, c, var, b)];
        s2 += S.hw2[S.idx2(c, var, b)];
      }
    }
    out[cb] = nEffFrom(s, s2);
  }
  return out;
}

}  // namespace fluxval
