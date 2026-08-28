// core_selftest.cxx -- ROOT-free checks of the statistics in Core.h.
//
// Builds with a plain compiler, no ROOT needed:
//     g++ -O2 -std=c++17 -Iinclude test/core_selftest.cxx -o selftest -pthread
//
// What it asserts, and why each one matters here:
//
//  1. SCALE INVARIANCE.  The two productions have deliberately different muon
//     rates.  Multiplying one sample's weights by 1000 must change nothing.
//     If this fails, every p-value below is a measurement of the generator
//     normalisation and not of the physics.
//  2. NULL CLOSURE.  Two independent draws from the same law must give
//     approximately uniform p-values.  This is the analogue of the A-vs-A test
//     in the existing Fisher code.
//  3. TAIL SENSITIVITY.  Two laws with identical mean AND identical variance
//     but different tails.  A mean-difference statistic (Fisher) sees nothing;
//     AD must see it.  This is the whole reason for the exercise.
//  4. FRACTIONS.  Injected fractions and charge ratios are recovered, and the
//     bootstrap error is of the right size.
//  5. PERMUTATION UNIT.  Correlated copies inside one event do not narrow the
//     null when the event is the permutation unit.
 
#include "muonflux/Core.h"
 
#include <cstdio>
#include <cstdlib>
#include <iostream>
 
using namespace mfc;
 
static int g_fail = 0;
 
static void check(bool ok, const std::string& what) {
  std::printf("  [%s] %s\n", ok ? " ok " : "FAIL", what.c_str());
  if (!ok) ++g_fail;
}
 
// ---------------------------------------------------------------------------
// Synthetic sample builders
// ---------------------------------------------------------------------------
 
enum class Law { Normal, MatchedTail, Shifted };
 
// Two-variable sample: var[0] is the variable under test, var[1] is a
// spectator used to check the region machinery.
static Sample makeSample(const std::string& lab, uint64_t seed, int nEvents, Law law,
                         double wScale = 1.0, double muPlusFrac = 0.5,
                         int copiesPerMuon = 1) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> gauss(0.0, 1.0);
  std::lognormal_distribution<double> lw(0.0, 0.8);
  std::poisson_distribution<int> mult(2.0);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
 
  Sample S;
  S.label = lab;
  S.nvar = 2;
  S.var.assign(2, {});
  for (int e = 0; e < nEvents; ++e) {
    const int n = mult(rng) + 1;
    for (int i = 0; i < n; ++i) {
      double x;
      switch (law) {
        case Law::Normal:
          x = gauss(rng);
          break;
        case Law::MatchedTail:
          // 90% N(0, 0.6455), 10% N(0, 2.5):  mean 0, variance 1, fat tails.
          x = (uni(rng) < 0.10) ? 2.5 * gauss(rng) : 0.6454972 * gauss(rng);
          break;
        default:
          x = gauss(rng) + 0.05;
          break;
      }
      const double w = wScale * lw(rng);
      const int q = (uni(rng) < muPlusFrac) ? +1 : -1;
      const double spec = gauss(rng);
      for (int c = 0; c < copiesPerMuon; ++c) {
        S.var[0].push_back(x);
        S.var[1].push_back(spec);
        S.w.push_back(w / copiesPerMuon);
        S.q.push_back(q);
        S.event.push_back(e);
        S.file.push_back(e % 8);
      }
    }
  }
  S.compactUnits();
  return S;
}
 
// Fisher-style standardised mean difference: what the current code measures.
static double meanSeparation(const Sample& A, const Sample& B, int j) {
  auto mv = [&](const Sample& S) {
    double sw = 0, m = 0, m2 = 0;
    for (std::size_t i = 0; i < S.rows(); ++i) {
      sw += S.w[i];
      m += S.w[i] * S.var[j][i];
      m2 += S.w[i] * S.var[j][i] * S.var[j][i];
    }
    m /= sw;
    return std::pair<double, double>(m, m2 / sw - m * m);
  };
  const auto a = mv(A), b = mv(B);
  const double sp = std::sqrt(0.5 * (a.second + b.second));
  return (b.first - a.first) / sp;
}
 
static ShapeResult run(const Sample& A, const Sample& B, int nperm, uint64_t seed,
                       int nbins = 128) {
  Sample a = A, b = B;
  a.normaliseWeights();
  b.normaliseWeights();
  const PooledVar P = buildPooled(a, b, 0, nbins);
  return shapeTest(P, "x", nperm, seed, 0, 0.0);
}
 
// ---------------------------------------------------------------------------
 
int main() {
  std::printf("\n=== 1. scale invariance ===\n");
  {
    const Sample A = makeSample("A", 1, 4000, Law::Normal, 1.0);
    const Sample B1 = makeSample("B", 2, 4000, Law::MatchedTail, 1.0);
    const Sample B2 = makeSample("B", 2, 4000, Law::MatchedTail, 1000.0);
    const ShapeResult r1 = run(A, B1, 60, 12345);
    const ShapeResult r2 = run(A, B2, 60, 12345);
    check(std::fabs(r1.obs.AD - r2.obs.AD) < 1e-9 * std::max(1.0, r1.obs.AD),
          "AD unchanged when sample B weights are scaled by 1000");
    check(std::fabs(r1.obs.CvM - r2.obs.CvM) < 1e-9 * std::max(1.0, r1.obs.CvM),
          "CvM unchanged under the same rescaling");
    check(r1.pAD == r2.pAD, "permutation p-value unchanged under rescaling");
    std::printf("        AD = %.6f (x1)   %.6f (x1000)\n", r1.obs.AD, r2.obs.AD);
  }
 
  std::printf("\n=== 2. null closure: same law, different seeds ===\n");
  {
    int nSmall = 0;
    const int nTrial = 12;
    double pmin = 1.0;
    for (int t = 0; t < nTrial; ++t) {
      const Sample A = makeSample("A", 100 + 2 * t, 3000, Law::Normal);
      const Sample B = makeSample("B", 101 + 2 * t, 3000, Law::Normal);
      const ShapeResult r = run(A, B, 100, 777 + t);
      pmin = std::min(pmin, r.pAD);
      if (r.pAD < 0.05) ++nSmall;
    }
    std::printf("        %d/%d trials with p_AD < 0.05 (expect ~1), min p = %.3f\n",
                nSmall, nTrial, pmin);
    check(nSmall <= 3, "null closure: p-values are not systematically small");
  }
 
  std::printf("\n=== 3. tail sensitivity at matched mean AND variance ===\n");
  {
    const Sample A = makeSample("A", 11, 6000, Law::Normal);
    const Sample B = makeSample("B", 12, 6000, Law::MatchedTail);
    const double d = meanSeparation(A, B, 0);
    const ShapeResult r = run(A, B, 200, 999);
    std::printf("        Fisher-style mean separation d = %+.4f  (blind)\n", d);
    std::printf("        AD  = %9.3f   p = %.4g   z = %+.1f\n", r.obs.AD, r.pAD, r.zAD);
    std::printf("        CvM = %9.3f   p = %.4g   z = %+.1f\n", r.obs.CvM, r.pCvM, r.zCvM);
    std::printf("        KS  = %9.3f   p = %.4g   z = %+.1f\n", r.obs.KS, r.pKS, r.zKS);
    check(std::fabs(d) < 0.06, "mean-difference statistic is blind to this difference");
    check(r.pAD <= 1.0 / 201.0 + 1e-12, "AD rejects at the permutation floor");
    check(r.zAD > r.zKS, "AD is more sensitive than KS on a tail difference");
  }
 
  std::printf("\n=== 4. fractions and charge ratio ===\n");
  {
    const Sample A = makeSample("A", 21, 5000, Law::Normal, 1.0, 0.50);
    const Sample B = makeSample("B", 22, 5000, Law::Normal, 500.0, 0.60);
    std::vector<Region> regs = {
        {"x > 2", 0, 2.0, 1e300, 0},
        {"mu+", 0, -1e300, 1e300, +1},
        {"mu-", 0, -1e300, 1e300, -1},
    };
    std::vector<Quantity> qs = {
        {"f(x > 2)", 1, 0},
        {"mu+ / mu-", 2, 3},
    };
    const auto res = compareQuantities(A, B, regs, qs, 300, 4242, 5, 200);
    for (const auto& r : res)
      std::printf("        %-12s A = %.4f +- %.4f   B = %.4f +- %.4f   pull %+.2f\n",
                  r.name.c_str(), r.a, r.ea, r.b, r.eb, r.pull);
    check(std::fabs(res[0].pull) < 3.0, "f(x>2) compatible: same law, 500x rate");
    check(res[1].b > res[1].a && std::fabs(res[1].pull) > 3.0,
          "injected charge-ratio change 1.00 -> 1.50 is detected");
    check(res[0].ea > 0.0 && res[1].ea > 0.0, "bootstrap errors are non-zero");
  }
 
  std::printf("\n=== 5. correlated copies inside one event ===\n");
  {
    // splitmult-like: each muon written 4 times with w/4, all inside one event.
    const Sample A = makeSample("A", 31, 2500, Law::Normal, 1.0, 0.5, 4);
    const Sample B = makeSample("B", 32, 2500, Law::Normal, 1.0, 0.5, 4);
    const ShapeResult r = run(A, B, 150, 31337);
    std::printf("        rows = %zu, events = %d, p_AD = %.3f\n", A.rows(), A.nEvents,
                r.pAD);
    check(r.pAD > 0.01,
          "event-level permutation is not fooled by within-event duplicates");
  }
 
  std::printf("\n=== 6. shape-ratio bootstrap band ===\n");
  {
    const Sample A = makeSample("A", 41, 4000, Law::Normal, 1.0);
    const Sample B = makeSample("B", 42, 4000, Law::Normal, 250.0);
    std::vector<double> vAll = A.var[0];
    vAll.insert(vAll.end(), B.var[0].begin(), B.var[0].end());
    std::vector<double> wAll = A.w;
    wAll.insert(wAll.end(), B.w.begin(), B.w.end());
    const std::vector<double> edges = equalWeightEdges(vAll, wAll, 20);
    const RatioBand band = shapeRatioBootstrap(A, B, 0, edges, 200, 5150);
    int nb = 0, npull = 0;
    for (std::size_t b = 0; b < band.ratio.size(); ++b) {
      if (!(band.err[b] > 0.0) || !std::isfinite(band.ratio[b])) continue;
      ++nb;
      if (std::fabs(band.ratio[b] - 1.0) / band.err[b] > 3.0) ++npull;
    }
    std::printf("        %d usable bins, %d beyond 3 sigma (expect ~0)\n", nb, npull);
    check(nb > 10, "band covers the range");
    check(npull <= 1, "shape ratio is compatible with 1 despite a 250x rate ratio");
  }
 
  std::printf("\n=== 7. streamed accumulator vs exact per-muon path ===\n");
  {
    // Same data, same units, same number of coarse bins.  The only difference
    // is that the streamed path bins onto a fine grid first.  If these two
    // disagree, the memory-efficient path is not measuring what the exact one
    // measures, and every large-scale result would be suspect.
    Sample A = makeSample("A", 51, 5000, Law::Normal);
    Sample B = makeSample("B", 52, 5000, Law::MatchedTail);
    // Force events == blocks so both paths use identical resampling units.
    const int nBlk = 60;
    for (auto* S : {&A, &B}) {
      for (auto& e : S->event) e = e % nBlk;
      S->compactUnits();
    }
    A.normaliseWeights();
    B.normaliseWeights();
 
    const ShapeResult exact =
        shapeTest(buildPooled(A, B, 0, 64), "x", 200, 4242, 0, 0.0);
 
    HistSpec sp;
    sp.nvar = 1;
    sp.nfine = 4096;
    double lo = 1e300, hi = -1e300;
    for (const Sample* S : {&A, &B})
      for (double x : S->var[0]) { lo = std::min(lo, x); hi = std::max(hi, x); }
    sp.lo = {lo - 1e-9};
    sp.hi = {hi + 1e-9};
 
    UnitHists HA, HB;
    HA.alloc(sp, nBlk);
    HB.alloc(sp, nBlk);
    for (int side = 0; side < 2; ++side) {
      const Sample& S = side ? B : A;
      UnitHists& H = side ? HB : HA;
      for (std::size_t i = 0; i < S.rows(); ++i) {
        const double v = S.var[0][i];
        H.fill(S.event[i], S.q[i], &v, S.w[i]);
      }
    }
    const std::vector<int> fe = equalWeightFineEdges(HA, HB, 0, -1, 64);
    const ShapeResult streamed =
        shapeTest(buildPooledFromHists(HA, HB, 0, -1, fe), "x", 200, 4242, 0, 0.0);
 
    std::printf("        exact    A^2 = %10.4f   p = %.4g\n", exact.obs.AD, exact.pAD);
    std::printf("        streamed A^2 = %10.4f   p = %.4g\n", streamed.obs.AD,
                streamed.pAD);
    const double rel =
        std::fabs(streamed.obs.AD - exact.obs.AD) / std::max(1e-12, exact.obs.AD);
    std::printf("        relative difference in A^2 = %.3f%%\n", 100.0 * rel);
    check(rel < 0.05, "streamed A^2 matches the exact value to better than 5%");
    check(streamed.pAD == exact.pAD, "same permutation p-value");
    check(HA.clamped == 0 && HB.clamped == 0, "nothing clamped at the range edges");
 
    // Scale invariance must survive the streamed path too.
    UnitHists HB2 = HB;
    for (auto& x : HB2.h) x *= 1000.0;
    for (auto& x : HB2.uSw) x *= 1000.0;
    for (auto& x : HB2.uSw2) x *= 1e6;
    const ShapeResult scaled =
        shapeTest(buildPooledFromHists(HA, HB2, 0, -1,
                                       equalWeightFineEdges(HA, HB2, 0, -1, 64)),
                  "x", 200, 4242, 0, 0.0);
    check(std::fabs(scaled.obs.AD - streamed.obs.AD) < 1e-9 * streamed.obs.AD,
          "streamed path is invariant to a 1000x weight rescaling");
  }
 
  std::printf("\n=== 8. splitUnits closure on the streamed path ===\n");
  {
    Sample A = makeSample("A", 61, 8000, Law::Normal);
    const int nBlk = 80;
    for (auto& e : A.event) e = e % nBlk;
    A.compactUnits();
    HistSpec sp;
    sp.nvar = 1;
    sp.nfine = 2048;
    double lo = 1e300, hi = -1e300;
    for (double x : A.var[0]) { lo = std::min(lo, x); hi = std::max(hi, x); }
    sp.lo = {lo - 1e-9};
    sp.hi = {hi + 1e-9};
    UnitHists H;
    H.alloc(sp, nBlk);
    for (std::size_t i = 0; i < A.rows(); ++i) {
      const double v = A.var[0][i];
      H.fill(A.event[i], A.q[i], &v, A.w[i]);
    }
    auto halves = splitUnits(H, 99);
    const ShapeResult r =
        shapeTest(buildPooledFromHists(halves.first, halves.second, 0, -1,
                                       equalWeightFineEdges(halves.first, halves.second,
                                                            0, -1, 48)),
                  "x", 200, 5150, 0, 0.0);
    std::printf("        units %d / %d, p_AD = %.3f\n", halves.first.nUnits,
                halves.second.nUnits, r.pAD);
    check(r.pAD > 0.01, "splitting one sample by units gives no false difference");
  }
 
  std::printf("\n=== 9. coarsening the fine grid is exact ===\n");
  {
    // Merging 4 adjacent fine bins must give bit-identical content to having
    // accumulated on the 4x coarser grid in the first place.  That is what
    // lets the binning-stability scan reuse one read.
    Sample A = makeSample("A", 71, 3000, Law::Normal);
    const int nBlk = 40;
    for (auto& e : A.event) e = e % nBlk;
    A.compactUnits();
    HistSpec fine, coarse;
    fine.nvar = coarse.nvar = 1;
    fine.nfine = 1024;
    coarse.nfine = 256;
    fine.lo = coarse.lo = {-6.0};
    fine.hi = coarse.hi = {6.0};
    UnitHists HF, HC;
    HF.alloc(fine, nBlk);
    HC.alloc(coarse, nBlk);
    for (std::size_t i = 0; i < A.rows(); ++i) {
      const double v = A.var[0][i];
      HF.fill(A.event[i], A.q[i], &v, A.w[i]);
      HC.fill(A.event[i], A.q[i], &v, A.w[i]);
    }
    const UnitHists HM = coarsenFine(HF, 4);
    double maxdiff = 0.0;
    for (std::size_t i = 0; i < HC.h.size(); ++i)
      maxdiff = std::max(maxdiff, std::fabs(HM.h[i] - HC.h[i]));
    std::printf("        max |merged - directly binned| = %.3g\n", maxdiff);
    check(HM.spec.nfine == HC.spec.nfine, "merged grid has the right size");
    check(maxdiff < 1e-12, "merging is exact, so the scan needs no extra read");
  }
 
  std::printf("\n%s  (%d failure%s)\n\n", g_fail ? "SELFTEST FAILED" : "SELFTEST PASSED",
              g_fail, g_fail == 1 ? "" : "s");
  return g_fail ? 1 : 0;
}
