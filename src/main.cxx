// flux_compare -- are the muons delivered by two productions compatible?
//
// The first sample is the REFERENCE, the second is the TEST.  Everything is
// reported as test relative to reference.
//
// Three tests, all invariant under the arbitrary per-production normalisation:
//   1. mu+ / mu-                 a relative amount, free of any rate scale
//   2. flux fractions in regions what share of the flux sits where
//   3. 1-D shape tests           weighted two-sample Anderson-Darling /
//                                Cramer-von Mises / Kolmogorov-Smirnov, with a
//                                permutation null over resampling units
//
// sum(w) is used only as a denominator.  It is never compared, because it is
// one of the things the productions vary on purpose.

#include "fluxval/Config.h"
#include "fluxval/Core.h"
#include "fluxval/Load.h"
#include "fluxval/Log.h"
#include "fluxval/Plot.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <vector>

namespace {

struct Opts {
  std::string fileRef, fileTest, labelRef, labelTest;
  std::string outDir = "outputs";
  std::string prefix = "flux";
  fluxval::LoadOpts load;
  int nperm = 500;
  int nboot = 400;
  int nbinsTest = 200;
  int nbinsPlot = 40;
  int nthreads = 0;
  unsigned long seed = 20260826;
  double adClip = 0.0;
  bool splitCharge = false;
  bool closure = false;
  int closureScan = 0;
  bool combinedPdf = false;
  bool binningCheck = true;
  bool png = false;
  std::vector<std::string> varNames;
};

std::string sfmt(const char* f, ...) {
  char b[2048];
  va_list ap;
  va_start(ap, f);
  std::vsnprintf(b, sizeof(b), f, ap);
  va_end(ap);
  return std::string(b);
}
std::string rule(int n) { return std::string(n, '-'); }

[[noreturn]] void usage(int code) {
  std::printf(R"(flux_compare REFERENCE TEST [options]

  REFERENCE, TEST   a .root file, a comma-separated list, a directory, or
                    @file-with-one-path-per-line.  All files of a sample are
                    read together as one dataset.

  --label-ref S / --label-test S    legend labels
  --outdir DIR      output directory (default outputs/)
  --prefix S        filename prefix (default flux)
  --combined-pdf    also write one multi-page PDF alongside the per-plot files
  --png             also write a .png next to every .pdf

  --tree NAME       default cbmsim
  --source S        plane | mctrack            (default plane)
  --branch NAME     scoring-plane branch (default PlaneHAPoint).  Any branch
                    exposing fX, fY, fPx, fPy, fPz, fPdgCode, fTrackID works.
  --pdg LIST        |PDG| codes to select (default 13 = muons).
                    e.g. --pdg 14 for muon neutrinos, --pdg 12,14,16 for all.
                    For a neutral species q*x is identically zero and is
                    dropped automatically; the sign split becomes nu / nubar.
  --pmin X / --pmax X   GeV, default 5 / 400
  --vars LIST       subset of log10_p,log10_pz,asinh_pT,log10_theta,x,y,qx,r
  --split-sign      also run every shape test per sign class separately

  --max-events N    read only the first N entries of each sample's file list.
                    Deterministic (a prefix, not a random subset), but it
                    disables multi-threading, so use it for quick checks only.
  --no-progress     suppress the progress bar (for batch logs)
  --max-files N     use only the first N files of each list

  --nperm N         permutations (default 500)
  --nboot N         bootstrap replicas (default 400)
  --nfine N         fine bins per variable in the accumulator (default 2048)
  --no-binning-check  skip the A^2 stability scan over coarser grids
  --nbins-test N    coarse bins for the EDF statistics (default 200)
  --nbins-plot N    coarse bins for the plots (default 40)
  --max-units N     cap on resampling units (default 256)
  --nblocks N       event blocks when there are too few files (default 200)
  --min-files N     use files as units at or above this count (default 5)
  --ad-clip X       drop bins with H(1-H) <= X
  --threads N       0 = all cores
  --seed N

  --closure         ignore TEST: split REFERENCE in half by unit and compare it
                    with itself.  Run this first; every p-value should be flat.
  --closure-scan N  repeat that N times with different seeds and emit a P-P
                    plot of the p-values against uniform.
  -h, --help
)");
  std::exit(code);
}

std::string shorten(const std::string& s, size_t n = 54) {
  return (s.size() <= n) ? s : "..." + s.substr(s.size() - (n - 3));
}
std::string deriveLabel(const std::string& spec, size_t nfiles) {
  std::string b = spec;
  if (!b.empty() && b[0] == '@') b = b.substr(1);
  if (nfiles > 1) b += "  (" + std::to_string(nfiles) + " files)";
  return shorten(b);
}

Opts parse(int argc, char** argv) {
  Opts o;
  std::vector<std::string> pos;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) usage(2);
      return argv[++i];
    };
    if (a == "-h" || a == "--help") usage(0);
    else if (a == "--label-ref") o.labelRef = next();
    else if (a == "--label-test") o.labelTest = next();
    else if (a == "--outdir") o.outDir = next();
    else if (a == "--prefix") o.prefix = next();
    else if (a == "--combined-pdf") o.combinedPdf = true;
    else if (a == "--no-binning-check") o.binningCheck = false;
    else if (a == "--png") o.png = true;
    else if (a == "--tree") o.load.tree = next();
    else if (a == "--source") o.load.source = next();
    else if (a == "--branch" || a == "--plane") o.load.branch = next();
    else if (a == "--pdg") {
      o.load.pdgAbs.clear();
      for (const std::string& v : fluxval::splitList(next()))
        o.load.pdgAbs.push_back(std::abs(std::atoi(v.c_str())));
    }
    else if (a == "--pmin") o.load.pmin = std::atof(next().c_str());
    else if (a == "--pmax") o.load.pmax = std::atof(next().c_str());
    else if (a == "--vars") o.varNames = fluxval::splitList(next());
    else if (a == "--split-sign" || a == "--split-charge") o.splitCharge = true;
    else if (a == "--max-events") o.load.maxEvents = std::atoll(next().c_str());
    else if (a == "--no-progress") o.load.progress = false;
    else if (a == "--max-files") o.load.maxFiles = std::atoi(next().c_str());
    else if (a == "--nperm") o.nperm = std::atoi(next().c_str());
    else if (a == "--nboot") o.nboot = std::atoi(next().c_str());
    else if (a == "--nfine") o.load.nfine = std::atoi(next().c_str());
    else if (a == "--nbins-test") o.nbinsTest = std::atoi(next().c_str());
    else if (a == "--nbins-plot") o.nbinsPlot = std::atoi(next().c_str());
    else if (a == "--max-units") o.load.maxUnits = std::atoi(next().c_str());
    else if (a == "--nblocks") o.load.nBlocks = std::atoi(next().c_str());
    else if (a == "--min-files") o.load.minFiles = std::atoi(next().c_str());
    else if (a == "--ad-clip") o.adClip = std::atof(next().c_str());
    else if (a == "--threads") o.nthreads = std::atoi(next().c_str());
    else if (a == "--seed") o.seed = std::strtoul(next().c_str(), nullptr, 10);
    else if (a == "--closure") o.closure = true;
    else if (a == "--closure-scan") o.closureScan = std::atoi(next().c_str());
    else if (!a.empty() && a[0] == '-') {
      std::printf("unknown option %s\n\n", a.c_str());
      usage(2);
    } else pos.push_back(a);
  }
  if (pos.empty()) usage(2);
  o.fileRef = pos[0];
  o.fileTest = (pos.size() > 1) ? pos[1] : pos[0];
  return o;
}

struct VarResult {
  int var = 0;
  std::string tag;
  fluxval::ShapeResult shape;
  std::vector<int> fineEdges;
  std::vector<double> edges;
  fluxval::RatioBand band;
  std::vector<double> neffRef, neffTest;
};

}  // namespace

int runMain(int argc, char** argv) {
  const Opts o = parse(argc, argv);
  TStopwatch clock;
  clock.Start();
  fluxval::setStyle();
  gSystem->mkdir(o.outDir.c_str(), kTRUE);
  const std::string base = o.outDir + "/" + o.prefix;

  // One PDF per plot is the default; the combined book is opt-in.
  const std::string bookPdf = base + "_all.pdf";
  bool bookOpen = false;
  int pageNo = 0;
  TCanvas* c = nullptr;
  auto emit = [&](const std::string& tag) {
    const std::string stem = sfmt("%s_%02d_%s", base.c_str(), pageNo, tag.c_str());
    c->Print((stem + ".pdf").c_str());
    if (o.png) c->Print((stem + ".png").c_str());
    if (o.combinedPdf) {
      if (!bookOpen) { c->Print((bookPdf + "[").c_str()); bookOpen = true; }
      c->Print(bookPdf.c_str());
    }
    ++pageNo;
  };

  const std::vector<std::string> filesRef = fluxval::splitList(o.fileRef);
  const std::vector<std::string> filesTest = fluxval::splitList(o.fileTest);

  fluxval::logf("================ flux compatibility ================");
  fluxval::logf("[cfg] reference : %zu file(s)", filesRef.size());
  fluxval::logf("[cfg] test      : %zu file(s)", o.closure ? 0 : filesTest.size());

  // The probe uses Range(), which ROOT rejects under implicit MT, so the range
  // pass runs single-threaded and MT is enabled only for the real pass.
  const fluxval::HistSpec spec =
      fluxval::probeRange(filesRef, o.closure ? std::vector<std::string>{} : filesTest,
                      o.load);
  for (int j = 0; j < fluxval::kNVar; ++j)
    fluxval::logQuiet(sfmt("[cfg] range %-12s [%+.4g, %+.4g]", fluxval::varDefs()[j].name,
                       spec.lo[j], spec.hi[j]));

  // RDataFrame rejects Range under implicit MT, so a limited read is
  // single-threaded.  Stated rather than silently slow.
  if (o.load.maxEvents > 0)
    fluxval::logf("[cfg] --max-events %lld: reading single-threaded (Range and "
                  "implicit MT are incompatible)", o.load.maxEvents);
  else
    ROOT::EnableImplicitMT(o.nthreads > 0 ? static_cast<unsigned>(o.nthreads) : 0u);
  const unsigned nSlots = std::max(1u, ROOT::GetThreadPoolSize());
  const double mbPerSlot = double(o.load.maxUnits) * 2.0 * int(fluxval::kNVar) * o.load.nfine * 8.0 /
                           (1024.0 * 1024.0);
  fluxval::logf("[cfg] accumulator %.0f MB/thread x %u threads = %.0f MB", mbPerSlot,
            nSlots, mbPerSlot * nSlots);

  fluxval::logf("[1/4] reading reference");
  TStopwatch st;
  st.Start();
  fluxval::UnitHists Ref = fluxval::loadHists(filesRef, o.load, spec, "reference");
  fluxval::logf("[1/4] done in %.0f s", st.RealTime());
  fluxval::UnitHists Test;
  std::string labRef = o.labelRef.empty() ? deriveLabel(o.fileRef, filesRef.size())
                                          : o.labelRef;
  std::string labTest = o.labelTest.empty()
                            ? deriveLabel(o.fileTest, filesTest.size())
                            : o.labelTest;

  if (o.closure || o.closureScan > 0) {
    auto h = fluxval::splitUnits(Ref, o.seed ^ 0xABCDEF);
    labTest = labRef + "  [closure half 2]";
    labRef += "  [closure half 1]";
    Test = h.second;
    Ref = h.first;
    fluxval::logf("[closure] reference split into %d / %d units", Ref.nUnits, Test.nUnits);
  } else {
    fluxval::logf("[2/4] reading test");
    st.Start();
    Test = fluxval::loadHists(filesTest, o.load, spec, "test");
    fluxval::logf("[2/4] done in %.0f s", st.RealTime());
  }

  fluxval::logf("[rate] sum(w): ref = %.6g   test = %.6g   ratio = %.4g", Ref.totalW(),
            Test.totalW(), Ref.totalW() > 0 ? Test.totalW() / Ref.totalW() : 0.0);
  fluxval::logf("[rate] NOT compared -- it is varied on purpose. Denominator only.");

  const fluxval::Species sp = fluxval::makeSpecies(o.load.pdgAbs);
  fluxval::logf("[cfg] species    %s (|PDG| =%s)%s", sp.name.c_str(),
            [&] {
              static std::string t;
              t.clear();
              for (int a : sp.pdgAbs) t += " " + std::to_string(a);
              return t.c_str();
            }(),
            sp.neutral ? "  -- neutral: q*x dropped, sign split is nu/nubar" : "");

  const auto& defs = fluxval::varDefs();
  std::vector<int> vars;
  if (o.varNames.empty()) {
    vars = fluxval::defaultVars(sp);
  } else {
    for (const std::string& nm : o.varNames) {
      int f = -1;
      for (int j = 0; j < fluxval::kNVar; ++j) if (nm == defs[j].name) f = j;
      if (f < 0) { std::printf("unknown variable '%s'\n", nm.c_str()); return 2; }
      if (sp.neutral && f == fluxval::kQX) {
        std::printf("qx is identically zero for a neutral species; drop it from "
                    "--vars\n");
        return 2;
      }
      vars.push_back(f);
    }
  }

  // ---------------- calibration scan ----------------
  if (o.closureScan > 0) {
    fluxval::logf("");
    fluxval::logf("---- closure scan: %d independent splits of the reference ----",
              o.closureScan);
    std::vector<std::vector<double>> pv(vars.size());
    std::vector<std::string> names;
    for (int j : vars) names.push_back(defs[j].name);
    for (int t = 0; t < o.closureScan; ++t) {
      auto h = fluxval::splitUnits(Ref, o.seed + 1000u * (t + 1));
      for (size_t vi = 0; vi < vars.size(); ++vi) {
        const std::vector<int> fe =
            fluxval::equalWeightFineEdges(h.first, h.second, vars[vi], -1, o.nbinsTest);
        pv[vi].push_back(fluxval::shapeTest(
                             fluxval::buildPooledFromHists(h.first, h.second, vars[vi], -1, fe),
                             names[vi], o.nperm, o.seed + 7777u * (t + 1) + vi,
                             o.nthreads, o.adClip)
                             .pAD);
      }
      std::printf("\r[scan] %d / %d", t + 1, o.closureScan);
      std::fflush(stdout);
    }
    std::printf("\r                    \r");
    fluxval::logf("%-14s %8s %8s %10s %11s", "variable", "mean p", "median", "frac<0.05",
              "KS vs unif");
    fluxval::logf("%s", rule(56).c_str());
    for (size_t vi = 0; vi < vars.size(); ++vi) {
      std::vector<double> p = pv[vi];
      std::sort(p.begin(), p.end());
      const int n = static_cast<int>(p.size());
      double m = 0.0, D = 0.0;
      for (double x : p) m += x;
      m /= n;
      int below = 0;
      for (double x : p) if (x < 0.05) ++below;
      for (int k = 0; k < n; ++k)
        D = std::max(D, std::max(std::fabs((k + 1.0) / n - p[k]),
                                 std::fabs(p[k] - double(k) / n)));
      fluxval::logf("%-14s %8.3f %8.3f %10.2f %11.3f", names[vi].c_str(), m, p[n / 2],
                double(below) / n, D);
    }
    fluxval::logf("%s", rule(56).c_str());
    fluxval::logf("Calibrated: mean p ~ 0.50, frac<0.05 ~ 0.05.");
    fluxval::logf("Above 0.50 = conservative (misses real differences).");
    fluxval::logf("Below 0.50 = anti-conservative (p-values not trustworthy).");

    c = new TCanvas("c", "scan", 950, 780);
    fluxval::drawPP(c, names, pv);
    emit("pp_calibration");
    std::vector<std::string> L(fluxval::logLines());
    for (size_t q = 0; q < L.size(); q += 34) {
      fluxval::drawTextPage(c, "closure scan -- summary",
                        {L.begin() + q, L.begin() + std::min(L.size(), q + 34)});
      emit("summary");
    }
    if (bookOpen) c->Print((bookPdf + "]").c_str());
    fluxval::writeLog(base + ".txt");
    std::printf("[out] %d plot file(s) and %s.txt in %s/\n", pageNo, o.prefix.c_str(),
                o.outDir.c_str());
    return 0;
  }

  // ---------------- tests 1 and 2 ----------------
  const std::vector<fluxval::Region> regs = fluxval::defaultRegions(sp);
  const std::vector<fluxval::Quantity> qs = fluxval::defaultQuantities(sp, regs);
  fluxval::logf("[3/4] fractions and bootstrap (%d replicas)", o.nboot);
  st.Start();
  const std::vector<fluxval::QuantityResult> qres = fluxval::compareQuantitySums(
      fluxval::unitSumsFromHists(Ref, regs), fluxval::unitSumsFromHists(Test, regs), qs,
      o.nboot, o.seed);

  fluxval::logf("[3/4] done in %.0f s", st.RealTime());
  std::vector<std::string> summary;
  auto say = [&](const std::string& s) {
    summary.push_back(s);
    fluxval::logLine(s);
  };
  auto fmtPM = [](double v, double e, bool ok) {
    if (!ok) return sfmt("%16s", "--");
    if (std::fabs(v) >= 1e-3 || v == 0.0) return sfmt("%7.5f+-%-7.5f", v, e);
    return sfmt("%8.2e+-%-6.1e", v, e);
  };

  say("");
  say(sfmt("---- tests 1 & 2: scale-free amounts (bootstrap unit: %s) ----",
           Ref.unitsAreFiles ? "files" : "event blocks"));
  say(sfmt("%-22s %-16s %-16s %10s %8s", "quantity", "reference", "test",
           "(T-R)/R", "pull"));
  say(rule(78));
  int nBad = 0, nDead = 0;
  for (const auto& r : qres) {
    const std::string rel =
        r.comparable ? sfmt("%+9.2f%%", 100.0 * r.relDiff) : sfmt("%10s", "--");
    const std::string pl =
        (r.comparable && r.ediff > 0.0) ? sfmt("%+8.2f", r.pull) : sfmt("%8s", "--");
    say(sfmt("%-22s %s %s %s %s%s", r.name.c_str(), fmtPM(r.a, r.ea, r.okA).c_str(),
             fmtPM(r.b, r.eb, r.okB).c_str(), rel.c_str(), pl.c_str(),
             (r.comparable && std::fabs(r.pull) > 3.0) ? "  <<<" : ""));
    if (r.comparable && std::fabs(r.pull) > 3.0) ++nBad;
    if (!r.comparable) ++nDead;
  }
  say(rule(78));
  say(sfmt("%d of %zu comparable quantities differ by more than 3 sigma", nBad,
           qres.size() - nDead));
  if (nDead)
    say(sfmt("%d marked '--': empty region, or known to worse than 50%%, so there is "
             "nothing to compare (not a null result).", nDead));

  // ---------------- test 3 ----------------
  std::vector<std::pair<std::string, int>> charges = {{"both", -1}};
  if (o.splitCharge) {
    charges.push_back({sp.posLabel, 0});
    charges.push_back({sp.negLabel, 1});
  }
  // Binning-stability scan.  A real difference does not care how finely the
  // substrate is binned; an artefact from a sharp feature straddling a bin
  // edge does.  Merging fine bins is exact, so all three resolutions come from
  // the single read already done -- no extra I/O.
  std::vector<fluxval::UnitHists> refGrid{Ref}, testGrid{Test};
  std::vector<int> gridFactor{1};
  if (o.binningCheck) {
    for (int f : {4, 16}) {
      if (Ref.spec.nfine % f != 0 || Ref.spec.nfine / f < 64) continue;
      refGrid.push_back(fluxval::coarsenFine(Ref, f));
      testGrid.push_back(fluxval::coarsenFine(Test, f));
      gridFactor.push_back(f);
    }
  }

  std::vector<VarResult> results;
  for (const auto& ch : charges) {
    say("");
    say(sfmt("---- test 3: 1-D shapes [%s], %d permutations, unit = %s ----",
             ch.first.c_str(), o.nperm, Ref.unitsAreFiles ? "file" : "event block"));
    std::string hdr = sfmt("%-14s %10s %9s %8s %10s %9s %9s", "variable", "A^2",
                           "p(AD)", "z(AD)", "omega^2", "p(CvM)", "p(KS)");
    for (size_t gi = 1; gi < refGrid.size(); ++gi)
      hdr += sfmt(" %9s", sfmt("A^2/%d", gridFactor[gi]).c_str());
    say(hdr);
    say(rule(static_cast<int>(hdr.size())));
    int nUnstable = 0;
    int done = 0;
    const int ntot = static_cast<int>(vars.size());
    for (int j : vars) {
      if (o.load.progress) {
        std::fprintf(stderr, "\r[4/4] shape tests [%s] %d/%d  %-14s", ch.first.c_str(),
                     done + 1, ntot, defs[j].name);
        std::fflush(stderr);
      }
      VarResult v;
      v.var = j;
      v.tag = ch.first;
      v.fineEdges = fluxval::equalWeightFineEdges(Ref, Test, j, ch.second, o.nbinsTest);
      v.shape = fluxval::shapeTest(
          fluxval::buildPooledFromHists(Ref, Test, j, ch.second, v.fineEdges), defs[j].name,
          o.nperm, o.seed + 100u * j, o.nthreads, o.adClip);
      const std::vector<int> pe =
          fluxval::equalWeightFineEdges(Ref, Test, j, ch.second, o.nbinsPlot);
      v.edges = fluxval::edgeValues(spec, j, pe);
      v.band = fluxval::shapeRatioFromHists(Ref, Test, j, ch.second, pe, o.nboot,
                                        o.seed + 200u * j);
      v.neffRef = fluxval::nEffPerBinFromHists(Ref, j, ch.second, pe);
      v.neffTest = fluxval::nEffPerBinFromHists(Test, j, ch.second, pe);
      // Same statistic on progressively coarser substrates, no permutations
      // needed -- only the value matters, not its p-value.
      std::string stab;
      double worst = 1.0;
      for (size_t gi = 1; gi < refGrid.size(); ++gi) {
        const std::vector<int> fe =
            fluxval::equalWeightFineEdges(refGrid[gi], testGrid[gi], j, ch.second,
                                      o.nbinsTest);
        const double a2 =
            fluxval::shapeTest(fluxval::buildPooledFromHists(refGrid[gi], testGrid[gi], j,
                                                     ch.second, fe),
                           defs[j].name, 0, o.seed, o.nthreads, o.adClip)
                .obs.AD;
        stab += sfmt(" %9.4g", a2);
        if (v.shape.obs.AD > 0.0)
          worst = std::max(worst, std::max(a2 / v.shape.obs.AD,
                                           v.shape.obs.AD / std::max(a2, 1e-12)));
      }
      const bool unstable = worst > 1.35;
      if (unstable) ++nUnstable;
      const bool fl = v.shape.pAD <= 1.0 / (o.nperm + 1.0) + 1e-12;
      say(sfmt("%-14s %10.4g %2s%-7.4g %+8.1f %10.4g %9.4g %9.4g%s%s", defs[j].name,
               v.shape.obs.AD, fl ? "<" : "", v.shape.pAD, v.shape.zAD, v.shape.obs.CvM,
               v.shape.pCvM, v.shape.pKS, stab.c_str(),
               unstable ? "  BINNING!" : ""));
      results.push_back(std::move(v));
      ++done;
    }
    if (o.load.progress) std::fprintf(stderr, "\r%70s\r", "");
    say(rule(static_cast<int>(hdr.size())));
    say("p with '<' is at the permutation floor 1/(nperm+1); raise --nperm.");
    say("z is NOT a significance: the AD null is strongly skewed.  Quote p.");
    if (!refGrid.empty() && refGrid.size() > 1) {
      say(sfmt("A^2/N columns repeat the statistic on a %dx and %dx coarser fine grid.",
               gridFactor.size() > 1 ? gridFactor[1] : 0,
               gridFactor.size() > 2 ? gridFactor[2] : 0));
      if (nUnstable)
        say(sfmt("%d variable(s) flagged BINNING!: A^2 moves by more than 35%% with the "
                 "grid, so the value is not converged.  Raise --nfine and re-check "
                 "before quoting it.", nUnstable));
      else
        say("All variables stable against the grid.");
    }
  }

  // ---------------- output ----------------
  c = new TCanvas("c", "flux", 950, 780);
  std::vector<std::string> header = {
      "reference : " + labRef, "test      : " + labTest, "",
      sfmt("source        %s, branch %s", o.load.source.c_str(),
           o.load.branch.c_str()),
      sfmt("species       %s", sp.name.c_str()),
      sfmt("momentum      %.1f - %.1f GeV", o.load.pmin, o.load.pmax),
      sfmt("reference     %lld muons, %d units, N_eff = %.0f", Ref.rows, Ref.nUnits,
           Ref.nEffTotal()),
      sfmt("test          %lld muons, %d units, N_eff = %.0f", Test.rows, Test.nUnits,
           Test.nEffTotal()),
      sfmt("sum(w) ratio  %.4g   -- recorded, never tested",
           Ref.totalW() > 0 ? Test.totalW() / Ref.totalW() : 0.0),
      sfmt("accumulator   %d fine bins, %d units", o.load.nfine, Ref.nUnits),
      sfmt("max-events    %lld  (0 = all)", o.load.maxEvents), ""};
  header.insert(header.end(), summary.begin(), summary.end());
  for (size_t p = 0; p < header.size(); p += 34) {
    fluxval::drawTextPage(c, p == 0 ? "flux compatibility -- summary"
                                : "summary (continued)",
                      {header.begin() + p,
                       header.begin() + std::min(header.size(), p + 34)});
    emit("summary");
  }

  fluxval::drawQuantityPulls(c, qres, labRef, labTest);
  emit("fractions");
  // One p-value page per charge selection: stacking all three in a single
  // frame gave three unlabelled blocks of identical variable names.
  for (const auto& ch : charges) {
    std::vector<fluxval::ShapeResult> all;
    for (const VarResult& v : results)
      if (v.tag == ch.first) all.push_back(v.shape);
    if (all.empty()) continue;
    fluxval::drawPValueOverview(c, all, ch.first, labRef, labTest);
    emit("pvalues_" + ch.first);
  }

  TFile fo((base + ".root").c_str(), "RECREATE");
  for (const VarResult& v : results) {
    const std::string sfx = std::string(defs[v.var].name) + "_" + v.tag;
    TH1D* hR = fluxval::densityHist(("ref_" + sfx).c_str(), v.edges, v.band.shapeA,
                                v.band.errA);
    TH1D* hT = fluxval::densityHist(("test_" + sfx).c_str(), v.edges, v.band.shapeB,
                                v.band.errB);
    TH1D* hRat = fluxval::ratioHist(("ratio_" + sfx).c_str(), v.edges, v.band.ratio,
                                v.band.err);
    const bool fl = v.shape.pAD <= 1.0 / (o.nperm + 1.0) + 1e-12;
    fluxval::drawWithRatio(c, hR, hT, hRat, defs[v.var].axis, defs[v.var].logyPlot, labRef,
                       labTest,
                       sfmt("%s [%s]   A^{2} = %.4g,  p %s %.3g", defs[v.var].name,
                            v.tag.c_str(), v.shape.obs.AD, fl ? "<" : "", v.shape.pAD));
    emit("shape_" + sfx);
    fo.cd();
    hR->Write();
    hT->Write();
    hRat->Write();
  }
  for (const VarResult& v : results) {
    fluxval::drawProfile(c, v.shape, defs[v.var].axis, v.tag, labRef, labTest);
    emit(sfmt("profile_%s_%s", defs[v.var].name, v.tag.c_str()));
  }
  for (const VarResult& v : results) {
    fluxval::drawNull(c, v.shape);
    emit(sfmt("null_%s_%s", defs[v.var].name, v.tag.c_str()));
  }
  for (const VarResult& v : results) {
    fluxval::drawNeff(c, v.edges, v.neffRef, v.neffTest,
                  std::string(defs[v.var].axis) + "   [" + v.tag + "]", labRef, labTest);
    emit(sfmt("neff_%s_%s", defs[v.var].name, v.tag.c_str()));
  }
  if (bookOpen) c->Print((bookPdf + "]").c_str());
  fo.Close();

  fluxval::logf("");
  fluxval::logf("[out] %d plot file(s) in %s/", pageNo, o.outDir.c_str());
  fluxval::logf("[out] %s.root", base.c_str());
  fluxval::logf("[time] %.1f s", clock.RealTime());
  fluxval::writeLog(base + ".txt");
  std::printf("[out] %s.txt\n", base.c_str());
  return 0;
}

int main(int argc, char** argv) {
  // A single unreadable path in a 150-line file list used to abort with an
  // uncaught exception and a core dump, which said nothing about which path.
  try {
    return runMain(argc, argv);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "\nERROR: %s\n", e.what());
    return 1;
  } catch (...) {
    std::fprintf(stderr, "\nERROR: unknown failure\n");
    return 1;
  }
}
