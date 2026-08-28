// muonflux_weight_audit -- is the split sample's weight bookkeeping sound?
//
// A "split" production writes each muon k times carrying w/k.  If that is done
// correctly it is pure variance reduction: the weighted distributions are
// unchanged and only the statistical error moves.  If the division is missing
// or wrong, every weighted comparison against an unsplit sample is biased, and
// no amount of care in the statistics downstream will save it.
//
// THE KEY INVARIANT
// -----------------
// Splitting happens inside one tree entry.  So the total weight carried by an
// entry, sum(w) over the muons in it, is INVARIANT under splitting:
//
//     k copies x (w/k)  =  w
//
// That makes the per-entry weight sum the sharpest available probe, and it
// needs no duplicate detection at all.  If the two samples disagree there,
// they disagree about how much flux an interaction produces -- which is a
// bookkeeping error, not physics.
//
// Three numbers should tell the same story:
//     k from multiplicity   = (muons/entry)_test / (muons/entry)_ref
//     k from mean weight    = (mean w)_ref / (mean w)_test
//     k from per-entry sums = 1  (the invariant above)
// The first two must agree with each other; the third must be 1.
//
// Usage:
//     muonflux_weight_audit REF.root TEST.root
//     muonflux_weight_audit @ref.txt @test.txt --max-events 0 --outdir audit
 
#include "muonflux/Load.h"
#include "muonflux/Log.h"
#include "muonflux/Plot.h"
 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TSystem.h>
 
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
 
namespace {
 
using RVecD = ROOT::RVec<double>;
using RVecI = ROOT::RVec<int>;
 
struct Audit {
  std::string label;
  long long entries = 0, muons = 0, emptyEntries = 0;
  double sw = 0.0, sw2 = 0.0;
  double swPerEntrySum = 0.0;  // sum over entries of (sum w in entry)
  TH1D* hLogW = nullptr;       // log10 of the per-muon weight
  TH1D* hLogSumW = nullptr;    // log10 of the per-entry weight sum
  TH1D* hMult = nullptr;       // muons per entry
 
  void book(const char* tag) {
    hLogW = new TH1D(Form("logw_%s", tag), "", 320, -12.0, 4.0);
    hLogSumW = new TH1D(Form("logsw_%s", tag), "", 360, -12.0, 6.0);
    hMult = new TH1D(Form("mult_%s", tag), "", 200, 0.0, 400.0);
    for (TH1D* h : {hLogW, hLogSumW, hMult}) h->SetDirectory(nullptr);
  }
  double muonsPerEntry() const { return entries ? double(muons) / entries : 0.0; }
  double meanWeight() const { return muons ? sw / muons : 0.0; }
  double meanEntryWeight() const { return entries ? swPerEntrySum / entries : 0.0; }
  double nEff() const { return sw2 > 0 ? sw * sw / sw2 : 0.0; }
};
 
/// Muon selection, kept deliberately identical to the main tool's.
inline bool keep(double px, double py, double pz, int pdg, double w, double pmin,
                 double pmax, const std::vector<int>& pdgAbs) {
  bool ok = false;
  for (int a : pdgAbs) ok = ok || (std::abs(pdg) == a);
  if (!ok || !std::isfinite(w)) return false;
  const double p = std::sqrt(px * px + py * py + pz * pz);
  return (p >= pmin) && (p <= pmax) && (pz > 0.0);
}
 
void run(const std::vector<std::string>& files, const mfl::LoadOpts& o,
         long long maxEntries, Audit& A) {
  ROOT::RDataFrame df(o.tree, files);
  // Single-threaded on purpose: Range is rejected under implicit MT, and a
  // diagnostic over a few hundred thousand entries takes seconds anyway.
  auto node = (maxEntries > 0) ? df.Range(0, static_cast<unsigned>(maxEntries))
                               : df.Range(0, 0);
 
  auto perEntry = [&](const double* w, const int* pdg, const double* px,
                      const double* py, const double* pz, size_t n) {
    double sEntry = 0.0;
    long long nsel = 0;
    for (size_t i = 0; i < n; ++i) {
      if (!keep(px[i], py[i], pz[i], pdg[i], w[i], o.pmin, o.pmax, o.pdgAbs))
        continue;
      A.sw += w[i];
      A.sw2 += w[i] * w[i];
      sEntry += w[i];
      ++nsel;
      if (w[i] > 0.0) A.hLogW->Fill(std::log10(w[i]));
    }
    ++A.entries;
    A.muons += nsel;
    A.hMult->Fill(std::min(399.0, double(nsel)));
    if (nsel == 0) ++A.emptyEntries;
    A.swPerEntrySum += sEntry;
    if (sEntry > 0.0) A.hLogSumW->Fill(std::log10(sEntry));
  };
 
  if (o.source == "plane") {
    const std::string P = o.branch;
    node.Foreach(
        [&](const RVecD& px, const RVecD& py, const RVecD& pz, const RVecI& pdg,
            const RVecI& trk, const RVecD& mcw) {
          std::vector<double> w(px.size(), 1.0);
          for (size_t i = 0; i < px.size(); ++i) {
            const int id = trk[i];
            w[i] = (id >= 0 && static_cast<size_t>(id) < mcw.size()) ? mcw[id] : 1.0;
          }
          perEntry(w.data(), pdg.data(), px.data(), py.data(), pz.data(), px.size());
        },
        {P + ".fPx", P + ".fPy", P + ".fPz", P + ".fPdgCode", P + ".fTrackID",
         "MCTrack.fW"});
  } else {
    node.Foreach(
        [&](const RVecD& px, const RVecD& py, const RVecD& pz, const RVecI& pdg,
            const RVecD& mcw) {
          perEntry(mcw.data(), pdg.data(), px.data(), py.data(), pz.data(), px.size());
        },
        {"MCTrack.fPx", "MCTrack.fPy", "MCTrack.fPz", "MCTrack.fPdgCode",
         "MCTrack.fW"});
  }
}
 
void overlay(TPad* p, TH1D* a, TH1D* b, const char* xt, const char* note,
             const std::string& labRef, const std::string& labTest) {
  p->cd();
  p->SetLogy();
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.10);
  p->SetBottomMargin(0.14);
  auto norm = [](TH1D* h) {
    if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
  };
  norm(a);
  norm(b);
  a->SetLineColor(mfl::colourA());
  b->SetLineColor(mfl::colourB());
  a->SetLineWidth(2);
  b->SetLineWidth(2);
  a->GetXaxis()->SetTitle(xt);
  a->GetYaxis()->SetTitle("fraction of entries");
  a->GetXaxis()->SetTitleSize(0.050);
  a->GetYaxis()->SetTitleSize(0.050);
  a->SetMaximum(std::max(a->GetMaximum(), b->GetMaximum()) * 30.0);
  a->Draw("HIST");
  b->Draw("HIST SAME");
  TLatex t;
  t.SetNDC();
  t.SetTextAlign(12);
  t.SetTextSize(0.042);
  t.DrawLatex(0.13, 0.955, note);
  t.SetTextSize(0.034);
  t.SetTextColor(mfl::colourA());
  t.DrawLatex(0.55, 0.88, ("ref:  " + labRef).c_str());
  t.SetTextColor(mfl::colourB());
  t.DrawLatex(0.55, 0.83, ("test: " + labTest).c_str());
}
 
[[noreturn]] void usage(int code) {
  std::printf(R"(muonflux_weight_audit REF TEST [options]
 
Checks that a split production's weights were divided correctly, using the
fact that sum(w) within one tree entry is invariant under splitting.
 
  --tree NAME       default cbmsim
  --source S        plane | mctrack        (default plane)
  --branch NAME     default PlaneHAPoint
  --pdg LIST        |PDG| codes (default 13)
  --pmin X/--pmax X GeV, default 5 / 400
  --max-events N    entries per sample, 0 = all (default 200000)
  --outdir DIR      default audit/
  --label-ref S / --label-test S
  -h, --help
)");
  std::exit(code);
}
 
}  // namespace
 
int runMain(int argc, char** argv) {
  mfl::LoadOpts o;
  long long maxEntries = 200000;
  std::string outDir = "audit", labRef, labTest;
  std::vector<std::string> pos;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) usage(2);
      return argv[++i];
    };
    if (a == "-h" || a == "--help") usage(0);
    else if (a == "--tree") o.tree = next();
    else if (a == "--source") o.source = next();
    else if (a == "--branch" || a == "--plane") o.branch = next();
    else if (a == "--pdg") {
      o.pdgAbs.clear();
      for (const std::string& v : mfl::splitList(next()))
        o.pdgAbs.push_back(std::abs(std::atoi(v.c_str())));
    }
    else if (a == "--pmin") o.pmin = std::atof(next().c_str());
    else if (a == "--pmax") o.pmax = std::atof(next().c_str());
    else if (a == "--max-events") maxEntries = std::atoll(next().c_str());
    else if (a == "--outdir") outDir = next();
    else if (a == "--label-ref") labRef = next();
    else if (a == "--label-test") labTest = next();
    else if (!a.empty() && a[0] == '-') usage(2);
    else pos.push_back(a);
  }
  if (pos.size() < 2) usage(2);
  if (labRef.empty()) labRef = "reference";
  if (labTest.empty()) labTest = "test";
 
  TStopwatch clock;
  clock.Start();
  mfl::setStyle();
  gSystem->mkdir(outDir.c_str(), kTRUE);
 
  Audit R, T;
  R.label = labRef;
  T.label = labTest;
  R.book("ref");
  T.book("test");
  run(mfl::splitList(pos[0]), o, maxEntries, R);
  run(mfl::splitList(pos[1]), o, maxEntries, T);
 
  const double kMult = R.muonsPerEntry() > 0 ? T.muonsPerEntry() / R.muonsPerEntry() : 0.0;
  const double kWeight = T.meanWeight() > 0 ? R.meanWeight() / T.meanWeight() : 0.0;
  const double entryRatio =
      R.meanEntryWeight() > 0 ? T.meanEntryWeight() / R.meanEntryWeight() : 0.0;
 
  mfl::logf("================ weight audit ================");
  mfl::logf("%-26s %16s %16s", "", labRef.c_str(), labTest.c_str());
  mfl::logf("%s", std::string(60, '-').c_str());
  mfl::logf("%-26s %16lld %16lld", "entries read", R.entries, T.entries);
  mfl::logf("%-26s %16lld %16lld", "muons selected", R.muons, T.muons);
  mfl::logf("%-26s %16lld %16lld", "entries with no muon", R.emptyEntries,
            T.emptyEntries);
  mfl::logf("%-26s %16.3f %16.3f", "muons / entry", R.muonsPerEntry(),
            T.muonsPerEntry());
  mfl::logf("%-26s %16.6g %16.6g", "mean weight / muon", R.meanWeight(),
            T.meanWeight());
  mfl::logf("%-26s %16.6g %16.6g", "mean sum(w) / entry", R.meanEntryWeight(),
            T.meanEntryWeight());
  mfl::logf("%-26s %16.0f %16.0f", "N_eff", R.nEff(), T.nEff());
  mfl::logf("%s", std::string(60, '-').c_str());
  mfl::logf("");
  mfl::logf("implied split factor k");
  mfl::logf("   from multiplicity        k = %.3f", kMult);
  mfl::logf("   from mean weight         k = %.3f", kWeight);
  const double agree =
      (kMult > 0 && kWeight > 0) ? std::fabs(std::log(kWeight / kMult)) : 99.0;
  mfl::logf("   agreement                %s (%.1f%% apart)",
            agree < 0.22 ? "CONSISTENT" : "INCONSISTENT",
            100.0 * std::fabs(kWeight / std::max(1e-12, kMult) - 1.0));
  mfl::logf("");
  mfl::logf("INVARIANT: mean sum(w) per entry, test / ref = %.4f", entryRatio);
  if (std::fabs(entryRatio - 1.0) < 0.10)
    mfl::logf("   -> consistent with 1: splitting preserved the weight per entry.");
  else if (kMult > 1.5 && std::fabs(entryRatio / kMult - 1.0) < 0.25)
    mfl::logf("   -> approximately k, NOT 1: the copies were never divided by k. "
              "Every weighted comparison against the unsplit sample is biased.");
  else
    mfl::logf("   -> neither 1 nor k.  The two samples disagree about the flux per "
              "interaction for some other reason; investigate before comparing.");
  mfl::logf("");
  mfl::logf("Caveat: this tests the weight bookkeeping in aggregate.  A per-muon");
  mfl::logf("misallocation that preserves the per-entry sum would pass it.");
 
  // ---- plots ----
  TCanvas* c = new TCanvas("c", "audit", 1000, 900);
  TPad* p1 = new TPad("p1", "", 0.0, 0.50, 1.0, 1.0);
  TPad* p2 = new TPad("p2", "", 0.0, 0.0, 0.5, 0.50);
  TPad* p3 = new TPad("p3", "", 0.5, 0.0, 1.0, 0.50);
  p1->Draw();
  p2->Draw();
  p3->Draw();
  overlay(p1, R.hLogSumW, T.hLogSumW, "log_{10}( #Sigma w per entry )",
          "THE INVARIANT: these must coincide if splitting is correct", labRef,
          labTest);
  overlay(p2, R.hLogW, T.hLogW, "log_{10}( w per muon )",
          Form("per-muon weight (offset by k = %.3g expected)", kMult), labRef,
          labTest);
  overlay(p3, R.hMult, T.hMult, "muons per entry",
          Form("multiplicity (ratio k = %.3g)", kMult), labRef, labTest);
  const std::string pdf = outDir + "/weight_audit.pdf";
  c->Print(pdf.c_str());
  mfl::writeLog(outDir + "/weight_audit.txt");
  mfl::logf("[out] %s and %s/weight_audit.txt   [%.1f s]", pdf.c_str(), outDir.c_str(),
            clock.RealTime());
  return 0;
}
 
int main(int argc, char** argv) {
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
