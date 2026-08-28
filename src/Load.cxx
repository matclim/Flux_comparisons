
#include "fluxval/Load.h"
#include "fluxval/Config.h"
#include "fluxval/Log.h"
 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
#include <ROOT/RDFHelpers.hxx>
#define MFL_HAVE_RDF_PROGRESS 1
#endif
#include <atomic>
#include <chrono>
#include <TFile.h>
#include <TROOT.h>
#include <TSystemDirectory.h>
#include <TTree.h>
 
#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
 
namespace fluxval {
 
using RVecD = ROOT::RVec<double>;
using RVecI = ROOT::RVec<int>;
 
namespace {
 
std::string trim(const std::string& s) {
  const size_t p = s.find_first_not_of(" \t\r\n");
  const size_t q = s.find_last_not_of(" \t\r\n");
  return (p == std::string::npos) ? std::string() : s.substr(p, q - p + 1);
}
 
/// splitmix64 -- a cheap, well-mixed hash.  Used to assign events to blocks
/// and to subsample, so both are deterministic and reproducible across runs
/// and across however many threads happen to be available.
inline uint64_t mix(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
 
/// Derived variables for one muon, written into `out` (kNVar entries).
/// Returns false if the muon fails the selection.
inline bool derive(double px, double py, double pz, double x, double y, int pdg,
                   double w, const LoadOpts& o, double* out, int& sgn) {
  bool ok = false;
  for (int a : o.pdgAbs) ok = ok || (std::abs(pdg) == a);
  if (!ok) return false;
  if (!std::isfinite(w)) return false;
  const double pt = std::hypot(px, py);
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (!(p >= o.pmin) || !(p <= o.pmax)) return false;
  if (o.requirePzPositive && !(pz > 0.0)) return false;
 
  sgn = signIndexOfPdg(pdg);
  const int q = chargeOfPdg(pdg);  // zero for neutrinos: q*x is then identically 0
  const double th = std::max(1e-4, std::atan2(pt, std::fabs(pz)) * 1000.0);
  out[kLogP] = std::log10(p);
  out[kLogPz] = std::log10(std::max(1e-6, pz));
  out[kAsinhPt] = std::asinh(pt / o.ptEps);
  out[kLogTheta] = std::log10(th);
  out[kX] = x;
  out[kY] = y;
  out[kQX] = q * x;
  out[kR] = std::hypot(x, y);
  return true;
}
 
/// Apply `fn(vals, q, fileIdx, entry)` to every selected muon in the dataset.
/// One RDataFrame over all files; DefinePerSample gives the file index, which
/// combined with rdfentry_ makes a unique event key whatever rdfentry_ counts
/// from across a chain boundary.
/// Fallback progress: no total is known, so report throughput rather than a
/// percentage.  Cheap enough to call on every entry.
struct Ticker {
  std::atomic<long long> entries{0};
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  long long lastPrint = 0;
  std::string label;
  void tick(long long muons) {
    const long long n = ++entries;
    if (n - lastPrint < 20000) return;
    lastPrint = n;
    const double dt = std::chrono::duration<double>(
                          std::chrono::steady_clock::now() - t0).count();
    std::fprintf(stderr, "\r[%s] %lld entries, %lld muons, %.0f entries/s   ",
                 label.c_str(), n, muons, dt > 0 ? n / dt : 0.0);
    std::fflush(stderr);
  }
};
 
template <class Fn>
void streamAll(const std::vector<std::string>& files, const LoadOpts& o,
               const std::string& label, bool progress, Fn&& fn) {
  ROOT::RDataFrame df(o.tree, files);
  std::vector<std::string> paths = files;
  // RNode erases the node type, so Range can be applied conditionally without
  // duplicating the whole event loop for the limited and unlimited cases.
  ROOT::RDF::RNode node = df.DefinePerSample(
      "_fileIdx", [paths](unsigned, const ROOT::RDF::RSampleInfo& id) {
        for (size_t i = 0; i < paths.size(); ++i)
          if (id.Contains(paths[i])) return static_cast<int>(i);
        return 0;
      });
  if (o.maxEvents > 0)
    node = node.Range(0u, static_cast<unsigned>(o.maxEvents));
#ifdef MFL_HAVE_RDF_PROGRESS
  // ROOT's own bar knows the total entry count, so it can show a percentage
  // and an ETA.  Only available from 6.30.
  if (progress) ROOT::RDF::Experimental::AddProgressBar(node);
  Ticker* tick = nullptr;
#else
  Ticker tickStorage;
  tickStorage.label = label;
  Ticker* tick = progress ? &tickStorage : nullptr;
#endif
  std::atomic<long long> muonCount{0};
 
  ROOT::RDF::RNode& sel = node;
 
  if (o.source == "plane") {
    const std::string P = o.branch;
    sel.ForeachSlot(
        [&](unsigned s, const RVecD& px, const RVecD& py, const RVecD& pz,
            const RVecD& xx, const RVecD& yy, const RVecI& pdg, const RVecI& trk,
            const RVecD& mcw, int fidx, ULong64_t ent) {
          const uint64_t key = mix((static_cast<uint64_t>(fidx) << 40) ^ ent);
          double vals[kNVar];
          int sgn = 0;
          long long n = 0;
          for (size_t i = 0; i < px.size(); ++i) {
            const int id = trk[i];
            const double w =
                (id >= 0 && static_cast<size_t>(id) < mcw.size()) ? mcw[id] : 1.0;
            if (derive(px[i], py[i], pz[i], xx[i], yy[i], pdg[i], w, o, vals, sgn)) {
              fn(s, vals, sgn, w, fidx, key);
              ++n;
            }
          }
          muonCount += n;
          if (tick) tick->tick(muonCount.load());
        },
        {P + ".fPx", P + ".fPy", P + ".fPz", P + ".fX", P + ".fY", P + ".fPdgCode",
         P + ".fTrackID", "MCTrack.fW", "_fileIdx", "rdfentry_"});
  } else {
    sel.ForeachSlot(
        [&](unsigned s, const RVecD& px, const RVecD& py, const RVecD& pz,
            const RVecD& xx, const RVecD& yy, const RVecI& pdg, const RVecD& mcw,
            int fidx, ULong64_t ent) {
          const uint64_t key = mix((static_cast<uint64_t>(fidx) << 40) ^ ent);
          double vals[kNVar];
          int sgn = 0;
          long long n = 0;
          for (size_t i = 0; i < px.size(); ++i)
            if (derive(px[i], py[i], pz[i], xx[i], yy[i], pdg[i], mcw[i], o, vals, sgn)) {
              fn(s, vals, sgn, mcw[i], fidx, key);
              ++n;
            }
          muonCount += n;
          if (tick) tick->tick(muonCount.load());
        },
        {"MCTrack.fPx", "MCTrack.fPy", "MCTrack.fPz", "MCTrack.fStartX",
         "MCTrack.fStartY", "MCTrack.fPdgCode", "MCTrack.fW", "_fileIdx",
         "rdfentry_"});
  }
  if (progress) std::fprintf(stderr, "\n");
}
 
}  // namespace
 
std::vector<std::string> splitList(const std::string& spec) {
  std::vector<std::string> out;
  if (!spec.empty() && spec[0] == '@') {
    std::ifstream in(spec.substr(1));
    if (!in) throw std::runtime_error("cannot open file list " + spec.substr(1));
    std::string line;
    while (std::getline(in, line)) {
      const std::string t = trim(line);
      if (!t.empty() && t[0] != '#') out.push_back(t);
    }
    return out;
  }
  // A bare directory means every .root inside it.
  if (spec.find(',') == std::string::npos && !spec.empty() &&
      spec.back() != 't') {  // cheap "does not end in .root"
    TSystemDirectory dir(spec.c_str(), spec.c_str());
    TList* l = dir.GetListOfFiles();
    if (l) {
      TIter next(l);
      while (TObject* ob = next()) {
        const std::string nm = ob->GetName();
        if (nm.size() > 5 && nm.substr(nm.size() - 5) == ".root")
          out.push_back(spec + "/" + nm);
      }
      std::sort(out.begin(), out.end());
      if (!out.empty()) return out;
    }
  }
  size_t a = 0;
  while (a <= spec.size()) {
    size_t b = spec.find(',', a);
    if (b == std::string::npos) b = spec.size();
    const std::string t = trim(spec.substr(a, b - a));
    if (!t.empty()) out.push_back(t);
    a = b + 1;
  }
  return out;
}
 
fluxval::HistSpec probeRange(const std::vector<std::string>& fa,
                         const std::vector<std::string>& fb, const LoadOpts& o) {
  std::vector<double> lo(kNVar, 1e300), hi(kNVar, -1e300);
  LoadOpts po = o;
  po.maxEvents = 0;
 
  auto probe = [&](const std::vector<std::string>& files) {
    if (files.empty()) return;
    // A handful of files is enough to fix a range; reading them all would
    // double the cost of the whole job for no benefit.
    std::vector<std::string> sub(files.begin(),
                                 files.begin() + std::min<size_t>(4, files.size()));
    ROOT::RDataFrame df(o.tree, sub);
    auto lim = df.Range(0, 20000);  // Range needs IMT off; the caller ensures it
    std::vector<std::vector<double>> slo, shi;
    const unsigned n = 1;
    slo.assign(n, std::vector<double>(kNVar, 1e300));
    shi.assign(n, std::vector<double>(kNVar, -1e300));
    auto upd = [&](const double* v) {
      for (int j = 0; j < kNVar; ++j) {
        slo[0][j] = std::min(slo[0][j], v[j]);
        shi[0][j] = std::max(shi[0][j], v[j]);
      }
    };
    if (o.source == "plane") {
      const std::string P = o.branch;
      lim.Foreach(
          [&](const RVecD& px, const RVecD& py, const RVecD& pz, const RVecD& xx,
              const RVecD& yy, const RVecI& pdg, const RVecI& trk, const RVecD& mcw) {
            double vals[kNVar];
            int sgn = 0;
            for (size_t i = 0; i < px.size(); ++i) {
              const int id = trk[i];
              const double w =
                  (id >= 0 && static_cast<size_t>(id) < mcw.size()) ? mcw[id] : 1.0;
              if (derive(px[i], py[i], pz[i], xx[i], yy[i], pdg[i], w, po, vals, sgn))
                upd(vals);
            }
          },
          {P + ".fPx", P + ".fPy", P + ".fPz", P + ".fX", P + ".fY", P + ".fPdgCode",
           P + ".fTrackID", "MCTrack.fW"});
    } else {
      lim.Foreach(
          [&](const RVecD& px, const RVecD& py, const RVecD& pz, const RVecD& xx,
              const RVecD& yy, const RVecI& pdg, const RVecD& mcw) {
            double vals[kNVar];
            int sgn = 0;
            for (size_t i = 0; i < px.size(); ++i)
              if (derive(px[i], py[i], pz[i], xx[i], yy[i], pdg[i], mcw[i], po, vals,
                         sgn))
                upd(vals);
          },
          {"MCTrack.fPx", "MCTrack.fPy", "MCTrack.fPz", "MCTrack.fStartX",
           "MCTrack.fStartY", "MCTrack.fPdgCode", "MCTrack.fW"});
    }
    for (int j = 0; j < kNVar; ++j) {
      lo[j] = std::min(lo[j], slo[0][j]);
      hi[j] = std::max(hi[j], shi[0][j]);
    }
  };
  probe(fa);
  probe(fb);
 
  fluxval::HistSpec sp;
  sp.nvar = kNVar;
  sp.nfine = o.nfine;
  sp.lo.resize(kNVar);
  sp.hi.resize(kNVar);
  for (int j = 0; j < kNVar; ++j) {
    if (!(lo[j] < hi[j])) { lo[j] = -1.0; hi[j] = 1.0; }
    // The probe sees a fraction of the data, so the true tails run past what
    // it found.  Pad generously; an over-wide range costs resolution, an
    // under-wide one piles the tail into the end bin, which is worse.
    const double pad = 0.35 * (hi[j] - lo[j]) + 1e-6;
    sp.lo[j] = lo[j] - pad;
    sp.hi[j] = hi[j] + pad;
  }
  // Momentum bounds are known exactly from the selection.
  sp.lo[kLogP] = std::log10(o.pmin) - 1e-6;
  sp.hi[kLogP] = std::log10(o.pmax) + 1e-6;
  sp.lo[kLogPz] = std::log10(o.pmin) - 1.0;
  sp.hi[kLogPz] = std::log10(o.pmax) + 1e-6;
  sp.lo[kAsinhPt] = 0.0;
  sp.hi[kAsinhPt] = std::asinh(o.pmax / o.ptEps) + 1e-6;
  sp.lo[kR] = 0.0;
  return sp;
}
 
fluxval::UnitHists loadHists(const std::vector<std::string>& files_in, const LoadOpts& o,
                         const fluxval::HistSpec& spec, const std::string& label) {
  if (files_in.empty()) throw std::runtime_error("no input files for " + label);
  std::vector<std::string> files = files_in;
  if (o.maxFiles > 0 && static_cast<int>(files.size()) > o.maxFiles)
    files.resize(o.maxFiles);
 
  const int nFiles = static_cast<int>(files.size());
  const bool byFile = nFiles >= o.minFiles;
  const int nUnits = byFile ? std::min(nFiles, o.maxUnits) : o.nBlocks;
 
  const unsigned nSlots = std::max(1u, ROOT::GetThreadPoolSize());
  std::vector<fluxval::UnitHists> per(nSlots);
  for (auto& u : per) u.alloc(spec, nUnits);
 
  logf("[io] %s: reading %d file(s)...", label.c_str(), nFiles);
  streamAll(files, o, label, o.progress,
            [&](unsigned s, const double* vals, int sgn, double w, int fidx,
                uint64_t key) {
    // Units: whole files when there are enough of them (files are the natural
    // correlated block in a production).  Files are grouped when there are
    // more than maxUnits, which keeps whole files inside one unit.
    const int unit = byFile ? (static_cast<long long>(fidx) * nUnits / nFiles)
                            : static_cast<int>(key % static_cast<uint64_t>(nUnits));
    per[s < nSlots ? s : 0].fillSign(unit, sgn, vals, w);
  });
 
  fluxval::UnitHists H;
  H.alloc(spec, nUnits);
  for (const auto& u : per) H.merge(u);
  H.label = label;
  H.unitsAreFiles = byFile;
 
  logf("[io] %-10s %12lld muons  %4d file(s)  %3d unit(s) (%s)  sum(w) = %.6g  "
       "N_eff = %.0f  muons/N_eff = %.0f",
       label.c_str(), H.rows, nFiles, nUnits, byFile ? "files" : "event blocks",
       H.totalW(), H.nEffTotal(), H.nEffTotal() > 0 ? H.rows / H.nEffTotal() : 0.0);
  if (H.clamped)
    logf("[io] %-10s WARNING %lld muons (%.3f%%) fell outside the histogram range and "
         "were piled into the end bins -- widen it with --nfine/--pmax or check the "
         "probe",
         label.c_str(), H.clamped, 100.0 * H.clamped / std::max(1LL, H.rows));
  if (H.rows == 0) throw std::runtime_error("sample " + label + " is empty");
  return H;
}
 
}  // namespace fluxval

