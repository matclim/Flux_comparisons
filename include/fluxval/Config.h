// Config.h -- the observables and the phase-space regions.
//
// Everything here is defined on the muon as it crosses the reference plane.
// No production-level or truth information is used anywhere: the comparison is
// between the muons the two productions deliver, not between the processes
// that made them.

#pragma once

#include "fluxval/Core.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace fluxval {

/// Variable indices.  Keep in sync with kVars and with fillRow() in Load.cxx.
/// The particle family under test.  Everything that depends on the species --
/// which PDG codes are accepted, whether q*x is meaningful, and what the two
/// sign classes are called -- is carried here rather than hard-coded, so that
/// muons and neutrinos are two configurations of one tool.
struct Species {
  std::vector<int> pdgAbs{13};
  bool neutral = false;
  std::string name = "muon";
  std::string posLabel = "mu+";   // sign index 0: PDG < 0, the antilepton
  std::string negLabel = "mu-";   // sign index 1: PDG > 0
};

/// Electric charge from a PDG code, for the species this tool handles.
/// Charged leptons: PDG 13 is mu-, so the charge is minus the sign of the code.
/// Neutrinos (12, 14, 16) are neutral.
inline int chargeOfPdg(int pdg) {
  const int a = std::abs(pdg);
  if (a == 12 || a == 14 || a == 16) return 0;
  return (pdg > 0) ? -1 : +1;
}

/// Sign class index: 0 for PDG < 0 (mu+, nubar), 1 for PDG > 0 (mu-, nu).
/// Defined on the PDG code rather than the charge so that it survives the
/// neutral case, where the charge carries no information.
inline int signIndexOfPdg(int pdg) { return (pdg < 0) ? 0 : 1; }

inline Species makeSpecies(const std::vector<int>& pdgAbs) {
  Species sp;
  sp.pdgAbs = pdgAbs.empty() ? std::vector<int>{13} : pdgAbs;
  sp.neutral = true;
  for (int a : sp.pdgAbs) sp.neutral = sp.neutral && (chargeOfPdg(a) == 0);
  bool allCharged = true;
  for (int a : sp.pdgAbs) allCharged = allCharged && (chargeOfPdg(a) != 0);
  if (!sp.neutral && !allCharged)
    throw std::runtime_error(
        "--pdg mixes charged and neutral species; the sign classes and the q*x "
        "observable would be inconsistent. Run them separately.");

  if (sp.neutral) {
    sp.name = (sp.pdgAbs.size() == 1)
                  ? (sp.pdgAbs[0] == 12 ? "nu_e"
                                        : (sp.pdgAbs[0] == 14 ? "nu_mu" : "nu_tau"))
                  : "neutrino";
    sp.posLabel = "nubar";
    sp.negLabel = "nu";
  } else {
    sp.name = "muon";
    sp.posLabel = "mu+";
    sp.negLabel = "mu-";
  }
  return sp;
}

enum VarIdx {
  kLogP = 0,   // log10(p / GeV)
  kLogPz,      // log10(pz / GeV)
  kAsinhPt,    // asinh(pT / 0.05 GeV)   -- linear near 0, log-like in the tail
  kLogTheta,   // log10(theta / mrad)
  kX,          // x [cm]
  kY,          // y [cm]
  kQX,         // q * x [cm]  -- the bending variable; unsigned x throws the sign away
  kR,          // sqrt(x^2 + y^2) [cm]
  kNVar
};

struct VarDef {
  const char* name;
  const char* axis;
  bool logyPlot;
};

inline const std::vector<VarDef>& varDefs() {
  static const std::vector<VarDef> v = {
      {"log10_p", "log_{10}(p / GeV)", true},
      {"log10_pz", "log_{10}(p_{z} / GeV)", true},
      {"asinh_pT", "asinh(p_{T} / 0.05 GeV)", true},
      {"log10_theta", "log_{10}(#theta / mrad)", true},
      {"x", "x  [cm]", true},
      {"y", "y  [cm]", true},
      {"qx", "q #upoint x  [cm]", true},
      {"r", "r = #sqrt{x^{2}+y^{2}}  [cm]", true},
  };
  return v;
}

/// Default regions.  All cuts are expressed in the transformed variable, so
/// the physical value is given in the name for readability.
///
/// These are the "relative amounts" of the flux: each is a fraction of the
/// sample's own total weight, hence immune to the deliberate rate difference.
/// q*x is identically zero for a neutral species: the binning degenerates and
/// A^2 is undefined, so it must be removed rather than merely ignored.
inline std::vector<int> defaultVars(const Species& sp) {
  std::vector<int> v;
  for (int j = 0; j < kNVar; ++j) {
    if (sp.neutral && j == kQX) continue;
    v.push_back(j);
  }
  return v;
}

inline std::vector<fluxval::Region> defaultRegions(const Species& sp) {
  auto L = [](double gev) { return std::log10(gev); };
  std::vector<fluxval::Region> r = {
      {"p > 50 GeV", kLogP, L(50.0), 1e300, 0},
      {"p > 100 GeV", kLogP, L(100.0), 1e300, 0},
      {"p > 200 GeV", kLogP, L(200.0), 1e300, 0},
      {"p > 300 GeV", kLogP, L(300.0), 1e300, 0},
      {"pT > 1 GeV", kAsinhPt, std::asinh(1.0 / 0.05), 1e300, 0},
      {"pT > 2 GeV", kAsinhPt, std::asinh(2.0 / 0.05), 1e300, 0},
      {"theta > 20 mrad", kLogTheta, L(20.0), 1e300, 0},
      {"theta > 50 mrad", kLogTheta, L(50.0), 1e300, 0},
      {"r > 100 cm", kR, 100.0, 1e300, 0},
      {"r > 200 cm", kR, 200.0, 1e300, 0},
      {sp.posLabel, kLogP, -1e300, 1e300, +1},
      {sp.negLabel, kLogP, -1e300, 1e300, -1},
      {sp.posLabel + ", p > 100", kLogP, L(100.0), 1e300, +1},
      {sp.negLabel + ", p > 100", kLogP, L(100.0), 1e300, -1},
  };
  if (!sp.neutral)
    r.insert(r.begin() + 10, {"qx > 0", kQX, 0.0, 1e300, 0});
  return r;
}

/// Quantities reported for each sample and compared between them.
/// Column 0 of the accumulator is the sample total, so {name, k, 0} is the
/// fraction of the flux in region k-1.
/// Quantities are resolved against the region list BY NAME.  Hand-maintained
/// numeric indices break silently the moment a region is inserted, which is a
/// poor property for shared code.
inline std::vector<fluxval::Quantity> defaultQuantities(
    const Species& sp, const std::vector<fluxval::Region>& regs) {
  auto idx = [&regs](const std::string& nm) {
    for (std::size_t i = 0; i < regs.size(); ++i)
      if (regs[i].name == nm) return static_cast<int>(i) + 1;  // column 0 = total
    throw std::runtime_error("no region named '" + nm + "'");
  };
  const std::string P = sp.posLabel, N = sp.negLabel;

  std::vector<fluxval::Quantity> q;
  for (const char* nm : {"p > 50 GeV", "p > 100 GeV", "p > 200 GeV", "p > 300 GeV",
                         "pT > 1 GeV", "pT > 2 GeV", "theta > 20 mrad",
                         "theta > 50 mrad", "r > 100 cm", "r > 200 cm"})
    q.push_back({std::string("f(") + nm + ")", idx(nm), 0});
  if (!sp.neutral) q.push_back({"f(qx > 0)", idx("qx > 0"), 0});
  q.push_back({P + " / " + N, idx(P), idx(N)});
  q.push_back({P + "/" + N + " , p>100", idx(P + ", p > 100"), idx(N + ", p > 100")});
  q.push_back({"f(p>100 | " + P + ")", idx(P + ", p > 100"), idx(P)});
  q.push_back({"f(p>100 | " + N + ")", idx(N + ", p > 100"), idx(N)});
  return q;
}

}  // namespace fluxval
