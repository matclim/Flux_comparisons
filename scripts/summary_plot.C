// summary_plot.C -- one picture of what differs and by how much, built from a
// run that has already happened.  No rerun, no ROOT file needed: everything is
// read back out of the tool's own text log.
//
//     root -l -b -q 'summary_plot.C("out9_full/muonflux.txt","out9_full/summary")'
//
// WHY THIS PLOT EXISTS
// --------------------
// At 2.3e9 muons every p-value sits at the permutation floor and every pull is
// enormous, so significance stops discriminating between results: a 0.5%
// difference and a 2% difference both come out "certain".  What still carries
// information is the SIZE of each difference, and the RANKING of A^2 across
// variables.  This macro plots exactly those two things:
//
//   page 1  importance of each deviation beside its significance, one row per
//           quantity.  A row near zero on the left but far right on the right
//           is significant and unimportant; a row far right on the left but
//           low on the right is important and unmeasured.  Both are common
//           here and both are invisible in a table of p-values.
//   page 2  A^2 per variable and per charge, log scale.  With every p at the
//           floor, only the relative magnitudes rank the observables, and A^2
//           scales with statistics so the ranking is meaningful.
 
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
 
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
 
namespace mfsum {
 
struct Frac {
  std::string name;
  double ref = 0, eref = 0, test = 0, etest = 0, rel = 0, pull = 0;
  double relErr() const {
    return (ref != 0.0) ? 100.0 * std::hypot(eref, etest) / std::fabs(ref) : 0.0;
  }
};
 
struct Shape {
  std::string var, charge;
  double a2 = 0, z = 0;
};
 
inline bool isNumber(const std::string& s) {
  if (s.empty()) return false;
  std::string t = s;
  if (t[0] == '<') t = t.substr(1);
  if (t.empty()) return false;
  char* end = nullptr;
  std::strtod(t.c_str(), &end);
  return end && *end == '\0' && end != t.c_str();
}
 
inline double toNum(const std::string& s) {
  std::string t = s;
  if (!t.empty() && t[0] == '<') t = t.substr(1);
  if (!t.empty() && t.back() == '%') t.pop_back();
  return std::atof(t.c_str());
}
 
inline bool splitPM(const std::string& s, double& v, double& e) {
  const size_t p = s.find("+-");
  if (p == std::string::npos) return false;
  const std::string a = s.substr(0, p), b = s.substr(p + 2);
  if (!isNumber(a) || !isNumber(b)) return false;
  v = std::atof(a.c_str());
  e = std::atof(b.c_str());
  return true;
}
 
inline std::vector<std::string> tokens(const std::string& line) {
  std::vector<std::string> t;
  std::istringstream is(line);
  std::string w;
  while (is >> w) t.push_back(w);
  return t;
}
 
struct Parsed {
  std::vector<Frac> fracs;
  std::vector<Shape> shapes;
  std::string labRef = "reference", labTest = "test";
};
 
// The log is console output, so it carries whatever mess the run left: stderr
// progress lines merged in by `tee`, rows whose label a carriage return ate,
// and "--" rows for empty regions.  Anything that does not parse cleanly is
// skipped rather than guessed at.
inline Parsed parseLog(const std::string& path) {
  Parsed out;
  std::ifstream in(path.c_str());
  if (!in) {
    ::Error("summary_plot", "cannot open %s", path.c_str());
    return out;
  }
  std::string line, charge = "both";
  while (std::getline(in, line)) {
    const size_t sh = line.find("1-D shapes [");
    if (sh != std::string::npos) {
      const size_t a = line.find('[', sh), b = line.find(']', a);
      if (a != std::string::npos && b != std::string::npos)
        charge = line.substr(a + 1, b - a - 1);
      continue;
    }
    if (line.compare(0, 12, "reference : ") == 0)
      out.labRef = line.substr(12);
    if (line.compare(0, 12, "test      : ") == 0)
      out.labTest = line.substr(12);
 
    std::vector<std::string> t = tokens(line);
    if (t.size() < 4) continue;
    if (t.back() == "<<<") t.pop_back();
 
    int pm = -1;
    for (size_t i = 0; i + 1 < t.size(); ++i)
      if (t[i].find("+-") != std::string::npos &&
          t[i + 1].find("+-") != std::string::npos) {
        pm = (int)i;
        break;
      }
    if (pm > 0 && (int)t.size() >= pm + 4) {
      Frac f;
      if (splitPM(t[pm], f.ref, f.eref) && splitPM(t[pm + 1], f.test, f.etest) &&
          t[pm + 2] != "--" && isNumber(t[pm + 3])) {
        f.rel = toNum(t[pm + 2]);
        f.pull = toNum(t[pm + 3]);
        for (int k = 0; k < pm; ++k) f.name += (k ? " " : "") + t[k];
        if (!f.name.empty()) out.fracs.push_back(f);
      }
      continue;
    }
 
    if (t.size() >= 7 && !isNumber(t[0]) && isNumber(t[1]) && isNumber(t[2]) &&
        isNumber(t[3]) && isNumber(t[4])) {
      if (t[0] == "variable" || t[0][0] == '-') continue;
      Shape s;
      s.var = t[0];
      s.charge = charge;
      s.a2 = toNum(t[1]);
      s.z = toNum(t[3]);
      if (s.a2 > 0) out.shapes.push_back(s);
    }
  }
  return out;
}
 
}  // namespace mfsum
 
// ---------------------------------------------------------------------------
 
void summary_plot(const char* logPath = "outputs/muonflux.txt",
                  const char* outPrefix = "outputs/summary",
                  const char* labRef = "", const char* labTest = "",
                  double sigLine = 3.0) {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPalette(kCool);
 
  mfsum::Parsed P = mfsum::parseLog(logPath);
  // The sample names are drawn onto the PDF summary page but never written to
  // the log, so they cannot be recovered from it; take them as arguments.
  if (labRef && *labRef) P.labRef = labRef;
  if (labTest && *labTest) P.labTest = labTest;
  if (P.fracs.empty() && P.shapes.empty()) {
    ::Error("summary_plot", "nothing parsed from %s", logPath);
    return;
  }
  ::Info("summary_plot", "parsed %zu fractions, %zu shape rows", P.fracs.size(),
         P.shapes.size());
 
  const int colBig = kOrange + 8, colSmall = kAzure + 2;
 
  // ============ page 1: effect size beside significance ============
  //
  // One row per quantity, so nothing overlaps.  Left panel: the importance of
  // the deviation, i.e. its size relative to the reference, with the bootstrap
  // uncertainty.  Right panel: its significance, on a log axis.  The two
  // panels share the row ordering, so a row near zero on the left and far
  // right on the right is a deviation that is significant and unimportant --
  // the dominant case at full statistics, and the thing a table of p-values
  // hides.
  //
  // Colour encodes significance only.  There is deliberately no per-sample
  // colour: every point is a comparison of both samples, not a property of one.
  {
    std::vector<mfsum::Frac> f = P.fracs;
    // Most significant at the top.
    for (size_t i = 0; i + 1 < f.size(); ++i)
      for (size_t j = i + 1; j < f.size(); ++j)
        if (std::fabs(f[j].pull) > std::fabs(f[i].pull)) std::swap(f[i], f[j]);
    const int n = (int)f.size();
    if (n == 0) return;
 
    TCanvas* c = new TCanvas("cv", "summary", 1100, 40 + 34 * n);
    TPad* pL = new TPad("pL", "", 0.00, 0.00, 0.66, 1.00);
    TPad* pR = new TPad("pR", "", 0.66, 0.00, 1.00, 1.00);
    pL->SetLeftMargin(0.34);
    pL->SetRightMargin(0.02);
    pL->SetTopMargin(0.13);
    pL->SetBottomMargin(0.13);
    pR->SetLeftMargin(0.02);
    pR->SetRightMargin(0.06);
    pR->SetTopMargin(0.13);
    pR->SetBottomMargin(0.13);
    pR->SetLogx();
    pL->Draw();
    pR->Draw();
 
    double xmax = 1.0, ymax = sigLine;
    for (int i = 0; i < n; ++i) {
      xmax = std::max(xmax, std::fabs(f[i].rel) + f[i].relErr());
      ymax = std::max(ymax, std::fabs(f[i].pull));
    }
    xmax *= 1.15;
 
    // ---- left: effect size ----
    pL->cd();
    TH2F* frL = new TH2F("frL", "", 100, -xmax, xmax, n, 0, n);
    frL->SetDirectory(0);
    for (int i = 0; i < n; ++i)
      frL->GetYaxis()->SetBinLabel(n - i, f[i].name.c_str());
    frL->GetYaxis()->SetLabelSize(0.030);
    frL->GetXaxis()->SetTitle("(test - ref) / ref   [%]");
    frL->GetXaxis()->SetTitleSize(0.036);
    frL->GetXaxis()->SetLabelSize(0.030);
    frL->Draw();
 
    TGraphErrors* gL = new TGraphErrors();
    for (int i = 0; i < n; ++i) {
      gL->SetPoint(i, f[i].rel, n - i - 0.5);
      gL->SetPointError(i, f[i].relErr(), 0.0);
    }
    gL->SetMarkerStyle(20);
    gL->SetMarkerSize(1.0);
    gL->SetLineWidth(2);
    gL->SetMarkerColor(kAzure + 2);
    gL->SetLineColor(kAzure + 2);
    gL->Draw("P SAME");
    TLine* zl = new TLine(0.0, 0.0, 0.0, n);
    zl->SetLineStyle(2);
    zl->SetLineColor(kGray + 2);
    zl->Draw();
 
    TLatex hL;
    hL.SetNDC();
    hL.SetTextAlign(12);
    hL.SetTextSize(0.038);
    hL.DrawLatex(0.02, 0.965, "Importance of the deviation");
    hL.SetTextSize(0.026);
    hL.SetTextColor(kGray + 3);
    hL.DrawLatex(0.34, 0.925, Form("ref: %s      test: %s", P.labRef.c_str(),
                                   P.labTest.c_str()));
 
    // ---- right: significance ----
    pR->cd();
    TH2F* frR = new TH2F("frR", "", 100, 0.03, ymax * 2.5, n, 0, n);
    frR->SetDirectory(0);
    frR->GetYaxis()->SetLabelSize(0.0);
    frR->GetXaxis()->SetTitle("Significance   [#sigma]");
    frR->GetXaxis()->SetTitleSize(0.036);
    frR->GetXaxis()->SetLabelSize(0.030);
    frR->Draw();
 
    TGraphErrors* gA = new TGraphErrors();
    TGraphErrors* gB = new TGraphErrors();
    int na = 0, nb = 0;
    for (int i = 0; i < n; ++i) {
      const double y = std::max(0.035, std::fabs(f[i].pull));
      if (std::fabs(f[i].pull) > sigLine) gB->SetPoint(nb++, y, n - i - 0.5);
      else gA->SetPoint(na++, y, n - i - 0.5);
    }
    gA->SetMarkerStyle(20);
    gA->SetMarkerSize(1.0);
    gA->SetMarkerColor(kAzure + 2);
    if (gA->GetN()) gA->Draw("P SAME");
    gB->SetMarkerStyle(20);
    gB->SetMarkerSize(1.0);
    gB->SetMarkerColor(kRed + 1);
    if (gB->GetN()) gB->Draw("P SAME");
 
    TLine* sl = new TLine(sigLine, 0.0, sigLine, n);
    sl->SetLineStyle(2);
    sl->SetLineColor(kRed + 1);
    sl->SetLineWidth(2);
    sl->Draw();
 
    TLatex hR;
    hR.SetNDC();
    hR.SetTextAlign(12);
    hR.SetTextSize(0.038);
    hR.DrawLatex(0.02, 0.965, "Significance of the deviation");
    hR.SetTextSize(0.026);
    hR.SetTextColor(kRed + 1);
    hR.DrawLatex(0.02, 0.925, Form("red: significance > %.0f#sigma", sigLine));
 
    c->Print(Form("%s_effectsize.pdf", outPrefix));
  }
 
  // ============ page 2: A^2 ranking ============
  if (!P.shapes.empty()) {
    std::vector<std::string> vars, charges;
    for (size_t i = 0; i < P.shapes.size(); ++i) {
      bool fv = false, fc = false;
      for (size_t k = 0; k < vars.size(); ++k) fv |= (vars[k] == P.shapes[i].var);
      for (size_t k = 0; k < charges.size(); ++k)
        fc |= (charges[k] == P.shapes[i].charge);
      if (!fv) vars.push_back(P.shapes[i].var);
      if (!fc) charges.push_back(P.shapes[i].charge);
    }
 
    TCanvas* c = new TCanvas("cr", "ranking", 950, 780);
    c->SetLeftMargin(0.22);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.14);
    c->SetBottomMargin(0.13);
    c->SetLogx();
 
    double amax = 1.0, amin = 1e30;
    for (size_t i = 0; i < P.shapes.size(); ++i) {
      amax = std::max(amax, P.shapes[i].a2);
      amin = std::min(amin, P.shapes[i].a2);
    }
    const int nv = (int)vars.size();
    TH2F* fr = new TH2F("fr2", "", 100, std::max(0.05, amin * 0.4), amax * 3.0, nv, 0,
                        nv);
    fr->SetDirectory(0);
    for (int i = 0; i < nv; ++i) fr->GetYaxis()->SetBinLabel(i + 1, vars[i].c_str());
    fr->GetYaxis()->SetLabelSize(0.038);
    fr->GetXaxis()->SetTitle("A^{2}   (larger = the distributions differ more)");
    fr->GetXaxis()->SetTitleSize(0.040);
    fr->Draw();
 
    const int mstyle[3] = {20, 22, 23};
    const int mcol[3] = {kAzure + 2, kOrange + 8, kSpring - 6};
    for (size_t ci = 0; ci < charges.size() && ci < 3; ++ci) {
      TGraphErrors* g = new TGraphErrors();
      int n = 0;
      for (int vi = 0; vi < nv; ++vi)
        for (size_t i = 0; i < P.shapes.size(); ++i)
          if (P.shapes[i].var == vars[vi] && P.shapes[i].charge == charges[ci]) {
            const double off = (charges.size() > 1)
                                   ? 0.28 + 0.44 * ci / (double)(charges.size() - 1)
                                   : 0.5;
            g->SetPoint(n++, P.shapes[i].a2, vi + off);
          }
      g->SetMarkerStyle(mstyle[ci]);
      g->SetMarkerSize(1.25);
      g->SetMarkerColor(mcol[ci]);
      g->SetLineColor(mcol[ci]);
      if (g->GetN()) g->Draw("P SAME");
 
      TText* lg = new TText(0.0, 0.0, charges[ci].c_str());
      lg->SetNDC();
      lg->SetX(0.70);
      lg->SetY(0.055 + 0.030 * (charges.size() - 1 - ci));
      lg->SetTextSize(0.026);
      lg->SetTextColor(mcol[ci]);
      lg->SetTextAlign(12);
      lg->Draw();
    }
 
    TLatex h;
    h.SetNDC();
    h.SetTextAlign(12);
    h.SetTextSize(0.034);
    h.DrawLatex(0.02, 0.975, "shape difference by variable and charge");
    h.SetTextSize(0.022);
    h.SetTextColor(kGray + 3);
    h.DrawLatex(0.02, 0.940,
                "log scale.  When every p-value is at the permutation floor, only "
                "these relative sizes rank the observables.");
    h.DrawLatex(0.02, 0.910,
                "theta and r follow from p geometrically; pT does not.");
    c->Print(Form("%s_a2ranking.pdf", outPrefix));
  }
 
  ::Info("summary_plot", "wrote %s_effectsize.pdf and %s_a2ranking.pdf", outPrefix,
         outPrefix);
}
 

