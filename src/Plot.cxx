#include "fluxval/Plot.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TString.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TText.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace fluxval {

namespace {
// Not const: setStyle() replaces these with colours drawn from the palette.
int kColA = kAzure + 2;
int kColB = kOrange + 8;

// Legends are rendered through TLatex, which turns "_" in a filename into a
// subscript.  TText does no markup interpretation, so labels survive intact.
void putLabel(double x, double y, const std::string& s, int col, double size = 0.030) {
  TText* t = new TText(x, y, s.c_str());
  t->SetNDC();
  t->SetTextColor(col);
  t->SetTextSize(size);
  t->SetTextAlign(12);
  t->Draw();
}
}  // namespace

void setStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetEndErrorSize(0);
  gStyle->SetPalette(kCool);

  // Draw the two sample colours from the ends of the active palette rather
  // than hard-coding them, so the overlays match anything else on the page.
  // 12% and 88% instead of 0% and 100%: the extreme ends of kCool are very
  // dark and very pale, and neither reads well as a line colour on white.
  const int n = TColor::GetNumberOfColors();
  if (n > 8) {
    kColA = TColor::GetColorPalette(static_cast<int>(0.12 * n));
    kColB = TColor::GetColorPalette(static_cast<int>(0.88 * n));
  }
}

int colourA() { return kColA; }
int colourB() { return kColB; }

TH1D* densityHist(const char* name, const std::vector<double>& edges,
                  const std::vector<double>& shape, const std::vector<double>& err) {
  const int nb = static_cast<int>(edges.size()) - 1;
  TH1D* h = new TH1D(name, "", nb, edges.data());
  h->SetDirectory(nullptr);
  for (int b = 0; b < nb; ++b) {
    // Bins hold equal pooled weight by construction, so the raw content is flat
    // and uninformative.  Divide by the bin width to plot a density, which is
    // the thing that actually has a shape.
    const double wdt = edges[b + 1] - edges[b];
    h->SetBinContent(b + 1, wdt > 0 ? shape[b] / wdt : 0.0);
    h->SetBinError(b + 1, wdt > 0 ? err[b] / wdt : 0.0);
  }
  return h;
}

TH1D* ratioHist(const char* name, const std::vector<double>& edges,
                const std::vector<double>& ratio, const std::vector<double>& err) {
  const int nb = static_cast<int>(edges.size()) - 1;
  TH1D* h = new TH1D(name, "", nb, edges.data());
  h->SetDirectory(nullptr);
  for (int b = 0; b < nb; ++b) {
    if (!std::isfinite(ratio[b])) continue;
    h->SetBinContent(b + 1, ratio[b]);
    h->SetBinError(b + 1, err[b]);
  }
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
  h->SetLineColor(kBlack);
  return h;
}

void drawWithRatio(TCanvas* c, TH1D* hA, TH1D* hB, TH1D* hR, const std::string& xtitle,
                   bool logy, const std::string& labA, const std::string& labB,
                   const std::string& annotation, double ratioSpan) {
  c->Clear();
  c->cd();
  TPad* p1 = new TPad("p1", "", 0, 0.32, 1, 1);
  TPad* p2 = new TPad("p2", "", 0, 0, 1, 0.32);
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(0.13);
  p1->SetTopMargin(0.07);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.32);
  p2->SetLeftMargin(0.13);
  p1->Draw();
  p2->Draw();

  p1->cd();
  if (logy) p1->SetLogy();
  hA->SetLineColor(kColA);
  hB->SetLineColor(kColB);
  hA->SetLineWidth(2);
  hB->SetLineWidth(2);
  hA->SetStats(0);
  hB->SetStats(0);
  hA->GetXaxis()->SetLabelSize(0);
  hA->GetYaxis()->SetTitle("normalised density");
  hA->GetYaxis()->SetTitleSize(0.050);
  hA->GetYaxis()->SetTitleOffset(1.25);

  double top = std::max(hA->GetMaximum(), hB->GetMaximum());
  double bot = 1e300;
  for (int b = 1; b <= hA->GetNbinsX(); ++b) {
    if (hA->GetBinContent(b) > 0) bot = std::min(bot, hA->GetBinContent(b));
    if (hB->GetBinContent(b) > 0) bot = std::min(bot, hB->GetBinContent(b));
  }
  if (!(bot < 1e299)) bot = 1e-8;
  hA->SetMaximum(logy ? top * 30.0 : top * 1.45);
  hA->SetMinimum(logy ? bot * 0.4 : 0.0);
  hA->Draw("HIST E");
  hB->Draw("HIST E SAME");
  putLabel(0.17, 0.885, labA, kColA);
  putLabel(0.17, 0.845, labB, kColB);
  if (!annotation.empty()) {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.037);
    tx.DrawLatex(0.13, 0.955, annotation.c_str());
  }

  p2->cd();
  hR->GetXaxis()->SetTitle(xtitle.c_str());
  hR->GetYaxis()->SetTitle("test / ref  (shape)");
  hR->GetYaxis()->SetNdivisions(505);
  hR->GetYaxis()->SetTitleSize(0.105);
  hR->GetYaxis()->SetTitleOffset(0.48);
  hR->GetYaxis()->SetLabelSize(0.095);
  hR->GetXaxis()->SetTitleSize(0.115);
  hR->GetXaxis()->SetTitleOffset(1.15);
  hR->GetXaxis()->SetLabelSize(0.095);

  double span = ratioSpan;
  if (span <= 0.0) {  // auto: cover the points and their errors, then round up
    double m = 0.0;
    for (int b = 1; b <= hR->GetNbinsX(); ++b) {
      if (hR->GetBinContent(b) == 0.0) continue;
      m = std::max(m, std::fabs(hR->GetBinContent(b) - 1.0) + hR->GetBinError(b));
    }
    span = std::min(2.0, std::max(0.15, 1.15 * m));
  }
  hR->SetMinimum(1.0 - span);
  hR->SetMaximum(1.0 + span);
  hR->Draw("E1");
  TLine* l = new TLine(hR->GetXaxis()->GetXmin(), 1.0, hR->GetXaxis()->GetXmax(), 1.0);
  l->SetLineStyle(2);
  l->SetLineColor(kGray + 2);
  l->Draw();
  c->Update();
}

void drawNull(TCanvas* c, const fluxval::ShapeResult& R) {
  c->Clear();
  c->cd();
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.09);
  if (R.nullAD.empty()) return;

  double hi = R.obs.AD;
  for (double x : R.nullAD) hi = std::max(hi, x);
  TH1D* h = new TH1D(("null_" + R.var).c_str(), "", 60, 0.0, hi * 1.15);
  h->SetDirectory(nullptr);
  for (double x : R.nullAD) h->Fill(x);
  h->GetXaxis()->SetTitle("Anderson-Darling A^{2}  (permutation null)");
  h->GetYaxis()->SetTitle("permutations");
  h->SetLineColor(kGray + 2);
  h->SetFillColorAlpha(kGray, 0.45);
  h->SetLineWidth(2);
  h->Draw("HIST");

  TLine* l = new TLine(R.obs.AD, 0.0, R.obs.AD, h->GetMaximum() * 0.92);
  l->SetLineColor(kRed + 1);
  l->SetLineWidth(3);
  l->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.036);
  tx.SetTextAlign(12);
  tx.DrawLatex(0.13, 0.955, Form("%s:  A^{2}_{obs} = %.3g,  p %s %.3g", R.var.c_str(),
                                 R.obs.AD,
                                 (R.pAD <= 1.0 / (R.nperm + 1.0) + 1e-12) ? "<" : "=",
                                 R.pAD));
  tx.SetTextSize(0.026);
  tx.SetTextAlign(32);  // right-aligned, so it can never run into the title
  tx.DrawLatex(0.73, 0.81, Form("N_{eff}^{ref} = %.0f", R.neffRef));
  tx.DrawLatex(0.73, 0.772, Form("N_{eff}^{test} = %.0f", R.neffTest));
  tx.DrawLatex(0.73, 0.734, Form("%d permutations, %d units", R.nperm, R.nUnits));
  tx.SetTextAlign(12);
  c->Update();
}

void drawNeff(TCanvas* c, const std::vector<double>& edges,
              const std::vector<double>& nA, const std::vector<double>& nB,
              const std::string& xtitle, const std::string& labA,
              const std::string& labB) {
  c->Clear();
  c->cd();
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.09);
  c->SetLogy();
  const int nb = static_cast<int>(edges.size()) - 1;
  TH1D* ha = new TH1D("neffA", "", nb, edges.data());
  TH1D* hb = new TH1D("neffB", "", nb, edges.data());
  ha->SetDirectory(nullptr);
  hb->SetDirectory(nullptr);
  for (int b = 0; b < nb; ++b) {
    ha->SetBinContent(b + 1, std::max(0.5, nA[b]));
    hb->SetBinContent(b + 1, std::max(0.5, nB[b]));
  }
  ha->SetLineColor(kColA);
  hb->SetLineColor(kColB);
  ha->SetLineWidth(2);
  hb->SetLineWidth(2);
  ha->GetXaxis()->SetTitle(xtitle.c_str());
  ha->GetYaxis()->SetTitle("N_{eff} per bin");
  ha->SetMinimum(0.5);
  ha->SetMaximum(std::max(ha->GetMaximum(), hb->GetMaximum()) * 30.0);
  ha->Draw("HIST");
  hb->Draw("HIST SAME");

  TLine* l = new TLine(edges.front(), 25.0, edges.back(), 25.0);
  l->SetLineStyle(2);
  l->SetLineColor(kRed + 1);
  l->Draw();
  putLabel(0.17, 0.875, labA, kColA);
  putLabel(0.17, 0.835, labB, kColB);
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.030);
  tx.SetTextColor(kRed + 1);
  tx.DrawLatex(0.62, 0.20, "N_{eff} = 25  (20% per bin)");
  tx.SetTextColor(kBlack);
  tx.SetTextSize(0.036);
  tx.DrawLatex(0.13, 0.955, "effective statistics per bin -- the ceiling on any test");
  c->Update();
}

void drawQuantityPulls(TCanvas* c, const std::vector<fluxval::QuantityResult>& q,
                       const std::string& labA, const std::string& labB) {
  c->Clear();
  c->cd();
  c->SetLeftMargin(0.30);
  c->SetTopMargin(0.13);
  c->SetRightMargin(0.16);
  c->SetBottomMargin(0.11);
  const int n = static_cast<int>(q.size());
  if (n == 0) return;

  // Range from the rows that carry information only.  A single unmeasurable
  // row spanning +-60% would otherwise squash every real result onto the zero
  // line, which is exactly what happened before.
  double xmax = 2.0;
  for (const auto& r : q) {
    if (!r.comparable) continue;
    const double e = (r.a != 0.0) ? r.ediff / std::fabs(r.a) : 0.0;
    xmax = std::max(xmax, 100.0 * (std::fabs(r.relDiff) + e) * 1.30);
  }

  TH2F* frame = new TH2F("pullframe", "", 100, -xmax, xmax, n, 0, n);
  frame->SetDirectory(nullptr);
  for (int i = 0; i < n; ++i) frame->GetYaxis()->SetBinLabel(i + 1, q[i].name.c_str());
  frame->GetYaxis()->SetLabelSize(0.034);
  frame->GetXaxis()->SetTitle("(test - ref) / ref   [%]");
  frame->GetXaxis()->SetTitleSize(0.042);
  frame->Draw();

  TGraphErrors* g = new TGraphErrors();
  int np = 0;
  for (int i = 0; i < n; ++i) {
    if (!q[i].comparable) continue;  // nothing to plot; the row is labelled below
    const double rel = 100.0 * q[i].relDiff;
    const double e = (q[i].a != 0.0) ? 100.0 * q[i].ediff / std::fabs(q[i].a) : 0.0;
    g->SetPoint(np, rel, i + 0.5);
    g->SetPointError(np, e, 0.0);
    ++np;
  }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetLineWidth(2);
  g->SetLineColor(kColB);
  g->SetMarkerColor(kColB);
  g->Draw("P SAME");

  TLine* z = new TLine(0.0, 0.0, 0.0, n);
  z->SetLineStyle(2);
  z->SetLineColor(kGray + 2);
  z->Draw();

  TLatex tx;  // user coordinates, so the labels track the rows exactly
  tx.SetTextSize(0.029);
  tx.SetTextAlign(12);
  for (int i = 0; i < n; ++i) {
    if (!q[i].comparable) {
      tx.SetTextColor(kGray + 1);
      tx.DrawLatex(xmax * 1.04, i + 0.5, "no stats");
      continue;
    }
    tx.SetTextColor(std::fabs(q[i].pull) > 3.0 ? kRed + 1 : kGray + 3);
    tx.DrawLatex(xmax * 1.04, i + 0.5, Form("%+.1f#sigma", q[i].pull));
  }
  // Header block only: nothing is drawn at the bottom, where it used to sit on
  // top of the axis title.
  TLatex hd;
  hd.SetNDC();
  hd.SetTextSize(0.031);
  hd.SetTextAlign(12);
  hd.DrawLatex(0.02, 0.965, "scale-free flux fractions: bootstrap errors");
  hd.SetTextSize(0.023);
  hd.SetTextColor(kColA);
  hd.DrawLatex(0.02, 0.050, ("ref:  " + labA).c_str());
  hd.SetTextColor(kColB);
  hd.DrawLatex(0.02, 0.020, ("test: " + labB).c_str());
  c->Update();
}

void drawPValueOverview(TCanvas* c, const std::vector<fluxval::ShapeResult>& res,
                        const std::string& caseTag, const std::string& labRef,
                        const std::string& labTest) {
  c->Clear();
  c->cd();
  c->SetLeftMargin(0.24);
  c->SetTopMargin(0.14);
  c->SetRightMargin(0.06);
  c->SetBottomMargin(0.12);
  const int n = static_cast<int>(res.size());
  if (n == 0) return;

  TH2F* frame = new TH2F("pframe", "", 100, 0.0, 1.0, n, 0, n);
  frame->SetDirectory(nullptr);
  for (int i = 0; i < n; ++i)
    frame->GetYaxis()->SetBinLabel(i + 1, res[i].var.c_str());
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetXaxis()->SetTitle("p(AD), permutation");
  frame->GetXaxis()->SetTitleSize(0.040);
  frame->GetXaxis()->SetTitleOffset(1.25);
  frame->Draw();

  TGraph* g = new TGraph(n);
  for (int i = 0; i < n; ++i) g->SetPoint(i, std::max(res[i].pAD, 0.004), i + 0.5);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.3);
  g->SetMarkerColor(kColB);
  g->Draw("P SAME");

  TLine* l = new TLine(0.05, 0.0, 0.05, n);
  l->SetLineStyle(2);
  l->SetLineColor(kRed + 1);
  l->SetLineWidth(2);
  l->Draw();

  // A p-value at the permutation floor is drawn at the left edge; mark it so
  // it is not mistaken for a resolved small number.
  TLatex fl;
  fl.SetTextSize(0.026);
  fl.SetTextAlign(12);
  fl.SetTextColor(kGray + 3);
  for (int i = 0; i < n; ++i)
    if (res[i].pAD <= 1.0 / (res[i].nperm + 1.0) + 1e-12)
      fl.DrawLatex(0.055, i + 0.5, "at floor");

  TLatex hd;
  hd.SetNDC();
  hd.SetTextAlign(12);
  hd.SetTextSize(0.035);
  hd.DrawLatex(0.02, 0.965, Form("shape-test p-values   [%s]", caseTag.c_str()));
  hd.SetTextSize(0.023);
  hd.SetTextColor(kColA);
  hd.DrawLatex(0.02, 0.050, ("ref:  " + labRef).c_str());
  hd.SetTextColor(kColB);
  hd.DrawLatex(0.02, 0.020, ("test: " + labTest).c_str());
  hd.SetTextColor(kRed + 1);
  hd.SetTextSize(0.022);
  hd.DrawLatex(0.24, 0.912, "dashed: p = 0.05");
  c->Update();
}

void drawProfile(TCanvas* c, const fluxval::ShapeResult& R, const std::string& xtitle,
                 const std::string& caseTag, const std::string& labRef,
                 const std::string& labTest) {
  c->Clear();
  c->cd();
  c->SetLeftMargin(0.14);
  c->SetTopMargin(0.14);
  c->SetRightMargin(0.13);
  c->SetBottomMargin(0.13);
  const int nb = R.nbins;
  if (nb < 2 || static_cast<int>(R.edges.size()) != nb + 1) return;

  double ymax = 0.0;
  for (int k = 0; k < nb; ++k)
    ymax = std::max(ymax, std::max(std::fabs(R.dCdf[k]),
                                   std::max(std::fabs(R.lo95[k]), std::fabs(R.hi95[k]))));
  if (!(ymax > 0.0)) {  // the two samples are bit-identical: nothing to profile
    TLatex t;
    t.SetNDC();
    t.SetTextAlign(22);
    t.SetTextSize(0.035);
    t.DrawLatex(0.5, 0.55, Form("%s [%s]", R.var.c_str(), caseTag.c_str()));
    t.SetTextSize(0.028);
    t.SetTextColor(kGray + 3);
    t.DrawLatex(0.5, 0.47, "reference and test are identical: F_{test} - F_{ref} = 0 "
                           "everywhere");
    c->Update();
    return;
  }
  ymax *= 1.35;

  TH2F* frame = new TH2F("profframe", "", 100, R.edges.front(), R.edges.back(), 100,
                         -ymax, ymax);
  frame->SetDirectory(nullptr);
  frame->GetXaxis()->SetTitle(xtitle.c_str());
  frame->GetYaxis()->SetTitle("F_{test} - F_{ref}");
  frame->GetXaxis()->SetTitleSize(0.042);
  frame->GetYaxis()->SetTitleSize(0.042);
  frame->GetYaxis()->SetTitleOffset(1.35);
  // Each axis takes the colour of the curve it describes, so there is never a
  // question of which scale a line is read against.
  frame->GetYaxis()->SetTitleColor(kColB);
  frame->GetYaxis()->SetLabelColor(kColB);
  frame->GetYaxis()->SetAxisColor(kColB);
  frame->Draw();

  // Permutation bands, widest first so the 68% sits on top of the 95%.
  auto bandGraph = [&](const std::vector<double>& lo, const std::vector<double>& hi,
                       int col, double alpha) {
    TGraph* g = new TGraph(2 * nb + 1);
    for (int k = 0; k < nb; ++k)
      g->SetPoint(k, 0.5 * (R.edges[k] + R.edges[k + 1]), hi[k]);
    for (int k = nb - 1; k >= 0; --k)
      g->SetPoint(2 * nb - 1 - k, 0.5 * (R.edges[k] + R.edges[k + 1]), lo[k]);
    g->SetPoint(2 * nb, 0.5 * (R.edges[0] + R.edges[1]), hi[0]);
    g->SetFillColorAlpha(col, alpha);
    g->SetLineWidth(0);
    g->Draw("F SAME");
  };
  bandGraph(R.lo95, R.hi95, kGray + 1, 0.35);
  bandGraph(R.lo68, R.hi68, kGray + 2, 0.45);

  TLine* z = new TLine(R.edges.front(), 0.0, R.edges.back(), 0.0);
  z->SetLineStyle(2);
  z->SetLineColor(kGray + 3);
  z->Draw();

  // Observed signed CDF difference.
  TGraph* obs = new TGraph(nb);
  for (int k = 0; k < nb; ++k)
    obs->SetPoint(k, 0.5 * (R.edges[k] + R.edges[k + 1]), R.dCdf[k]);
  obs->SetLineColor(kColB);
  obs->SetLineWidth(3);
  obs->Draw("L SAME");

  // Cumulative A^2, mapped onto the same pad and read off the right axis.
  TGraph* cum = new TGraph(nb);
  for (int k = 0; k < nb; ++k)
    cum->SetPoint(k, 0.5 * (R.edges[k] + R.edges[k + 1]),
                  -ymax + R.cumA2[k] * 2.0 * ymax);
  cum->SetLineColor(kColA);
  cum->SetLineWidth(3);
  cum->SetLineStyle(7);
  cum->Draw("L SAME");

  TGaxis* ax = new TGaxis(R.edges.back(), -ymax, R.edges.back(), ymax, 0.0, 1.0, 510,
                          "+L");
  ax->SetTitle("cumulative share of A^{2}");
  ax->SetTitleSize(0.038);
  ax->SetLabelSize(0.034);
  ax->SetTitleOffset(1.30);
  ax->SetLineColor(kColA);
  ax->SetLabelColor(kColA);
  ax->SetTitleColor(kColA);
  ax->Draw();

  TLatex hd;
  hd.SetNDC();
  hd.SetTextAlign(12);
  hd.SetTextSize(0.034);
  hd.DrawLatex(0.02, 0.975,
               Form("%s [%s]   CDF difference and A^{2} contribution", R.var.c_str(),
                    caseTag.c_str()));
  hd.SetTextSize(0.023);
  hd.SetTextColor(kColB);
  hd.DrawLatex(0.02, 0.940, "solid: F_{test} - F_{ref}   (above 0 = test distribution "
                            "sits further left)");
  hd.SetTextColor(kColA);
  hd.DrawLatex(0.02, 0.908, "dashed: running share of A^{2}   (fastest rise at "
                            "largest disagreement bins)");
  hd.SetTextColor(kGray + 3);
  hd.SetTextSize(0.021);
  hd.DrawLatex(0.55, 0.20, "grey: 68% / 95% permutation bands");
  hd.DrawLatex(0.55, 0.168, "pointwise -- read runs of bins, not single ones");
  c->Update();
}

void drawPP(TCanvas* c, const std::vector<std::string>& varNames,
            const std::vector<std::vector<double>>& pvals) {
  c->Clear();
  c->cd();
  c->SetLeftMargin(0.14);
  c->SetTopMargin(0.09);
  c->SetRightMargin(0.06);
  c->SetBottomMargin(0.13);

  TH2F* frame = new TH2F("ppframe", "", 100, 0.0, 1.0, 100, 0.0, 1.0);
  frame->SetDirectory(nullptr);
  frame->GetXaxis()->SetTitle("expected  k / (n+1)");
  frame->GetYaxis()->SetTitle("observed  sorted p(AD)");
  frame->GetXaxis()->SetTitleSize(0.042);
  frame->GetYaxis()->SetTitleSize(0.042);
  frame->Draw();

  TLine* d = new TLine(0, 0, 1, 1);
  d->SetLineStyle(2);
  d->SetLineColor(kGray + 2);
  d->SetLineWidth(2);
  d->Draw();

  const int nv = static_cast<int>(varNames.size());
  const int ncol = TColor::GetNumberOfColors();
  for (int v = 0; v < nv; ++v) {
    std::vector<double> p = pvals[v];
    if (p.empty()) continue;
    std::sort(p.begin(), p.end());
    const int n = static_cast<int>(p.size());
    TGraph* g = new TGraph(n);
    for (int k = 0; k < n; ++k) g->SetPoint(k, (k + 1.0) / (n + 1.0), p[k]);
    const int col = (ncol > 8)
        ? TColor::GetColorPalette(static_cast<int>((0.10 + 0.78 * v / std::max(1, nv - 1)) * ncol))
        : (kAzure + v);
    g->SetMarkerStyle(20 + (v % 4));
    g->SetMarkerSize(0.9);
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->Draw("PL SAME");
    TLatex lg;
    lg.SetNDC();
    lg.SetTextSize(0.026);
    lg.SetTextColor(col);
    lg.DrawLatex(0.62, 0.36 - 0.033 * v, varNames[v].c_str());
  }

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.033);
  tx.DrawLatex(0.14, 0.955, "P-P plot of closure p-values against uniform");
  tx.SetTextSize(0.025);
  tx.SetTextColor(kGray + 3);
  tx.DrawLatex(0.16, 0.90, "above the diagonal = conservative (misses real differences)");
  tx.DrawLatex(0.16, 0.865, "below the diagonal = anti-conservative (invents significance)");
  c->Update();
}

void drawTextPage(TCanvas* c, const std::string& title,
                  const std::vector<std::string>& lines) {
  c->Clear();
  c->cd();
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.035);
  tx.DrawLatex(0.06, 0.955, title.c_str());
  tx.SetTextFont(82);  // fixed pitch, so the columns line up
  tx.SetTextSize(0.0215);
  double y = 0.905;
  for (const std::string& s : lines) {
    if (y < 0.03) break;
    tx.DrawLatex(0.06, y, s.c_str());
    y -= 0.0255;
  }
  c->Update();
}

}  // namespace fluxval
