// Plot.h -- overlay + ratio panels, permutation nulls, pull summaries.

#pragma once

#include "fluxval/Core.h"

#include <string>
#include <vector>

class TCanvas;
class TH1D;

namespace fluxval {

void setStyle();

/// The two sample colours, taken from the active palette by setStyle().
int colourA();
int colourB();

/// Build a density histogram (unit area, then divided by bin width) on the
/// variable-width quantile binning, with bootstrap errors.
TH1D* densityHist(const char* name, const std::vector<double>& edges,
                  const std::vector<double>& shape, const std::vector<double>& err);

/// B/A shape ratio as a histogram, errors from the bootstrap band.
TH1D* ratioHist(const char* name, const std::vector<double>& edges,
                const std::vector<double>& ratio, const std::vector<double>& err);

/// The standard two-pad layout: overlay on top, ratio underneath.
void drawWithRatio(TCanvas* c, TH1D* hA, TH1D* hB, TH1D* hR, const std::string& xtitle,
                   bool logy, const std::string& labA, const std::string& labB,
                   const std::string& annotation, double ratioSpan = 0.0);

/// Permutation null of the AD statistic with the observed value marked.
void drawNull(TCanvas* c, const fluxval::ShapeResult& R);

/// Effective statistics per bin -- the ceiling on what any test can resolve.
void drawNeff(TCanvas* c, const std::vector<double>& edges,
              const std::vector<double>& nA, const std::vector<double>& nB,
              const std::string& xtitle, const std::string& labA,
              const std::string& labB);

/// Relative difference (B-A)/A of every scale-free quantity, with pulls.
void drawQuantityPulls(TCanvas* c, const std::vector<fluxval::QuantityResult>& q,
                       const std::string& labA, const std::string& labB);

/// p(AD) for every variable of ONE charge selection, 0.05 marked.
void drawPValueOverview(TCanvas* c, const std::vector<fluxval::ShapeResult>& res,
                        const std::string& caseTag, const std::string& labRef,
                        const std::string& labTest);

/// Discrepancy profile: signed CDF difference with permutation bands, and the
/// cumulative contribution to A^2, overlaid on one pad with a right-hand axis.
void drawProfile(TCanvas* c, const fluxval::ShapeResult& R, const std::string& xtitle,
                 const std::string& caseTag, const std::string& labRef,
                 const std::string& labTest);

/// P-P plot: sorted p-values against k/(n+1).  Under the null they lie on the
/// diagonal; systematically above means the test is conservative.
void drawPP(TCanvas* c, const std::vector<std::string>& varNames,
            const std::vector<std::vector<double>>& pvals);

/// A text page carrying the numerical summary into the PDF.
void drawTextPage(TCanvas* c, const std::string& title,
                  const std::vector<std::string>& lines);

}  // namespace fluxval
