// Load.h -- stream FairShip trees into the binned accumulator.
//
// All files of a sample are handed to a single RDataFrame and read as one
// dataset with implicit MT, rather than opened one after another: with a few
// hundred files that is the difference between minutes and hours.  Nothing is
// kept per muon, so memory does not depend on how much data goes in.
 
#pragma once
 
#include "fluxval/Core.h"
 
#include <string>
#include <vector>
 
namespace fluxval {
 
struct LoadOpts {
  std::string tree = "cbmsim";
  std::string source = "plane";
  // Any scoring-plane branch exposing fX, fY, fPx, fPy, fPz, fPdgCode and
  // fTrackID.  The member names are assumed identical across such branches.
  std::string branch = "PlaneHAPoint";
  std::vector<int> pdgAbs{13};  // accepted |PDG| codes
  double pmin = 5.0;
  double pmax = 400.0;
  double ptEps = 0.05;
  bool requirePzPositive = true;
 
  // Accumulator geometry.  Memory is nUnits * 2 * nvar * nfine * 8 bytes per
  // thread; the driver prints the estimate before reading anything.
  int nfine = 2048;    // fine bins per variable
  int maxUnits = 256;  // cap on resampling units (files are grouped above it)
  int nBlocks = 200;   // random event blocks when there are too few files
  int minFiles = 5;    // use files as units at or above this count
 
  bool progress = true;  // show a progress bar while reading
 
  // Input reduction.  maxEvents takes the FIRST N entries of the sample's file
  // list, in list order, via RDataFrame::Range.  A prefix is used rather than
  // a hashed pseudo-random subset because a prefix is trivially reproducible:
  // two loads of the same list necessarily see the same entries.  A hashed
  // subset did not provide that -- comparing one file against itself under it
  // yielded two different samples -- so it was removed.  Range is rejected
  // under implicit MT, so setting this disables threading for the read.
  long long maxEvents = 0;  // 0 = no limit
  int maxFiles = 0;
};
 
/// "a.root,b.root", "@list.txt", or a directory (all *.root inside it).
std::vector<std::string> splitList(const std::string& spec);
 
/// Read a limited number of entries from the front of both file lists to fix
/// the histogram ranges.  Both samples must share one grid, or pooling them
/// would be meaningless.
fluxval::HistSpec probeRange(const std::vector<std::string>& a,
                         const std::vector<std::string>& b, const LoadOpts& o);
 
/// One pass over every file, filling the accumulator.
fluxval::UnitHists loadHists(const std::vector<std::string>& files, const LoadOpts& o,
                         const fluxval::HistSpec& spec, const std::string& label);
 
}  // namespace fluxval

