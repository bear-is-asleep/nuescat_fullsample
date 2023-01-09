#include "TString.h"

const FidVol fvndNuEScat{  +1.5, +190.,  // x
                      -190., +190.,  // y
                       +20., +480.}; // z

template <typename T>
struct Plot {
  TString name = "";
  T variable; 
  Binning binning;
  TString axes_labels;
  TString label;
  std::vector<float> legBox = {.1, .1, .4, .3};
};

struct CutDef {
  TString name = "";
  TString label = "";
  SpillCut cut = kNoSpillCut;
  int colour = kBlack;
};

struct TrueCategory {
  TString name = "";
  SpillCut cut = kNoSpillCut;
  int colour = kBlack;
  TString label = "";
};