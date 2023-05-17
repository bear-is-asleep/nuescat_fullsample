#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

//Systs
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

//#include "SystEnsembles.h"

using namespace ana;

#include "Constants.h"
#include "Structs.h"
#include "utils.h"
#include "TrueEventCategories.h"
//#include "NuEScatterRecoVars.h"
//#include "NuEScatterCuts.h"
#include "NuEScatterTrueVars.h"
#include "plotStyle.C"

#include <string>
#include <fstream>
#include <iostream>
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TFile.h"

using namespace std;

const string state_fname = "test_flux.root";
const bool do_systematics = true;

void test_flux_var(){
  selectionstyle();
  setlocale(LC_NUMERIC, "");

  const std::string inputNameNu = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_flat_caf_sbnd";
  const double gPOT = 10e20;
  const string surName = "test_flux_var";
  const TString saveDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/tests/"+get_date()+"_"+surName;
  const TString stateDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/tests/"+get_date()+"_"+surName;

  SpectrumLoader loader(inputNameNu);

  //Get bins for both spectrum
  const Binning binsEnergy = Binning::Simple(80, 0, 4);
  HistAxis ax("True E_{#nu} (GeV)",binsEnergy,kTrueNuESlice);

  std::vector<Spectrum*> sNu;
  const std::vector<std::string>& flux_systs = GetSBNBoosterFluxWeightNames();
  for (auto const& name : flux_systs){
    //std::cout<<name<<endl;
  }
  vector<Var> weis(1000,kUnweighted); //Initialize weights
  vector<Var> weis_flux(1000,kUnweighted); //Initialize weights
  for (unsigned j=0; j<1000; j++){
    for (unsigned i =0; i<flux_systs.size(); i++){
      weis_flux[j] = weis_flux[j]*GetUniverseWeight(flux_systs[i]+"_Flux", j); //Get weight from spill
    }
    //weis_flux_reweighted[j] = nue_reweight[j]*weis_flux[j]; //Reweight universe by nu+e constraint
    weis[j] = weis[j]*weis_flux[j]; //Get total flux weight and evenly weight it with the genie weight
    //weis[j] = weis[j]*kFluxWeight_0; //Reweight by cross section and number of targets to get flux 
  }

  Spectrum *s_Nu0_CV = new Spectrum("ENu_CV", binsEnergy, loader, kTrueNuESlice, kNoCut, kNoShift, kUnweighted);
  Spectrum *s_Nu0_uni1 = new Spectrum("ENu_uni1", binsEnergy, loader, kTrueNuESlice, kNoCut, kNoShift, weis[0]);

  loader.Go();
  //Save state
  gSystem->Exec("mkdir -p " + stateDir);
  TFile *file_all = new TFile(stateDir+"/state_all.root", "RECREATE");

  file_all->cd();
  TH1D *hist_nom = s_Nu0_CV->ToTH1(gPOT,kRed);
  TH1D *hist_uni = s_Nu0_uni1->ToTH1(gPOT,kRed);

  hist_nom->Write();
  hist_uni->Write();
  file_all->Close();

  




}