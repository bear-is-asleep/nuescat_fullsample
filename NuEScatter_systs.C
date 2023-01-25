#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

//Cov mat
#include "sbnana/CAFAna/Core/Utilities.h"
//#include "/sbnd/app/users/brindenc/mysbn/srcs/sbnana/sbnana/CAFAna/Core/Utilities.h"

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
#include "TrueEventCategories.h"
#include "NuEScatterRecoVars.h"
#include "NuEScatterCuts.h"
//#include "TrueEventCategories_simple.h"
//#include "NuEScatterRecoVars_simple.h"
//#include "NuEScatterCuts_simple.h"
#include "utils.h"
#include "Reducer.h"
#include "plotStyle.C"

#include <string>
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TFile.h"

using namespace std;

const string state_fname = "NuEScatter_state_all.root";
const bool do_systematics = true;

void NuEScatter_ana(bool reload = true)
{

  selectionstyle();
  setlocale(LC_NUMERIC, "");

  const std::string inputNameNu = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_flat_caf_sbnd";


  const double gPOT = 10e20;
  const bool save = true;
  const string surName = "systs";
  const TString saveDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/"+get_date()+"_"+surName;
  const TString stateDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/"+get_date()+"_"+surName;

  SpectrumLoader loaderNu(inputNameNu);

  std::vector<Spectrum*> sNu;

  const std::vector<std::string>& flux_systs = GetSBNBoosterFluxWeightNames();
  vector<Var> weis; //Initialize weights
  for (unsigned j=0; j<1000; j++) {
    weis.push_back(kUnweighted);
  }
  const Var kTrueE = SIMPLEVAR(truth.E);
  //const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  //const Binning binsEnergy = Binning::Simple(6, 0, 4);
  const vector<double>& edges = {0,0.8,1.2,1.9,2.6,3.5,10};
  const Binning binsEnergy = Binning::Custom(edges);
  HistAxis ax("True Energy (GeV)",binsEnergy,kTrueE);
  for (unsigned j=0; j<1000; j++){
    weis[j] = weis[j]*GetUniverseWeight("multisim_Genie", j); //Get weight from 
  }

  //Flux syst
  vector<Var> weis_flux(1000,kUnweighted); //Initialize weights
  vector<vector<Var>> weis_flux_all(flux_systs.size(),weis_flux); 
  for (unsigned i =0; i<flux_systs.size(); i++){
    for (unsigned j=0; j<1000; j++){
      weis_flux[j] = weis_flux[j]*GetUniverseWeight(flux_systs[i], j); //Get weight from 
      weis_flux_all[i][j] = weis_flux_all[i][j]*GetUniverseWeight(flux_systs[i], j);
    }
  }
  std::vector<TString> syst_names;// = {"GENIE_multisim"};
  std::vector<TString> syst_labels;// = {"GENIE"};
  std::vector<int> syst_colors;// = {kCyan};
  
  int counter = 0;
  for (auto const& flux_syst : flux_systs){
    syst_names.push_back(flux_syst);
    syst_labels.push_back(flux_syst.substr(0,9));
    syst_colors.push_back(colors[counter]);
    counter+=1;
  }
  syst_names.insert(syst_names.end(),{"Flux_multisim","GENIE_multisim"});
  syst_labels.insert(syst_labels.end(),{"Flux","GENIE"});
  syst_colors.insert(syst_colors.end(),{colors[counter],kCyan});

  std::vector<SystEnsemble> systs;
  if (do_systematics){
    for (unsigned i=0; i<syst_names.size(); i++){
      if (syst_names[i] == "GENIE_multisim"){
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueE,
        new EnsembleSpectrum(loaderNu,ax,kEthetaSelection,kNoCut,weis)};
        systs.emplace_back(syst);
      }
      else if (syst_names[i] == "Flux_multisim"){
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueE,
        new EnsembleSpectrum(loaderNu,ax,kEthetaSelection,kNoCut,weis_flux)};
        systs.emplace_back(syst);
      }
      else { //Important to put all of the flux weights in alignment with syst_name indexing
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueE,
        new EnsembleSpectrum(loaderNu,ax,kEthetaSelection,kNoCut,weis_flux_all[i])};
        systs.emplace_back(syst);
      }
    }
  }
  //std::cout.setstate(std::ios_base::failbit);
  loaderNu.Go();
  std::cout.clear();
  //loaderIntime.Go();

  TFile fout(state_fname.c_str(),"RECREATE");

  gSystem->Exec("mkdir -p " + stateDir);

  TCanvas *canvas = new TCanvas("c");
  TLegend *legend = new TLegend(.69,.57,.89,.95);
  double max_val = -9999;
  for (auto const& syst : systs){
    TH1D *hist_nom = syst.syst->Nominal().ToTH1(gPOT, syst.color);
    for(unsigned int i = 0; i < syst.syst->NUniverses(); ++i){
      TH1D *hist = syst.syst->Universe(i).ToTH1(gPOT, syst.color-10);
      if (hist->GetMaximum() > max_val){max_val = hist->GetMaximum();}
      hist->SetLineColorAlpha(syst.color, 0.1);
      hist->Draw("hist same");
    }
    hist_nom->Draw("hist same");

    legend->AddEntry(hist_nom, syst.label, "lf");
  }
	legend->SetTextSize(0.03);
  legend->SetLineColorAlpha(kGray,0.9);
  legend->SetLineWidth(2);
  legend->SetFillColorAlpha(kWhite, 0.7);
  legend->Draw();
  gSystem->Exec("mkdir -p "+saveDir+"/systs");
  gPad->Print(saveDir + "/systs/all.png");

  for (auto const& syst: systs){
    new TCanvas;
    TGraphAsymmErrors* band = syst.syst->ErrorBand(gPOT,"error_band_"+syst.name,stateDir+"/");
    DrawErrorBand(syst.syst->Nominal().ToTH1(gPOT, syst.color), band);
    gPad->Print(saveDir + "/systs/error_band_"+syst.name+".png");
  }

}