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
#include "utils.h"
#include "TrueEventCategories.h"
#include "NuEScatterRecoVars.h"
#include "NuEScatterTrueVars.h"
#include "NuEScatterCuts.h"
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

const bool do_systematics = true;
const Cut kReweightSelection = kIsNuMuNCQEOnArgonVar && kIsTrueAVVar;
//Use CRUMBs to select best slice to reject cosmics

unsigned int nuni = 100; //number of universes to check

void NuEScatter_reweight(bool save = true)
{

  selectionstyle();
  setlocale(LC_NUMERIC, "");

  const std::string inputNameNu = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_flat_caf_sbnd";


  const double gPOT = 10e20;
  const string surName = "reweight_all_flux_var";
  const TString saveDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/"+get_date()+"_"+surName;
  const TString stateDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/"+get_date()+"_"+surName;
  const TString weightsName = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut/weights.txt";

  SpectrumLoader loaderNu(inputNameNu);

  std::vector<Spectrum*> sNu;

  const std::vector<std::string>& flux_systs = GetSBNBoosterFluxWeightNames();
  vector<Var> weis(nuni,kUnweighted); //Initialize weights
  vector<double> nue_reweight(nuni,1.);
  //Read in reweighted values
  ifstream reweight_file(weightsName);
  int weight,i = 0;
  while(reweight_file >> weight && i <= nuni){
    nue_reweight[i++] = weight; //assign weight to array
  }
  //const Var kTrueNuESlice = SIMPLEVAR(truth.E); //Neutrino energy
  const Binning binsEnergy = Binning::Simple(100, 0, 5);
  HistAxis ax("True E_{#nu} (GeV)",binsEnergy,kTrueNuESlice);

  //Flux syst
  vector<Var> weis_flux(nuni,kUnweighted); //Initialize weights
  vector<Var> weis_flux_reweighted(nuni,kUnweighted); //Initialize weights
  vector<Var> weis_tot(nuni,kUnweighted); //Initialize weights
  vector<Var> weis_tot_reweighted(nuni,kUnweighted); //Initialize weights
  vector<vector<Var>> weis_flux_all(flux_systs.size(),weis_flux); 
  for (unsigned j=0; j<nuni; j++){
    weis_tot[j] = weis_tot[j]*GetUniverseWeight("multisim_Genie", j); //Get weight from genie
    for (unsigned i =0; i<flux_systs.size(); i++){
      weis_flux[j] = weis_flux[j]*GetUniverseWeight(flux_systs[i], j); //Get weight from 
      weis_flux_all[i][j] = weis_flux_all[i][j]*GetUniverseWeight(flux_systs[i], j);
    }
    //weis_flux_reweighted[j] = nue_reweight[j]*weis_flux[j]; //Reweight universe by nu+e constraint
    //weis_flux_reweighted[j] = weis_flux[j];
    weis_tot[j] = weis_tot[j]*weis_flux[j]; //Get total flux weight and evenly weight it with the genie weight
    //weis_tot_reweighted[j] = weis_tot[j]*weis_flux_reweighted[j]; //Reweight universe by nu+e constraint
  }
  std::vector<TString> syst_names;// = {"GENIE_multisim"};
  std::vector<TString> syst_labels;// = {"GENIE"};
  std::vector<int> syst_colors;// = {kCyan};

  //Add specific flux vals
  int counter = 0;
  for (auto const& flux_syst : flux_systs){
    syst_names.push_back(flux_syst);
    syst_labels.push_back(flux_syst.substr(0,9));
    syst_colors.push_back(colors[counter]);
    counter+=1;
  }
  
  //Add combined flux count
  syst_names.insert(syst_names.end(),{"total_multisim","total_reweighted","total_flux","total_flux_unweighted"});
  syst_labels.insert(syst_labels.end(),{"Total","Total Reweighted","Flux","Flux Unweighted"});
  syst_colors.insert(syst_colors.end(),{kBlack,kCyan,kRed,kBlue});

  std::vector<SystEnsemble> systs;
  if (do_systematics){
    for (unsigned i=0; i<syst_names.size(); i++){
      //Assign different weights to each
      if (syst_names[i] == "total_multisim"){
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueNuESlice,
        new EnsembleSpectrum(loaderNu,ax,kNoSpillCut,kReweightSelection,weis_tot)};
        systs.emplace_back(syst);
      }
      else if (syst_names[i] == "total_reweighted"){
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueNuESlice,
        new EnsembleSpectrum(loaderNu,ax,kNoSpillCut,kReweightSelection,weis_tot,kFluxWeightVar)};
        systs.emplace_back(syst);
      }
      else if (syst_names[i] == "total_flux"){
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueNuESlice,
        new EnsembleSpectrum(loaderNu,ax,kNoSpillCut,kReweightSelection,weis_flux,kFluxWeightVar)};
        systs.emplace_back(syst);
      }
      else if (syst_names[i] == "total_flux_unweighted"){
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueNuESlice,
        new EnsembleSpectrum(loaderNu,ax,kNoSpillCut,kReweightSelection,weis_flux)};
        systs.emplace_back(syst);
      }
      else{
        SystEnsemble syst = {syst_names[i],syst_colors[i],syst_labels[i],kTrueNuESlice,
        new EnsembleSpectrum(loaderNu,ax,kNoSpillCut,kReweightSelection,weis_flux_all[i],kFluxWeightVar)};
        systs.emplace_back(syst);
      }
    }
  }
  //std::cout.setstate(std::ios_base::failbit);
  loaderNu.Go();
  //std::cout.clear();

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

  //This is to save the universe sample values to construct the covariance matrix
  if (save){

    // Fill hist with all systs.
    for (auto const& syst: systs){
      const TString state_fname = "state_"+syst.name;
      // Create a new TFile to write to
      // TFile *file_nominal = new TFile(stateDir+"/"+state_fname+"_nom.root", "RECREATE");
      // file_nominal->cd();

      // //Save to nominal file
      // TH1D *hist_nominal = syst.syst->Nominal().ToTH1(gPOT, syst.color);
      // hist_nominal->Write();
      // file_nominal->Close();

      TFile *file_all = new TFile(stateDir+"/"+state_fname+".root", "RECREATE");
      file_all->cd();
      //Save to all file
      TH1D *hist_all = syst.syst->Nominal().ToTH1(gPOT, syst.color);
      hist_all->Write();
      for (unsigned iUni = 0; iUni < syst.syst->NUniverses(); ++iUni){
        syst.syst->Universe(iUni).ToTH1(gPOT, syst.color)->Write();
      }
      file_all->Close();
    }
  }
  
  gSystem->Exec("cp NuEScatter_reweight.C NuEScatterRecoVars.h NuEScatterCuts.h TrueEventCategories.h Constants.h Structs.h " + stateDir);
}