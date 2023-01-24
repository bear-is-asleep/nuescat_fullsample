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

void NuEScatter_ana(bool reload = true)
{

  selectionstyle();
  setlocale(LC_NUMERIC, "");

  //gROOT->SetBatch(kTRUE);
  
  //Use my sample with nu e scat
  const std::string inputNameNuE = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_58_02/CAFnue1_full.root";
  const std::string inputNameNuECC = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_58_02/CAFnuecc1_full.root";
  const std::string inputNameNu = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_flat_caf_sbnd";
  const std::string inputNameNu_noflat = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_caf_sbnd";


  const double gPOT = 10e20;
  const bool save = true;
  const string surName = "fullsample_originalcuts_wsysts";
  const TString saveDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/"+get_date()+"_"+surName;
  const TString stateDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/"+get_date()+"_"+surName;
  //const TString date = get_date();

  SpectrumLoader loaderNu(inputNameNu);
  //SpectrumLoader loaderIntime(inputNameIntime);

  std::vector<Spectrum*> sNu;
  //std::vector<EnsembleSpectrum*> sNuSyst;
  //sIntime;

  SpillCut previousCuts = kNoSpillCut;

  const std::vector<std::string>& flux_systs = GetSBNBoosterFluxWeightNames();
  //systs = {"multisim_Genie"};
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
  //std::cout<<"Initialize weights"<<endl;
  vector<Var> weis_flux(1000,kUnweighted); //Initialize weights
  vector<vector<Var>> weis_flux_all(flux_systs.size(),weis_flux); 
  for (unsigned i =0; i<flux_systs.size(); i++){
    for (unsigned j=0; j<1000; j++){
      weis_flux[j] = weis_flux[j]*GetUniverseWeight(flux_systs[i], j); //Get weight from 
      weis_flux_all[i][j] = weis_flux_all[i][j]*GetUniverseWeight(flux_systs[i], j);
    }
  }
  //cout<<"Made weights"<<endl;
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
  //syst_labels.push_back("Flux");
  //syst_colors.push_back(counter+2);



  //cout<<"Made vectors for weights"<<endl;
  // std::vector<TString> syst_names = {"CCNuE","nuescat","all"};
  // std::vector<TString> syst_labels = {"CC #nu_{e}","#nu + e","All"};
  // std::vector<int> syst_colors = {kCyan,kRed,kOrange};

  for (auto const& cut : nuescatter_cuts){
    //std::cout<<cut.name<<std::endl;
    for(auto const& plot : recoPlots){
      //std::cout<<"--"<<plot.name<<std::endl;
      for(auto const& category : nuescat_sel_categories){
        string name = *category.label.Data() + "_" + *plot.name.Data();
	      sNu.emplace_back(new Spectrum("nu_" + name, plot.binning, loaderNu, 
					    plot.variable, 
					    category.cut && previousCuts && cut.cut));
	      // sIntime.emplace_back(new Spectrum("intime_" + name, plot.binning, 
				// 		loaderIntime, plot.variable,  
				// 		category.cut && previousCuts && cut.cut));
        
      }
    }
    previousCuts = previousCuts && cut.cut;
    if (cut.name == "Etheta2 Cut"){
      //std::cout<<"reducer started";
      //reducer(previousCuts,inputNameNu_noflat);
      //std::cout<<"reducer finished"<<std::endl;
      
    }
  }
  //std::vector<SpillCut> syst_cuts = {previousCuts && kCCNuE, previousCuts && kNuEScat, previousCuts};
  //std::vector<SystEnsemble> systs = setup_systematics(loaderNu,kTrueE,syst_names,syst_cuts,syst_colors,syst_labels,weis,ax);
  std::vector<SystEnsemble> systs;
  //std::vector<SystEnsemble> systs_flux;
  std::cout<<"Making ensemble"<<endl;

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

    // SystEnsemble syst_flux = {syst_names[i],syst_colors[i],syst_labels[i],kTrueE,
    //   new EnsembleSpectrum(loaderNu,ax,syst_cuts[i],kNoCut,weis_flux)};
    // systs_flux.emplace_back(syst_flux);
    //syst = {};
  }
  //std::vector<std::vector<SystEnsemble>> systematics = {systs,systs_flux};

  //if(reload || TFile(state_fname.c_str()).IsZombie()){
  std::cout.setstate(std::ios_base::failbit);
  loaderNu.Go();
  std::cout.clear();
  //loaderIntime.Go();

  TFile fout(state_fname.c_str(),"RECREATE");

  gSystem->Exec("mkdir -p " + stateDir);
  //gSystem->Exec("mkdir -p " + stateDir + "/nuescat");
  //std::cout<<"Saving states"<<std::endl;
  // sNuSyst.SaveTo(fout.mkdir("all"));
  // sNueScatSyst.SaveTo(fout.mkdir("nuescat"));
  // sNueScatSyst.SaveTo(fout.mkdir("nuecc"));
  //sNcpi0Syst.SaveTo(fout.mkdir("nuescat"));
  //std::cout<<"Saved states"<<std::endl;

  const double nomInt           = 5e12;
  const double targetLive       = gPOT / nomInt;
  const double nuPOT            = sNu[0]->POT();
  const double nuLive           = sNu[0]->Livetime();
  const double nuScaledLive     = nuLive * gPOT / nuPOT;
  const double intimeScaledLive = targetLive - nuScaledLive;
  const int nuEvents            = nuLive;
  const double cosmicPOT        = nuEvents * nomInt / ( 1 - (1.f / 21.8f));

  unsigned i = 0, j = 0;

  double initialSignal = 0.;

  for(auto const& cut : nuescatter_cuts)
    {

    gSystem->Exec("mkdir -p " + saveDir + "/" + cut.label);
    for(auto const& plot : recoPlots)
	    {

        TCanvas *canvas = new TCanvas("c " + plot.name + " " + cut.name, 
              "c " + plot.name + " " + cut.name);
        canvas->cd();
        TLegend *legend = new TLegend(plot.legBox[0], plot.legBox[1], plot.legBox[2], plot.legBox[3]);
        THStack *stack = new THStack("stack_" + plot.name + "_" + cut.name, plot.axes_labels);
	  
	      double selectedSig = 0., selected = 0.;

        for(auto const& category : nuescat_sel_categories)
          {
            TH1D *hist = sNu[i]->ToTH1(gPOT);
            //hist->Add(sIntime[i]->ToTH1(intimeScaledLive, kLivetime));
            double integral = sNu[i]->Integral(gPOT);
            //integral += sIntime[i]->Integral(intimeScaledLive, 0, kLivetime);

            if(j == 0 && category.label == "NuEScat")
        initialSignal = integral;
            
            if(category.label == "NuEScat")
        selectedSig = integral;

            selected += integral;

            hist->SetLineColor(category.colour);
            hist->SetFillColorAlpha(category.colour, 0.3);

            stack->Add(hist);
            legend->AddEntry(hist, Form("%s (%'.0f)", category.name.Data(), integral), "lf");
            ++i;
          }
        if(plot.name == "FV Cut")
          stack->SetMaximum(1.8 * stack->GetMaximum());
          stack->Draw("hist");

	  const std::vector<std::string> binlabels = plot.binning.Labels();
	  if(!binlabels.empty())
	    {
	      for(int bin = 1; bin <= stack->GetXaxis()->GetNbins(); ++bin)
		stack->GetXaxis()->SetBinLabel(bin, binlabels[bin-1].c_str());
	    }

	  legend->SetLineColorAlpha(kGray,0.9);
    legend->SetLineWidth(2);
	  legend->SetTextSize(0.04);
    legend->SetFillColorAlpha(kWhite, 0.7);
	  legend->Draw();

	  const double eff = selectedSig / initialSignal * 100;
	  const double pur = selectedSig / selected * 100;

	  double box[4] = {.62, .36, .82, .48};
	  
	  if(plot.name == "FV Cut")
	    { box[0] = .25; box[1] = .68; box[2] = .45; box[3] = .80; }
	  else if(plot.name == "CRUMBS Score" || plot.name == "Leading Shower PDG" || plot.name == "SubLeading Shower PDG")
	    { box[0] = .25; box[1] = .43; box[2] = .47; box[3] = .55; }
	  else if(plot.name == "Longest Track PDG")
 	    { box[0] = .35; box[1] = .43; box[2] = .57; box[3] = .55; }
	  TPaveText *pvt = new TPaveText(box[0],box[1],box[2],box[3], "blNDC");
	  pvt->SetTextAlign(12);
	  pvt->SetTextSize(0.04);
	  pvt->AddText(Form("Eff = %.1f%%", eff));
	  pvt->AddText(Form("Pur = %.1f%%", pur));
    pvt->SetFillColorAlpha(kWhite, 0.5);
    pvt->SetLineColorAlpha(kGray,0.9);
    pvt->SetLineWidth(2);
    pvt->SetBorderSize(1);
	  pvt->Draw();	  

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + cut.label + "/" + plot.label + "_" + cut.label + ".png");
	      //canvas->SaveAs(saveDir + "/" + cut.label + "/" + plot.label + "_" + cut.label + ".pdf");
	    }
	  delete canvas, legend, stack, pvt;
	}
      ++j;
    }
  
  
  // //Get cov. matrix
  // unsigned ii = 0;
  // TH1D *hist_uni;
  // TH1D *hist_reco;
  // for (auto const& cut: nuescatter_cuts){
  //   for(auto const& plot : recoPlots){
  //     for(auto const& category : nuescat_sel_categories){
  //       if (cut.name == "one_shw" && category.name == "Signal"){
  //         if (plot.name == "reco_E"){
  //           hist_reco = sNu[ii]->ToTH1(gPOT);
  //         }
  //       }
  //     ++ii;}
  //   }
  // }
  // for (auto const& syst: systs){
  //   if (syst.name == "nuescat"){
  //     hist_uni = syst.syst->Universe(0).ToTH1(gPOT,syst.color);
  //   }
  // }
  // const std::vector<TH1D*> binSets = {hist_reco,hist_uni};
  // std::unique_ptr<TMatrixD> covmat_eng = CalcCovMx(binSets);

  // TFile fcovmat("covmat_recoE.root", "RECREATE");
  // fcovmat.WriteObject(&covmat_eng, "covmat_reco_uni"); 
  // fcovmat.Close();
  
  //TFile fin(state_fname.c_str());
  //EnsembleSpectrum* sNuSyst = LoadFrom<EnsembleSpectrum>(fin.GetDirectory("nue_scat")).release();

  TCanvas *canvas = new TCanvas("c");
  TLegend *legend = new TLegend(.69,.57,.89,.95);
  double max_val = -9999;
  for (auto const& syst : systs){
    TH1D *hist_nom = syst.syst->Nominal().ToTH1(gPOT, syst.color);
    for(unsigned int i = 0; i < syst.syst->NUniverses(); ++i){
      TH1D *hist = syst.syst->Universe(i).ToTH1(gPOT, syst.color-10);
      if (hist->GetMaximum() > max_val){max_val = hist->GetMaximum();}
      hist->SetLineColorAlpha(syst.color, 0.1);
      //hist->SetMaximum(40000);
      hist->Draw("hist same");
    }
    //Redraw nominal over top
    //hist_nom->SetMaximum(40000);
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

  //delete legend,canvas;

  // new TCanvas;
  // TLegend *legend_flux = new TLegend(.59,.67,.89,.85);
  // max_val = -9999;
  // for (auto const& syst : systs_flux){
  //   TH1D *hist_nom = syst.syst->Nominal().ToTH1(gPOT, syst.color);
  //   for(unsigned int i = 0; i < syst.syst->NUniverses(); ++i){
  //     TH1D *hist = syst.syst->Universe(i).ToTH1(gPOT, syst.color-10);
  //     if (hist->GetMaximum() > max_val){max_val = hist->GetMaximum();}
  //     hist->SetLineColorAlpha(syst.color, 0.1);
  //     hist->SetMaximum(500);
  //     hist->Draw("hist same");
  //   }
  //   //Redraw nominal over top
  //   hist_nom->SetMaximum(500);
  //   hist_nom->Draw("hist same");

  //   legend_flux->AddEntry(hist_nom, syst.label, "lf");
  // }
  // legend_flux->SetLineColorAlpha(0,0);
	// legend_flux->SetTextSize(0.04);
  // legend_flux->Draw();
  // gSystem->Exec("mkdir -p "+saveDir+"/systs");
  //gPad->Print(saveDir + "/systs/all_flux.png");
  
  //delete legend,canvas;

  
  //TCanvas *canvas = new TCanvas("c");
  //TLegend *legend = new TLegend(.59,.67,.89,.85);

  // for (auto const& syst: systs_flux){
  //   new TCanvas;
  //   //std::cout<<"In loop"<<std::endl;
  //   TGraphAsymmErrors* band = syst.syst->ErrorBand(gPOT,"error_band_flux_"+syst.name,stateDir+"/");
  //   DrawErrorBand(syst.syst->Nominal().ToTH1(gPOT, syst.color), band);
  //   //std::cout<<"Drew error band"<<std::endl;
  //   gPad->Print(saveDir + "/systs/error_band_flux_"+syst.name+".png");
  // }

  for (auto const& syst: systs){
    new TCanvas;
    //std::cout<<"In loop"<<std::endl;
    TGraphAsymmErrors* band = syst.syst->ErrorBand(gPOT,"error_band_"+syst.name,stateDir+"/");
    DrawErrorBand(syst.syst->Nominal().ToTH1(gPOT, syst.color), band);
    //std::cout<<"Drew error band"<<std::endl;
    gPad->Print(saveDir + "/systs/error_band_"+syst.name+".png");
  }

  // sNuSyst.Nominal().ToTH1(gPOT, kRed)->Draw("hist");
  // //sNcpi0Syst.Nominal().ToTH1(gPOT, kRed)->Draw("hist");
  // for(unsigned int i = 0; i < sNuSyst.NUniverses(); ++i){
  //   sNuSyst.Universe(i).ToTH1(gPOT, kRed-10)->Draw("hist same");
  // }
  // //Redraw nominal over top
  // sNuSyst.Nominal().ToTH1(gPOT, kRed)->Draw("hist same");
  // gSystem->Exec("mkdir -p "+saveDir+"/systs");
  // gPad->Print(saveDir + "/systs/all.png");

  // new TCanvas;

  // sNueScatSyst.Nominal().ToTH1(gPOT, kOrange)->Draw("hist");
  // for(unsigned int i = 0; i < sNueScatSyst.NUniverses(); ++i){
  //   sNueScatSyst.Universe(i).ToTH1(gPOT, kOrange-10)->Draw("hist same");
  // }
  // //Redraw nominal over top
  // sNueScatSyst.Nominal().ToTH1(gPOT, kOrange)->Draw("hist same");
  // legend->AddEntry(sNueScatSyst.Nominal().ToTH1(gPOT, kOrange), Form("%s (%'.0f)", , integral), "lf");


  // gPad->Print(saveDir + "/systs/nuescat.png");

  // new TCanvas;

  // TGraphAsymmErrors* band = sNuSyst.ErrorBand(gPOT,"error_band_all",stateDir+"/");
  // DrawErrorBand(sNuSyst.Nominal().ToTH1(gPOT, kRed), band);

  // gPad->Print(saveDir + "/systs/all_band.png");

  // new TCanvas;

  // TGraphAsymmErrors* band_nue = sNueScatSyst.ErrorBand(gPOT,"error_band_nuescat",stateDir+"/");
  // DrawErrorBand(sNuSyst.Nominal().ToTH1(gPOT, kOrange), band_nue);

  // gPad->Print(saveDir + "/systs/nuescat_band.png");

}