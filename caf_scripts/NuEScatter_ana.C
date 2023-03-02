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
#include "NuEScatter_debug.h"
#include "TrueEventCategories.h"
#include "NuEScatterRecoVars.h"
#include "NuEScatterCuts.h"
#include "utils.h"
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

const std::vector<CutDef> cuts = truth_int_type_cuts;
const std::vector<Plot<SpillVar>> plots = recoPlots_all;
const std::vector<TrueCategory> categories = no_cosmic_sel;

void NuEScatter_ana(bool reload = true)
{

  selectionstyle();
  setlocale(LC_NUMERIC, "");
  
  //Use my sample with nu e scat
  const std::string inputNameNuE = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_58_02/CAFnue1_full.root";
  const std::string inputNameNuE_new = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_54_00/CAFnue_full.root";
  const std::string inputNameNuECC = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_58_02/CAFnuecc1_full.root";
  const std::string inputNameNu = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_flat_caf_sbnd";
  const std::string inputNameNu_noflat = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_caf_sbnd";


  const double gPOT = 10e20;
  const bool save = true;
  const string surName = "fullsample_inttype_cuts_tmatch_slice";
  const TString saveDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/"+get_date()+"_"+surName;
  const TString stateDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/"+get_date()+"_"+surName;

  SpectrumLoader loaderNu(inputNameNu);
  //SpectrumLoader loaderIntime(inputNameIntime);

  std::vector<Spectrum*> sNu;

  SpillCut previousCuts = kNoSpillCut;

  for (auto const& cut : cuts){
    for(auto const& plot : plots){
      for(auto const& category : categories){
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
  }
  //std::cout.setstate(std::ios_base::failbit);
  loaderNu.Go();
  std::cout.clear();
  //loaderIntime.Go();

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

  for(auto const& cut : cuts)
    {

    gSystem->Exec("mkdir -p " + saveDir + "/" + cut.label);
    for(auto const& plot : plots)
	    {

        TCanvas *canvas = new TCanvas("c " + plot.name + " " + cut.name, 
              "c " + plot.name + " " + cut.name);
        canvas->cd();
        TLegend *legend = new TLegend(plot.legBox[0], plot.legBox[1], plot.legBox[2], plot.legBox[3]);
        THStack *stack = new THStack("stack_" + plot.name + "_" + cut.name, plot.axes_labels);
	  
	      double selectedSig = 0., selected = 0.;

        for(auto const& category : categories)
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
	  //pvt->Draw();

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + cut.label + "/" + plot.label + "_" + cut.label + ".png");
	      stack->SetMinimum(1);
        gPad->SetLogy(1); 
        canvas->SaveAs(saveDir + "/" + cut.label + "/" + plot.label + "_" + cut.label + "_logy.png");
	    }
	  delete canvas, legend, stack, pvt;
	}
      ++j;
    }
  gSystem->Exec("cp NuEScatter_events.C NuEScatterRecoVars.h NuEScatterCuts.h TrueEventCategories.h Constants.h Structs.h " + saveDir);
}