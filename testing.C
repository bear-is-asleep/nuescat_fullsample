// #include "sbnana/CAFAna/Core/Spectrum.h"
// #include "sbnana/CAFAna/Core/SpectrumLoader.h"
// //#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
// #include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

//using namespace ana;

// #include "Constants.h"
// #include "Structs.h"
// #include "TrueEventCategories.h"
// #include "NuEScatterRecoVars.h"
// #include "NuEScatterCuts.h"
// #include "utils.h"

#include <string>
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TSystem.h"

//#include "sbnana/CAFAna/Core/Utilities.h"
//#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace std;

TH2D CalcCovMat(std::vector<TH1D*> hists){
  int nbins = hists[0]->GetNbinsX();
  assert(nbins == hists[1]->GetNbinsX());

  TH2D covmat = new TH2D()
  for (unsigned hist_id = 0; hist_id < hists.size(); hist_id++){
    for (unsigned i = 0; i < nbins; i++){

    }
  }
}

void testing(const std::string inputName = "/pnfs/sbnd/persistent/sbndpro/mcp/mc/workshop/SBNWorkshop0421/prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc/v09_19_00_01/caf/flat_caf_0-9f00feff-e742-419d-9856-9fe7428b93a9.root")
{
  TH1* h = new TH1D(
    /* name */ "h2",
    /* title */ "Hist with constant bin width",
    /* X-dimension */ 100, 0.0, 4.0);
  h->Fill(2,3);
  h->Fill(1,4);
  h->Fill(0.5,3);
  h->Fill(3.5,6);
  h->Draw();

  cout<<h->GetMinimumBin()<<endl;

  // const std::vector<TH1D*> binSets = {tarr,tarr};
  // std::unique_ptr<TMatrixD> covmat_eng = CalcCovMx(binSets);

  // TFile fcovmat("testing_covmat.root", "RECREATE");
  // fcovmat.WriteObject(&covmat_eng, "covmat_reco_uni"); 
  // fcovmat.Close();





}