#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"


//Systs
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

using namespace ana;

#include "Constants.h"
#include "Structs.h"
#include "TrueEventCategories.h"
#include "NuEScatterTruthVars.h"
#include "NuEScatterRecoVars.h"
#include "NuEScatterCuts.h"
#include "utils.h"
#include "Reducer.h"
#include "plotStyle.C"

#include <string>
#include <iostream>
#include <fstream>
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TFile.h"

using namespace std;

//Use my sample with nu e scat
const string inputNameNuE = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_58_02/CAFnue1_full.root";
const string inputNameNuECC = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_58_02/CAFnuecc1_full.root";
const string inputNameNu = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_flat_caf_sbnd";
const string inputNameNu_noflat = "defname: official_MCP2022A_prodoverlay_corsika_cosmics_proton_genie_rockbox_sce_reco2_concat_caf_sbnd";

const TString stateDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/"+get_date();

int kEventType(const caf::SRSpillProxy* sp){
  if (kNuEScat(sp)){
    return 0; //"NuEScat";
  }
  else if (kNCPiZero(sp)){
    return 1; //"NCpi0";
  }
  else if (kNC(sp) && ! kNCPiZero(sp)){
    return 2 //"NC";
  }
  else if (kCCNuMu(sp)){
    return 3 //"CCNuMu";
  }
  else if (kCCNuE(sp)){
    return 4 //"CCNuE";
  }
  else if (kDirt(sp)){
    return 5 //"Dirt";
  }
  else{
    return 6 //"Other";
  }
}

void NuEScatter_events()
{
  vector<int> run;
  vector<int> subrun;
  vector<int> evt;

  vector<int> nshw;
  vector<int> ntrk;
  vector<int> nstub;
  vector<int> nslc;

  vector<int> evt_type;

  vector<double> Etheta;

  vector<int> truenshw;
  vector<int> truentrk;

  vector<double> vtxx;
  vector<double> vtxy;
  vector<double> vtxz;


  //Make TTree
  TFile* file = new TFile(stateDir+"/cut_events.root", "RECREATE");
  TTree *tree = new TTree("rectree","CAF Tree");

  tree->Branch("run",&run);
  tree->Branch("subrun",&subrun);
  tree->Branch("evt",&evt);
  tree->Branch("nshw",&nshw);
  tree->Branch("ntrk",&ntrk);
  tree->Branch("truenshw",&nshw);
  tree->Branch("truentrk",&ntrk);
  tree->Branch("nstub",&nstub);
  tree->Branch("nslc",&nslc);
  tree->Branch("Etheta",&Etheta);
  tree->Branch("evt_type",&evt_type);

  tree->Branch("vtx.x",&vtxx);
  tree->Branch("vtx.y",&vtxy);
  tree->Branch("vtx.z",&vtxz);



  const double gPOT = 10e20;
  SpectrumLoader loader(inputNameNu);

  gSystem->Exec("mkdir -p " + stateDir);
  ofstream out(stateDir+"/selected_events.txt");
  out <<"Run\t Subrun\t Event\t SpillID\n";
  const SpillVar dummy_var([&](const caf::SRSpillProxy* sp){
    //for(const auto& slc: sp->slc) {
      if(kFullSelection(sp)) {
        out << sp->hdr.run << "\t" << sp->hdr.subrun
            << "\t" << sp->hdr.evt << "\t" << kBestSlcID(sp) <<"\n";
            //<< "\tVertex: (" 
            //<< slc.vertex.x << ", " << slc.vertex.y << ", " << slc.vertex.z << ")\n";
        run.push_back(sp->hdr.run);
        subrun.push_back(sp->hdr.subrun);
        evt.push_back(sp->hdr.evt);
        nshw.push_back(kNShowers(sp));
        ntrk.push_back(kNTracks(sp));
        nstub.push_back(kNStubs(sp));
        nslc.push_back(sp->nslc);
        Etheta.push_back(kEtheta2Var(sp));
        evt_type.push_back(kEventType(sp));
        truenshw.push_back(kTruthShw(sp));
        truentrk.push_back(kTruthTrk(sp));
        vtxx.push_back(kNuX(sp));
        vtxy.push_back(kNuY(sp));
        vtxz.push_back(kNuZ(sp));
        return 1;
      }
    //} 
    return 0;
  });

  const Binning dummy_bins = Binning::Simple(2, 0, 2);
  Spectrum dummy_spec("Dummy Label", dummy_bins, loader, dummy_var, kNoSpillCut);

  loader.Go();
  tree->Fill();
  file->Write();

}
