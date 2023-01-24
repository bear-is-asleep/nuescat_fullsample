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
    return 2; //"NC";
  }
  else if (kCCNuMu(sp)){
    return 3; //"CCNuMu";
  }
  else if (kCCNuE(sp)){
    return 4; //"CCNuE";
  }
  else if (kDirt(sp)){
    return 5; //"Dirt";
  }
  else{
    return 6; //"Other";
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
  vector<int> nele;
  vector<int> nrele;
  vector<int> nrph;

  //vector<double> shw_conversion_gap;

  vector<int> evt_type;

  vector<double> Etheta;

  vector<int> truenshw;
  vector<int> truentrk;

  vector<double> vtxx;
  vector<double> vtxy;
  vector<double> vtxz;

  vector<double> true_spill_eng;
  vector<double> reco_eng;
  vector<double>  reco_theta;
  vector<double> true_theta;

  vector<double> lshw_eng;
  vector<double> lshw_dedx;
  vector<double> lshw_cnvgap;
  vector<double> lshw_dens;
  vector<double> lshw_len;
  vector<double> lshw_openangle;

  vector<double> slshw_eng;
  vector<double> slshw_dedx;
  vector<double> slshw_cnvgap;
  vector<double> slshw_dens;
  vector<double> slshw_len;
  vector<double> slshw_openangle;

  vector<int> genie_inttype;
  vector<int> genie_mode;


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
  tree->Branch("nele",&nele);
  tree->Branch("razzle.electrons",&nrele);
  tree->Branch("razzle.photons",&nrph);
  tree->Branch("Etheta",&Etheta);
  tree->Branch("evt_type",&evt_type);
  //tree->Branch("shw.conversion_gap",&shw_conversion_gap);

  tree->Branch("vtx.x",&vtxx);
  tree->Branch("vtx.y",&vtxy);
  tree->Branch("vtx.z",&vtxz);

  tree->Branch("true_spill_eng",&true_spill_eng);
  tree->Branch("reco_eng",&reco_eng);
  tree->Branch("reco_theta",&reco_theta);
  tree->Branch("true_theta",&true_theta);

  tree->Branch("lshw.eng",&lshw_eng);
  tree->Branch("lshw.dedx",&lshw_dedx);
  tree->Branch("lshw.cnvgap",&lshw_cnvgap);
  tree->Branch("lshw.dens",&lshw_dens);
  tree->Branch("lshw.len",&lshw_len);
  tree->Branch("lshw.openangle",&lshw_openangle);

  tree->Branch("slshw.eng",&slshw_eng);
  tree->Branch("slshw.dedx",&slshw_dedx);
  tree->Branch("slshw.cnvgap",&slshw_cnvgap);
  tree->Branch("slshw.dens",&slshw_dens);
  tree->Branch("slshw.len",&slshw_len);
  tree->Branch("slshw.openangle",&slshw_openangle);

  tree->Branch("genie_inttype",&genie_inttype);
  tree->Branch("genie_mode",&genie_mode);




  const double gPOT = 10e20;
  SpectrumLoader loader(inputNameNu);

  gSystem->Exec("mkdir -p " + stateDir);
  ofstream out(stateDir+"/selected_events.txt");
  out <<"Run\t Subrun\t Event\t SpillID\n";
  const SpillVar dummy_var([&](const caf::SRSpillProxy* sp){
    //for(const auto& slc: sp->slc) {
      if(kEthetaSelection(sp)) {
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
        nele.push_back(kNElectrons(sp));
        nrele.push_back(kNRazzleElectrons(sp));
        nrph.push_back(kNRazzlePhotons(sp));
        Etheta.push_back(kEtheta2Var(sp));
        evt_type.push_back(kEventType(sp));
        truenshw.push_back(kTruthShw(sp));
        truentrk.push_back(kTruthTrk(sp));
        vtxx.push_back(kNuX(sp));
        vtxy.push_back(kNuY(sp));
        vtxz.push_back(kNuZ(sp));

        true_spill_eng.push_back(kTrueSliceE(sp));
        reco_eng.push_back(kRecoE(sp));
        reco_theta.push_back(kRecoSmallestTheta(sp));
        true_theta.push_back(kTrueSmallestTheta(sp));

        lshw_eng.push_back(kLeadingShwEnergy(sp));
        lshw_dedx.push_back(kLeadingShwdEdx(sp));
        lshw_cnvgap.push_back(kLeadingShwCnvGap(sp));
        lshw_dens.push_back(kLeadingShwDensity(sp));
        lshw_len.push_back(kLeadingShwLen(sp));
        lshw_openangle.push_back(kLeadingShwOpenAngle(sp));

        slshw_eng.push_back(kSubLeadingShwEnergy(sp));
        slshw_dedx.push_back(kSubLeadingShwdEdx(sp));
        slshw_cnvgap.push_back(kSubLeadingShwCnvGap(sp));
        slshw_dens.push_back(kSubLeadingShwDensity(sp));
        slshw_len.push_back(kSubLeadingShwLen(sp));
        slshw_openangle.push_back(kSubLeadingShwOpenAngle(sp));

        genie_inttype.push_back(sp->slc[kBestSlcID(sp)].truth.genie_inttype);
        genie_mode.push_back(sp->slc[kBestSlcID(sp)].truth.genie_mode);
        return 1;
      }
    //} 
    return 0;
  });

  const Binning dummy_bins = Binning::Simple(2, 0, 2);
  Spectrum dummy_spec("Dummy Label", dummy_bins, loader, dummy_var, kNoSpillCut);
  //Spectrum dummy_spec_univ("Dummy weighted",binsEnergy,loader,kTrueE,kNoCut,weis[0]);
  //EnsembleSpectrum dummy_spec_univ(loader,ax,kPreSelection,kNoCut,weis);

  loader.Go();
  tree->Fill();
  file->Write();
  gSystem->Exec("mv " + stateDir+"/cut_events.root"+" "+stateDir);

}
