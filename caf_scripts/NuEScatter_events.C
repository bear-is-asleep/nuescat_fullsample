#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Core/Var.h"

//Systs
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

using namespace ana;

#include "Constants.h"
#include "Structs.h"
#include "utils.h"
#include "TrueEventCategories.h"
#include "NuEScatterRecoVars.h"
#include "NuEScatterCuts.h"
#include "plotStyle.C"

#include <string>
#include <iostream>
#include <fstream>
#include "TVector.h"
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
const std::string inputNameNuE_new = "/sbnd/data/users/brindenc/analyze_sbnd/nue/v09_54_00/CAFnue_full.root";

//const string surName = "fullsample_softetheta_recoslc";
const string surName = "pure_nue_truefv_recoslc";
const TString stateDir = "/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/"+get_date()+"_"+surName;

int kEventType(const caf::SRSpillProxy* sp){
  if (kNuEScat(sp) && kTrueAV(sp)){
    return 0; //"NuEScat";
  }
  else if (kNCPiZero(sp) && kTrueAV(sp)){
    return 1; //"NCpi0";
  }
  else if (kNC(sp) && ! kNCPiZero(sp) && kTrueAV(sp)){
    return 2; //"NC";
  }
  else if (kCCNuMu(sp) && kTrueAV(sp)){
    return 3; //"CCNuMu";
  }
  else if (kCCNuE(sp) && kTrueAV(sp)){
    return 4; //"CCNuE";
  }
  else if (kDirt(sp)){
    return 5; //"Dirt";
  }
  else if (kCosmicSpill(sp)){
    return 6; //Cosmic;
  }
  else{
    return 7; //"Other";
  }
}

int kReserveSpace = 2000000; //2,000,000

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

  vector<int> evt_type;

  vector<double> Etheta;
  vector<double> true_Etheta;

  vector<int> truenshw;
  vector<int> truentrk;

  vector<double> true_vtxx;
  vector<double> true_vtxy;
  vector<double> true_vtxz;

  vector<double> reco_vtxx;
  vector<double> reco_vtxy;
  vector<double> reco_vtxz;

  vector<double> true_slice_eng;
  vector<double> reco_eng;
  vector<double>  reco_theta;
  vector<double> true_theta;
  vector<double> crumbs_score;

  vector<double> lshw_eng;
  vector<double> lshw_dedx;
  vector<double> lshw_cnvgap;
  vector<double> lshw_dens;
  vector<double> lshw_len;
  vector<double> lshw_openangle;
  vector<double> lshw_razzle_electron;
  vector<double> lshw_razzle_photon;
  vector<double> lshw_razzle_other;
  vector<double> lshw_start_x;
  vector<double> lshw_start_y;
  vector<double> lshw_start_z;
  vector<double> lshw_true_pdg;
  vector<double> lshw_true_angle;
  vector<double> lshw_true_len;
  vector<double> lshw_angle;

  vector<double> slshw_eng;
  vector<double> slshw_dedx;
  vector<double> slshw_cnvgap;
  vector<double> slshw_dens;
  vector<double> slshw_len;
  vector<double> slshw_openangle;
  vector<double> slshw_razzle_electron;
  vector<double> slshw_razzle_photon;
  vector<double> slshw_razzle_other;
  vector<double> slshw_start_x;
  vector<double> slshw_start_y;
  vector<double> slshw_start_z;
  vector<double> slshw_true_pdg;
  vector<double> slshw_true_angle;
  vector<double> slshw_true_len;
  vector<double> slshw_angle;

  vector <double> ltrk_eng;
  vector <double> ltrk_len;
  vector <double> ltrk_npts;
  vector <double> ltrk_pfptrkscore;
  vector <double> ltrk_muonscore;
  vector <double> ltrk_pionscore;
  vector <double> ltrk_protonscore;
  vector<double> ltrk_start_x;
  vector<double> ltrk_start_y;
  vector<double> ltrk_start_z;
  vector<double> ltrk_true_pdg;
  vector<double> ltrk_true_angle;
  vector<double> ltrk_true_len;
  vector<double> ltrk_angle;
  vector<double> ltrk_crttrk_time;
  vector<double> ltrk_crttrk_angle;
  vector<double> ltrk_dedx;

  vector <double> sltrk_eng;
  vector <double> sltrk_len;
  vector <double> sltrk_npts;
  vector <double> sltrk_pfptrkscore;
  vector <double> sltrk_muonscore;
  vector <double> sltrk_pionscore;
  vector <double> sltrk_protonscore;
  vector<double> sltrk_start_x;
  vector<double> sltrk_start_y;
  vector<double> sltrk_start_z;
  vector<double> sltrk_true_pdg;
  vector<double> sltrk_true_angle;
  vector<double> sltrk_true_len;
  vector<double> sltrk_angle;
  vector<double> sltrk_crttrk_time;
  vector<double> sltrk_crttrk_angle;
  vector<double> sltrk_dedx;

  vector<int> genie_inttype;
  vector<int> genie_mode;

  vector<double> fmatch_time;
  vector<double> fmatch_score;


  //REserve space
  run.reserve(kReserveSpace);
  subrun.reserve(kReserveSpace);
  evt.reserve(kReserveSpace);

  nshw.reserve(kReserveSpace);
  ntrk.reserve(kReserveSpace);
  nstub.reserve(kReserveSpace);
  nslc.reserve(kReserveSpace);
  nele.reserve(kReserveSpace);
  nrele.reserve(kReserveSpace);
  nrph.reserve(kReserveSpace);

  evt_type.reserve(kReserveSpace);

  Etheta.reserve(kReserveSpace);
  true_Etheta.reserve(kReserveSpace);

  truenshw.reserve(kReserveSpace);
  truentrk.reserve(kReserveSpace);

  true_vtxx.reserve(kReserveSpace);
  true_vtxy.reserve(kReserveSpace);
  true_vtxz.reserve(kReserveSpace);

  reco_vtxx.reserve(kReserveSpace);
  reco_vtxy.reserve(kReserveSpace);
  reco_vtxz.reserve(kReserveSpace);
  crumbs_score.reserve(kReserveSpace);

  true_slice_eng.reserve(kReserveSpace);
  reco_eng.reserve(kReserveSpace);
  reco_theta.reserve(kReserveSpace);
  true_theta.reserve(kReserveSpace);

  lshw_eng.reserve(kReserveSpace);
  lshw_dedx.reserve(kReserveSpace);
  lshw_cnvgap.reserve(kReserveSpace);
  lshw_dens.reserve(kReserveSpace);
  lshw_len.reserve(kReserveSpace);
  lshw_openangle.reserve(kReserveSpace);
  lshw_razzle_electron.reserve(kReserveSpace);
  lshw_razzle_photon.reserve(kReserveSpace);
  lshw_razzle_other.reserve(kReserveSpace);
  lshw_start_x.reserve(kReserveSpace);
  lshw_start_y.reserve(kReserveSpace);
  lshw_start_z.reserve(kReserveSpace);
  lshw_true_pdg.reserve(kReserveSpace);
  lshw_true_angle.reserve(kReserveSpace);
  lshw_true_len.reserve(kReserveSpace);
  lshw_angle.reserve(kReserveSpace);

  slshw_eng.reserve(kReserveSpace);
  slshw_dedx.reserve(kReserveSpace);
  slshw_cnvgap.reserve(kReserveSpace);
  slshw_dens.reserve(kReserveSpace);
  slshw_len.reserve(kReserveSpace);
  slshw_openangle.reserve(kReserveSpace);
  slshw_razzle_electron.reserve(kReserveSpace);
  slshw_razzle_photon.reserve(kReserveSpace);
  slshw_razzle_other.reserve(kReserveSpace);
  slshw_start_x.reserve(kReserveSpace);
  slshw_start_y.reserve(kReserveSpace);
  slshw_start_z.reserve(kReserveSpace);
  slshw_true_pdg.reserve(kReserveSpace);
  slshw_true_angle.reserve(kReserveSpace);
  slshw_true_len.reserve(kReserveSpace);
  slshw_angle.reserve(kReserveSpace);

  ltrk_eng.reserve(kReserveSpace);
  ltrk_len.reserve(kReserveSpace);
  ltrk_npts.reserve(kReserveSpace);
  ltrk_pfptrkscore.reserve(kReserveSpace);
  ltrk_muonscore.reserve(kReserveSpace);
  ltrk_pionscore.reserve(kReserveSpace);
  ltrk_protonscore.reserve(kReserveSpace);
  ltrk_start_x.reserve(kReserveSpace);
  ltrk_start_y.reserve(kReserveSpace);
  ltrk_start_z.reserve(kReserveSpace);
  ltrk_true_pdg.reserve(kReserveSpace);
  ltrk_true_angle.reserve(kReserveSpace);
  ltrk_true_len.reserve(kReserveSpace);
  ltrk_angle.reserve(kReserveSpace);
  ltrk_crttrk_time.reserve(kReserveSpace);
  ltrk_crttrk_angle.reserve(kReserveSpace);
  ltrk_dedx.reserve(kReserveSpace);

  sltrk_eng.reserve(kReserveSpace);
  sltrk_len.reserve(kReserveSpace);
  sltrk_npts.reserve(kReserveSpace);
  sltrk_pfptrkscore.reserve(kReserveSpace);
  sltrk_muonscore.reserve(kReserveSpace);
  sltrk_pionscore.reserve(kReserveSpace);
  sltrk_protonscore.reserve(kReserveSpace);
  sltrk_start_x.reserve(kReserveSpace);
  sltrk_start_y.reserve(kReserveSpace);
  sltrk_start_z.reserve(kReserveSpace);
  sltrk_true_pdg.reserve(kReserveSpace);
  sltrk_true_angle.reserve(kReserveSpace);
  sltrk_true_len.reserve(kReserveSpace);
  sltrk_angle.reserve(kReserveSpace);
  sltrk_crttrk_time.reserve(kReserveSpace);
  sltrk_crttrk_angle.reserve(kReserveSpace);
  sltrk_dedx.reserve(kReserveSpace);

  genie_inttype.reserve(kReserveSpace);
  genie_mode.reserve(kReserveSpace);
  fmatch_time.reserve(kReserveSpace);
  fmatch_score.reserve(kReserveSpace);


  //Make TTree
  TFile* file = new TFile(stateDir+"/cut_events_basic.root", "RECREATE");
  TTree *tree = new TTree("rectree","CAF Tree");

  TFile* file2 = new TFile(stateDir+"/cut_events_shw.root", "RECREATE");
  TTree *tree2 = new TTree("rectree","CAF Tree");

  TFile* file3 = new TFile(stateDir+"/cut_events_trk.root", "RECREATE");
  TTree *tree3 = new TTree("rectree","CAF Tree");

  tree->Branch("run",&run);
  tree->Branch("subrun",&subrun);
  tree->Branch("evt",&evt);

  tree2->Branch("run",&run);
  tree2->Branch("subrun",&subrun);
  tree2->Branch("evt",&evt);

  tree3->Branch("run",&run);
  tree3->Branch("subrun",&subrun);
  tree3->Branch("evt",&evt);

  tree->Branch("nshw",&nshw);
  tree->Branch("ntrk",&ntrk);
  tree->Branch("truenshw",&truenshw);
  tree->Branch("truentrk",&truentrk);
  tree->Branch("nstub",&nstub);
  tree->Branch("nslc",&nslc);
  tree->Branch("nele",&nele);
  //tree->Branch("razzle.electrons",&nrele);
  //tree->Branch("razzle.photons",&nrph);
  tree->Branch("Etheta",&Etheta);
  tree->Branch("true_Etheta",&true_Etheta);
  tree->Branch("evt_type",&evt_type);
  //tree->Branch("shw.conversion_gap",&shw_conversion_gap);

  tree->Branch("true_vtx.x",&true_vtxx);
  tree->Branch("true_vtx.y",&true_vtxy);
  tree->Branch("true_vtx.z",&true_vtxz);

  tree->Branch("reco_vtx.x",&reco_vtxx);
  tree->Branch("reco_vtx.y",&reco_vtxy);
  tree->Branch("reco_vtx.z",&reco_vtxz);

  tree->Branch("true_slice_eng",&true_slice_eng);
  tree->Branch("reco_eng",&reco_eng);
  tree->Branch("reco_theta",&reco_theta);
  tree->Branch("true_theta",&true_theta);
  tree->Branch("fmatch.time",&fmatch_time);
  tree->Branch("fmatch.score",&fmatch_score);
  tree->Branch("crumbs.score",&crumbs_score);

  tree2->Branch("lshw.eng",&lshw_eng);
  tree2->Branch("lshw.dedx",&lshw_dedx);
  tree2->Branch("lshw.cnvgap",&lshw_cnvgap);
  tree2->Branch("lshw.dens",&lshw_dens);
  tree2->Branch("lshw.len",&lshw_len);
  tree2->Branch("lshw.openangle",&lshw_openangle);
  tree2->Branch("lshw.electron",&lshw_razzle_electron);
  tree2->Branch("lshw.photon",&lshw_razzle_photon);
  tree2->Branch("lshw.other",&lshw_razzle_other);
  tree2->Branch("lshw.start.x",&lshw_start_x);
  tree2->Branch("lshw.start.y",&lshw_start_y);
  tree2->Branch("lshw.start.z",&lshw_start_z);
  tree2->Branch("lshw.true.pdg",&lshw_true_pdg);
  tree2->Branch("lshw.true.angle",&lshw_true_angle);
  tree2->Branch("lshw.true.len",&lshw_true_len);
  tree2->Branch("lshw.angle",&lshw_angle);

  tree2->Branch("slshw.eng",&slshw_eng);
  tree2->Branch("slshw.dedx",&slshw_dedx);
  tree2->Branch("slshw.cnvgap",&slshw_cnvgap);
  tree2->Branch("slshw.dens",&slshw_dens);
  tree2->Branch("slshw.len",&slshw_len);
  tree2->Branch("slshw.openangle",&slshw_openangle);
  tree2->Branch("slshw.electron",&slshw_razzle_electron);
  tree2->Branch("slshw.photon",&slshw_razzle_photon);
  tree2->Branch("slshw.other",&slshw_razzle_other);
  tree2->Branch("slshw.start.x",&slshw_start_x);
  tree2->Branch("slshw.start.y",&slshw_start_y);
  tree2->Branch("slshw.start.z",&slshw_start_z);
  tree2->Branch("slshw.true.pdg",&slshw_true_pdg);
  tree2->Branch("slshw.true.angle",&slshw_true_angle);
  tree2->Branch("slshw.true.len",&slshw_true_len);
  tree2->Branch("slshw.angle",&slshw_angle);
  

  tree->Branch("genie_inttype",&genie_inttype);
  tree->Branch("genie_mode",&genie_mode);

  tree3->Branch("ltrk.eng",&ltrk_eng);
  tree3->Branch("ltrk.len",&ltrk_len);
  tree3->Branch("ltrk.npts",&ltrk_npts);
  tree3->Branch("ltrk.pfptrkscore",&ltrk_pfptrkscore);
  tree3->Branch("ltrk.muonscore",&ltrk_muonscore);
  tree3->Branch("ltrk.pionscore",&ltrk_pionscore);
  tree3->Branch("ltrk.protonscore",&ltrk_protonscore);
  tree3->Branch("ltrk.start.x",&ltrk_start_x);
  tree3->Branch("ltrk.start.y",&ltrk_start_y);
  tree3->Branch("ltrk.start.z",&ltrk_start_z);
  tree3->Branch("ltrk.true.pdg",&ltrk_true_pdg);
  tree3->Branch("ltrk.true.angle",&ltrk_true_angle);
  tree3->Branch("ltrk.true.len",&ltrk_true_len);
  tree3->Branch("ltrk.angle",&ltrk_angle);
  tree3->Branch("ltrk.crttrk.time",&ltrk_crttrk_time);
  tree3->Branch("ltrk.crttrk.angle",&ltrk_crttrk_angle);
  tree3->Branch("ltrk.dedx",&ltrk_dedx);

  tree3->Branch("sltrk.eng",&sltrk_eng);
  tree3->Branch("sltrk.len",&sltrk_len);
  tree3->Branch("sltrk.npts",&sltrk_npts);
  tree3->Branch("sltrk.pfptrkscore",&sltrk_pfptrkscore);
  tree3->Branch("sltrk.muonscore",&sltrk_muonscore);
  tree3->Branch("sltrk.pionscore",&sltrk_pionscore);
  tree3->Branch("sltrk.protonscore",&sltrk_protonscore);
  tree3->Branch("sltrk.start.x",&sltrk_start_x);
  tree3->Branch("sltrk.start.y",&sltrk_start_y);
  tree3->Branch("sltrk.start.z",&sltrk_start_z);
  tree3->Branch("sltrk.true.pdg",&sltrk_true_pdg);
  tree3->Branch("sltrk.true.angle",&sltrk_true_angle);
  tree3->Branch("sltrk.true.len",&sltrk_true_len);
  tree3->Branch("sltrk.angle",&sltrk_angle);
  tree3->Branch("sltrk.crttrk.time",&sltrk_crttrk_time);
  tree3->Branch("sltrk.crttrk.angle",&sltrk_crttrk_angle);
  tree3->Branch("sltrk.dedx",&sltrk_dedx);

  SpectrumLoader loader(inputNameNuE_new);

  gSystem->Exec("mkdir -p " + stateDir);
  ofstream out(stateDir+"/selected_events.txt");
  //out <<"Run\t Subrun\t Event\t SpillID\n";
  const SpillVar dummy_var([&](const caf::SRSpillProxy* sp){
    //for(const auto& slc: sp->slc) {
      if(true) {
        out << sp->hdr.run << "\t" << sp->hdr.subrun << "\t" << sp->hdr.evt <<"\n";
            // << "\t" << kBestSlcID(sp) <<"\n";
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
        true_Etheta.push_back(kTruthEtheta2Var(sp));
        evt_type.push_back(kEventType(sp));
        truenshw.push_back(kTruthShw(sp));
        truentrk.push_back(kTruthTrk(sp));
        true_vtxx.push_back(kNuX(sp));
        true_vtxy.push_back(kNuY(sp));
        true_vtxz.push_back(kNuZ(sp));

        reco_vtxx.push_back(kRecoVtxX(sp));
        reco_vtxy.push_back(kRecoVtxY(sp));
        reco_vtxz.push_back(kRecoVtxZ(sp));

        true_slice_eng.push_back(kTrueSliceEnergy(sp));
        reco_eng.push_back(kRecoE(sp));
        reco_theta.push_back(kRecoSmallestTheta(sp));
        true_theta.push_back(kTrueSmallestTheta(sp));
        crumbs_score.push_back(kCRUMBSScore(sp));

        lshw_eng.push_back(kLeadingShwEnergy(sp));
        lshw_dedx.push_back(kLeadingShwdEdx(sp));
        lshw_cnvgap.push_back(kLeadingShwCnvGap(sp));
        lshw_dens.push_back(kLeadingShwDensity(sp));
        lshw_len.push_back(kLeadingShwLen(sp));
        lshw_openangle.push_back(kLeadingShwOpenAngle(sp));
        lshw_razzle_electron.push_back(kLeadingShwRazzleElectronScore(sp));
        lshw_razzle_photon.push_back(kLeadingShwRazzlePhotonScore(sp));
        lshw_razzle_other.push_back(kLeadingShwRazzleOtherScore(sp));
        lshw_start_x.push_back(kLeadingShwStartX(sp));
        lshw_start_y.push_back(kLeadingShwStartY(sp));
        lshw_start_z.push_back(kLeadingShwStartZ(sp));
        lshw_true_pdg.push_back(kLeadingShwTruePdg(sp));
        lshw_true_angle.push_back(kLeadingShwTrueAngle(sp));
        lshw_true_len.push_back(kLeadingShwTrueLength(sp));
        lshw_angle.push_back(kLeadingShwAngle(sp));

        slshw_eng.push_back(kSubLeadingShwEnergy(sp));
        slshw_dedx.push_back(kSubLeadingShwdEdx(sp));
        slshw_cnvgap.push_back(kSubLeadingShwCnvGap(sp));
        slshw_dens.push_back(kSubLeadingShwDensity(sp));
        slshw_len.push_back(kSubLeadingShwLen(sp));
        slshw_openangle.push_back(kSubLeadingShwOpenAngle(sp));
        slshw_razzle_electron.push_back(kSubLeadingShwRazzleElectronScore(sp));
        slshw_razzle_photon.push_back(kSubLeadingShwRazzlePhotonScore(sp));
        slshw_razzle_other.push_back(kSubLeadingShwRazzleOtherScore(sp));
        slshw_start_x.push_back(kSubLeadingShwStartX(sp));
        slshw_start_y.push_back(kSubLeadingShwStartY(sp));
        slshw_start_z.push_back(kSubLeadingShwStartZ(sp));
        slshw_true_pdg.push_back(kSubLeadingShwTruePdg(sp));
        slshw_true_angle.push_back(kSubLeadingShwTrueAngle(sp));
        slshw_true_len.push_back(kSubLeadingShwTrueLength(sp));
        slshw_angle.push_back(kSubLeadingShwAngle(sp));

        ltrk_eng.push_back(kLeadingTrkEnergy(sp));
        ltrk_len.push_back(kLeadingTrkLen(sp));
        ltrk_npts.push_back(kLeadingTrkNPts(sp));
        ltrk_pfptrkscore.push_back(kLeadingTrkPFPTrkScore(sp));
        ltrk_muonscore.push_back(kLeadingTrkDazzleMuonScore(sp));
        ltrk_pionscore.push_back(kLeadingTrkDazzlePionScore(sp));
        ltrk_protonscore.push_back(kLeadingTrkDazzleProtonScore(sp));
        ltrk_start_x.push_back(kLeadingTrkStartX(sp));
        ltrk_start_y.push_back(kLeadingTrkStartY(sp));
        ltrk_start_z.push_back(kLeadingTrkStartZ(sp));
        ltrk_true_pdg.push_back(kLeadingTrkTruePdg(sp));
        ltrk_true_angle.push_back(kLeadingTrkTrueAngle(sp));
        ltrk_true_len.push_back(kLeadingTrkTrueLength(sp));
        ltrk_angle.push_back(kLeadingTrkAngle(sp));
        ltrk_crttrk_time.push_back(kLeadingTrkCRTTrkTime(sp));
        ltrk_crttrk_angle.push_back(kLeadingTrkCRTTrkAngle(sp));
        ltrk_dedx.push_back(kLeadingTrkdEdx(sp));

        sltrk_eng.push_back(kSubLeadingTrkEnergy(sp));
        sltrk_len.push_back(kSubLeadingTrkLen(sp));
        sltrk_npts.push_back(kSubLeadingTrkNPts(sp));
        sltrk_pfptrkscore.push_back(kSubLeadingTrkPFPTrkScore(sp));
        sltrk_muonscore.push_back(kSubLeadingTrkDazzleMuonScore(sp));
        sltrk_pionscore.push_back(kSubLeadingTrkDazzlePionScore(sp));
        sltrk_protonscore.push_back(kSubLeadingTrkDazzleProtonScore(sp));
        sltrk_start_x.push_back(kSubLeadingTrkStartX(sp));
        sltrk_start_y.push_back(kSubLeadingTrkStartY(sp));
        sltrk_start_z.push_back(kSubLeadingTrkStartZ(sp));
        sltrk_true_pdg.push_back(kSubLeadingTrkTruePdg(sp));
        sltrk_true_angle.push_back(kSubLeadingTrkTrueAngle(sp));
        sltrk_true_len.push_back(kSubLeadingTrkTrueLength(sp));
        sltrk_angle.push_back(kSubLeadingTrkAngle(sp));
        sltrk_crttrk_time.push_back(kSubLeadingTrkCRTTrkTime(sp));
        sltrk_crttrk_angle.push_back(kSubLeadingTrkCRTTrkAngle(sp));
        sltrk_dedx.push_back(kSubLeadingTrkdEdx(sp));

        genie_inttype.push_back(kGenieType(sp));
        genie_mode.push_back(kGenieMode(sp));
        fmatch_time.push_back(kFMatchTime(sp));
        fmatch_score.push_back(kFMatchScore(sp));
        return 1;
      }
    //} 
    return 0;
  });

  const Binning dummy_bins = Binning::Simple(2, 0, 2);
  Spectrum dummy_spec("Dummy Label", dummy_bins, loader, dummy_var, kNoSpillCut);
  //Spectrum dummy_spec_univ("Dummy weighted",binsEnergy,loader,kTrueE,kNoCut,weis[0]);
  //EnsembleSpectrum dummy_spec_univ(loader,ax,kPreSelection,kNoCut,weis);
  //std::cout.setstate(std::ios_base::failbit);
  loader.Go();
  std::cout.clear();
  
  tree->Fill();
  file->Write();
  file->Close();

  tree2->Fill();
  file2->Write();

  tree3->Fill();
  file3->Write();

  
  //gSystem->Exec("mv " + stateDir+"/cut_events.root"+" "+stateDir);
  gSystem->Exec("cp NuEScatter_events.C NuEScatterRecoVars.h NuEScatterCuts.h TrueEventCategories.h Constants.h Structs.h " + stateDir);
}
