#include "TVector3.h"

const SpillCut kHasSlc([](const caf::SRSpillProxy* sp) {
    // if (sp->nslc == 0){
    //   std::cout<<sp->hdr.run<<","<<sp->hdr.subrun<<","<<sp->hdr.evt<<std::endl;
    //   std::cout<<sp->mc.nu.position.x<<","<<sp->mc.nu.position.y<<","<<sp->mc.nu.position.z<<","<<std::endl;
    // }
    return sp->nslc != 0;
  });

const SpillCut kHasNuSlc([](const caf::SRSpillProxy* sp) {
    for(auto const& slc : sp->slc)
      if(!slc.is_clear_cosmic) return true;

    return false;
  });

const SpillCut kHasNuFVSlc([](const caf::SRSpillProxy* sp) {
    for(auto const& slc : sp->slc)
      {
	if(slc.is_clear_cosmic) continue;

	if(PtInVolAbsX(slc.vertex, fvndNuEScat)) return true;
      }
    return false;
  });

const SpillCut kHasCRUMBSSlc([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];

    if(isnan(slc.crumbs_result.score)) return false;

    return slc.crumbs_result.score > 0;
  });

const SpillCut kIsFV([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];
    
    return PtInVolAbsX(slc.vertex, fvndNuEScat);
  });

const SpillCut kIsTrueFV([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];
    //std::cout<<"run,subrun,evt : "<<sp->hdr.run<<","<<sp->hdr.subrun<<","<<sp->hdr.evt<<std::endl;
    //print_all_prim_matched(sp,kBestSlcID(sp));
    //print_all_prim_info_slc(sp,kBestSlcID(sp));
    return PtInVolAbsX(slc.truth.position, fvndNuEScat);
  });

const SpillCut kIsTrueAV([](const caf::SRSpillProxy* sp) {
  auto const& slc = sp->slc[kBestSlcID(sp)];
  //std::cout<<"run,subrun,evt : "<<sp->hdr.run<<","<<sp->hdr.subrun<<","<<sp->hdr.evt<<std::endl;
  //print_all_prim_matched(sp,kBestSlcID(sp));
  //print_all_prim_info_slc(sp,kBestSlcID(sp));
  return PtInVolAbsX(slc.truth.position, avnd);
});

const SpillCut kEtheta2([](const caf::SRSpillProxy* sp) {
  // --Return true if a single object has Etheta < Ethetacut
  // --Future want to redesign this so that we sum over fragmented reco objects

  auto const& slc = sp->slc[kBestSlcID(sp)];
  TVector3 nu_direction = kNuDir(&slc);
  //Showers
  for (auto const& shw : slc.reco.shw){
    TVector3 shw_dir(shw.dir.x,shw.dir.y,shw.dir.z);
    double shw_theta = acos(shw_dir.Dot(nu_direction));
    //pdg = abs(shw.razzle.pdg); //Use razzle pdg
    //m=pdgmass.at(pdg);//mass using table
    double Eng = shw.bestplane_energy; //ke
    if (Eng*shw_theta*shw_theta < kEtheta2Cut && Eng != 0){
      return true;
    }
  }
  //Tracks
  for (auto const& trk : slc.reco.trk){
    TVector3 trk_dir(trk.dir.x,trk.dir.y,trk.dir.z);
    double trk_theta = acos(trk_dir.Dot(nu_direction));
    //pdg = abs(trk.dazzle.pdg); //Use dazzle pdg
    //m=pdgmass.at(pdg);//mass using table
    double Eng = trk.calo[trk.bestplane].ke*1e-3; //ke
    if (Eng*trk_theta*trk_theta < kEtheta2Cut && Eng != 0){ 
      return true;
    }
  }
  return false;
});

const SpillCut kHasNoTrks([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return slc.reco.ntrk == 0;
  });

const SpillCut kHasOneShw([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return slc.reco.nshw == 1;
  });

const SpillCut kTooManyRecoObjects([](const caf::SRSpillProxy* sp) {
  auto const& slc = sp->slc[kBestSlcID(sp)];

  if (slc.reco.nshw + slc.reco.ntrk <= 5){
    if (slc.reco.ntrk <= 3 && slc.reco.nshw <= 4){ //additional condition on track cut based on pure nue study - loses 3 events
      return true;
    }
  }
  return false;

});

const SpillCut kSoftRecoAngleCut([](const caf::SRSpillProxy* sp) {
  return kRecoSmallestTheta(sp) < 0.6; //Small angle requirement 
});

const SpillCut kSoftEthetaCut([](const caf::SRSpillProxy* sp) {
  return kEtheta2Var(sp) < 0.05; //Small Etheta requirement
});


const SpillCut kRazzleCut([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];

    for (auto const& shw : slc.reco.shw){
      if (shw.razzle.electronScore <= kRazzleCutVal){
        return false;
      }
    }
    return true;
  });

const SpillCut kTruthIDNuEScat([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];
    return slc.truth.genie_inttype == 1098;
});

const SpillCut kTruthEleCut([](const caf::SRSpillProxy* sp) {
    auto const& slc = sp->slc[kBestSlcID(sp)];
    int ne = 0;
    for (auto const& prim : slc.truth.prim){
      if (abs(prim.pdg) == 12 || abs(prim.pdg) == 14) continue; //Don't count these
      if (abs(prim.pdg) == 11) {ne++;}
    }
    return ne==1;
  });

const SpillCut kTruthShwCut([](const caf::SRSpillProxy *sp) {
  int nshw = 0;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  for (auto const& prim : slc.truth.prim){
    int pdg = prim.pdg;
    if (abs(pdg) == 11 || pdg == 22){++nshw;}
    if (pdg == 111){nshw+=2;}  
  }
  return nshw==1;
});

const SpillCut kTruthTrkCut([](const caf::SRSpillProxy *sp) {
  int ntrk = 0;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  for (auto const& prim : slc.truth.prim){
    int pdg = prim.pdg;
    if (abs(pdg) == 2212){++ntrk;}
    else if (abs(pdg) == 13){++ntrk;}
    else if (abs(pdg) == 211){++ntrk;}   
  }
  return ntrk==0;
});

const SpillCut kTruthEtheta2([](const caf::SRSpillProxy* sp) {
  // --Return true if a single object has Etheta < Ethetacut

  return abs(kTruthEtheta2Var(sp)) < kEtheta2Cut; //Use abs incase the mask value is wrong
});

const SpillCut kTruthSlcGenieType([](const caf::SRSpillProxy* sp) {
  auto const& slc = sp->slc[kBestSlcID(sp)];
  return slc.truth.genie_inttype == 1098;
});

std::vector<CutDef> truth_int_type_cuts = { { "No Cut", "no_cut", kNoSpillCut },
  { "Has Slc", "has_slc", kHasSlc },
  {"Genie Int Type Slice","int_type",kTruthSlcGenieType},
  { "Is FV", "is_fv", kIsTrueFV },
};


std::vector<CutDef> truth_cuts = { { "No Cut", "no_cut", kNoSpillCut },
  { "Has Slc", "has_slc", kHasSlc },
  { "Is FV", "is_fv", kIsTrueFV },
  {"No Tracks","no_trk",kTruthTrkCut},
  {"One Shower","one_shw",kTruthShwCut},
  {"Etheta2 Truth Cut","Etheta2",kTruthEtheta2},
  {"One Electron","one_ele",kTruthEleCut},
  //{ "Is NuEScat", "is_nue_scat", kTruthIDNuEScat}
};

std::vector<CutDef> original_cuts = { { "No Cut", "no_cut", kNoSpillCut },
  { "Has Slc", "has_slc", kHasSlc },
  { "Has Nu Slc", "has_nu_slc", kHasNuSlc },
  { "Has Nu FV Slc", "has_nu_fv_slc", kHasNuFVSlc },
  //{ "Has CRUMBS Slc", "has_crumbs_slc", kHasCRUMBSSlc },
  { "Is FV", "is_fv", kIsFV },
  {"Has One Shower","one_shw",kHasOneShw},
  {"Has No Tracks","no_trks",kHasNoTrks},
  {"Electron Razzle","erazzle",kRazzleCut},
  {"Etheta2 Cut","Etheta2",kEtheta2},
};

std::vector<CutDef> nuescatter_cuts = { { "No Cut", "no_cut", kNoSpillCut },
  { "Has Slc", "has_slc", kHasSlc },
  { "Has Nu Slc", "has_nu_slc", kHasNuSlc },
  { "Has Nu FV Slc", "has_nu_fv_slc", kHasNuFVSlc },
  { "Has CRUMBS Slc", "has_crumbs_slc", kHasCRUMBSSlc },
  { "Is FV", "is_fv", kIsFV },
  {"Has One Shower","one_shw",kHasOneShw},
  {"Has No Tracks","no_trks",kHasNoTrks},
  {"Electron Razzle","erazzle",kRazzleCut},
  {"Etheta2 Cut","Etheta2",kEtheta2},
};

std::vector<CutDef> Etheta_cuts = { { "No Cut", "no_cut", kNoSpillCut },
  { "Has Slc", "has_slc", kHasSlc },
  { "Has Nu Slc", "has_nu_slc", kHasNuSlc },
  { "Has Nu FV Slc", "has_nu_fv_slc", kHasNuFVSlc },
  { "Has CRUMBS Slc", "has_crumbs_slc", kHasCRUMBSSlc },
  { "Is FV", "is_fv", kIsFV },
  {"Etheta2 Cut","Etheta2",kEtheta2},
};

std::vector<CutDef> no_cuts = { { "No Cut", "no_cut", kNoSpillCut }
};

std::vector<CutDef> preselection_cuts = { { "No Cut", "no_cut", kNoSpillCut },
  { "Has Slc", "has_slc", kHasSlc },
  { "Has Nu Slc", "has_nu_slc", kHasNuSlc },
  { "Has Nu FV Slc", "has_nu_fv_slc", kHasNuFVSlc },
};


const SpillCut kPreSelection = kHasSlc && kHasNuSlc && kHasNuFVSlc;
const SpillCut kCosmicRej    = kPreSelection && kHasCRUMBSSlc && kIsFV;
const SpillCut kNuFluxSelection = kHasSlc && kHasNuSlc && kTrueAV && !kCosmicSpill;
const SpillCut kEthetaSelection = kPreSelection && kCosmicRej && kEtheta2;
const SpillCut kSoftSelection = kPreSelection && kIsFV && kTooManyRecoObjects && kSoftEthetaCut;
const SpillCut kFullSelection = kEthetaSelection && kHasNoTrks && kHasOneShw && kRazzleCut;
const SpillCut kTrueSelection = kHasSlc && kIsTrueFV && kTruthTrkCut && kTruthShwCut && kTruthEleCut && kTruthEtheta2;
