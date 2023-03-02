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

const SpillCut kEtheta2([](const caf::SRSpillProxy* sp) {
  // --Return true if a single object has Etheta < Ethetacut

  double theta = 0;
  double ke = 0;
  double m = 0;
  double Eng = 0;
  int pdg = 0;
  int bestplane = -1;
  double Etheta2 = 999;
  //std::cout<<"nslc: "<<sp->nslc<<std::endl;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  int cnt = 0;
  //Showers
  for (auto const& shw : slc.reco.shw){
    cnt++;
    theta = acos(shw.dir.z); //Longitudinal angle
    bestplane = shw.bestplane; //Best plane
    pdg = abs(shw.razzle.pdg); //Use razzle pdg
    if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          m=0;
        }
    else{
      m = pdgmass.at(pdg); //mass
    }
    ke = shw.bestplane_energy; //ke
    Eng = sqrt(m*m+ke*ke); //Energy
    Etheta2 = Eng*theta*theta;
    if (Etheta2 < kEtheta2Cut && Eng != 0){ 
      if (kNuEScat(sp)){
        // std::cout<<"passed shw "<<shw.dir.z<<","<<theta<<","<<Eng<<","<<Etheta2<<","<<cnt<<std::endl;
        // std::cout<<"n shw "<< slc.reco.nshw<< std::endl;
        // std::cout<<"n trk "<< slc.reco.ntrk<< std::endl;
        // std::cout<<"n stub "<< slc.reco.nstub<< std::endl;
      }
      // std::cout<<"passed "<<shw.dir.z<<","<<theta<<","<<Eng<<","<<Etheta2<<std::endl;
      // std::cout<<"nue scat?  "<< kNuEScat(sp)<< std::endl;
      // std::cout<<"n shw "<< slc.reco.nshw<< std::endl;
      // std::cout<<"n trk "<< slc.reco.ntrk<< std::endl;
      // std::cout<<"n stub "<< slc.reco.nstub<< std::endl;
      return true;
    }
    else if(kNuEScat(sp)){
      TVector3 mom(shw.truth.p.genp.x,shw.truth.p.genp.y,shw.truth.p.genp.z);
      double true_angle = acos(mom.Z()/mom.Mag());
      // std::cout<<mom.X()<<","<<mom.Y()<<","<<mom.Z()<<std::endl;
      // std::cout<<"failed shw "<<shw.dir.z<<","<<theta<<",true angle "<<true_angle<<Eng<<","<<Etheta2<<","<<cnt<<std::endl;
      // std::cout<<"n shw "<< slc.reco.nshw<< std::endl;
      // std::cout<<"n trk "<< slc.reco.ntrk<< std::endl;
      // std::cout<<"n stub "<< slc.reco.nstub<< std::endl;
    }
  }
  //Tracks
  for (auto const& trk : slc.reco.trk){
    cnt++;
    theta = acos(trk.dir.z); //Longitudinal angle
    bestplane = trk.bestplane; //Best plane
    pdg = abs(trk.dazzle.pdg); //Use dazzle pdg
    ke = trk.calo[bestplane].ke*1e-3; //ke
    m = pdgmass.at(pdg); //mass
    Eng = sqrt(m*m+ke*ke); //Energy
    Etheta2 = Eng*theta*theta;
    if (Etheta2 < kEtheta2Cut && Eng != 0){ 
      if (kNuEScat(sp)){
        // std::cout<<"passed trk "<<trk.dir.z<<","<<theta<<","<<Eng<<","<<Etheta2<<","<<cnt<<std::endl;

        // std::cout<<"n shw "<< slc.reco.nshw<< std::endl;
        // std::cout<<"n trk "<< slc.reco.ntrk<< std::endl;
        // std::cout<<"n stub "<< slc.reco.nstub<< std::endl;
      }
      // std::cout<<"passed "<<trk.dir.z<<","<<theta<<","<<Eng<<","<<Etheta2<<std::endl;
      // std::cout<<"nue scat?  "<< kNuEScat(sp)<< std::endl;
      // std::cout<<"n shw "<< slc.reco.nshw<< std::endl;
      // std::cout<<"n trk "<< slc.reco.ntrk<< std::endl;
      // std::cout<<"n stub "<< slc.reco.nstub<< std::endl;
      return true;
    }
    else if(kNuEScat(sp)){
      TVector3 mom(trk.truth.p.genp.x,trk.truth.p.genp.y,trk.truth.p.genp.z);
      double true_angle = acos(mom.Z()/mom.Mag());
      // std::cout<<mom.X()<<","<<mom.Y()<<","<<mom.Z()<<std::endl;
      // std::cout<<"failed trk "<<trk.dir.z<<","<<theta<<",true angle "<<true_angle<<Eng<<","<<Etheta2<<","<<cnt<<std::endl;
      // std::cout<<"n shw "<< slc.reco.nshw<< std::endl;
      // std::cout<<"n trk "<< slc.reco.ntrk<< std::endl;
      // std::cout<<"n stub "<< slc.reco.nstub<< std::endl;
    }
  }
  //std::cout<<Etheta2<<std::endl;
  //std::cout<<"cut val"<<kEtheta2Cut<<std::endl;
  //Stubs -- not yet implemented
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
  return kRecoSmallestTheta(sp) < 0.65; //Small angle requirement 
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

  double theta = 0;
  double ke = 0;
  double m = 0;
  double Eng = 0;
  int pdg = 0;
  int bestplane = -1;
  double Etheta2 = 999;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  for (auto const& prim : slc.truth.prim){
    TVector3 pvec(prim.genp.x,prim.genp.y,prim.genp.z);
    double p = pvec.Mag();
    theta = acos(prim.genp.z/p); //Longitudinal angle
    pdg = abs(prim.pdg); //Use razzle pdg
    Eng = prim.genE; //Energy
    Etheta2 = Eng*theta*theta;
    if (Etheta2 < kEtheta2Cut && Eng != 0){ 
      return true;
    }
  }
  return false;
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
  {"Etheta2 Cut","Etheta2",kEtheta2},
  {"Has One Shower","one_shw",kHasOneShw},
  {"Has No Tracks","no_trks",kHasNoTrks},
  {"Electron Razzle","erazzle",kRazzleCut},
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
const SpillCut kEthetaSelection = kPreSelection && kCosmicRej && kEtheta2;
const SpillCut kSoftSelection = kPreSelection && kIsFV && kTooManyRecoObjects && kSoftRecoAngleCut;
const SpillCut kFullSelection = kEthetaSelection && kHasNoTrks && kHasOneShw && kRazzleCut;
const SpillCut kTrueSelection = kHasSlc && kIsTrueFV && kTruthTrkCut && kTruthShwCut && kTruthEleCut && kTruthEtheta2;
