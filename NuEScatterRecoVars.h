#include <math.h>

const SpillVar kBestSlcID([](const caf::SRSpillProxy* sp) -> unsigned {
    unsigned i = 0;
    unsigned returnID = 0;
    double Etheta2 = 999;
    double theta = 0;
    double ke = 0;
    double m = 0;
    double Eng = 0;
    int pdg = 0;
    int bestplane = -1;
    for(auto const& slc : sp->slc)
    {
      if(slc.is_clear_cosmic || isnan(slc.crumbs_result.score) || !PtInVolAbsX(slc.vertex, fvndNuEScat)) 
        { ++i; continue; }
      
      //Showers
      for (auto const& shw : slc.reco.shw){
        pdg = abs(shw.razzle.pdg); //Use razzle pdg
        if (pdg == 5 || pdg == 0){ //Why is this sometimes -5? - default value in rzzle id
          {continue; }
        }
        theta = acos(shw.dir.z); //Longitudinal angle
        bestplane = shw.bestplane; //Best plane
        ke = shw.bestplane_energy; //ke
        m = pdgmass.at(pdg); //mass
        Eng = sqrt(m*m+ke*ke); //Energy
        //std::cout<<pdg<<","<<ke<<std::endl;
        if (Eng*theta*theta < Etheta2){
          Etheta2 = Eng*theta*theta;
          returnID = i;
        }
      }
      //Tracks
      for (auto const& trk : slc.reco.trk){
        pdg = abs(trk.dazzle.pdg); //Use dazzle pdg
        if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          {continue; }
        }
        theta = acos(trk.dir.z); //Longitudinal angle
        bestplane = trk.bestplane; //Best plane
        ke = trk.calo[bestplane].ke*1e-3; //ke
        m = pdgmass.at(pdg); //mass
        Eng = sqrt(m*m+ke*ke); //Energy
        //std::cout<<pdg<<","<<ke<<std::endl;
        if (Eng*theta*theta < Etheta2){ 
          // Second condition is to make sure if we have a shower with lower threshold, we'll keep that one 
          Etheta2 = Eng*theta*theta;
          returnID = i;
        }
      }
      //std::cout<<Etheta2<<std::endl;
      //Stubs - not yet implemented
      // for (int j = 0; j++; j<=slc.reco.nstub){
      //   theta = acos(slc.reco.stub.dir.z); //Longitudinal angle
      //   bestplane = slc.reco.trk.bestplane; //Best plane
      //   pdg = abs(slc.reco.shw.dazzle.pdg); //Use dazzle pdg
      //   ke = slc.reco.tk.calo[bestplane].ke; //ke
      //   m = pdgmass.at(pdg); //mass
      //   Eng = sqrt(m*m+ke*ke); //Energy
      //   if (Eng*theta*theta < Etheta2 && Etheta2 > Etheta2Cut){ 
      //     // Second condition is to make sure if we have a shower with lower threshold, we'll keep that one 
      //     Etheta2 = Eng*theta*theta;
      //     returnID = i;
      //   }
      // }

	    ++i;
    }
    return returnID;
  });

// const SpillVar kTrueE([](const caf::SRSpillProxy *sp) ->double {
//     return sp->slc[kBestSlcID(sp)].truth.E;
//   });

const SpillVar kFlashTrigVar([](const caf::SRSpillProxy *sp) ->unsigned {
    return sp->pass_flashtrig ? 1 : 0;
  });

const SpillVar kNSlices([](const caf::SRSpillProxy* sp) -> unsigned {
    return sp->nslc;
  });

const SpillVar kNNuSlices([](const caf::SRSpillProxy* sp) -> unsigned {
    unsigned i = 0;

    for(auto const& slc : sp->slc)
      {
	if(!slc.is_clear_cosmic) ++i;
      }

    return i;
  });

const SpillVar kNFVSlices([](const caf::SRSpillProxy* sp) -> unsigned {
    unsigned i = 0;

    for(auto const& slc : sp->slc)
      {
	if(slc.is_clear_cosmic) continue;

	if(PtInVolAbsX(slc.vertex, fvndNuEScat)) ++i;
      }

    return i;
  });

const SpillVar kCRUMBSScore([](const caf::SRSpillProxy* sp) -> double {
    if(sp->nslc==0) return -1.5;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    if(slc.is_clear_cosmic || isnan(slc.crumbs_result.score)) return -1.5;
    return slc.crumbs_result.score;
  });

const SpillVar kIsFVVar([](const caf::SRSpillProxy* sp) -> double {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return PtInVolAbsX(slc.vertex, fvndNuEScat) ? 1 : 0;
  });

const SpillVar kNTracks([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return slc.reco.ntrk;
  });

const SpillVar kNRazzleElectrons([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    unsigned i = 0;

    for(auto const& shw : slc.reco.shw)
      {
	if(shw.razzle.pdg == 11) ++i;
      }

    return i;
  });

const SpillVar kNRazzlePhotons([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    unsigned i = 0;

    for(auto const& shw : slc.reco.shw)
      {
	if(shw.razzle.pdg == 22) ++i;
      }

    return i;
  });

const SpillVar kNStubs([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return slc.reco.nstub;
  });

const SpillVar kNShowers([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return slc.reco.nshw;
  });

const SpillVar kNPrims([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return slc.truth.nprim;
  });

const SpillVar kNElectrons([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];
    int ne = 0;
    for (auto const& prim : slc.truth.prim){
      if (abs(prim.pdg) == 12 || abs(prim.pdg) == 14) continue; //Don't count these
      if (abs(prim.pdg) == 11) {ne++;}
    }
    return ne;
  });

const SpillVar kEtheta2Var([](const caf::SRSpillProxy* sp) -> double {

    if(sp->nslc==0) return 999;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    double Etheta2 = 999.;
    double theta = 999.;
    double ke = 0;
    double m = 0;
    double Eng = 0;
    int pdg = 0;
    int bestplane = -1;

    //Showers
    for (auto const& shw : slc.reco.shw){
      theta = acos(shw.dir.z); //Longitudinal angle
      bestplane = shw.bestplane; //Best plane
      pdg = abs(shw.razzle.pdg); //Use razzle pdg
      if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          {continue; }
      }
      ke = shw.bestplane_energy; //ke
      m = pdgmass.at(pdg); //mass
      Eng = sqrt(m*m+ke*ke); //Energy
      if (Eng*theta*theta < Etheta2 && Eng != 0){
        Etheta2 = Eng*theta*theta;
      }
    }
    //Tracks
    for (auto const& trk : slc.reco.trk){
      theta = acos(trk.dir.z); //Longitudinal angle
      bestplane = trk.bestplane; //Best plane
      pdg = abs(trk.dazzle.pdg); //Use dazzle pdg
      if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          {continue; }
      }
      ke = trk.calo[bestplane].ke*1e-3; //ke
      m = pdgmass.at(pdg); //mass
      Eng = sqrt(m*m+ke*ke); //Energy
      if (Eng*theta*theta < Etheta2 && Eng != 0){ 
        // Second condition is to make sure if we have a shower with lower threshold, we'll keep that one 
        Etheta2 = Eng*theta*theta;
      }
    }
    //std::cout<<theta<<","<<Eng<<","<<Etheta2<<std::endl;
    //Stubs - not yet implemented
    return Etheta2;
    //Return lowest Etheta value 
  });

const SpillVar kLeadingShwID([](const caf::SRSpillProxy* sp) -> int {
    if(kNShowers(sp) == 0) return -1;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    std::vector<std::pair<int, double> > list;
    unsigned i = 0;

    for(auto const& shw : slc.reco.shw)
      {
	list.push_back(std::pair<int, double>(i, shw.bestplane_energy));
	++i;
      }
    
    std::sort(list.begin(), list.end(), [](const std::pair<int, double> &x,
					   const std::pair<int, double> &y)
	      { return x.second > y.second; });
    
    return list[0].first;
  });

const SpillVar kSubLeadingShwID([](const caf::SRSpillProxy* sp) -> int {
    if(kNShowers(sp) < 2) return -1;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    std::vector<std::pair<int, double> > list;
    unsigned i = 0;

    for(auto const& shw : slc.reco.shw)
      {
	list.push_back(std::pair<int, double>(i, shw.bestplane_energy));
	++i;
      }
    
    std::sort(list.begin(), list.end(), [](const std::pair<int, double> &x,
					   const std::pair<int, double> &y)
	      { return x.second > y.second; });
    
    return list[1].first;
  });

const SpillVar kLeadingShwEnergy([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.bestplane_energy;
  });

const SpillVar kLeadingShwdEdx([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.bestplane_dEdx;
  });
const SpillVar kLeadingShwCnvGap([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.conversion_gap;
  });
const SpillVar kLeadingShwDensity([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.density;
  });
const SpillVar kLeadingShwLen([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.len;
  });

const SpillVar kLeadingShwOpenAngle([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return TMath::RadToDeg() * shw.open_angle;
  });

const SpillVar kLeadingShwDensityGrad([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.selVars.densityGradient;
  });

const SpillVar kLeadingShwDensityGradPower([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.selVars.densityGradientPower;
  });

const SpillVar kLeadingShwTrkLength([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.selVars.trackLength;
  });

const SpillVar kLeadingShwTrkWidth([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];
    
    return shw.selVars.trackWidth;
  });
const SpillVar kLeadingShwRazzleElectronScore([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -.5;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];

    return shw.razzle.electronScore;
  });

const SpillVar kLeadingShwRazzlePhotonScore([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -.5;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];

    return shw.razzle.photonScore;
  });

const SpillVar kLeadingShwRazzleOtherScore([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -.5;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)];

    return shw.razzle.otherScore;
  });
const SpillVar kSubLeadingShwEnergy([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.bestplane_energy;
  });

const SpillVar kSubLeadingShwdEdx([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.bestplane_dEdx;
  });

const SpillVar kSubLeadingShwCnvGap([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.conversion_gap;
  });
const SpillVar kSubLeadingShwDensity([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.density;
  });

const SpillVar kSubLeadingShwLen([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.len;
  });

const SpillVar kSubLeadingShwOpenAngle([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return TMath::RadToDeg() * shw.open_angle;
  });

const SpillVar kSubLeadingShwCosmicDist([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.cosmicDist;
  });

const SpillVar kSubLeadingShwDensityGrad([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.selVars.densityGradient;
  });

const SpillVar kSubLeadingShwDensityGradPower([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.selVars.densityGradientPower;
  });
const SpillVar kSubLeadingShwTrkLength([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.selVars.trackLength;
  });

const SpillVar kSubLeadingShwTrkWidth([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -1.;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];
    
    return shw.selVars.trackWidth;
  });
const SpillVar kSubLeadingShwRazzleElectronScore([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -.5;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];

    return shw.razzle.electronScore;
  });

const SpillVar kSubLeadingShwRazzlePhotonScore([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -.5;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];

    return shw.razzle.photonScore;
  });

const SpillVar kSubLeadingShwRazzleOtherScore([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -.5;
    auto const& shw = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)];

    return shw.razzle.otherScore;
  });

const SpillVar kRecoE([](const caf::SRSpillProxy* sp) -> double {
    if(sp->nslc==0) return -1;
    auto const& slc = sp->slc[kBestSlcID(sp)];
    double Eng = 0;

    //Showers
    for (auto const& shw : slc.reco.shw){
      int bestplane = shw.bestplane; //Best plane
      int pdg = abs(shw.razzle.pdg); //Use razzle pdg
      if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          {continue; }
      }
      double ke = shw.bestplane_energy; //ke
      double m = pdgmass.at(pdg); //mass
      Eng += sqrt(m*m+ke*ke); //Energy
    }
    //Tracks
    for (auto const& trk : slc.reco.trk){
      int bestplane = trk.bestplane; //Best plane
      int pdg = abs(trk.dazzle.pdg); //Use dazzle pdg
      if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          {continue; }
      }
      double ke = trk.calo[bestplane].ke*1e-3; //ke
      double m = pdgmass.at(pdg); //mass
      Eng += sqrt(m*m+ke*ke); //Energy
    }
    //Implement stubs in future
    return Eng;
  });

const SpillVar kTrueSpillE([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return -1;
  auto const& slc = sp->slc[kBestSlcID(sp)];
    double Eng = 0;

    //Showers
    for (auto const& shw : slc.reco.shw){
      int pdg = abs(shw.truth.p.pdg); //Use razzle pdg
      if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          {continue; }
      }
      Eng += shw.truth.p.genE;
    }
    //Tracks
    for (auto const& trk : slc.reco.trk){
      int pdg = abs(trk.dazzle.pdg); //Use dazzle pdg
      if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
          {continue; }
      }
      Eng += trk.truth.p.genE;
    }
  return Eng;
});

const SpillVar kTrueNuE([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return -1;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  return slc.truth.E;
});

const SpillVar kRecoSmallestTheta([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return 999;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  double E = kRecoE(sp);
  double theta = 999;

  //Showers
  for (auto const& shw : slc.reco.shw){
    int pdg = abs(shw.razzle.pdg); //Use razzle pdg
    if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
        {continue;}
    }
    theta = acos(shw.dir.z); //Longitudinal angle
    if (acos(shw.dir.z) < theta){
      theta = acos(shw.dir.z);
    }
  }
  //Tracks
  for (auto const& trk : slc.reco.trk){
    int pdg = abs(trk.dazzle.pdg); //Use dazzle pdg
    if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
        {continue; }
    }
    if (acos(trk.dir.z)<theta){ 
      theta = acos(trk.dir.z);
    }
  }
  return theta;
});

const SpillVar kTrueSmallestTheta([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return 999;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  double theta = 999;
  for (auto const& prim : slc.truth.prim){
    if (abs(prim.pdg) == 12 || abs(prim.pdg) == 14) continue; //Don't count these
    double px = prim.genp.x;
    double py = prim.genp.y;
    double pz = prim.genp.z;
    double p = sqrt(px*px+py*py+pz*pz);
    if (acos(pz/p)<theta){
      theta = acos(pz/p);
    }
  }
  return theta; //Lowest angle
});

const SpillVar kNuAngle([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return 999;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  double xNu = slc.truth.position.x;
  double yNu = slc.truth.position.y;

  double dx = xNu-kPrismCentroid[0];
  double dy = yNu-kPrismCentroid[1];
  double dz = kDistanceFromBNB;

  return atan(sqrt(dx*dx+dy*dy)/dz)*180/PI;
});

const SpillVar kNuX([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return 9999;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  if (isnan(slc.truth.position.x)){return 9999;}
  return slc.truth.position.x;
});

const SpillVar kNuY([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return 9999;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  if (isnan(slc.truth.position.y)){return 9999;}
  return slc.truth.position.y;
});

const SpillVar kNuZ([](const caf::SRSpillProxy* sp) -> double {
  if(sp->nslc==0) return 9999;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  if (isnan(slc.truth.position.z)){return 9999;}
  return slc.truth.position.z;
});

const SpillVar kTruthShw([](const caf::SRSpillProxy *sp) {
  //if(!kSignal(sp)) return 999999;
  if(sp->nslc==0) return 9999;
  int nshw = 0;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  for (auto const& prim : slc.truth.prim){
    int pdg = prim.pdg;
    if (abs(pdg) == 11 || pdg == 22){++nshw;}
    if (pdg == 111){++nshw;}  
  }
  return nshw;
});

const SpillVar kTruthTrk([](const caf::SRSpillProxy *sp) {
  //if(!kSignal(sp)) return 999999;
  if(sp->nslc==0) return 9999;
  int ntrk = 0;
  auto const& slc = sp->slc[kBestSlcID(sp)];
  for (auto const& prim : slc.truth.prim){
    int pdg = prim.pdg;
    if (abs(pdg) == 2212){++ntrk;}
    else if (abs(pdg) == 13){++ntrk;}
    else if (abs(pdg) == 211){++ntrk;}   
  }
  return ntrk;
});

const SpillVar kNDazzleMuons([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    unsigned i = 0;

    for(auto const& trk : slc.reco.trk)
      {
	if(trk.dazzle.pdg == 13) ++i;
      }

    return i;
  });

const SpillVar kNDazzlePions([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    unsigned i = 0;

    for(auto const& trk : slc.reco.trk)
      {
	if(trk.dazzle.pdg == 211) ++i;
      }

    return i;
  });

const SpillVar kNDazzleProtons([](const caf::SRSpillProxy* sp) -> unsigned {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    unsigned i = 0;

    for(auto const& trk : slc.reco.trk)
      {
	if(trk.dazzle.pdg == 2212) ++i;
      }

    return i;
  });
const SpillVar kSubLeadingShwPDGPlot([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) < 2) return -.5;

    const int pdg = sp->slc[kBestSlcID(sp)].reco.shw[kSubLeadingShwID(sp)].truth.p.pdg;

    switch(pdg) {
    case 11:
      return 0.5;
    case -11:
      return 1.5;
    case 13:
      return 2.5;
    case -13:
      return 3.5;
    case 211:
      return 4.5;
    case -211:
      return 5.5;
    case 2212:
      return 6.5;
    case 22:
      return 7.5;
    default:
      return 8.5;
    }
  });
const SpillVar kLeadingShwPDGPlot([](const caf::SRSpillProxy* sp) -> double {
    if(kNShowers(sp) == 0) return -.5;

    const int pdg = sp->slc[kBestSlcID(sp)].reco.shw[kLeadingShwID(sp)].truth.p.pdg;

    switch(pdg) {
    case 11:
      return 0.5;
    case -11:
      return 1.5;
    case 13:
      return 2.5;
    case -13:
      return 3.5;
    case 211:
      return 4.5;
    case -211:
      return 5.5;
    case 2212:
      return 6.5;
    case 22:
      return 7.5;
    default:
      return 8.5;
    }
  });


const Binning binsPDGShw = Binning::Simple(10,-1,9, {"No Shw", "e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "#pi^{+}", "#pi^{-}", "p", "#gamma", "Other"});

std::vector<Plot<SpillVar>> recoPlots = {//{ "E#theta^{2}", kEtheta2Var, Binning::Simple(15,-10,10),";E#theta^{2};Events", "E_theta",{.59,.57,.89,.85}},
            { "E#theta^{2}", kEtheta2Var, Binning::Simple(10,0.,5.),";E#theta^{2};Events", "E_theta",{.59,.57,.89,.85}},
            { "#theta_{#nu}", kNuAngle, Binning::Simple(10,0.,2.),";#theta_{#nu};Events", "theta_nu",{.59,.57,.89,.85}},
            { "#theta_{#nu} (Zoom)", kNuAngle, Binning::Simple(10,0.,0.2),";#theta_{#nu};Events", "theta_nu_zoom",{.59,.57,.89,.85}},
            { "E#theta^{2} (Zoom)", kEtheta2Var, Binning::Simple(15,0.,0.05),";E#theta^{2};Events", "E_theta_zoom",{.59,.57,.89,.85}},
            { "E#theta^{2} (Zoom x 2)", kEtheta2Var, Binning::Simple(15,0.,0.005),";E#theta^{2};Events", "E_theta_zoom_zoom",{.59,.57,.89,.85}},
            { "N Slices", kNSlices, Binning::Simple(30,0,30), ";nSlices;Events", "n_slices", {.59,.57,.89,.85} },
					  { "N Nu Slices", kNNuSlices, Binning::Simple(10,0,10), ";nNuSlices;Events", "n_nu_slices", {.59,.57,.89,.85} },
					  { "N FV Slices", kNFVSlices, Binning::Simple(10,0,10), ";nFVSlices;Events", "n_fv_slices", {.59,.57,.89,.85} },
					  { "CRUMBS Score", kCRUMBSScore, Binning::Simple(42,-1.5,.6), ";CRUMBS Score;Events", "crumbs_score", {.22,.57,.52,.85} }, 
					  { "FV Cut", kIsFVVar, Binning::Simple(2,0,2), ";FV Cut;Events", "fv_cut", {.59,.57,.89,.85} }, 
					  { "N Showers", kNShowers, Binning::Simple(10,0,10), ";nShowers;Events", "n_showers", {.59,.57,.89,.85} }, 
					  { "N Tracks", kNTracks, Binning::Simple(10,0,10), ";nTracks;Events", "n_tracks", {.59,.57,.89,.85} }, 
					  { "N Stubs", kNStubs, Binning::Simple(10,0,10), ";nStubs;Events", "n_stubs", {.59,.57,.89,.85} }, 
            { "e Razzle", kNRazzleElectrons, Binning::Simple(5,0,5),";Razzle Electrons;Events", "e_razzle",{.59,.57,.89,.85}},
            { "ph Razzle", kNRazzlePhotons, Binning::Simple(5,0,5),";Razzle Photons;Events", "ph_razzle",{.59,.57,.89,.85}},
            { "Reco E", kRecoE, Binning::Simple(30,0,3),";Reco E [GeV];Events", "reco_E",{.59,.57,.89,.85}},
            { "True E_{#nu}", kTrueNuE, Binning::Simple(30,0,3),";True E_{#nu}[GeV];Events", "true_Enu",{.59,.57,.89,.85}},
            { "True E", kTrueSpillE, Binning::Simple(30,0,3),";True E [GeV];Events", "true_E",{.59,.57,.89,.85}},
            { "Reco #theta", kRecoSmallestTheta, Binning::Simple(15,0,1),";Reco #theta ;Events", "reco_theta",{.59,.57,.89,.85}},
            { "Reco #theta (Zoom)", kRecoSmallestTheta, Binning::Simple(15,0,0.1),";Reco #theta;Events", "reco_theta_zoom",{.59,.57,.89,.85}},
            { "True #theta", kTrueSmallestTheta, Binning::Simple(15,0,1),";True #theta ;Events", "true_theta",{.59,.57,.89,.85}},
            { "True #theta (Zoom)", kTrueSmallestTheta, Binning::Simple(15,0,0.1),";True #theta ;Events", "true_theta_zoom",{.59,.57,.89,.85}},
            { "N Prim", kNPrims, Binning::Simple(20,0,20), ";nPrim;Events", "n_prims", {.59,.57,.89,.85} }, 
            { "N Ele", kNElectrons, Binning::Simple(5,0,5), ";nEle;Events", "n_eles", {.59,.57,.89,.85} }, 
            { "Leading Shower Energy", kLeadingShwEnergy, Binning::Simple(40,0,1), ";Leading Shower Energy (GeV);Events", "leading_shw_energy", {.59,.57,.89,.85} },
					  { "Leading Shower dEdx", kLeadingShwdEdx, Binning::Simple(40,0,10), ";Leading Shower dE/dx (MeV/cm);Events", "leading_shw_dedx", {.59,.57,.89,.85} },
					  { "Leading Shower Conversion Gap", kLeadingShwCnvGap, Binning::Simple(40,0,20), ";Leading Shower Conversion Gap (cm);Events", "leading_shw_cnv_gap", {.59,.57,.89,.85} },
					  { "Leading Shower Density", kLeadingShwDensity, Binning::Simple(40,0,20), ";Leading Shower Density (MeV/cm);Events", "leading_shw_density", {.59,.57,.89,.85} },
					  { "Leading Shower Length", kLeadingShwLen, Binning::Simple(40,0,200), ";Leading Shower Length (cm);Events", "leading_shw_len", {.59,.57,.89,.85} },
					  { "Leading Shower Opening Angle", kLeadingShwOpenAngle, Binning::Simple(36,0,90), ";Leading Shower Opening Angle (#circ);Events", "leading_shw_open_angle", {.59,.57,.89,.85} },
					  //{ "Leading Shower Cosmic Dist", kLeadingShwCosmicDist, Binning::Simple(40,0,400), ";Leading Shower Cosmic Dist (cm);Events", "leading_shw_cosmic_dist", {.59,.57,.89,.85} },
					  { "Leading Shower Density Gradient", kLeadingShwDensityGrad, Binning::Simple(40,0,4), ";Leading Shower Density Gradient (MeV/cm);Events", "leading_shw_density_grad", {.59,.57,.89,.85} },
					  { "Leading Shower Density Gradient Power", kLeadingShwDensityGradPower, Binning::Simple(4,0,4), ";Leading Shower Density Gradient Power;Events", "leading_shw_density_grad_power", {.59,.57,.89,.85} },
					  { "Leading Shower Track Stub Length", kLeadingShwTrkLength, Binning::Simple(24,0,24), ";Leading Shower Track Stub Length (cm);Events", "leading_shw_trk_stub_length", {.59,.57,.89,.85} },
					  { "Leading Shower Track Stub Width", kLeadingShwTrkWidth, Binning::Simple(24,0,24), ";Leading Shower Track Stub Width (cm);Events", "leading_shw_trk_stub_width", {.59,.57,.89,.85} },
					  { "Leading Shower PDG", kLeadingShwPDGPlot, binsPDGShw, ";Leading Shower PDG;Events", "leading_shw_pdg",  {.22,.57,.52,.85} },
					  { "Leading Shower Razzle Electron Score", kLeadingShwRazzleElectronScore, Binning::Simple(40,0,1), ";Leading Shower Razzle Electron Score;Events", "leading_shw_razzle_electron_score",  {.42,.57,.72,.85} },
					  { "Leading Shower Razzle Photon Score", kLeadingShwRazzlePhotonScore, Binning::Simple(40,0,1), ";Leading Shower Razzle Photon Score;Events", "leading_shw_razzle_photon_score",  {.42,.57,.72,.85} },
					  { "Leading Shower Razzle Other Score", kLeadingShwRazzleOtherScore, Binning::Simple(40,0,1), ";Leading Shower Razzle Other Score;Events", "leading_shw_razzle_other_score",  {.42,.57,.72,.85} },
					  { "SubLeading Shower Energy", kSubLeadingShwEnergy, Binning::Simple(40,0,1), ";SubLeading Shower Energy (GeV);Events", "subleading_shw_energy", {.59,.57,.89,.85} },
					  { "SubLeading Shower dEdx", kSubLeadingShwdEdx, Binning::Simple(40,0,10), ";SubLeading Shower dE/dx (MeV/cm);Events", "subleading_shw_dedx", {.59,.57,.89,.85} },
					  { "SubLeading Shower Conversion Gap", kSubLeadingShwCnvGap, Binning::Simple(40,0,20), ";SubLeading Shower Conversion Gap (cm);Events", "subleading_shw_cnv_gap", {.59,.57,.89,.85} },
					  { "SubLeading Shower Density", kSubLeadingShwDensity, Binning::Simple(40,0,20), ";SubLeading Shower Density (MeV/cm);Events", "subleading_shw_density", {.59,.57,.89,.85} },
					  { "SubLeading Shower Length", kSubLeadingShwLen, Binning::Simple(40,0,200), ";SubLeading Shower Length (cm);Events", "subleading_shw_len", {.59,.57,.89,.85} },
					  { "SubLeading Shower Opening Angle", kSubLeadingShwOpenAngle, Binning::Simple(36,0,90), ";SubLeading Shower Opening Angle (#circ);Events", "subleading_shw_open_angle", {.59,.57,.89,.85} },
					  { "SubLeading Shower Cosmic Dist", kSubLeadingShwCosmicDist, Binning::Simple(40,0,400), ";SubLeading Shower Cosmic Dist (cm);Events", "subleading_shw_cosmic_dist", {.59,.57,.89,.85} },
					  { "SubLeading Shower Density Gradient", kSubLeadingShwDensityGrad, Binning::Simple(40,0,4), ";SubLeading Shower Density Gradient (MeV/cm);Events", "subleading_shw_density_grad", {.59,.57,.89,.85} },
					  { "SubLeading Shower Density Gradient Power", kSubLeadingShwDensityGradPower, Binning::Simple(4,0,4), ";SubLeading Shower Density Gradient Power;Events", "subleading_shw_density_grad_power", {.59,.57,.89,.85} },
					  { "SubLeading Shower Track Stub Length", kSubLeadingShwTrkLength, Binning::Simple(24,0,24), ";SubLeading Shower Track Stub Length (cm);Events", "subleading_shw_trk_stub_length", {.59,.57,.89,.85} },
					  { "SubLeading Shower Track Stub Width", kSubLeadingShwTrkWidth, Binning::Simple(24,0,24), ";SubLeading Shower Track Stub Width (cm);Events", "subleading_shw_trk_stub_width", {.59,.57,.89,.85} },
					  { "SubLeading Shower PDG", kSubLeadingShwPDGPlot, binsPDGShw, ";SubLeading Shower PDG;Events", "subleading_shw_pdg", {.22,.57,.52,.85} },
					  { "SubLeading Shower Razzle Electron Score", kSubLeadingShwRazzleElectronScore, Binning::Simple(40,0,1), ";SubLeading Shower Razzle Electron Score;Events", "subleading_shw_razzle_electron_score",  {.42,.57,.72,.85} },
					  { "SubLeading Shower Razzle Photon Score", kSubLeadingShwRazzlePhotonScore, Binning::Simple(40,0,1), ";SubLeading Shower Razzle Photon Score;Events", "subleading_shw_razzle_photon_score",  {.42,.57,.72,.85} },
					  { "SubLeading Shower Razzle Other Score", kSubLeadingShwRazzleOtherScore, Binning::Simple(40,0,1), ";SubLeading Shower Razzle Other Score;Events", "subleading_shw_razzle_other_score",  {.42,.57,.72,.85} },
};

// std::vector<Plot2D<SpillVar>> recoPlots2d = {
//   {"Nu Prism",kNuX,kNuY,Binning::Simple(-300,300,60),Binning::Simple(-300,300,60),";x_{#nu};y_{#nu}","nu_pos",{.59,.57,.89,.85}},
// };