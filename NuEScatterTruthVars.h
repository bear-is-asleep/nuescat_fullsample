const SpillVar kPiZeroID([](const caf::SRSpillProxy* sp) -> unsigned {
    //if(!kSignal(sp)) return 999999;

    for(int i = 0; i < sp->mc.nu[kBestNuID(sp)].nprim; ++i) {
      auto const& prim = sp->mc.nu[kBestNuID(sp)].prim[i];
      if(prim.pdg == 111) return i;
    }

    return 999999;
  });

// const SpillVar kNuEScatID([](const caf::SRSpillProxy* sp) -> unsigned {
//     if(!kSignal(sp)) return 999999;
//     if(sp->mc.nu[kBestNuID(sp)].nprim > 1) return 999999;
//     auto const& prim = sp->mc.nu[kBestNuID(sp)].prim[0];
//     if(abs(prim.pdg) == 11) return i; //Is an electron

//     return 999999;
//   });

const SpillVar kNuECCID([](const caf::SRSpillProxy* sp) -> unsigned {
    //if(!kSignal(sp)) return 999999;
    if(sp->mc.nu[kBestNuID(sp)].isnc) return 999999;
    for(int i = 0; i < sp->mc.nu[kBestNuID(sp)].nprim; ++i) {
      auto const& prim = sp->mc.nu[kBestNuID(sp)].prim[i];
      if(prim.pdg == 11) return i;
    }

    return 999999;
  });

const SpillCut kTwoGamma([](const caf::SRSpillProxy* sp) {
    //if(!kSignal(sp)) return false;

    for(int i = 0; i < sp->mc.nu[kBestNuID(sp)].nprim; ++i) {
      auto const& prim = sp->mc.nu[kBestNuID(sp)].prim[i];
      if(prim.pdg == 111 && prim.daughters.size() == 2) return true;
    }
    
    return false;
  });

const SpillVar kPiZeroLeadingPhotonID([](const caf::SRSpillProxy* sp) -> unsigned {
    if(!kTwoGamma(sp)) return false;

    const unsigned pion_id = kPiZeroID(sp);

    auto const& pion = sp->mc.nu[kBestNuID(sp)].prim[pion_id];
    auto const& daughters = pion.daughters;

    if(daughters.size() != 2)
      return 999999;

    int photon_one_id = daughters[0];
    int photon_two_id = daughters[1];

    double en_one(-999999.), en_two(-999999.);

    unsigned one_id(999999), two_id(999999);

    for(int i = 0; i < sp->ntrue_particles; ++i)
      {
        auto const& part = sp->true_particles[i];
        if(part.G4ID == photon_one_id)
          {
            en_one = part.genE;
            one_id = i;
          }
        if(part.G4ID == photon_two_id)
          {
            en_two = part.genE;
            two_id = i;
          }
      }

    return en_one > en_two ? one_id : two_id;
  });

const SpillVar kPiZeroSubLeadingPhotonID([](const caf::SRSpillProxy* sp) -> unsigned {
    if(!kTwoGamma(sp)) return false;

    const unsigned pion_id = kPiZeroID(sp);

    auto const& pion = sp->mc.nu[kBestNuID(sp)].prim[pion_id];
    auto const& daughters = pion.daughters;

    if(daughters.size() != 2)
      return 999999;

    int photon_one_id = daughters[0];
    int photon_two_id = daughters[1];

    double en_one(-999999.), en_two(-999999.);

    unsigned one_id(999999), two_id(999999);

    for(int i = 0; i < sp->ntrue_particles; ++i)
      {
        auto const& part = sp->true_particles[i];
        if(part.G4ID == photon_one_id)
          {
            en_one = part.genE;
            one_id = i;
          }
        if(part.G4ID == photon_two_id)
          {
            en_two = part.genE;
            two_id = i;
          }
      }

    return en_one > en_two ? two_id : one_id;
  });

const SpillVar kPiZeroThresholdRecoStatus([](const caf::SRSpillProxy *sp) {
    if(!kTwoGamma(sp)) return -1.;
    
    int photon_one_id = sp->true_particles[kPiZeroLeadingPhotonID(sp)].G4ID;
    int photon_two_id = sp->true_particles[kPiZeroSubLeadingPhotonID(sp)].G4ID;

    bool reco_one = false, reco_two = false;

    for(auto const& shw : sp->reco.shw)
      {
	if(shw.truth.bestmatch.G4ID == photon_one_id 
	   && shw.truth.bestmatch.hit_completeness > 0.5
	   && shw.truth.bestmatch.hit_purity > 0.5)
	    reco_one = true;
	else if(shw.truth.bestmatch.G4ID == photon_two_id 
		&& shw.truth.bestmatch.hit_completeness > 0.5
		&& shw.truth.bestmatch.hit_purity > 0.5)
	    reco_two = true;
      }

    if(reco_one && reco_two)
      return 1.;
    else if(reco_one)
      return 2.;
    else if(reco_two)
      return 3.;
    else 
      return 4.;
  });

// const SpillVar kTruthShw([](const caf::SRSpillProxy *sp) {
//   //if(!kSignal(sp)) return 999999;
//   int nshw = 0;
//   auto const& slc = sp->slc[kBestSlcID(sp)];
//   for (auto const& pdg : slc.truth.pdg){
//     if (abs(pdg) == 11):{++nshw;}
//     else (if pdg == 22):{++nshw;}
//     else (if pdg == 111):{nshw+=2;}   
//   }
//   return nshw;
// });

// const SpillVar kTruthTrk([](const caf::SRSpillProxy *sp) {
//   //if(!kSignal(sp)) return 999999;
//   int ntrk = 0;
//   auto const& slc = sp->slc[kBestSlcID(sp)];
//   for (auto const& pdg : slc.truth.pdg){
//     if (abs(pdg) == 2212):{++ntrk;}
//     else if (abs(pdg) == 13):{++ntrk;}
//     else if (abs(pdg) == 211):{++ntrk;}   
//   }
//   return ntrk;
// });

//Not yet implemented
// const SpillVar kTruthStub([](const caf::SRSpillProxy *sp) {
//   if(!kSignal(sp)) return 999999;
//   int ntrk = 0;
//   auto const& slc = sp->slc[kBestSlcID(sp)];
//   for (auto const& pdg : sp->truth.pdg){
//     if abs(pdg) == 2212:{++ntrk;}
//     else if abs(pdg) == 13:{++ntrk;}
//     else if abs(pdg) == 211:{++ntrk;}   
//   }
//   return ntrk;
// });


