const SpillCut kCosmicSpill([](const caf::SRSpillProxy* sp) {
    return sp->mc.nnu == 0;
  });

const SpillCut kMultiNu([](const caf::SRSpillProxy* sp) {
    return sp->mc.nnu > 1;
  });

// Choose the neutrino that deposited the most energy in the TPC, this neutrino will
// be used to determine the event type (e.g. CCNuMu) and be the target of selection.
const SpillVar kBestNuID([](const caf::SRSpillProxy* sp) {
    int id = std::numeric_limits<int>::max(), counter = 0;
    double most_en = -std::numeric_limits<double>::max();

    for(auto const& nu : sp->mc.nu)
      {
        double visE = nu.plane[0][0].visE + nu.plane[0][1].visE + nu.plane[0][2].visE;

        if(visE > most_en)
          {
            most_en = visE;
            id      = counter;
          }
        ++counter;
      }
    return id;
  });

//Chose neutrino event with least primaries - optimal for selecting nu e scat?
const SpillVar kBestNuIDLeast([](const caf::SRSpillProxy* sp) {
    int id = std::numeric_limits<int>::max(), counter = 0;
    int nprim = std::numeric_limits<int>::max();

    for(auto const& nu : sp->mc.nu)
      {

        if(nu.nprim < nprim)
          {
            nprim = nu.nprim;
            id      = counter;
          }
        ++counter;
      }
    return id;
  });

// Chose the slice that has the highest match to our "best neutrino" this is the
// slice we are aiming to select.
const SpillVar kTrueBestSlice([](const caf::SRSpillProxy* sp) -> unsigned {
    if(kCosmicSpill(sp))
      return 999999;

    auto const& nuid = sp->mc.nu[kBestNuID(sp)].index;

    float best_eff(-1.f);
    unsigned id = 0, returnid = 999999;

    for(auto const &slc : sp->slc)
      {
        if(slc.tmatch.index == nuid && slc.tmatch.eff > best_eff)
          {
            best_eff = slc.tmatch.eff;
            returnid = id;
          }
        ++id;
      }

    return returnid;
  });

// Determine if this event contains multiple interactions in the AV. Note this is 
// different to the kMultiNu cut as kMultiNu does not exclude dirt neutrinos.
const SpillCut kPileUp([](const caf::SRSpillProxy* sp) {
    if(sp->mc.nnu < 2) return false;

    int nnuint = 0;

    for(auto const& nu : sp->mc.nu)
      {
        if(PtInVolAbsX(nu.position, avnd))
          ++nnuint;
      }

    return nnuint > 1;
  });

const SpillCut kTrueFV([](const caf::SRSpillProxy *sp) {
    return !kCosmicSpill(sp) && PtInVolAbsX(sp->mc.nu[kBestNuID(sp)].position, fvndNuEScat);
  });

const SpillCut kDirt([](const caf::SRSpillProxy *sp) {
    return !kCosmicSpill(sp) && !PtInVolAbsX(sp->mc.nu[kBestNuID(sp)].position, avnd);
  });

const SpillCut kNC([](const caf::SRSpillProxy* sp) {
    return !kDirt(sp) && !kCosmicSpill(sp) && sp->mc.nu[kBestNuID(sp)].isnc &&
      sp->mc.nu[kBestNuID(sp)].hitnuc != 0;
  });

const SpillCut kCC([](const caf::SRSpillProxy* sp) {
    return !kDirt(sp) && !kCosmicSpill(sp) && sp->mc.nu[kBestNuID(sp)].iscc &&
      sp->mc.nu[kBestNuID(sp)].hitnuc != 0;
  });

const SpillCut kCCNuMu([](const caf::SRSpillProxy* sp) {
    return kCC(sp) && std::abs(sp->mc.nu[kBestNuID(sp)].pdg) == 14;
  });

const SpillCut kCCNuE([](const caf::SRSpillProxy* sp) {
    return kCC(sp) && std::abs(sp->mc.nu[kBestNuID(sp)].pdg) == 12 &&
      sp->mc.nu[kBestNuID(sp)].hitnuc != 0;
  });
    
const SpillVar kNPiPlus([](const caf::SRSpillProxy* sp) -> unsigned {
    if(kCosmicSpill(sp)) return std::numeric_limits<unsigned>::max();

    unsigned int npiplus = 0;

    for(auto const& prim : sp->mc.nu[kBestNuID(sp)].prim)
      if(prim.pdg == 211) npiplus++;

    return npiplus;
  });

const SpillVar kNPiMinus([](const caf::SRSpillProxy* sp) -> unsigned {
    if(kCosmicSpill(sp)) return std::numeric_limits<unsigned>::max();

    unsigned int npiminus = 0;

    for(auto const& prim : sp->mc.nu[kBestNuID(sp)].prim)
      if(prim.pdg == -211) npiminus++;

    return npiminus;
  });

const SpillVar kNPiZero([](const caf::SRSpillProxy* sp) -> unsigned {
    if(kCosmicSpill(sp)) return std::numeric_limits<unsigned>::max();

    unsigned int npizero = 0;

    for(auto const& prim : sp->mc.nu[kBestNuID(sp)].prim)
      if(prim.pdg == 111) npizero++;

    return npizero;
  });

const SpillVar kNElectron([](const caf::SRSpillProxy* sp) -> unsigned {
    if(kCosmicSpill(sp)) return std::numeric_limits<unsigned>::max();

    unsigned int nele = 0;

    for(auto const& prim : sp->mc.nu[kBestNuID(sp)].prim)
      if(prim.pdg == 11) nele++;

    return nele;
  });

const SpillVar kNChargedPi = kNPiPlus + kNPiMinus;
const SpillVar kNPi = kNChargedPi + kNPiZero;

const SpillVar kNProton([](const caf::SRSpillProxy* sp) -> unsigned {
    if(kCosmicSpill(sp)) return std::numeric_limits<unsigned>::max();

    unsigned int nproton = 0;

    for(auto const& prim : sp->mc.nu[kBestNuID(sp)].prim)
      if(prim.pdg == 2212) nproton++;

    return nproton;
  });

const SpillCut kNCPiZero([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) == 0 && kNPiZero(sp) >= 1;
  });
const SpillVar kNPrim([](const caf::SRSpillProxy* sp) {
  return sp->mc.nu[kBestNuID(sp)].nprim;
});

const SpillVar kNuEScat([](const caf::SRSpillProxy* sp) -> unsigned {
    if(kDirt(sp) || kCosmicSpill(sp) || kCCNuE(sp) || kCCNuMu(sp)){return false;}
    if (kNElectron(sp) == 1 && kNPrim(sp) == 2 && sp->mc.nu[kBestNuID(sp)].hitnuc == 0){
      std::cout<<"kNuEScat"<<kNPrim(sp)<<std::endl;
      std::cout<<"Nele"<<kNElectron(sp)<<std::endl;
      std::cout<<"Hit Nuc: "<<sp->mc.nu[kBestNuID(sp)].hitnuc<<std::endl;
      //Loop over primaries. What are these?
      for (int i =0; i<sp->mc.nu[kBestNuID(sp)].nprim; i++){
        std::cout<<"Prim: "<<i<<" pdg: "<<sp->mc.nu[kBestNuID(sp)].prim[i].pdg<<std::endl;
      }
      std::cout<<"run,subrun,event "<<sp->hdr.run<<","<<sp->hdr.subrun<<","<<sp->hdr.evt<<std::endl;
      //return true;
    }
    return kNElectron(sp) == 1 && kNPrim(sp) == 2 && sp->mc.nu[kBestNuID(sp)].hitnuc == 0; //&& sp->mc.nu[kBestNuID(sp)].nprim==2;

  });

const SpillCut kSignal([](const caf::SRSpillProxy* sp) {
    return kTrueFV(sp) && kNuEScat(sp);
  });

const SpillCut kNCNPiZero([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) == 0 && kNPiZero(sp) > 1;
  });

const SpillCut kNCPiZeroProton([](const caf::SRSpillProxy* sp) {
    return kNCPiZero(sp) && kNProton(sp) == 1;
  });

const SpillCut kNCPiZeroNoProton([](const caf::SRSpillProxy* sp) {
    return kNCPiZero(sp) && kNProton(sp) == 0;
  });

const SpillCut kNCPiPlusMinus([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) == 1 && kNPiZero(sp) == 0;
  });

const SpillCut kNCNPiPlusMinus([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) > 1 && kNPiZero(sp) == 0;
  });

const SpillCut kNCMultiPi([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) > 0 && kNPiZero(sp) > 0;
  });

const SpillCut kNCZeroPi([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNPi(sp) == 0 && kNProton(sp) > 0;
  });

const SpillCut kInvNC([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNPi(sp) + kNProton(sp) == 0;
  });

const SpillCut kCCNuMuPiZero([](const caf::SRSpillProxy* sp) {
    return kCCNuMu(sp) && kNChargedPi(sp) == 0 && kNPiZero(sp) == 1;
  });

const SpillCut kCCNuMuNPiZero([](const caf::SRSpillProxy* sp) {
    return kCCNuMu(sp) && kNChargedPi(sp) == 0 && kNPiZero(sp) > 1;
  });

const SpillCut kCCNuMuPiPlusMinus([](const caf::SRSpillProxy* sp) {
    return kCCNuMu(sp) && kNChargedPi(sp) == 1 && kNPiZero(sp) == 0;
  });

const SpillCut kCCNuMuNPiPlusMinus([](const caf::SRSpillProxy* sp) {
    return kCCNuMu(sp) && kNChargedPi(sp) > 1 && kNPiZero(sp) == 0;
  });

const SpillCut kCCNuMuMultiPi([](const caf::SRSpillProxy* sp) {
    return kCCNuMu(sp) && kNChargedPi(sp) > 0 && kNPiZero(sp) > 0;
  });

const SpillCut kCCNuMuZeroPi([](const caf::SRSpillProxy* sp) {
    return kCCNuMu(sp) && kNPi(sp) == 0;
  });

std::vector<TrueCategory> nuescat_sel_categories = {
  {"#nu + e",kSignal,kOrange+2,"Signal"},
  {"NC N#pi^{0}", kNCPiZero, kMagenta+2, "NCpi0"},
  {"Other NC", kNC && !kNCPiZero, kYellow+2, "NC"},
  {"CC #nu_{#mu}", kCCNuMu, kGreen+2, "CCNuMu"},
  {"CC #nu_{e}", kCCNuE, kCyan+2, "CCNuE"},
  {"Dirt", kDirt, kOrange+3, "Dirt"},
  {"Cosmic", kCosmicSpill, kRed+1, "Cosmic"},
  //{"Other", !kNC && !kCCNuMu && !kCCNuE && !kDirt && !kCosmicSpill, kBlack, "Other"}
};