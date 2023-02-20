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

const SpillCut kDirtNuEScat([](const caf::SRSpillProxy *sp) {
    return !kCosmicSpill(sp) && !PtInVolAbsX(sp->mc.nu[kBestNuID(sp)].position, avnd) && 
    sp->mc.nu[kBestNuID(sp)].genie_inttype == 1098;
  });

const SpillCut kNuEScat([](const caf::SRSpillProxy* sp) {
    return !kDirt(sp) && sp->mc.nu[kBestNuID(sp)].genie_inttype == 1098; //Interaction number for nu e scat
  });

const SpillCut kNC([](const caf::SRSpillProxy* sp) {
    return !kDirt(sp) && !kCosmicSpill(sp) && sp->mc.nu[kBestNuID(sp)].isnc;
      //sp->mc.nu[kBestNuID(sp)].hitnuc != 0;
  });

const SpillCut kCC([](const caf::SRSpillProxy* sp) {
    return !kDirt(sp) && !kCosmicSpill(sp) && sp->mc.nu[kBestNuID(sp)].iscc;
  });

const SpillCut kCCNuMu([](const caf::SRSpillProxy* sp) {
    return kCC(sp) && std::abs(sp->mc.nu[kBestNuID(sp)].pdg) == 14 && !kNuEScat(sp);
  });

const SpillCut kCCNuE([](const caf::SRSpillProxy* sp) {
    return kCC(sp) && std::abs(sp->mc.nu[kBestNuID(sp)].pdg) == 12 && !kNuEScat(sp);
  });
const SpillCut kNuEEvent([](const caf::SRSpillProxy* sp) {
    return std::abs(sp->mc.nu[kBestNuID(sp)].pdg) == 12 && !kNuEScat(sp);
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

const SpillVar kNChargedHadron = kNChargedPi + kNProton; //Include kaons later?

const SpillCut kNCPiZero([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) == 0 && kNPiZero(sp) >= 1; 
  });
const SpillVar kNPrim([](const caf::SRSpillProxy* sp) {
  return sp->mc.nu[kBestNuID(sp)].nprim;
});

const SpillVar kCCQENuE([](const caf::SRSpillProxy* sp) {
    return kCCNuE(sp) && sp->mc.nu[kBestNuID(sp)].genie_inttype == 1001;
  });

// const SpillCut kCCNuEZeroHadron([](const caf::SRSpillProxy* sp){
//   return kCCNuE(sp) &&
// });

const SpillCut kSignal([](const caf::SRSpillProxy* sp) {
    return kNuEScat(sp) && kTrueFV(sp);
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

const SpillCut kOther([](const caf::SRSpillProxy* sp) {
  if(kNC(sp) || kCCNuMu(sp) || kCCNuE(sp) || kDirt(sp) || kCosmicSpill(sp) || kNuEScat(sp)){return false;}
  else{
    std::cout<<"kOther----------------"<<std::endl;
    std::cout<<"run info: "<<sp->hdr.run<<", "<<sp->hdr.subrun<<", "<<sp->hdr.evt<<std::endl;
    std::cout<<"nprim: "<< kNPrim(sp) <<std::endl;
    std::cout<<"nele: "<< kNElectron(sp) <<std::endl;
    std::cout<<"hitnuc: "<< sp->mc.nu[kBestNuID(sp)].hitnuc <<std::endl;
    std::cout<<"charged hadr: "<<kNChargedHadron(sp)<<std::endl;
    for(auto const& prim : sp->mc.nu[kBestNuID(sp)].prim)
      std::cout<<"--prim pdg: "<<prim.pdg;
    return true;
  }
});

const SpillCut kOtherFV([](const caf::SRSpillProxy* sp) {
  return kOther(sp) && kTrueFV(sp);
});

std::vector<TrueCategory> full_sel_categories = {
  {"#nu + e",kNuEScat && kTrueFV,kBlue,"NuEScat"},
  {"NC N#pi^{0}", kNCPiZero && kTrueFV, kMagenta+2, "NCpi0"},
  {"Other NC", kNC && !kNCPiZero && !kNuEScat && kTrueFV, kYellow+2, "NC"},
  {"CC #nu_{#mu}", kCCNuMu && kTrueFV, kRed+2, "CCNuMu"},
  {"CC #nu_{e}", kCCNuE && kTrueFV, kTeal+2, "CCNuE"},
  {"Dirt", kDirt, kOrange+3, "Dirt"},
  {"Cosmic", kCosmicSpill, kRed+1, "Cosmic"},
};

std::vector<TrueCategory> nuecc_sel = {
  {"CC #nu_{e}", kCCNuE && kTrueFV, kTeal+2, "CCNuE"},
};

std::vector<TrueCategory> nue_sel = {
  {"#nu + e",kNuEScat && kTrueFV,kBlue,"NuEScat"},
};

std::vector<TrueCategory> nue_nuedirt_sel = {
  {"#nu + e (FV)",kNuEScat && kTrueFV,kBlue,"NuEScat"},
  {"#nu + e",kNuEScat && !kTrueFV,kOrange+3,"NuEScat_dirt"},
};

std::vector<TrueCategory> nuecc_nue_sel = {
  {"#nu + e",kNuEScat && kTrueFV,kBlue,"NuEScat"},
  {"CC #nu_{e}", kCCNuE && kTrueFV, kTeal+2, "CCNuE"},
};

std::vector<TrueCategory> no_cosmic_sel = {
  {"#nu + e",kNuEScat && kTrueFV,kBlue,"NuEScat"},
  {"NC N#pi^{0}", kNCPiZero && kTrueFV, kMagenta+2, "NCpi0"},
  {"Other NC", kNC && !kNCPiZero && !kNuEScat && kTrueFV, kYellow+2, "NC"},
  {"CC #nu_{#mu}", kCCNuMu && kTrueFV, kRed+2, "CCNuMu"},
  {"CC #nu_{e}", kCCNuE && kTrueFV, kTeal+2, "CCNuE"},
  {"Dirt", kDirt, kOrange+3, "Dirt"},
};



//{"#nu + e",kSignal,kOrange+2,"Signal"},
//{"#nu + e",kNuEScat && kTrueFV,kBlue,"NuEScat"},
//{"#nu + e",kNuEScat,kBlue+3,"NuEScat"},
// {"NC N#pi^{0}", kNCPiZero && kTrueFV, kMagenta+2, "NCpi0"},
// {"Other NC", kNC && !kNCPiZero && kTrueFV, kYellow+2, "NC"},
// {"CC #nu_{#mu}", kCCNuMu && kTrueFV, kRed+2, "CCNuMu"},
//{"CC #nu_{e}", kCCNuE && kTrueFV, kTeal+2, "CCNuE"},
//{"#nu_{e}",kNuEEvent, kTeal,"NuE"}
// {"Dirt", kDirt, kOrange+3, "Dirt"},
//{"Cosmic", kCosmicSpill, kRed+1, "Cosmic"},
//{"Other", kOther, kBlack, "Other"},
//{"Other FV", kOtherFV, kOrange+2, "OtherFV"},
