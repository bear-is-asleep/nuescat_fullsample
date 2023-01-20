const SpillCut kCosmicSpill([](const caf::SRSpillProxy* sp) {
    return sp->mc.nnu == 0;
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
    return kCC(sp) && std::abs(sp->mc.nu[kBestNuID(sp)].pdg) == 14 &&
      sp->mc.nu[kBestNuID(sp)].hitnuc != 0;
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

const SpillVar kNChargedPi = kNPiPlus + kNPiMinus;

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

const SpillCut kNCPiZero([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) == 0 && kNPiZero(sp) >= 1;
  });
const SpillVar kNPrim([](const caf::SRSpillProxy* sp) {
  return sp->mc.nu[kBestNuID(sp)].nprim;
});

const SpillCut kNuEScat([](const caf::SRSpillProxy* sp) {
    if(kDirt(sp) || kCosmicSpill(sp) || kCCNuE(sp) || kCCNuMu(sp)){return false;}
    return kNElectron(sp) == 1 && kNPrim(sp) == 2 && sp->mc.nu[kBestNuID(sp)].hitnuc == 0;
  });

const SpillCut kNCNPiZero([](const caf::SRSpillProxy* sp) {
    return kNC(sp) && kNChargedPi(sp) == 0 && kNPiZero(sp) > 1;
  });

const SpillCut kOther([](const caf::SRSpillProxy* sp) {
  if(kNC(sp) || kCCNuMu(sp) || kCCNuE(sp) || kDirt(sp) || kCosmicSpill(sp) || kNuEScat(sp)){return false;}
  else{
    return true;
  }
});

const SpillCut kOtherFV([](const caf::SRSpillProxy* sp) {
  return kOther(sp) && kTrueFV(sp);
});

std::vector<TrueCategory> nuescat_sel_categories = {
  //{"#nu + e",kSignal,kOrange+2,"Signal"},
  {"#nu + e FV",kNuEScat && kTrueFV,kBlue,"NuEScatFV"},
  {"#nu + e",kNuEScat,kBlue+3,"NuEScat"},
  //{"NC N#pi^{0}", kNCPiZero && kTrueFV, kMagenta+2, "NCpi0"},
  //{"Other NC", kNC && !kNCPiZero && kTrueFV, kYellow+2, "NC"},
  //{"CC #nu_{#mu}", kCCNuMu && kTrueFV, kRed+2, "CCNuMu"},
  //{"CC #nu_{e}", kCCNuE && kTrueFV, kTeal+2, "CCNuE"},
  //{"Dirt", kDirt && kTrueFV, kOrange+3, "Dirt"},
  //{"Cosmic", kCosmicSpill, kRed+1, "Cosmic"},
  {"Other", kOther, kBlack, "Other"},
  {"Other FV", kOtherFV, kOrange+2, "OtherFV"},
};