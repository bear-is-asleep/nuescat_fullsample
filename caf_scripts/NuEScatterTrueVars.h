const Var kTrueNuESlice([](const caf::SRSliceProxy* sr) -> double {
  if (isnan(sr->truth.E)) return -9999;
  return sr->truth.E;
});

const Var kTrueVisEVar([](const caf::SRSliceProxy* sr) -> double {
  double visE = 0;
  for (auto const& prim : sr->truth.prim){
    if (abs(prim.pdg) == 12 || abs(prim.pdg) == 14){continue;} //skip neutrinos
    visE+=(prim.plane[0][0].visE + prim.plane[0][1].visE + prim.plane[0][2].visE)/3;
  }
  return visE; //average of 3 planes
});

const Var kTrueElectronEVar([](const caf::SRSliceProxy* sr) -> double {
  double E = 0;
  for (auto const& prim : sr->truth.prim){
    if (prim.pdg == 11){return prim.genE;} //skip neutrinos
  }
  return -1; //average of 3 planes
});

const Var kTrueEVar([](const caf::SRSliceProxy* sr) -> double {
  double E = 0;
  for (auto const& prim : sr->truth.prim){
    if (abs(prim.pdg) == 12 || abs(prim.pdg) == 14){continue;} //skip neutrinos
    E+=prim.genE;
  }
  return E; //average of 3 planes
});

const SpillVar kTrueNuE([](const caf::SRSpillProxy* sr) -> double {
  auto const& nu = sr->mc.nu[kBestNuID(sr)];
  return nu.E;
});

const SpillVar kTrueNuQ2([](const caf::SRSpillProxy* sr) -> double {
  auto const& nu = sr->mc.nu[kBestNuID(sr)];
  return nu.Q2;
});

const SpillVar kTrueNuY([](const caf::SRSpillProxy* sr) -> double {
  auto const& nu = sr->mc.nu[kBestNuID(sr)];
  return nu.inelasticityY;
});

const SpillVar kTrueVisE([](const caf::SRSpillProxy* sr) -> double {
  auto const& nu = sr->mc.nu[kBestNuID(sr)];
  return (nu.plane[0][0].visE + nu.plane[0][1].visE + nu.plane[0][2].visE)/3; //average of 3 planes
});

const SpillVar kTrueNuXSec([](const caf::SRSpillProxy* sr) -> double {
  auto const& nu = sr->mc.nu[kBestNuID(sr)];
  return nu.xsec;
});

const SpillVar kTrueNuE_0([](const caf::SRSpillProxy* sr) -> double {
  return sr->mc.nu[0].E;
});

const SpillVar kTrueNuE_1([](const caf::SRSpillProxy* sr) -> double {
  if (sr->mc.nnu<=1){return -9999;} //Return dummy value
  return sr->mc.nu[1].E;
});

const SpillVar kTrueNuXSec_0([](const caf::SRSpillProxy* sr) -> double {
  return sr->mc.nu[0].xsec;
});

const SpillVar kTrueNuXSec_1([](const caf::SRSpillProxy* sr) -> double {
  if (sr->mc.nnu<=1){return -9999;} //Return dummy value
  return sr->mc.nu[1].xsec;
});

const Var kTrueNuXSecVar([](const caf::SRSliceProxy* sr) -> double {
  return sr->truth.xsec;
});

const Var kFluxWeightVar([](const caf::SRSliceProxy* sr) -> double {
  return (GeV2perm2/(kTrueNuXSecVar(sr)*nt)*1e-6);
});

const SpillVar kFluxWeight_0([](const caf::SRSpillProxy *sp) -> double{
  return (GeV2perm2/(kTrueNuXSec_0(sp)*nt)*1e-6); //per cm^2
});

const SpillVar kFluxWeight_1([](const caf::SRSpillProxy *sp) -> double{
  return GeV2perm2/(kTrueNuXSec_1(sp)*nt)*1e-6; //per cm^2
});

const SpillCut kIsNCQEOnArgon_0([](const caf::SRSpillProxy* sr){
  auto const& nu = sr->mc.nu[0];
  return !nu.iscc && nu.pdg == 14 &&
      nu.genie_mode == caf::kQE && nu.genie_inttype == caf::kNCQE &&
      nu.targetPDG == 1000180400 /* Argon 40 */ && 
      !nu.ischarm &&
      nu.hitnuc == 2112;
});

const Cut kIsNuMuNCQEOnArgonVar([](const caf::SRSliceProxy* sr,const int pdg = 14){
  auto const& truth = sr->truth;
  return !truth.iscc && truth.pdg == pdg &&
      truth.genie_mode == caf::kQE && truth.genie_inttype == caf::kNCQE &&
      truth.targetPDG == 1000180400 /* Argon 40 */ && 
      !truth.ischarm &&
      truth.hitnuc == 2112;
});

const Cut kIsTrueAVVar([](const caf::SRSliceProxy* sp) {
  return PtInVolAbsX(sp->truth.position, avnd);
});

const SpillCut kIsTrueAV_0([](const caf::SRSpillProxy* sp) {
  if (sp->mc.nnu==0){return false;} //Return dummy value
  auto const& nu = sp->mc.nu[0];
  return PtInVolAbsX(nu.position, avnd);
});

const SpillCut kIsTrueAV_1([](const caf::SRSpillProxy* sp) {
  if (sp->mc.nnu<=1){return false;} //Return dummy value
  auto const& nu = sp->mc.nu[1];
  return PtInVolAbsX(nu.position, avnd);
});

