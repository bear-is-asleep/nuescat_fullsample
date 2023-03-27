const Var kTrueNuESlice([](const caf::SRSliceProxy* sr) -> double {
  if (isnan(sr->truth.E)) return -9999;
  return sr->truth.E;
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

const SpillVar kFluxWeight_0([](const caf::SRSpillProxy *sp) -> double{
  return kTrueNuXSec_0(sp)*nt;
});

SpillVar GetUniverseWeight_0(const std::string& psetName, double x)
{
    return SpillVar(UniverseWeight(psetName, x));
}