#include <math.h>

const SpillVar kBestSlcID_Etheta([](const caf::SRSpillProxy* sp) -> unsigned {
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

	    ++i;
    }
    return returnID;
  });


const SpillVar kIsFVVar([](const caf::SRSpillProxy* sp) -> double {
    if(sp->nslc==0) return 0;
    auto const& slc = sp->slc[kBestSlcID(sp)];

    return PtInVolAbsX(slc.vertex, fvndNuEScat) ? 1 : 0;
  });

std::vector<Plot<SpillVar>> recoPlots = {
					  { "FV Cut", kIsFVVar, Binning::Simple(2,0,2), ";FV Cut;Events", "fv_cut", {.59,.57,.89,.85} }, 
					  
};

// std::vector<Plot2D<SpillVar>> recoPlots2d = {
//   {"Nu Prism",kNuX,kNuY,Binning::Simple(-300,300,60),Binning::Simple(-300,300,60),";x_{#nu};y_{#nu}","nu_pos",{.59,.57,.89,.85}},
// };