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
        if (pdg == 5 || pdg == 0){ //Why is this sometimes -5?
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
        if (Eng*theta*theta < Etheta2 && Etheta2 > kEtheta2Cut){ 
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



std::vector<Plot<SpillVar>> recoPlots = {//{ "E#theta^{2}", kEtheta2Var, Binning::Simple(15,-10,10),";E#theta^{2};Events", "E_theta",{.59,.57,.89,.85}},
            { "E#theta^{2}", kEtheta2Var, Binning::Simple(10,0.,5.),";E#theta^{2};Events", "E_theta",{.59,.57,.89,.85}},
            { "E#theta^{2} (Zoom)", kEtheta2Var, Binning::Simple(15,0.,kEtheta2Cut),";E#theta^{2};Events", "E_theta_zoom",{.59,.57,.89,.85}},
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
            { "ph Razzle", kNRazzlePhotons, Binning::Simple(5,0,5),";Razzle Photons;Events", "ph_razzle",{.59,.57,.89,.85}}
};
