int print_all_prim_info_slc(const caf::SRSpillProxy* sp,unsigned slcID){
  if(sp->nslc==0) return 0;
  auto const& slc = sp->slc[slcID];
  std::cout<<"----- primary slice info (slc) -----"<<std::endl;
  for (auto const& prim : slc.truth.prim){
    std::cout<<"--pdg : "<<prim.pdg<<", "<<std::endl;
  }
  return 1;

}

int print_all_prim_info_nu(const caf::SRSpillProxy* sp,unsigned nuID){
  auto const& nu = sp->mc.nu[nuID];
  std::cout<<"----- primary slice info (nu) -----"<<std::endl;
  for (auto const& prim :nu.prim){
    std::cout<<"--pdg : "<<prim.pdg<<", "<<std::endl;
  }
  std::cout<<"-genie_inttype : "<<nu.genie_inttype<<std::endl;
  return 1;
}

int print_all_prim_matched(const caf::SRSpillProxy* sp,unsigned slcID){
  auto const& slc = sp->slc[slcID];
  std::cout<<"----- matched pdg info -----"<<std::endl;
  for (auto const& shw : slc.reco.shw){
    std::cout<<"Match shw pdg : "<<shw.truth.p.pdg<<std::endl;
  }
  for (auto const& trk : slc.reco.trk){
    std::cout<<"Match trk pdg : "<<trk.truth.p.pdg<<std::endl;
  }
  return 1;
}