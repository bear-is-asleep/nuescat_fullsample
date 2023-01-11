#include "sbnana/CAFAna/Core/FileReducer.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

using namespace ana;

//Save reduced file with cuts applied for further investigation
void reducer(SpillCut cut, 
  std::string fin, 
  std::string fout = "reduced.root"){
    FileReducer reducer(fin,fout);

    reducer.AddSpillCut(cut); //Should be combination of cuts

    // And when we do keep them, remove their true particle list
    //reducer.AddReductionStep(ClearTrueParticles);

    reducer.Go();
  }

