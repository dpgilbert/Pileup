#ifndef ScanChain_h
#define ScanChain_h

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "Math/LorentzVector.h"
#include "../../CORE/IsoTrackVeto.h"
#include "../../CORE/IsolationTools.h"
#include "../../CORE/CMS3.h"
#include "../../CORE/Tools/utils.h"
#include "../../CORE/MuonSelections.h"
#include "../../CORE/MCSelections.h"
#include <vector>

class PileupCorrection {
 public:
  PileupCorrection() {};
  ~PileupCorrection() {};
  
  void ScanChain(TTree *, const char* outname = "TestSample", int MaxEvents = -1);

 private:
  
};

#endif
