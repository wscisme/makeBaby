#ifndef ScanChain_h
#define ScanChain_h

// C++ includes
#include <string>
#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"

#include "/home/users/aaivazis/CORE/ssSelections.h"
#include "/home/users/aaivazis/CORE/muonSelections.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class babyMaker {

 public:

  babyMaker() {};
  ~babyMaker() {
    delete BabyFile_;
    delete BabyTree_;
  };

  void ScanChain(TChain* chain, std::string baby_name = "testSample", unsigned int numEvents = 0, float customScale = -1);

  void MakeBabyNtuple(const char *);
  void InitBabyNtuple();
  void FillBabyNtuple();
  void CloseBabyNtuple();

 private:

  TFile *BabyFile_;
  TTree *BabyTree_;

  //baby ntuple variables
  float _pf_sumet;
  float _pf_tcmet;
  float _pf_tcmetPhi;


  float _met;
  float _metPhi;

  //std::vector<LorentzVector> jets_p4;
  std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > _jets_p4;
  std::vector<float> _jets_p4Correction;

  int type;

  LorentzVector _ll_p4;
  LorentzVector _lt_p4;
  LorentzVector _total_p4;
  LorentzVector _lr_p4;
  LorentzVector _lr_vtx_p4;

  int _ll_id;
  int _lt_id;
  int _ll_charge;
  int _lt_charge;
  int _ll_index;
  int _lt_index;

  int _lr_id;
  int _lr_charge;
  int _lr_index;

  int _isRealData;
  int _nvtxs;
  int _njets;
  int _nbTag;

  bool _isZmetLep;
  bool _isHLT1Lep;
  bool _isHLT2Lep;
  bool _isStopLep;
  bool _isSamesignLep;

  bool _trackingProblemVeto;
  bool _tauVeto;
  bool _isoTrackVeto;

  float _scale_1fb;
  std::vector<float> _btagDiscriminant;
  
  int _numEvents;

  int _eventNumber;
  int _runNumber;
  int _lumiBlock;

  string file;

};

#endif

void babyMaker::MakeBabyNtuple(const char *BabyFilename){

  //
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  BabyFile_ = new TFile(Form("%s", BabyFilename), "RECREATE");
  BabyFile_->cd();
  BabyTree_ = new TTree("tree", "A Baby Ntuple");

  BabyTree_->Branch("pf_sumet", &_pf_sumet );
  BabyTree_->Branch("pf_tcmet", &_pf_tcmet );
  BabyTree_->Branch("pf_tcmetPhi", &_pf_tcmetPhi );

  BabyTree_->Branch("met", &_met );
  BabyTree_->Branch("metPhi", &_metPhi );
  BabyTree_->Branch("jets_p4", &_jets_p4 );
  BabyTree_->Branch("jets_p4Correction", &_jets_p4Correction );
  BabyTree_->Branch("type", &type);

  BabyTree_->Branch("ll_p4", &_ll_p4);
  BabyTree_->Branch("lt_p4", &_lt_p4);
  BabyTree_->Branch("lr_p4", &_lr_p4);
  BabyTree_->Branch("total_p4", &_total_p4);
  BabyTree_->Branch("lr_vtx_p4", &_lr_vtx_p4);

  BabyTree_->Branch("ll_id", &_ll_id);
  BabyTree_->Branch("lt_id", &_lt_id);
  BabyTree_->Branch("lr_id", &_lr_id);

  BabyTree_->Branch("isSamesignLep", &_isSamesignLep);
  BabyTree_->Branch("isStopLep", &_isStopLep);
  BabyTree_->Branch("isZmetLep", &_isZmetLep);
  BabyTree_->Branch("isHLT1Lep", &_isHLT1Lep);
  BabyTree_->Branch("isHLT2Lep", &_isHLT2Lep);

  BabyTree_->Branch("trackingProblemVeto", &_trackingProblemVeto);
  BabyTree_->Branch("tauVeto", &_tauVeto);
  BabyTree_->Branch("isoTrackVeto", &_isoTrackVeto);

  BabyTree_->Branch("ll_charge" , &_ll_charge);
  BabyTree_->Branch("lt_charge",  &_lt_charge);
  BabyTree_->Branch("lr_charge",  &_lr_charge);
  BabyTree_->Branch("ll_index",   &_ll_index);
  BabyTree_->Branch("lt_index",   &_lt_index);
  BabyTree_->Branch("lr_index",   &_lr_index);
  
  BabyTree_->Branch("eventNumber", &_eventNumber);
  BabyTree_->Branch("runNumber",   &_runNumber);
  BabyTree_->Branch("lumiBlock",   &_lumiBlock);
  
  BabyTree_->Branch("scale_1fb", &_scale_1fb);
  BabyTree_->Branch("btagDiscriminant", &_btagDiscriminant);

  BabyTree_->Branch("isRealData", &_isRealData );
  BabyTree_->Branch("nvtxs", &_nvtxs );
  BabyTree_->Branch("njets", &_njets);
  BabyTree_->Branch("nbTag", &_nbTag);

  BabyTree_->Branch("numEvents", &_numEvents);

  BabyTree_->Branch("file", &file);

  return;
}

void babyMaker::InitBabyNtuple () {

  type = -1;
  _pf_sumet = -999.0;
  _pf_tcmet = -999.0;
  _pf_tcmetPhi = -999.0;
  _met = -999.0;
  _metPhi = -999.0;
  _ll_id = -1;
  _lt_id = -1;
  _lr_id = -1;
  _ll_charge = -999;
  _lt_charge = -999;
  _lr_charge = -999;
  _ll_index = -1;
  _lt_index = -1;
  _lr_index = -1;
  _scale_1fb = 1;
  _isRealData = -999;
  _nvtxs = -1;
  _njets = -1;
  _nbTag = -1;
  _eventNumber = -1;
  _runNumber = -1;
  _lumiBlock = 1;
  _numEvents = 0;

  return;
}

void babyMaker::FillBabyNtuple(){
  BabyTree_->Fill();
  return;
}

void babyMaker::CloseBabyNtuple(){
  BabyFile_->cd();
  BabyTree_->Write();
  BabyFile_->Close();
  return;
}

