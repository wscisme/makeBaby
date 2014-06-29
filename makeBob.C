// C++
#include <iostream>
#include <vector>
#include <set>
#include <fstream>

// ROOT
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

// CMS2
#include "CMS2.h"
#include "electronSelections.h"
#include "muonSelections.h"
#include "ssSelections.h"
#include "trackSelections.h"
#include "eventSelections.h"
#include "susySelections.h"

#include "/home/users/jgran/CMSSW_5_3_2_patch4_V05-03-23/src/CMS2/NtupleMacros/Tools/goodrun.cc"
#include "/home/users/aaivazis/susy/babymaker/Include.C"

// header
#include "makeBaby.h"
#include "/home/users/sicheng/build/NSBABY.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;
//using namespace tas;
using namespace ROOT::Math::VectorUtil;
using namespace ROOT::Math;
using namespace samesign;
using namespace nsbaby;


void babyMaker::ScanChain(TChain* chain, std::string baby_name, unsigned int numEvent, float customScale){

  MakeBabyNtuple( Form("babies/%s.root", baby_name.c_str()) );

  // File Loop
  int nDuplicates = 0;
  int nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  int fileCounter = 0;
  int goodCounter = 0;
  int pt20Counter = 0; 
  int osCounter = 0; 
  int etaCounter = 0;
  int muIdCounter = 0;
  int elIdCounter = 0; 
  int llCounter = 0;
  int ltCounter = 0;
  int muCounter = 0;
  int elCounter = 0;
  int noGoodVtx = 0;
  int less_jets = 0;

  const float bTagDiscriminant = 0.244;

  // set good run list
  set_goodrun_file("/home/users/jgran/analysis/sswh/fakes/json/final_19p49fb_cms2.txt");

  //File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    fileCounter++;
    
    if (fileCounter > numEvent && numEvent !=0) {
        break;
    }
    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("tree");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    baby.Init(tree);

    // event Loop
    unsigned int nEventsTree = tree->GetEntriesFast();

    for(unsigned int event = 0; event < nEventsTree; ++event) {

      // get event content
      tree->LoadTree(event);
      baby.GetEntry(event);
      ++nEventsTotal;

       // select good 
/*
      if(evt_isRealData() && !goodrun(evt_run(), evt_lumiBlock())) continue;

      if(evt_isRealData()){
        DorkyEventIdentifier id = {evt_run(), evt_event(), evt_lumiBlock()};
        if (is_duplicate(id)){
	nDuplicates++;
          continue;
        }
      }
*/  
    // Progress
      CMS2::progress( nEventsTotal, nEventsChain );

      // count number of hypotheses
      int index = 0;
      //float _maxPt = 0.0;

      // count number of hypotheses

      // To avoid no good vertex warning
      // if ( firstGoodVertex() < 0) { noGoodVtx++; }//continue;}

      int n_mus = 0;
      int n_els = 0;      
      // apply cuts to hypotheses
/*
      for (unsigned int i = 0; i< els_p4().size(); i++){
      	//SELECTION
      	if (els_p4().at(i).pt() < 20)	        {	pt20Counter++;  continue; }
	if (fabs(els_p4().at(i).eta()) > 2.4)   {	etaCounter++;   continue; }
	if ( ! isNumeratorLepton(11, i) )       {	elIdCounter++;  continue; }
      	// electron id			        
        // if ( !passElectronSelection_ZMet2012_v3_Iso(i)) continue;
      	// if ( electronIsoValuePF2012_FastJetEffArea_v3(i) > 0.15)   continue;
      	// if ( ! ETHSelection::goodElectron(i) ) continue;

      	n_els++;
      	index = i;
      }    
      for (unsigned int i = 0; i< mus_p4().size(); i++){
	//SELECTION
	if (mus_p4().at(i).pt() < 20)		{	pt20Counter++;  continue; } 
	if (fabs(mus_p4().at(i).eta()) > 2.4)	{	etaCounter++;   continue; }
	if ( ! isNumeratorLepton(13, i) )	{	muIdCounter++;  continue; }

	// if ( !muonId(i, ZMet2012_v1)) continue;
     	// if ( muonIsoValuePF2012_deltaBeta(i) > 0.2 ) continue;
     	// if ( ! ETHSelection::goodMuon(i) ) continue;

	n_mus++;
	index = i;
      }
      
      if ((n_els + n_mus) != 1 ) continue;
*/
      if( jets_p4().size() < 2) {less_jets++; continue;} // pre-selection check
      

      int n_bTag =0;
      int n_jets = 0;

      for(unsigned int c = 0; c < jets_p4().size(); c++) {
      	float _jetPt = jets_p4().at(c).pt() * jets_p4Correction().at(c);
      	// jet pt times correction to jets
      	if(_jetPt < 30) continue;       // disgard those with small pt

      	if(fabs(jets_p4().at(c).eta()) > 2.5) continue;   // Jet eta check

      	float dr_lr = DeltaR(jets_p4().at(c), lr_p4() );
      	// delta R is the distance in the eta-phi plane
      	if(dr_lr < 0.4) continue;      // disgard small delta R 

      	float _bTag = btagDiscriminant().at(c);
      	if ( _bTag > bTagDiscriminant)  n_bTag++;	 // mark as the has qualified b-jet, _bTag should be 0~1

      	n_jets++;					 // as qualified jet multiplicity
      }//pfjets_p4().size()
      
      if(n_jets < 2 || n_bTag != 0) {less_jets++; continue;}


      // create the ntuple
      InitBabyNtuple();
      // set the branch values

      _pf_sumet		 = pf_sumet(); 
      _pf_tcmet		 = pf_tcmet(); 
      _pf_tcmetPhi	 = pf_tcmetPhi();
      _lr_id	         = lr_id();
      _njets		 = n_jets;
      _nbTag		 = n_bTag;
      _met		 = met();
      _metPhi		 = metPhi();
      _lr_index		 = lr_index();
      _jets_p4		 = jets_p4();
      _jets_p4Correction = jets_p4Correction();
      _btagDiscriminant  = btagDiscriminant();
      _lr_p4		 = lr_p4();
      _lr_charge	 = lr_charge();
      _lr_vtx_p4	 = lr_vtx_p4();
      _eventNumber 	 = eventNumber();
      _runNumber   	 = runNumber();
      _lumiBlock   	 = lumiBlock();
      _scale_1fb   	 = scale_1fb();
      _isRealData  	 = isRealData();
      _nvtxs       	 = nvtxs();
      _isZmetLep 	 = isZmetLep();
      _isHLT1Lep 	 = isHLT1Lep();
      _isHLT2Lep 	 = isHLT2Lep();
      _isStopLep 	 = isStopLep();
      _isSamesignLep 	 = isSamesignLep();
      _tauVeto		 = tauVeto();
      _isoTrackVeto	 = isoTrackVeto();
      _trackingProblemVeto = trackingProblemVeto();

      file = Form("%s", currentFile->GetTitle());

      FillBabyNtuple();


    } // end of loop over events in file

    delete tree;
    f.Close();

  } // end of loop over files

  if (nEventsChain != nEventsTotal){
      cout << "WARNING: The number of events added is not equal to the total number of events!" << endl;
  }

  cout << nDuplicates << " duplicate events were skipped." << endl;

  CloseBabyNtuple();

  ofstream myfile;
  myfile.open(Form("%s.txt", baby_name.c_str() ));
  if(myfile.is_open())
  {

    myfile << "total events :" << nEventsTotal << endl;
    myfile << "good  events :" << goodCounter << endl;
    myfile << "Duplicates :" << nDuplicates << endl;
    myfile << "no Good Vertex :" << noGoodVtx << endl;

    myfile << "# of ll leptons :" << llCounter << endl;
    myfile << "# of lt leptons :" << ltCounter << endl;

    myfile << "# of electrons :" << elCounter << endl;
    myfile << "# of muons:" << muCounter << endl;

    myfile << "# of leptons cut at 20 GeV :" << pt20Counter << endl;


    myfile << "# of events cut at OS :" << osCounter << endl;
    myfile << "# of events cut at ETA :" << etaCounter << endl;

    myfile << "# of events cut at muon ID/ISO :" << muIdCounter << endl;
    myfile << "# of events cut at electron ID/ISO :" << elIdCounter << endl;
  }

  return;
}

