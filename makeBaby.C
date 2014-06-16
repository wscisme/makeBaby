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

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;
using namespace tas;
using namespace ROOT::Math::VectorUtil;
using namespace samesign;


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
  int samesign = 0;
  int zmet = 0;

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
    TTree *tree = (TTree*)f.Get("Events");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    cms2.Init(tree);

    // event Loop
    unsigned int nEventsTree = tree->GetEntriesFast();

    for(unsigned int event = 0; event < nEventsTree; ++event) {

      // get event content
      tree->LoadTree(event);
      cms2.GetEntry(event);
      ++nEventsTotal;

       // select good 
      if(evt_isRealData() && !goodrun(evt_run(), evt_lumiBlock())) continue;

      if(evt_isRealData()){
        DorkyEventIdentifier id = {evt_run(), evt_event(), evt_lumiBlock()};
        if (is_duplicate(id)){
	nDuplicates++;
          continue;
        }
      }
  
    // Progress
      CMS2::progress( nEventsTotal, nEventsChain );

      // count number of hypotheses
      int index = 0;
      //float _maxPt = 0.0;

      // count number of hypotheses

      // To avoid no good vertex warning
      if ( firstGoodVertex() < 0) { noGoodVtx++; }//continue;}

      /////// Event Cleaning //////////

      /// Rho
      if ( evt_ww_rho_vor() > 40 )		continue;

      /// Filters ///
      if(  evt_cscTightHaloId()    )		continue;
      if( !evt_hbheFilter()        )		continue;
      if( !filt_hcalLaser()        )		continue;
      if( !filt_ecalTP()           )		continue;
      if( !filt_trackingFailure()  )		continue;
      if( !filt_eeBadSc()          )		continue;
      if( !passHBHEFilter()	   )		continue;

      /// Tracking Problem Veto
      // TCut TrackingProblemVeto("Sum$(( abs(pfjets->Eta()) > 0.9 && abs(pfjets->Eta()) < 1.9 && (pfjets_chm - pfjets_neu) > 40)) == 0");
      // if any one jet meet the problem of the above, veto the whole event.
      bool trkProbVeto = false;

      bool tpv = TrackingProblemVeto == 1; 

      bool isZmet = false;
      bool isHLT1 = false;
      bool isHLT2 = false;

      int n_mus = 0;
      int n_els = 0;      
      // apply cuts to hypotheses
      for (unsigned int i = 0; i< els_p4().size(); i++){
      	//SELECTION
      	if (els_p4().at(i).pt() < 20)	        {	pt20Counter++;  continue; }
	if (fabs(els_p4().at(i).eta()) > 2.4)   {	etaCounter++;   continue; }
	if ( ! isNumeratorLepton(11, i) )       {	elIdCounter++;  continue; }

	isZmet = passElectronSelection_ZMet2012_v3_Iso(i);
	isHLT1 = els_HLT_Ele27_WP80()[i];
	isHLT2 = els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22()[i];

	// if ( isNumeratorLepton(11, i) && ( !passElectronSelection_ZMet2012_v3_Iso(i) )) {samesign++;}  
	// if ( (! isNumeratorLepton(11, i)) && ( passElectronSelection_ZMet2012_v3_Iso(i) )) {zmet++;}  
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
	if (fabs(mus_p4().at(i).eta()) > 2.1)	{	etaCounter++;   continue; }
	if ( ! isNumeratorLepton(13, i) )	{	muIdCounter++;  continue; }

	isZmet = muonId(i, ZMet2012_v1);
	isHLT1 = mus_HLT_IsoMu24_eta2p1()[i];
	isHLT2 = mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1()[i];

	// if ( ( isNumeratorLepton(13, i)) && (! muonId(i, ZMet2012_v1)) )	{	samesign++;  }
	// if ( (!isNumeratorLepton(13, i)) && muonId(i, ZMet2012_v1) )	        {	zmet++;  }
	// if ( !muonId(i, ZMet2012_v1)) continue;
     	// if ( muonIsoValuePF2012_deltaBeta(i) > 0.2 ) continue;
     	// if ( ! ETHSelection::goodMuon(i) ) continue;

	n_mus++;
	index = i;
      }
      
      if ((n_els + n_mus) != 1 ) continue;


      if( pfjets_p4().size() < 2) {less_jets++; continue;} // pre-selection check
      
      LorentzVector lrp4;

      int n_bTag =0;
      int n_jets = 0;

      if(n_els > 0)	  lrp4 = els_p4().at(index);
      if(n_mus > 0)	  lrp4 = mus_p4().at(index);

      for(unsigned int c = 0; c < pfjets_p4().size(); c++) {
      	float _jetPt = pfjets_p4().at(c).pt() * pfjets_corL1FastL2L3().at(c);
      	// jet pt times correction to jets
      	if(_jetPt < 30) continue;       // disgard those with small pt

      	if(fabs(pfjets_p4().at(c).eta()) > 2.5) continue; // Jet eta check

      	float dr_lr = DeltaR(pfjets_p4().at(c), lrp4 );
      	// delta R is the distance in the eta-phi plane
      	if(dr_lr < 0.4) continue;      // disgard small delta R 

      	float _bTag = pfjets_combinedSecondaryVertexBJetTag().at(c);
      	if ( _bTag > bTagDiscriminant)  n_bTag++;	 // mark as the has qualified b-jet, _bTag should be 0~1

	if (( abs(pfjets->Eta()) > 0.9 && abs(pfjets->Eta()) < 1.9 &&
	      (pfjets_chargedMultiplicity() - pfjets_neutralMultiplicity()) > 40) == 0) 
	  trkProbVeto = true;
	
      	n_jets++;					 // as qualified jet multiplicity
      }//pfjets_p4().size()
      
      if(n_jets < 2 || n_bTag == 0) {less_jets++; continue;}
      // if(trkProbVeto )     continue;

      goodCounter++;

      // create the ntuple
      InitBabyNtuple();
      // set the branch values

      pf_sumet = evt_pfsumet(); //float
      pf_tcmet = evt_pf_tcmet();  //float
      pf_tcmetPhi = evt_pf_tcmetPhi(); // float

      if(n_els > 0)	  lr_id = els_charge().at(index) * -11;
      if(n_mus > 0)	  lr_id = mus_charge().at(index) * -13;

      // njets = n_jets;
      // nbTag = n_bTag;

      met = evt_pfmet_type1cor();
      metPhi = evt_pfmetPhi_type1cor();
      
      lr_index = index;
      isZmetLep = isZmet;
      isHLT1Lep = isHLT1;
      isHLT2Lep = isHLT2;

      jets_p4 = pfjets_p4();
      jets_p4Correction = pfjets_corL1FastL2L3();
      btagDiscriminant = pfjets_combinedSecondaryVertexBJetTag();

      if(abs(lr_id) == 11)	  lr_p4 = els_p4().at(lr_index);
      if(abs(lr_id) == 13)	  lr_p4 = mus_p4().at(lr_index);

      if(abs(lr_id) == 11)	  lr_p4 = els_p4().at(lr_index);
      if(abs(lr_id) == 13)	  lr_p4 = mus_p4().at(lr_index);

      if(abs(lr_id) == 11)	  lr_charge = els_charge().at(lr_index);
      if(abs(lr_id) == 13)	  lr_charge = mus_charge().at(lr_index);
      
      if(abs(lr_id) == 11)	  lr_vtx_p4 = els_vertex_p4().at(lr_index);
      if(abs(lr_id) == 13)	  lr_vtx_p4 = mus_vertex_p4().at(lr_index);
      

      file = Form("%s", currentFile->GetTitle());

      eventNumber = evt_event();
      runNumber = evt_run();
      lumiBlock = evt_lumiBlock();
      scale_1fb = evt_scale1fb();
      isRealData = evt_isRealData();
      nvtxs = evt_nvtxs();
      trackingProblemVeto = trkProbVeto;

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
    myfile << "# of cuts of samesign vs zmet : " << samesign << "   " << zmet << endl;

    myfile << "# of events cut at muon ID/ISO :" << muIdCounter << endl;
    myfile << "# of events cut at electron ID/ISO :" << elIdCounter << endl;
  }

  return;
}
