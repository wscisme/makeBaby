doBaby_test(int numEvents = 0, std::string file_name = "file_name")
{
  gSystem->AddIncludePath(Form("-I%s/CORE", gSystem->Getenv("HOME")));
  gSystem->Load(Form("%s/CORE/libCMS2NtupleMacrosCORE.so", gSystem->Getenv("HOME")));
  //gSystem->Load("/home/users/cgeorge/CORE/libCMS2NtupleMacrosCORE.so");
  // gROOT->ProcessLine(".L ETHLeptonSelections.cc+");
  gROOT->ProcessLine(".L IsoTrackVeto.cc+");

  gROOT->ProcessLine(".L makeBaby.C+");

  babyMaker *baby = new babyMaker();
  
  TChain *sl = new TChain("Events"); 

  //sl->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-28/merged_ntuple_1.root");
  //sl->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleElectron_Run2012B-13Jul2012-v1_AOD/merged/merged_ntuple_1.root");
  sl->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/merged_ntuple_1.root");
  // sl->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_HT-200To250_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/V05-03-28/merged_ntuple_1.root");
  // sl->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_HT-250To300_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-28/*.root");
  // sl->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_HT-300To400_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-28/*.root");
  // sl->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-28/*.root");

  baby->ScanChain(sl, "babytest", numEvents); 
 
}
