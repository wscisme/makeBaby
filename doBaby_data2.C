doBaby_data2(int numEvents = 0, std::string file_name = "file_name")
{
  gSystem->AddIncludePath(Form("-I%s/CORE", gSystem->Getenv("HOME")));
  gSystem->Load(Form("%s/CORE/libCMS2NtupleMacrosCORE.so", gSystem->Getenv("HOME")));
  //gSystem->Load("/home/users/cgeorge/CORE/libCMS2NtupleMacrosCORE.so");
  //gROOT->ProcessLine(".L ETHLeptonSelections.cc+");
  gROOT->ProcessLine(".L IsoTrackVeto.cc+");
  gROOT->ProcessLine(".L makeBaby.C+");

  babyMaker *baby = new babyMaker();

  TChain *data = new TChain("Events");

  //data->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleElectron_Run2012A-13Jul2012-v1_AOD/merged/*.root");
  //data->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root");
  data->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleElectron_Run2012B-13Jul2012-v1_AOD/merged/*.root");
  // data->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/*.root");
  // data->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/*.root");
  // data->Add("/hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-24/SingleMu_Run2012B-13Jul2012-v1_AOD/merged/*.root");

  baby->ScanChain(data, "data_sl_8_2", numEvents); 
}
