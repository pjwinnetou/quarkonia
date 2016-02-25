#include <ctime>

#include "cutsAndBin.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"
#include <TLorentzVector.h>
#include "TriggerManipulation.h" 
#include "commonUtility.h"
static const long MAXTREESIZE = 10000000000;

TString getDayAndTime();
bool isTrackMatched(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) ;
void onia2ySkim( int nevt = 2000, int fileID = kPPDATA, int trigId=kL1DoubleMu0, int epSelection = kEPOppositeHF, bool saveTracks=false, TString skimVersion="unIdentified"  ) {

  using namespace std;
    
  bool isMC = false; 
  if ( (fileID == kPPMC) || (fileID == kPAMC) || (fileID == kAAMC) )
    isMC = true; 
  TFile *f1;
  if (fileID == kPPDATA)
    f1 = new TFile("Onia5TeV/ppData/OniaTree_DoubleMu_Run2015E-PromptReco-v1_Run_262157_262328.root"); 
  else if  (fileID == kAADATA) 
    f1 = new TFile("Onia5TeV/PbPbData/OniaTree_HIOniaL1DoubleMu0ABCD_HIRun2015-PromptReco-v1_Run_262548_263757.root");
  else if  (fileID == kAADATAPeri) 
    f1 = new TFile("Onia5TeV/PbPbData/OniaTree_HIOniaPeripheral30100_HIRun2015-PromptReco-v1_Run_262548_263757.root"); 
  else if  (fileID == kAADATACentL3) 
    f1 = new TFile("Onia5TeV/PbPbData/OniaTree_HIOniaCentral30L2L3_HIRun2015-PromptReco-v1_Run_262548_263757.root");  // Jan 29th
  else if  (fileID == kAAMC)
    f1 = new TFile("Onia5TeV/PbPbMC/OniaTree_Ups1S2SMM_5p02TeV_TuneCUETP8M1_ptUps1S2S06_noCUT.root");
  else if  (fileID == kPPMC)
    f1 = new TFile("Onia5TeV/ppMC/OniaTree_Ups1S2SMM_pp5p02TeV_TuneCUETP8M1_noCUT_v2.root");  
  
  TTree *mytree = (TTree*)f1->Get("hionia/myTree");
  if ( mytree == 0 ){ 
    cout << " it's not the oniaAndFriends tree made by TFileServie.. Taking myTree instead of hionia/myTree" << endl;
    mytree = (TTree*)f1->Get("myTree");
  }
  TTree *trkTree;
  if ( (fileID == kPPDATA) || (fileID == kPADATA) || (fileID == kPPMC) || (fileID == kPAMC) )
    trkTree = (TTree*)f1->Get("ppTrack/trackTree");
  else 
    trkTree = (TTree*)f1->Get("anaTrack/trackTree");
  if ( trkTree == 0 )  {
    cout << " There is no track collection.  Proceeding w/o it...." << endl;
    saveTracks = false ;
  }
  
  
  // *==*==*==*==*==*==* Trigger selection *==*==*==*==*==*==* //
  TString trigName = getTrig(trigId);
  cout << endl << endl << "*==*==*==*==*==*==* Trigger selection  *==*==*==*==*==*==*" << endl;
  cout << " trigger selection : " << trigName << endl;
  TH1F* hStats = (TH1F*)f1->Get("hionia/hStats");
  if ( hStats == 0 ){
    cout << " it's not the oniaAndFriends tree made by TFileServie.. hionia/hStats -> hStats myTree instead of hionia/myTree" << endl;
    hStats = (TH1F*)f1->Get("hStats");
  }
  

  hltIndex hltBits = getTrigIndex(trigId, hStats);
  cout << endl << endl << "*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*" << endl;
  
  
  // *==*==*==*==*==*==* Event Plane  *==*==*==*==*==*==* //
  TString epName = getEPSel(epSelection) ;
  cout << " Event Plane : " << epName << endl;
  
  
  // *==*==*==*==*==*==* Output file  *==*==*==*==*==*==* //
  TFile* newfile;
  if (fileID == kPPDATA)
    newfile = new TFile(Form("skimmedFiles/yskimPP_L1DoubleMu0PD_Trig-%s_%s_%s.root",trigName.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  else if (fileID == kAADATA) 
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_L1DoubleMu0PD_Trig-%s_EP-%s_%s_%s.root",trigName.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  else if  (fileID == kAADATACentL3) 
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_CentralPD_Trig-%s_EP-%s_%s_%s.root",trigName.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  else if  (fileID == kAADATAPeri) 
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_PeripheralPD_Trig-%s_EP-%s_%s_%s.root",trigName.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  else if (fileID == kPPMC)
    newfile = new TFile(Form("skimmedFiles/yskimPP_MC_Trig-%s_EP-%s_%s_%s.root",trigName.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  else if (fileID == kAAMC)
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_MC_Trig-%s_EP-%s_%s_%s.root",trigName.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 

   
  // import the tree to the RooDataSet
  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  ULong64_t       Reco_QQ_trig[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[200];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_4mom;   //!
  TBranch        *b_Reco_QQ_mumi_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_QQ_mupl_highPurity[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_highPurity[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_highPurity;   //!
  TBranch        *b_Reco_QQ_mumi_highPurity;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
  mytree->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);


  
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  //  mytree->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  //  mytree->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Int_t           Reco_QQ_mupl_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_QQ_mupl_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_normChi2_global[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
  TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
  mytree->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_QQ_mupl_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nMuValHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
  TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_QQ_mupl_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_StationsMatched[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
  TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
  mytree->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_QQ_mupl_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dxy[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_dxy;   //!
  TBranch        *b_Reco_QQ_mumi_dxy;   //!
  TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
  TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  mytree->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_QQ_mupl_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dz[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_dz;   //!
  TBranch        *b_Reco_QQ_mumi_dz;   //!
  TBranch        *b_Reco_QQ_mupl_dzErr;   //!
  TBranch        *b_Reco_QQ_mumi_dzErr;   //!
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  mytree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  mytree->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_QQ_mupl_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_QQ_mupl_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_mu_TMOneStaTight[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
  TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
  TBranch        *b_Reco_mu_TMOneStaTight;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
  mytree->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_QQ_mupl_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nPixWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_QQ_mupl_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
  Int_t           Reco_QQ_mumi_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
  Float_t         Reco_QQ_mupl_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
  Float_t         Reco_QQ_mumi_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);

  /////////////////////////////////////////
  ////// Gen QQ 
  /////////////////////////////////////////
  Int_t           Gen_QQ_size;
  Int_t           Gen_QQ_type[200];   //[Gen_QQ_size]
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_QQ_mupl_4mom;
  TClonesArray    *Gen_QQ_mumi_4mom;
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_QQ_type;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_QQ_mupl_4mom;   //!
  TBranch        *b_Gen_QQ_mumi_4mom;   //!
  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;
  if (isMC) { 
    mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    mytree->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
    mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    mytree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
    mytree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
  }
  
  
  


  
  const int NMAXTRK = 10000;
  Int_t           nTrk;
  Float_t         trkPt[NMAXTRK];   //[nTrk]
  Float_t         trkPtError[NMAXTRK];   //[nTrk]
  UChar_t         trkNHit[NMAXTRK];   //[nTrk]
  Float_t         trkEta[NMAXTRK];   //[nTrk]
  Float_t         trkPhi[NMAXTRK];   //[nTrk]
  Int_t           trkCharge[NMAXTRK];   //[nTrk]
  Bool_t          highPurity[NMAXTRK];   //[nTrk]
  Bool_t          tight[NMAXTRK];   //[nTrk]
  Bool_t          loose[NMAXTRK];   //[nTrk]
  Float_t         trkChi2[NMAXTRK];   //[nTrk]
  UChar_t         trkNdof[NMAXTRK];   //[nTrk]
  Float_t         trkDxy1[NMAXTRK];   //[nTrk]
  Float_t         trkDxyError1[NMAXTRK];   //[nTrk]
  Float_t         trkDz1[NMAXTRK];   //[nTrk]
  Float_t         trkDzError1[NMAXTRK];   //[nTrk]
  Bool_t          trkFake[NMAXTRK];   //[nTrk]
  UChar_t         trkAlgo[NMAXTRK];   //[nTrk]
  Float_t         trkMVA[NMAXTRK];   //[nTrk]
  
  TBranch        *b_nTrk;   //!
  TBranch        *b_trkPt;   //!
  TBranch        *b_trkPtError;   //!
  TBranch        *b_trkNHit;   //!
  TBranch        *b_trkEta;   //!
  TBranch        *b_trkPhi;   //!
  TBranch        *b_trkCharge;   //!
  TBranch        *b_highPurity;   //!
  TBranch        *b_tight;   //!
  TBranch        *b_loose;   //!
  TBranch        *b_trkChi2;   //!
  TBranch        *b_trkNdof;   //!
  TBranch        *b_trkDxy1;   //!
  TBranch        *b_trkDxyError1;   //!
  TBranch        *b_trkDz1;   //!
  TBranch        *b_trkDzError1;   //!
  TBranch        *b_trkFake;   //!
  TBranch        *b_trkAlgo;   //!
  TBranch        *b_trkMVA;   //!
  if (saveTracks) { 
    trkTree->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
    trkTree->SetBranchAddress("trkPt", trkPt, &b_trkPt);
    trkTree->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
    trkTree->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
    trkTree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
    trkTree->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
    trkTree->SetBranchAddress("trkCharge", trkCharge, &b_trkCharge);
    trkTree->SetBranchAddress("highPurity", highPurity, &b_highPurity);
    trkTree->SetBranchAddress("tight", tight, &b_tight);
    trkTree->SetBranchAddress("loose", loose, &b_loose);
    trkTree->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
    trkTree->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
    trkTree->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
    trkTree->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
    trkTree->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
    trkTree->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
    trkTree->SetBranchAddress("trkFake", trkFake, &b_trkFake);
    trkTree->SetBranchAddress("trkAlgo", trkAlgo, &b_trkAlgo);
    trkTree->SetBranchAddress("trkMVA", trkMVA, &b_trkMVA);
  }

  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  dimuon tree 
  ////////////////////////////////////////////////////////////////////////
  DiMuon dm;
  TTree *mmTree = new TTree("mm","diMuonPairs");
  mmTree->SetMaxTreeSize(MAXTREESIZE);
  mmTree->Branch("mm",&dm.run,branchString.Data());
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* dphiEp2Var   = new RooRealVar("dphiEp2","Delta Phi from 2nd order event plane", -100,100,"");

  //  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *ep2Var, *pt1Var, *eta1Var, *pt2Var, *eta2Var);
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var);
  if ( (fileID == kAAMC) || (fileID == kAADATA) || (fileID == kAADATAPeri) || (fileID == kAADATACentL3) )
    { argSet->add(*ep2Var); argSet->add(*dphiEp2Var); }
  if ( (fileID == kAAMC) || (fileID == kAADATA) || (fileID == kAADATAPeri) || (fileID == kAADATACentL3) || (fileID == kPAMC) || (fileID == kPADATA) )
    argSet->add(*cBinVar);

  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);

  ////////////////////////////////////////////////////////////////////////
  //////////////////  Gen QQ tree        
  ////////////////////////////////////////////////////////////////////////
  DiMuon dmGen;
  TTree *mmGenTree;
  if (isMC)  {
    mmGenTree = new TTree("mmGen","Gen Di-muon Pairs");
    mmGenTree->SetMaxTreeSize(MAXTREESIZE);
    mmGenTree->Branch("mmGen",&dmGen.run,branchString.Data());
  }
    
  ////////////////////////////////////////////////////////////////////////
  //////////////////  Track tree        
  ////////////////////////////////////////////////////////////////////////
  TTree *trkOutTree = new TTree("track","trackTree");
  static const int MAXTRK  = 10000;   // This must be  enough.
  int nTrkOut;
  float trkOutPt[MAXTRK];
  float trkOutEta[MAXTRK];
  float trkOutPhi[MAXTRK];
  float trkOutDphi[MAXTRK];
  float trkOutDeta[MAXTRK];
  float trkOutDr[MAXTRK];
  float trkPiTripleMass[MAXTRK];
  float trkKaTripleMass[MAXTRK];
  float trkOutCorr[MAXTRK];
  trkOutTree->SetMaxTreeSize(MAXTREESIZE);
  trkOutTree->Branch("nTrack",&nTrkOut,"nTrack/I");
  trkOutTree->Branch("pt",trkOutPt,"pt[nTrack]/F");
  trkOutTree->Branch("eta",trkOutEta,"eta[nTrack]/F");
  trkOutTree->Branch("phi",trkOutPhi,"phi[nTrack]/F");
  trkOutTree->Branch("dphi",trkOutDphi,"dphi[nTrack]/F");
  trkOutTree->Branch("deta",trkOutDeta,"deta[nTrack]/F");
  trkOutTree->Branch("dr",trkOutDr,"dr[nTrack]/F");
  trkOutTree->Branch("corr",trkOutCorr,"corr[nTrack]/F");
  trkOutTree->Branch("massPiTriple",trkPiTripleMass,"massPiTriple[nTrack]/F");
  trkOutTree->Branch("massKaTriple",trkKaTripleMass,"massKaTriple[nTrack]/F");
  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  Event selection tree 
  ////////////////////////////////////////////////////////////////////////
  TH1D* hEvent = new TH1D("hFilter","",20 ,0.5, 21.5);

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  TLorentzVector* JP_Gen = new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;
  
  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  for(int iev=0; iev<nevt ; ++iev){
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;
    
  
    hEvent->GetXaxis()->SetBinLabel(1,"Events total");     hEvent->Fill(1);  
    ///////////////////////////
    ///// Call the values /////
    ///////////////////////////
    mytree->GetEntry(iev);
    if (saveTracks)  trkTree->GetEntry(iev);
        
    ///////////////////////////
    ///// Vertex cut //////////
    ///////////////////////////
    if (TMath::Abs(zVtx) > 15.) continue;
    hEvent->GetXaxis()->SetBinLabel(2,"Events z vertex < 15cm");    hEvent->Fill(2);
    
    //     if ( !((trigId == kL1DoubleMu0) && (fileID==kAADATA)) && (  isTrigMatched(hltBits,HLTriggers) == false ) )  // Andre fixed the L1_DoubleMu0ABCD problem on Feb 10-ish! 
    if ( isTrigMatched(hltBits,HLTriggers) == false ) 
      continue; // trigger selection
    
    
    hEvent->GetXaxis()->SetBinLabel(3,Form("Events Trig %s",trigName.Data()));    hEvent->Fill(3);
    
    
    dm.clear();
    dm.run = runNb;
    dm.lumi = LS ;
    dm.event = eventNb ;
    dm.vz = zVtx;
    dm.cBin = Centrality ;
    
    
    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) {
      hEvent->GetXaxis()->SetBinLabel(4,"Di-muons Total");      hEvent->Fill(4);
      //struct Condition Jpsi_Reco; 
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_QQ_mupl_4mom->At(irqq);
      mumi_Reco = (TLorentzVector*) Reco_QQ_mumi_4mom->At(irqq);
      //       if ( !  ( acceptance( mupl_Reco->Pt(), mupl_Reco->Eta(), mupl_Reco->P()) && acceptance( mumi_Reco->Pt(), mumi_Reco->Eta(), mumi_Reco->P())) )
      
      if ( !  ( acceptance( mupl_Reco->Pt(), mupl_Reco->Eta() ) && acceptance( mumi_Reco->Pt(), mumi_Reco->Eta()) )   )
	continue;      
      hEvent->GetXaxis()->SetBinLabel(5,"Di-muons Accep");      hEvent->Fill(5);
      
      // vertex probablitiy cut:
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
	continue;
      hEvent->GetXaxis()->SetBinLabel(6,"Di-muons Vtx prob.");      hEvent->Fill(6);
      
      
      if (  isTrigMatched(hltBits,Reco_QQ_trig[irqq]) == false )
	continue; // trigger selection
      hEvent->GetXaxis()->SetBinLabel(7,"Di-muons trig");      hEvent->Fill(7);
      
      // muon id cut :  tight cut  pixel cut is missing!!
      // Reference : https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Soft_Muon
      bool muplSoft = ( (Reco_QQ_mupl_TMOneStaTight[irqq]==true) &&
			(Reco_QQ_mupl_nTrkWMea[irqq] > 5) &&
			(Reco_QQ_mupl_nPixWMea[irqq] > 0) &&
			(Reco_QQ_mupl_dxy[irqq]<0.3) &&
			(Reco_QQ_mupl_dz[irqq]<20.)   			//			 &&  (Reco_QQ_mupl_highPurity[irqq]==true) 
			) ; 
      
      bool mumiSoft = ( (Reco_QQ_mumi_TMOneStaTight[irqq]==true) &&
			(Reco_QQ_mumi_nTrkWMea[irqq] > 5) &&
			(Reco_QQ_mumi_nPixWMea[irqq] > 0) &&
			(Reco_QQ_mumi_dxy[irqq]<0.3) &&
			(Reco_QQ_mumi_dz[irqq]<20.)  // && (Reco_QQ_mumi_highPurity[irqq]==true)
			) ; 
      
      bool muplHighPtCut = ( (Reco_QQ_mupl_nMuValHits[irqq]>0) &&
			     (Reco_QQ_mupl_StationsMatched[irqq]>1) &&
			     (Reco_QQ_mupl_ptErr_global[irqq]/mupl_Reco->Pt() < 0.3 ) && 
			     (Reco_QQ_mupl_dxy[irqq]<0.2) &&
			     (Reco_QQ_mupl_dz[irqq] <0.5) &&
			     (Reco_QQ_mupl_nPixValHits[irqq] > 0 ) &&
			     (Reco_QQ_mupl_nTrkWMea[irqq] > 5) 
			     );
      bool mumiHighPtCut = ( (Reco_QQ_mumi_nMuValHits[irqq]>0) &&
			     (Reco_QQ_mumi_StationsMatched[irqq]>1) &&
			     (Reco_QQ_mumi_ptErr_global[irqq]/mumi_Reco->Pt() < 0.3 ) && 
			     (Reco_QQ_mumi_dxy[irqq]<0.2) &&
			     (Reco_QQ_mumi_dz[irqq] <0.5) &&
			     (Reco_QQ_mumi_nPixValHits[irqq] > 0 ) &&
			     (Reco_QQ_mumi_nTrkWMea[irqq] > 5) 
			     );
      
      
      if ( !(muplSoft && mumiSoft) ) 
	continue;   
      hEvent->GetXaxis()->SetBinLabel(8,"Di-muons mu ID");      hEvent->Fill(8);
      if ( Reco_QQ_sign[irqq] != 0 ) 
	continue;
      hEvent->GetXaxis()->SetBinLabel(9,"Di-muoons charge sign");      hEvent->Fill(9);
      
      dm.softFlag = 0;
      if (muplSoft && mumiSoft) 
	dm.softFlag = 1;
      
      dm.highPtFlag = 0;
      if ( muplHighPtCut && mumiHighPtCut ) 
	dm.highPtFlag = 1;
      
      dm.mass   = JP_Reco->M();
      dm.pt     = JP_Reco->Pt();
      dm.phi    = JP_Reco->Phi();
      dm.y      = JP_Reco->Rapidity();
      dm.eta      = JP_Reco->Eta();
      dm.pt1  = mupl_Reco->Pt();
      dm.eta1 = mupl_Reco->Eta();
      dm.phi1 = mupl_Reco->Phi();
      dm.pt2  = mumi_Reco->Pt();
      dm.eta2 = mumi_Reco->Eta();
      dm.phi2 = mumi_Reco->Phi();
      
      if      ( epSelection == kEPl2HF )            dm.ep2 = rpAng[8];  
      else if ( epSelection == kEPOppositeHF )      {
	if ( dm.y > 0 )     dm.ep2 = rpAng[6];  // [6] = HFm2;
	else                 dm.ep2 = rpAng[7];  // [7] = HFm2;
      }
      else if ( epSelection == kEPSameSideHF )      {
	if ( dm.y > 0 )     dm.ep2 = rpAng[7];  
	else                 dm.ep2 = rpAng[6]; 
      }
      
      dm.dphiEp2 = getDPHI( dm.phi, dm.ep2);
      
      dm.oniaIndex = irqq;
      mmTree->Fill();
      
      massVar->setVal( (double)dm.mass ) ;
      ptVar->setVal(   (double)dm.pt   ) ;
      yVar->setVal(    (double)dm.y    ) ;
      dphiEp2Var->setVal(   (double)dm.dphiEp2  ) ;
      pt1Var->setVal(  (double)dm.pt1  ) ;
      eta1Var->setVal( (double)dm.eta1 ) ;
      pt2Var->setVal(  (double)dm.pt2  ) ;
      eta2Var->setVal( (double)dm.eta2 ) ;
      ep2Var->setVal( (double)dm.ep2 ) ;
      cBinVar->setVal( (double)dm.cBin ) ;
      dataSet->add( *argSet);

      /////////////////////////////////////////////////////
      // Track tree.  This block must be in the onia loop
      /////////////////////////////////////////////////////
      nTrkOut = 0;
      for ( int ipart=0  ;  ipart<nTrk  ;  ipart++) { 
	if ( (fileID == kPPDATA)||(fileID == kPPMC) ) {
	  //	 cout << "trkNHit[ipart] =  " << (int)trkNHit[ipart] << endl;
	  if ( highPurity[ipart] == false)    continue;
	  if(( trkMVA[ipart]<0.5 &&  trkMVA[ipart]!=-99) || (int)trkNHit[ipart]<8) continue;
	  if( fabs( trkDz1[ipart]/ trkDzError1[ipart])>3 ||  fabs(trkDxy1[ipart]/ trkDxyError1[ipart])>3) continue;
	  if( trkPtError[ipart]/ trkPt[ipart]>0.3) continue;
	}
	else 
	  continue;
	////////////////    Matching with muons?
	if ( isTrackMatched(dm.pt1, dm.eta1, dm.phi1,  trkPt[ipart], trkEta[ipart], trkPhi[ipart]) )
	  continue;
	if ( isTrackMatched(dm.pt2, dm.eta2, dm.phi2,  trkPt[ipart], trkEta[ipart], trkPhi[ipart]) )
	  continue;
	
	trkOutPt[nTrkOut] = trkPt[ipart];
	trkOutEta[nTrkOut] = trkEta[ipart];
	trkOutPhi[nTrkOut] = trkPhi[ipart];
	trkOutDphi[nTrkOut] = getDPHI(trkPhi[ipart], dm.phi);
	trkOutDeta[nTrkOut] = trkEta[ipart] - dm.eta;
	trkOutDr[nTrkOut] = sqrt ( trkOutDphi[nTrkOut] * trkOutDphi[nTrkOut] + trkOutDeta[nTrkOut]*trkOutDeta[nTrkOut] ) ;
	
	// triplet mass :
	float selectMass=pdgMass.Y1S; // tentative
	TLorentzVector theQQ, thePion, theKaon, piTriple, kaTriple;
	theQQ.SetPtEtaPhiM( dm.pt, dm.eta, dm.phi, selectMass );
	thePion.SetPtEtaPhiM( trkPt[ipart], trkEta[ipart], trkPhi[ipart], pdgMass.PiPlus );
	theKaon.SetPtEtaPhiM( trkPt[ipart], trkEta[ipart], trkPhi[ipart], pdgMass.KaPlus );
	piTriple = thePion + theQQ ;
	kaTriple = theKaon + theQQ ;
	trkPiTripleMass[ipart] = piTriple.M();
	trkKaTripleMass[ipart] = kaTriple.M();
	nTrkOut++;
      }
      
      if (saveTracks) trkOutTree->Fill();
    } // end of dimuon loop


    /////////////////////////////////////////////////////
    /////////  gen QQ loop first
    /////////////////////////////////////////////////////
    if (isMC) { 
      dmGen.clear();
      dmGen.run = runNb;
      dmGen.lumi = LS ;
      dmGen.event = eventNb ;
      dmGen.vz = zVtx;
      dmGen.cBin = Centrality ;
      for (Int_t irqq=0; irqq<Gen_QQ_size; ++irqq) {
	JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(irqq);
	mupl_Gen = (TLorentzVector*) Gen_QQ_mupl_4mom->At(irqq);
	mumi_Gen = (TLorentzVector*) Gen_QQ_mumi_4mom->At(irqq);
	dmGen.mass   = JP_Gen->M();
	dmGen.pt     = JP_Gen->Pt();
	dmGen.phi    = JP_Gen->Phi();
	dmGen.y      = JP_Gen->Rapidity();
	dmGen.eta      = JP_Gen->Eta();
	dmGen.pt1  = mupl_Gen->Pt();
	dmGen.eta1 = mupl_Gen->Eta();
	dmGen.phi1 = mupl_Gen->Phi();
	dmGen.pt2  = mumi_Gen->Pt();
	dmGen.eta2 = mumi_Gen->Eta();
	dmGen.phi2 = mumi_Gen->Phi();
	
	dmGen.oniaIndex = irqq;
	mmGenTree->Fill();
      } //end of Gen QQ tree
    } //end of isMC condition for Gen QQ tree
    
  } //end of event loop
  
  dataSet->Write();
  mmTree->Write();  // Don't need to call Write() for trees
  hEvent->Write();
  if (isMC) mmGenTree->Write();
  //   if (saveTracks) trkOutTree->Write();
  newfile->Close();
} 



/*
      Jpsi_Reco.theMass = JP_Reco->M();
      Jpsi_Reco.theRapidity = JP_Reco->Rapidity(); // y_{lab}	
      Jpsi_Reco.thePt = JP_Reco->Pt();
      Jpsi_Reco.thePhi = JP_Reco->Phi();
      
      Jpsi_Reco.mupl_p = sqrt( (mupl_Reco->Px())*(mupl_Reco->Px()) + (mupl_Reco->Py())*(mupl_Reco->Py()) + (mupl_Reco->Pz())*(mupl_Reco->Pz()) );
      Jpsi_Reco.mumi_p = sqrt( (mumi_Reco->Px())*(mumi_Reco->Px()) + (mumi_Reco->Py())*(mumi_Reco->Py()) + (mumi_Reco->Pz())*(mumi_Reco->Pz()) );
      Jpsi_Reco.mupl_pt = mupl_Reco->Pt();
      Jpsi_Reco.mumi_pt = mumi_Reco->Pt();
      Jpsi_Reco.mupl_eta = mupl_Reco->Eta();
      Jpsi_Reco.mumi_eta = mumi_Reco->Eta();
      */    

/*  ////// Gen_QQ_size loop
    for (Int_t igqq=0; igqq<Gen_QQ_size; ++igqq) {
    //struct Condition Jpsi_Gen; 
    mupl_Gen = (TLorentzVector*) Gen_QQ_mupl_4mom->At(igqq);
    mumi_Gen = (TLorentzVector*) Gen_QQ_mumi_4mom->At(igqq);
    JP_Gen_tmp_qq = (TLorentzVector*) Gen_QQ_4mom->At(igqq); // Gen Jpsi (for filling NoCut)
    *JP_Gen_tmp = *mupl_Gen +  *mumi_Gen; // Gen dimuon pairs (for filling NoCut)
	  
	  // variables only used for Reco
	  Jpsi_Gen.HLTriggers = 0;
	  Jpsi_Gen.Reco_QQ_trig  = 0;
	  Jpsi_Gen.theSign = 0; //already +- pair
	  
	  Jpsi_Gen.theCentrality = 97.5; // for pp!
	  Jpsi_Gen.theType = Gen_QQ_type[igqq]; // PR or NP
	  Jpsi_Gen.mupl_p = sqrt( (mupl_Gen->Px())*(mupl_Gen->Px()) + (mupl_Gen->Py())*(mupl_Gen->Py()) + (mupl_Gen->Pz())*(mupl_Gen->Pz()) );
	  Jpsi_Gen.mumi_p = sqrt( (mumi_Gen->Px())*(mumi_Gen->Px()) + (mumi_Gen->Py())*(mumi_Gen->Py()) + (mumi_Gen->Pz())*(mumi_Gen->Pz()) );
	  Jpsi_Gen.mupl_pt = mupl_Gen->Pt();
	  Jpsi_Gen.mumi_pt = mumi_Gen->Pt();
	  Jpsi_Gen.mupl_eta = mupl_Gen->Eta();
	  Jpsi_Gen.mumi_eta = mumi_Gen->Eta();
	  
	  // --- cut01 : GEN - No cut (only DimuCut = +-pair from 443)
	  h2D_NoCut_Gen_pt_y->Fill(JP_Gen_tmp->Rapidity(),JP_Gen_tmp->Pt()); // Gen dimuon
	  h2D_NoCut_GenJpsi_pt_y->Fill(JP_Gen_tmp_qq->Rapidity(),JP_Gen_tmp_qq->Pt()); // Gen Jpsi
	  
	  // --- cut02 : GEN for denominator
	  bool yn_gen = false;
	  if ( kineCut(Jpsi_Gen.mupl_pt, Jpsi_Gen.mupl_eta, Jpsi_Gen.mupl_p) 
			&& kineCut(Jpsi_Gen.mumi_pt, Jpsi_Gen.mumi_eta, Jpsi_Gen.mumi_p)) {
			*JP_Gen = *mupl_Gen +  *mumi_Gen; // used for actual denominator ( GEN dimuon pairs)
			Jpsi_Gen.theMass = JP_Gen->M();
			Jpsi_Gen.theRapidity = JP_Gen->Rapidity();	
			Jpsi_Gen.thePt = JP_Gen->Pt();
			Jpsi_Gen.thePhi = JP_Gen->Phi();
			if ( massCut1(Jpsi_Gen.theMass)  
			&& minpt<=Jpsi_Gen.thePt && Jpsi_Gen.thePt <maxpt 
			&& minylab<=Jpsi_Gen.theRapidity && Jpsi_Gen.theRapidity <maxylab) {
			yn_gen = true;
			}
			}
			if (yn_gen == true){
			h2D_Den_pt_y->Fill(Jpsi_Gen.theRapidity,Jpsi_Gen.thePt,zWeight01*zWeight02);
			} // end of yn_gen
			} //end of Gen_QQ_size loop
      */



TString getDayAndTime() { 
  time_t currentTime;
  struct tm *localTime;
  
  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time
  
  int Day    = localTime->tm_mday;
  int Month  = localTime->tm_mon + 1;
  int Year   = localTime->tm_year + 1900;
  int Hour   = localTime->tm_hour;
  int Min    = localTime->tm_min;
  //  int Sec    = localTime->tm_sec;
  return Form("%d%d%d%d%d",Year,Month,Day,Hour,Min);
}


bool isTrackMatched(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
  if (  (pt1-pt2)/pt1 > 0.02 ) 
    return false;
  if ( fabs( getDPHI(phi1, phi2) ) > 0.1 )
    return false;
  if ( fabs( eta1 - eta2 ) > 0.1 )
    return false;

  return true;
}
