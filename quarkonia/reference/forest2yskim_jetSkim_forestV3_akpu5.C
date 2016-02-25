///////////////////////////////////////////////////////////////////                                
// forest2yskim.C                                                //                                                 
// Creator : Yongsun Kim (MIT), jazzitup@mit.edu                 //                                                 
// Function : Transform hiForest files into yskim file           //
// yskims for MinBias1, Minbias2 and photon jet skims            //
///////////////////////////////////////////////////////////////////         
///////////////////////////////////////////////////////////////////
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <TMath.h>
#include "../../hiForestV3/hiForest.h"
#include "../CutAndBinCollection2012.h"
#include <time.h>
#include <TRandom3.h>

using namespace std;

static const int MAXTRK  = 10000;   // This is very enough.
static const int MAXMTRK = 30000;   // Again this is very enough for 10 mixing  
static const int MAXMJET = 2000;
static const long MAXTREESIZE = 10000000000;


vector<jetKinem> nullVec;




void forest2yskim_jetSkim_forestV3_akpu5(TString inputFile_="forestFiles/pA/pA_photonSkimForest_v85_skimPhotonPt50_eta1.5.root",
				   std::string MinbiasFname = "skim_trackJet_minbiasTrackJet_mc.root",
				   float cutphotonPt  = 40, 
				   std::string outname = "testPhotonSkim.root",
				   sampleType colli=kPADATA,
				   bool doMix = false,
				   bool doJetResCorrection = 0,
				   float addJetEnergyRes = 0,
				   float addFlatJetEnergyRes = 0,
                                   float jetEnergyScale = 1,
				   bool useGenJetColl = 0
				   )
{ 
  bool isMC=true;
  if ((colli==kPPDATA)||(colli==kPADATA)||(colli==kHIDATA))
    isMC=false;
  
  int seconds = time(NULL);         cout << " time = " <<seconds%10000<< endl;
  TRandom3 rand(seconds%10000);
  TString datafname  = "";
  float cutphotonEta = 1.44;
  float preCutPhotonEt = 40;
   
  const int nMaxPho = 100;
  
  
  HiForest *c;
  if((colli==kPADATA)||(colli==kPAMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPPb, isMC );
  }
  else if  ((colli==kPPDATA)||(colli==kPPMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPP, isMC );
  }
  else if  ((colli==kHIDATA)||(colli==kHIMC)) {
    c = new HiForest(inputFile_.Data(), "forest", cPbPb, isMC );
    c->GetEnergyScaleTable("../photonEnergyScaleTable_lowPt_v6.root");
  }
  else {
    cout << " Error!  No such collision type" << endl;
    return;
  }
  
  c->InitTree();
  
  
  // now open new root file
  TFile* newfile_data = new TFile(outname.data(),"recreate");
  
  TTree* newtreePhoton;
  float newPt[nMaxPho];  // <<= temporary space                                                           
  int order[nMaxPho];
  float corrPt[nMaxPho];
  newtreePhoton = c->photonTree->CloneTree(0);
  newtreePhoton->SetName("yPhotonTree");
  newtreePhoton->SetMaxTreeSize(MAXTREESIZE);
  newtreePhoton->Branch("order",  order, "order[nPhotons]/I");
  newtreePhoton->Branch("corrPt", corrPt,"corrPt[nPhotons]/F");

  TTree* treeFullJet;
  if (   (colli==kPPDATA) || (colli==kPPMC) ) {
    treeFullJet = c->ak5jetTree->CloneTree(0);
    cout << "pp collision.  Using ak5PF Jet Algo" << endl<<endl;
  }
  else {
    treeFullJet = c->akPu5jetTree->CloneTree(0);
    cout << "pPb or PbPb collision. Using akPu5PF Jet Algo" << endl<<endl;
  }
  treeFullJet->SetName("fullJet");
  treeFullJet->SetMaxTreeSize(MAXTREESIZE);
   
  // jet tree! 
  int nJet;    
  float jetPt[500];
  float jetEta[500];
  float jetPhi[500];
  float jetDphi[500];
  int jetSubid[500];
  float jetRefPt[500];
  float jetRefEta[500];
  float jetRefPhi[500];
  float jetRefDphi[500];
  float jetRefPartonPt[500];
  int  jetRefPartonFlv[500];

  TTree *newtreeJet = new TTree("yJet","jets");
  newtreeJet->SetMaxTreeSize(MAXTREESIZE);
  newtreeJet->Branch("nJet",&nJet,"nJet/I");
  newtreeJet->Branch("pt",jetPt,"pt[nJet]/F");
  newtreeJet->Branch("eta",jetEta,"eta[nJet]/F");
  newtreeJet->Branch("phi",jetPhi,"phi[nJet]/F");
  newtreeJet->Branch("dphi",jetDphi,"dphi[nJet]/F");
  if ( isMC )  {
    newtreeJet->Branch("subid",jetSubid,"subid[nJet]/I");
    newtreeJet->Branch("refPt",jetRefPt,"refPt[nJet]/F");
    newtreeJet->Branch("refEta",jetRefEta,"refEta[nJet]/F");
    newtreeJet->Branch("refPhi",jetRefPhi,"refPhi[nJet]/F");
    newtreeJet->Branch("refDphi",jetRefDphi,"refDphi[nJet]/F");
    newtreeJet->Branch("refPartonPt",jetRefPartonPt,"refPartonPt[nJet]/F");
    newtreeJet->Branch("refPartonFlv",jetRefPartonFlv,"refPartonFlv[nJet]/I");
  }



  int nMjet;
  float mJetPt[MAXMJET];
  float mJetEta[MAXMJET];
  float mJetPhi[MAXMJET];
  float mJetDphi[MAXMJET];
  
  TTree * tmixJet = new TTree("mJet","Jet from minbias events");
  tmixJet->SetMaxTreeSize(MAXTREESIZE);
  tmixJet->Branch("nJet",&nMjet,"nJet/I");
  tmixJet->Branch("pt",mJetPt,"pt[nJet]/F");
  tmixJet->Branch("eta",mJetEta,"eta[nJet]/F");
  tmixJet->Branch("phi",mJetPhi,"phi[nJet]/F");
  tmixJet->Branch("dphi", mJetDphi, "dphi[nJet]/F");
  
  
  
  // Imb = Input MinBias events 
  EvtSel          evtImb;
  Int_t           nJetImb;
  Float_t         jetPtImb[100];
  Float_t         jetEtaImb[100];
  Float_t         jetPhiImb[100];
  TBranch        *b_evt;
  TBranch        *b_nJetImb;
  TBranch        *b_jetPtImb;
  TBranch        *b_jetEtaImb;
  TBranch        *b_jetPhiImb; 
  
  
  int nCentBins =  nCentBinSkim;
  if ((colli==kPADATA)||(colli==kPAMC)) {
    nCentBins = nCentBinSkimPA;
  } 

  TChain   *tjmb[100][nVtxBin+1];
  int nMB[100][nVtxBin+1] ; //= 199109;                                                                                          
  int mbItr[100][nVtxBin+1];
  if ( doMix ) {
    cout <<"  Tree initialization for MinBias mixing" << endl;
    for( int icent = 0 ; icent< nCentBins ; icent++) {
      for( int ivz = 1 ; ivz<=nVtxBin ; ivz++) {
	
	tjmb[icent][ivz] = new TChain(Form("trkAndJets_first_cBin2icent%d_ivz%d",icent,ivz));
	tjmb[icent][ivz]->Add(MinbiasFname.data());
	tjmb[icent][ivz]->SetBranchAddress("evt", &evtImb,&b_evt);
	tjmb[icent][ivz]->SetBranchAddress("nJet",   &nJetImb,   &b_nJetImb);
	tjmb[icent][ivz]->SetBranchAddress("jetPt",  &jetPtImb,  &b_jetPtImb);
	tjmb[icent][ivz]->SetBranchAddress("jetEta", &jetEtaImb, &b_jetEtaImb);
	tjmb[icent][ivz]->SetBranchAddress("jetPhi", &jetPhiImb, &b_jetPhiImb);	
	
	/*  tjmb[icent][ivz]->SetBranchAddress("nTrk", &nTrkImb,&b_nTrkImb);
	    tjmb[icent][ivz]->SetBranchAddress("trkPt", trkPtImb, &b_trkPtImb);
	    tjmb[icent][ivz]->SetBranchAddress("trkEta", trkEtaImb, &b_trkEtaImb);
	    tjmb[icent][ivz]->SetBranchAddress("trkPhi", trkPhiImb, &b_trkPhiImb);
	    tjmb[icent][ivz]->SetBranchAddress("trkWeight", trkWeightImb, &b_trkWeightImb);
	    if ( isMC) {
	    tjmb[icent][ivz]->SetBranchAddress("nGp", &nGpImb, &b_nGpImb);
	    tjmb[icent][ivz]->SetBranchAddress("gpPt", gpPtImb, &b_gpPtImb);
	    tjmb[icent][ivz]->SetBranchAddress("gpEta", gpEtaImb, &b_gpEtaImb);
	    tjmb[icent][ivz]->SetBranchAddress("gpPhi", gpPhiImb, &b_gpPhiImb);
	    tjmb[icent][ivz]->SetBranchAddress("gpChg", gpChgImb, &b_gpChgImb);
	    }      
	*/
	nMB[icent][ivz] = tjmb[icent][ivz]->GetEntries();
	cout << "number of evetns in (icent = " << icent << ", ivtxZ = "<< ivz << ")  = " << nMB[icent][ivz] << endl;
	int primeSeed = rand.Integer(37324);
	mbItr[icent][ivz] = primeSeed%(nMB[icent][ivz]);
	cout <<" initial itr = " << mbItr[icent][ivz] << endl;
      }
    }
  }
  else 
    cout << endl << endl << "  Mixing process is skipped" << endl << endl << endl ;

  
  // reweighting factor should go here
  int eTot(0), eSel(0);
  EvtSel evt;
  GammaJet gj;
  Isolation isol;
  TTree *tgj;
  tgj = new TTree("tgj","gamma jet tree");
  tgj->SetMaxTreeSize(MAXTREESIZE);
  tgj->Branch("evt",&evt.run,"run/I:evt:cBin:pBin:trig/O:offlSel:noiseFilt:anaEvtSel:vz/F:reweight/F:hf4Pos:hf4Neg:hf4Sum");
  tgj->Branch("lpho",&gj.photonEt,"photonEt/F:photonRawEt:photonEta:photonPhi:hovere:r9:sigmaIetaIeta:sumIso:genIso:genPhotonEt:genMomId/I:lJetPt/F:lJetEta:lJetPhi:lJetDphi:lJetSubid/I");
  tgj->Branch("isolation",&isol.cc1,"cc1:cc2:cc3:cc4:cc5:cr1:cr2:cr3:cr4:cr5:ct1PtCut20:ct2PtCut20:ct3PtCut20:ct4PtCut20:ct5PtCut20:ecalIso:hcalIso:trackIso");  // ecalIso,hcalIso,trackIso are the pp style isolation                                       
  

  float vzCut = vtxCutPhotonAna;
  TH1F* hvz = new TH1F("hvz","",nVtxBin,-vzCut,vzCut);
  // event plane hitogram                                                                                              
  TH1F* hEvtPlnBin = new TH1F("hEvtPlnBin", "", nPlnBin, -PI/2., PI/2.);
  // jet algos                                                                                                         
  Jets* theJet;
  if (   (colli==kPPDATA) || (colli==kPPMC) ) {
    theJet = &(c->ak5PF) ; 
    cout << "pp collision.  Using ak5PF Jet Algo" << endl<<endl;
  }
  else { 
    theJet = &(c->akPu5PF) ;   
    cout << "pa, aa collision. Using akPu5PF Jet Algo" << endl<<endl;
  }
  
  // Loop starts.
  int nentries = c->GetEntries();
  cout << "number of entries = " << nentries << endl;
  for (Long64_t jentry = 0 ; jentry < nentries; jentry++) {
    eTot++;
    if (jentry% 2000 == 0)  {
      cout <<jentry<<" / "<<nentries<<" "<<setprecision(2)<<(double)jentry/nentries*100<<endl;
    }
    c->GetEntry(jentry);
    evt.clear();
    evt.run   = c->evt.run;
    evt.evt = c->evt.evt;
    evt.hf4Pos = c->evt.hiHFplusEta4;
    evt.hf4Neg = c->evt.hiHFminusEta4;
    evt.hf4Sum = evt.hf4Pos + evt.hf4Neg;
    evt.cBin = -99; 
    evt.pBin   = -99 ; 
    if ((colli==kHIDATA)||(colli==kHIMC))   {
      evt.cBin = c->evt.hiBin;
      evt.pBin   = hEvtPlnBin->FindBin( c->evt.hiEvtPlanes[theEvtPlNumber] ) ;
    }
    else if ((colli==kPADATA)||(colli==kPAMC))   {
      evt.cBin =  getHfBin(evt.hf4Sum);
      if (  ((evt.cBin) < 0) || (evt.cBin) > 18 )
        cout << " Check the pA centrality..  cbin = " << evt.cBin << endl;
    }
    
    evt.reweight = 1;
    evt.trig = 0;
    evt.offlSel = (c->skim.pcollisionEventSelection > 0);
    evt.noiseFilt = (c->skim.pHBHENoiseFilter > 0);
    evt.anaEvtSel = c->selectEvent() && evt.trig;
    evt.vz = c->evt.vz;
    
    
    if ( ( (colli==kHIDATA)||(colli==kHIMC)||(colli==kPADATA)||(colli==kPAMC) ) && ( c->selectEvent() == 0 ))
      continue;
    if ( ( (colli==kPADATA)||(colli==kPPDATA) ) && ( c->skim.pVertexFilterCutGplus ==0 ) ) // No Pile up events
      continue;
    
    eSel++;      // OK.  This event is a collisional and no-noise event.  

    int cBin = evt.cBin;
    int vzBin = hvz->FindBin(evt.vz)  ;
    hvz->Fill(evt.vz) ;
    if ( (vzBin<1) || ( vzBin > nVtxBin) )
      continue;
    
    /// correct the photon energy and make order 

    for (int j=0;j< c->photon.nPhotons;j++) {
	
    if (  ( c->photon.pt[j] > preCutPhotonEt ) && ( fabs( c->photon.eta[j] ) < cutphotonEta ) ) { 
      newPt[j] = c->getCorrEt(j);
    }
    else 
      newPt[j] = c->photon.pt[j] - 10000;
    
    if ( (c->isSpike(j)) || (c->photon.hadronicOverEm[j]>0.2) ||  (c->photon.isEle[j]) )
      newPt[j] = newPt[j] - 20000;
    if (c->photon.seedTime[j] ==0 )   // clustering bug
      newPt[j] = newPt[j] - 30000;
    
    corrPt[j] = newPt[j];
    }
    
    TMath::Sort(c->photon.nPhotons, newPt, order);
    
     
    // Select the leading photon
    gj.clear();
    int leadingIndex=-1;
    for (int j=0;j<c->photon.nPhotons;j++) {
      if ( c->photon.pt[j]  < preCutPhotonEt ) continue; 
      if ( fabs(c->photon.eta[j]) > cutphotonEta ) continue;
      if (c->isSpike(j)) continue;       
      if (!(c->isLoosePhoton(j))) continue;
	
      // sort using corrected photon pt        
      float theCorrPt= corrPt[j];
      if ( theCorrPt > gj.photonEt) {
	gj.photonEt = theCorrPt;
	leadingIndex = j;
      }
    }
    if ( (gj.photonEt < cutphotonPt) ) 
      continue;
    
    
    /// Save leading photons 
    if (leadingIndex!=-1) {
      gj.photonRawEt=c->photon.pt[leadingIndex];
      gj.photonEta=c->photon.eta[leadingIndex];
      gj.photonPhi=c->photon.phi[leadingIndex];
      gj.hovere=c->photon.hadronicOverEm[leadingIndex];
      gj.r9=c->photon.r9[leadingIndex];
      gj.sigmaIetaIeta=c->photon.sigmaIetaIeta[leadingIndex];
      gj.sumIsol = (c->photon.cr4[leadingIndex]+c->photon.cc4[leadingIndex]+c->photon.ct4PtCut20[leadingIndex]) / 0.9;
      gj.genIso   = c->photon.genCalIsoDR04[leadingIndex];
      gj.genPhotonEt = c->photon.genMatchedPt[leadingIndex];
      gj.genMomId = c->photon.genMomId[leadingIndex];
      isol.Set(c,leadingIndex);
    }
    else {
      gj.clear();
    }
    if ( (colli==kPADATA) && ( evt.run > 211256 ) )  {
      gj.photonEta = - gj.photonEta;
      //      cout << " reflect eta" << endl;                                                                         
    }
    


    ///////////////////// Skimmed Jet tree ///////////////////////////////////
    nJet = 0 ;

    int jetEntries = 0;
    if (useGenJetColl )    jetEntries = theJet->ngen;
    else                   jetEntries = theJet->nref;
    
    for (int ij=0; ij< jetEntries ; ij++) {  
      if (  useGenJetColl )   {
        jetPt[nJet] = theJet->genpt[ij];
        jetEta[nJet] = theJet->geneta[ij];
        jetPhi[nJet] = theJet->genphi[ij];
      }
      else  {
        jetPt[nJet] = theJet->jtpt[ij];
        jetEta[nJet] = theJet->jteta[ij];
        jetPhi[nJet] = theJet->jtphi[ij];
      }
      
      // smear the jet pT 
      float smeared = jetPt[nJet] * rand.Gaus(1,addJetEnergyRes/jetPt[nJet])   *  rand.Gaus(1, addFlatJetEnergyRes) ;
      // then multiply jet energy sclae
      smeared = smeared * jetEnergyScale;
      
      // resCorrection 
      float resCorrection =1. ;
      if (  (doJetResCorrection) && ((colli==kHIDATA)||(colli==kHIMC)) ) { // do the residual correction
	if ( evt.cBin  < 12 ) {  // central
	  if ( jetPt[nJet] < 100 ) 
	    resCorrection  =   1.066 - jetPt[nJet]*0.001117 ; 
	  else 
	    resCorrection  =   0.975 ;
	}
	else                  {  // peripheral 
	  if ( jetPt[nJet] < 100 )
            resCorrection  =  1.03440 - jetPt[nJet]*0.000601 ;
          else
            resCorrection  =   0.98 ;
	}
      }
      
      jetPt[nJet] = smeared/resCorrection;
     
      if ( jetPt[nJet] < cutjetPtSkim)
	 continue;
      if ( fabs( jetEta[nJet] ) > cutjetEtaSkim )    
	continue;
      if ( getDR( jetEta[nJet], jetPhi[nJet], gj.photonEta, gj.photonPhi) < 0.3 )
	continue;
      
      if ( (colli==kPADATA) && ( evt.run > 211256 ) )  {
	jetEta[nJet] = -jetEta[nJet];
	//	cout << " reflect eta" << endl;
      }
      
      if (jetPt[nJet] >0) 
	jetDphi[nJet] = getAbsDphi( jetPhi[nJet], gj.photonPhi) ;
      else
	jetDphi[nJet] = -1;
      
      
      if (  useGenJetColl )   {
	jetSubid[nJet] = -9999;
        jetRefPt[nJet] = -9999;
        jetRefEta[nJet] = -9999;
        jetRefPhi[nJet] = -9999;
        jetRefPt[nJet] =  -9999;
	jetRefPartonPt[nJet] = -9999;
        jetRefPartonFlv[nJet] = -9999;
      }
      else {
	jetSubid[nJet] = theJet->subid[ij];
	jetRefPt[nJet] = theJet->refpt[ij];
	jetRefEta[nJet] = theJet->refeta[ij];
	jetRefPhi[nJet] = theJet->refphi[ij];
	if (jetRefPt[nJet] >0)
	  jetRefDphi[nJet] = getAbsDphi( jetRefPhi[nJet] , gj.photonPhi) ;
	else
	  jetRefDphi[nJet] = -1;
	jetRefPartonPt[nJet] = theJet->refparton_pt[ij];
	jetRefPartonFlv[nJet] = theJet->refparton_flavor[ij];
      }

      nJet++ ; 
    }

    //////// Leading jet kinematics in dphi>7pi/8
    float maxJpt = 0;
    int jetLeadingIndex = -1;
    
    for (int ij=0; ij< nJet ; ij++) {
      if ( jetDphi[ij] < awayRange )  // const float awayRange= PI * 7./8.;
	continue; 
      if ( fabs( jetEta[ij] ) > cutjetEta )  // double cutjetEta = 1.6;
	 continue;
      if ( jetPt[ij] > maxJpt) {
	maxJpt = jetPt[ij] ; 
	jetLeadingIndex = ij;
      }
    }
    if ( jetLeadingIndex > -1 ) {
      gj.lJetPt = jetPt[jetLeadingIndex];
      gj.lJetEta = jetEta[jetLeadingIndex];
      gj.lJetPhi = jetPhi[jetLeadingIndex];
      gj.lJetDphi =  jetDphi[jetLeadingIndex];
      gj.lJetSubid=  jetSubid[jetLeadingIndex];
    }
    else {
      gj.lJetPt = -1;
      gj.lJetEta = 999;
      gj.lJetPhi = 999;
      gj.lJetDphi = 0;
      gj.lJetSubid=  -99;
    }
    
    
    //////////////////////////////////////////// New collection of tracks                                                 
    // Forget about the tracks at the moment
    /*    nTrk = 0;
	  for (int it=0; it < c->track.nTrk; it++ ) {
	  if ( c->track.highPurity[it] == 0 )
	  continue;   // only high purity high pt tracks and all pixel track          
	  if ( c->track.trkPt[it] < cuttrkPt )
	  continue;
	  if ( fabs( c->track.trkEta[it]  ) > cuttrkEtaSkim )
	  continue;
	  trkPt[nTrk]  = c->track.trkPt[it];
	  trkEta[nTrk] = c->track.trkEta[it];
	  trkPhi[nTrk] = c->track.trkPhi[it];
	  if  (trkPt[nTrk] > 0 )
	  trkDphi[nTrk] = getAbsDphi ( trkPhi[nTrk], gj.photonPhi) ;
	  else
	  trkDphi[nTrk] = -1;
	  trkWeight[nTrk] = c->getTrackCorrection(it);
	  nTrk++;
	  }
    */
    
    /// Gen Particle 
    // Forget about the gen particle information 
    /* nGp = 0;
       if ( isMC ) {
       for ( int ig = 0 ; ig < c->genparticle.mult ; ig++) {
       if ( c->genparticle.pt[ig] < cuttrkPt )  /// not possibly though..                                 
       continue;
       if ( fabs( c->genparticle.eta[ig] ) > cuttrkEtaSkim )
       continue;
       if ( c->genparticle.chg[ig] == 0 )
       continue;
	 
       gpChg[nGp] =  c->genparticle.chg[ig] ;
	 gpPt[nGp]  = c->genparticle.pt[ig];
	 gpEta[nGp] = c->genparticle.eta[ig];
	 gpPhi[nGp] = c->genparticle.phi[ig];
	 gpSube[nGp] = c->genparticle.sube[ig];
	 if ( gpPt[nGp]>0 )
	 gpDphi[nGp] = getAbsDphi (gpPhi[nGp], gj.photonPhi) ; 
	 else
	 gpDphi[nGp] = -1 ;
	 nGp++;
	 }
	 }
    */
    
     // now let's mix minbias tracks
     
     //    nMtrk = 0;
     // nMgp  = 0;
    
    int nMixing = nMixing1;
    nMjet = 0;
    bool noSuchEvent = false;
    int iMix=0;
    int loopCounter=0;
    
    if ( !doMix ) 
      iMix = nMixing+1;   // Mixing step will be skipped 
    
    while (iMix<nMixing)  {
      loopCounter++;
      if ( loopCounter > nMB[cBin][vzBin]+1) {
	iMix = 999999 ;
	noSuchEvent = true;
	cout << " no such event!! :  icent = " << cBin << ",  vzBin = " << vzBin << ",  pBin = " << evt.pBin << endl;
	continue;
      }
      
      mbItr[cBin][vzBin] = mbItr[cBin][vzBin] + 1;
      if ( mbItr[cBin][vzBin] == nMB[cBin][vzBin] )
	mbItr[cBin][vzBin] =  mbItr[cBin][vzBin] - nMB[cBin][vzBin];
      
      /// Load the minbias tracks!!
      tjmb[cBin][vzBin]->GetEntry(mbItr[cBin][vzBin]);
      
      // Event plane is out of control for both pA and PbPb Sept 2013
      //    if  ( evt.pBin != evtImb.pBin ) 
      //   continue;
      
      // ok found the event!! ///////////
      loopCounter =0;  // Re-initiate loopCounter
      /*  forget about the tracks!   
	  for (int it = 0 ; it < nTrkImb ; it++) {
	  if ( trkPtImb[it] < cuttrkPt )
	  continue;
	  if ( fabs( trkEtaImb[it] ) > cuttrkEtaSkim )
	  continue;
	  
	  mTrkPt[nMtrk]    = trkPtImb[it];
	  mTrkEta[nMtrk]   = trkEtaImb[it];
	  mTrkPhi[nMtrk]   = trkPhiImb[it];
	  mTrkWeight[nMtrk]= trkWeightImb[it];
	  if  ( mTrkPt[nMtrk]>0 )  
	  mTrkDphi[nMtrk]  =   getAbsDphi(mTrkPhi[nMtrk], gj.photonPhi) ;
	  else 
	  mTrkDphi[nMtrk]=-1;
	  
	  nMtrk++;
	  }
      */
      // Jet mixing 
      for (int it = 0 ; it < nJetImb ; it++) {
	
	// smear the jet pT
	float smeared = jetPtImb[it] * rand.Gaus(1,addJetEnergyRes/jetPtImb[it]) *  rand.Gaus(1, addFlatJetEnergyRes) ; 
	smeared = jetEnergyScale * smeared ;  // multiply jet energy scale 
	float resCorrection =1. ;
	if (  (doJetResCorrection) && ((colli==kHIDATA)||(colli==kHIMC)) ) { // do the residual correction   
	  if ( evt.cBin  < 12 ) {  // central                                                                            
	    if ( jetPtImb[it] < 100 )
	      resCorrection  =   1.066 - jetPtImb[it]*0.001117 ;
	    else
	      resCorrection  =   0.975 ;
	  }
	  else                  { // peripheral       
	    if ( jetPtImb[it] < 100 )
	      resCorrection  =  1.03440 - jetPtImb[it]*0.000601 ;
	    else
	      resCorrection  =   0.98 ;
	  }
	}
	float smearedCorrected  = smeared / resCorrection; // residual correction

	if ( smearedCorrected < cutjetPtSkim ) 
	  continue;
	if ( fabs( jetEtaImb[it] ) > cutjetEtaSkim ) 
	  continue;
	
	mJetPt[nMjet]    = smearedCorrected;
	mJetEta[nMjet]   = jetEtaImb[it];
	mJetPhi[nMjet]   = jetPhiImb[it];
	if  ( mJetPt[nMjet]>0 )
	  mJetDphi[nMjet]  =   getAbsDphi(mJetPhi[nMjet], gj.photonPhi) ;
	else
	  mJetDphi[nMjet]=-1;
	
	nMjet++; // < == Important! 
      }
      
      /*      

	 if ( isMC) { 
	 for (int it = 0 ; it < nGpImb ; it++) {
	 if ( gpPtImb[it] < cuttrkPt )
	 continue;
	 if ( fabs(gpEtaImb[it]) > cuttrkEtaSkim )
	 continue;
	 if ( gpChgImb[it] == 0 ) 
	 continue;
	 mGpPt[nMgp]  =   gpPtImb[it];
	 mGpEta[nMgp]  =   gpEtaImb[it];
	 mGpPhi[nMgp]  =   gpPhiImb[it];
	 mGpChg[nMgp]  =   gpChgImb[it];
	 if  (mGpPt[nMgp] > 0 ) 
	 mGpDphi[nMgp] = getAbsDphi ( mGpPhi[nMgp], gj.photonPhi );
	 else 
	 mGpDphi[nMgp]=-1;
	 nMgp++;
	 }
	 }
      */
      
      iMix++;
    }
    if ( noSuchEvent )
      continue;
    
    tgj->Fill();
    newtreeJet->Fill();
    tmixJet->Fill();     
    newtreePhoton->Fill();
    //    treeFullJet->Fill();
    /*
      newtreeTrack->Fill();
      tmixTrk->Fill();   // No tracks saved at the moment 
      if ( isMC ) {
      newtreeGp->Fill();
      tmixGp->Fill();
      } 
    */ 
    
    
  }
  
  
  
  newfile_data->Write();
  //   newfile_data->Close();   // <<=== If there is close() function. writing stucks in the middle of looping.. I don't know why!!
  cout << " Done! "<< endl;
  cout << "    " << eSel<<" out of total "<<eTot<<" events were analyzed."<<endl;
}

