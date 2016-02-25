#ifndef CutAndBinCollection_C
#define CutAndBinCollection_C

#include <TF1.h>
#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>



struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlus, KaPlus; };
ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.139570, 0.49367 };


int kPPDATA = 0 ;
int kPADATA = 1 ;
int kAADATA = 2 ; // L1 doubleMu 0
int kPPMC = 3 ;
int kPAMC = 4 ;
int kAAMC = 5 ;
int kAADATAPeri = 6 ;
int kAADATACentL3 = 7 ;

TString getCollID( int collid ) {
  if ( collid == kPPDATA ) return "PP_DATA";
  else if ( collid == kPADATA ) return "PA_DATA";
  else if ( collid == kAADATA ) return "AA_DATA";
  else if ( collid == kPPMC ) return "PP_MC";
  else if ( collid == kPAMC ) return "PA_MC";
  else if ( collid == kAAMC ) return "AA_MC";
  else if ( collid == kAADATAPeri ) return "AA_DATA_PeriL1";
  else if ( collid == kAADATACentL3 ) return "AA_DATA_CentL3";
  else return "none";
}

int kEPl2HF = 0;
int kEPOppositeHF = 1;
int kEPSameSideHF = 2;


TString getEPSel( int eventPln) {
  if ( eventPln == kEPl2HF)  return "BothHFs";
  else if ( eventPln == kEPOppositeHF ) return "OppositeHF" ;
  else if ( eventPln == kEPSameSideHF ) return "SameSideHF" ;
  else return "none";
}


int kSoftMuCut = 0;
int kHighPtMuCut = 0;




class DiMuon {
 public:
 DiMuon() :
  run(0),   lumi(0), event(0), cBin(0), ep2(0), dphiEp2(0),
    vz(-99),  mass(-1), pt(-1), y(999), phi(999), eta(999),
    pt1(-1), eta1(-1), phi1(-1),        
    pt2(-1), eta2(-1), phi2(-1),        
    oniaIndex(-1), softFlag(0), highPtFlag(0)
    {}
  
  int run;
  int lumi;
  int event;
  int cBin;
  float ep2;
  float dphiEp2;
  float vz;
  float mass;
  float pt;
  float y;
  float phi;    
  float eta;
  float pt1; 
  float eta1;
  float phi1;
  float pt2;
  float eta2;
  float phi2;    
  int oniaIndex;
  int softFlag;
  int highPtFlag;

  void clear() {
    run = -99;  lumi=-99; event=-99; cBin=-99; ep2=-99, dphiEp2=-99; 
    vz=-99;     mass = -99; pt=-99; y=-99; phi=-99; eta=-99;      
    pt1=-99; eta1=-99; phi1=-99; pt2=-99; eta2=-99; phi2=-99;
    oniaIndex=-1; softFlag=-1; highPtFlag=-1; 
  }

};
TString branchString = "run/I:lumi:event:cBin:ep2/F:dphiEp2:vz:mass:pt:y:phi:eta:pt1:eta1:phi1:pt2:eta2:phi2:oniaIndex/I:softFlag:highPtFlag";

class  PSet3SingleCB { 
 public:
 PSet3SingleCB() :
  n1s(0), n2s(0), n3s(0),
    alpha1s(0), alpha2s(0), alpha3s(0),
    sigma1s(0), sigma2s(0), sigma3s(0)
    {}
  
  float n1s;
  float n2s;
  float n3s;
  float alpha1s;
  float alpha2s;
  float alpha3s;
  float sigma1s;
  float sigma2s;
  float sigma3s;

  void setNAlphaSigma(float n1_, float n2_, float n3_, float a1_, float a2_, float a3_, float s1_, float s2_, float s3_) {
    n1s = n1_ ;    n2s = n2_ ;    n3s = n3_ ;
    alpha1s = a1_ ;    alpha2s = a2_ ;    alpha3s = a3_ ;
    sigma1s = s1_ ;    sigma2s = s2_ ;    sigma3s = s3_ ;
  } 
  void reset() {
    n1s = 0 ;    n2s = 0 ;    n3s = 0 ;
    alpha1s = 0 ;    alpha2s = 0 ;    alpha3s = 0 ;
    sigma1s = 0 ;    sigma2s = 0 ;    sigma3s = 0 ;
  } 

};



int const nPtBin =4 ;
float ptBin[nPtBin+1] =  {3,6.5,12,20,30};
int const nYBin =2 ;
float yBin[nYBin+1] =  {0, 1.6, 2.4};

TString getKineLabel(int collId, float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut, int cLow, int cHigh, float dphiEp2Low, float dphiEp2High) {
  TString kineLabel = Form("%s_pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f",getCollID(collId).Data(), ptLow,ptHigh, yLow, yHigh, (float)muPtCut) ;
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATAPeri ) || ( collId == kAADATACentL3) )
    kineLabel = kineLabel+ Form("_centrality%d-%d_dphiEp_%.2fPI_%.2fPI",(int)cLow, (int)cHigh, (float)dphiEp2Low, (float)dphiEp2High ) ;
  return kineLabel;
}

#endif
