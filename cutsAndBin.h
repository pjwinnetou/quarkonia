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
int kPPMCUps1S = 8 ;
int kPPMCUps2S = 9 ;
int kPPMCUps3S = 10 ;
int kAAMCUps1S = 11 ;
int kAAMCUps2S = 12 ;
int kAAMCUps3S = 13 ;
int kPPAADATASIMUL = 20 ; // 2 and 0 simultaneous fit
int kPPAADATAPeriSIMUL = 60 ; // 6 and 0 simultaneous fit

TString getCollID( int collid ) {
  if ( collid == kPPDATA ) return "PP_DATA";
  else if ( collid == kPADATA ) return "PA_DATA";
  else if ( collid == kAADATA ) return "AA_DATA";
  else if ( collid == kPPMC ) return "PP_MC";
  else if ( collid == kPAMC ) return "PA_MC";
  else if ( collid == kAAMC ) return "AA_MC";
  else if ( collid == kAADATAPeri ) return "AA_DATA_PeriL1";
  else if ( collid == kAADATACentL3 ) return "AA_DATA_CentL3";
  else if ( collid == kPPMCUps1S ) return "PP_MC_Ups1S";
  else if ( collid == kPPMCUps2S ) return "PP_MC_Ups2S";
  else if ( collid == kPPMCUps3S ) return "PP_MC_Ups3S";
  else if ( collid == kAAMCUps1S ) return "AA_MC_Ups1S";
  else if ( collid == kAAMCUps2S ) return "AA_MC_Ups2S";
  else if ( collid == kAAMCUps3S ) return "AA_MC_Ups3S";
  else if ( collid == kPPAADATASIMUL ) return "PP_AA_DATA_SIMUL";
  else if ( collid == kPPAADATAPeriSIMUL ) return "PP_AA_DATA_PeriL1_SIMUL";

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
    sigma1s(0), sigma2s(0), sigma3s(0),
    N1s_1(0), N2s_1(0), N3s_1(0),
    N1s_2(0), N2s_2(0), N3s_2(0),
    Alpha1s_1(0), Alpha2s_1(0), Alpha3s_1(0),
    Alpha1s_2(0), Alpha2s_2(0), Alpha3s_2(0),
    Sigma1s_1(0), Sigma2s_1(0), Sigma3s_1(0),
    Sigma1s_2(0), Sigma2s_2(0), Sigma3s_2(0),
    MCN(0), MCAlpha(0), MCSigma1S(0), MCSigma2S(0), MCSigma3S(0), MCM0(0), MCf(0), MCX(0),
    bkg_mu(0), bkg_sigma(0), bkg_lambda(0),
    bkg_mu_res(0), bkg_sigma_res(0), bkg_lambda_res(0), bkg_mass_res(0),
    // Only For Systematics
    bkg_mu1(0), bkg_sigma1(0), bkg_lambda1(0), bkg_mu2(0), bkg_sigma2(0), bkg_lambda2(0), rBkg2over1(0), // double ErrFunction
    ch3_k1(0), ch3_k2(0), ch3_k3(0),
    ch4_k1(0), ch4_k2(0), ch4_k3(0), ch4_k4(0),
    nSignal1s(0), nSignal2s(0), nSignal3s(0), nBkg(0),
    bkg4_mu(0), bkg4_sigma(0), bkg4_lambda(0),  bkg4_lambda2(0), rBkg42over1(0) //bkg4 = exp + err*exp
    
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
  
  float N1s_1;
  float N1s_2;
  float Alpha1s_1;
  float Alpha1s_2;
  float Sigma1s_1;
  float Sigma1s_2;
  float N2s_1;
  float N2s_2;
  float Alpha2s_1;
  float Alpha2s_2;
  float Sigma2s_1;
  float Sigma2s_2;
  float N3s_1;
  float N3s_2;
  float Alpha3s_1;
  float Alpha3s_2;
  float Sigma3s_1;
  float Sigma3s_2;

  float MCN, MCAlpha, MCSigma1S, MCSigma2S, MCSigma3S, MCM0, MCf, MCX;
  float bkg_mu, bkg_sigma, bkg_lambda;
  
  float bkg_mu_res, bkg_sigma_res, bkg_lambda_res, bkg_mass_res;

  float bkg_mu1, bkg_sigma1, bkg_lambda1, bkg_mu2, bkg_sigma2, bkg_lambda2, rBkg2over1; // double ErrFunction
  float ch3_k1, ch3_k2, ch3_k3 ; 
  float ch4_k1, ch4_k2, ch4_k3, ch4_k4 ; 

  float nSignal1s, nSignal2s, nSignal3s, nBkg;

  float bkg4_mu, bkg4_sigma, bkg4_lambda, bkg4_lambda2, rBkg42over1;

  void setParBkg(float bkg_mu_, float bkg_sigma_, float bkg_lambda_)
  {
    bkg_mu = bkg_mu_; bkg_sigma = bkg_sigma_; bkg_lambda = bkg_lambda_;
  }
  
  void setParBkgRes(float bkg_mu_res_, float bkg_sigma_res_, float bkg_lambda_res_, float bkg_mass_res_)
  {
    bkg_mu_res = bkg_mu_res_; bkg_sigma_res = bkg_sigma_res_; bkg_lambda_res = bkg_lambda_res_; bkg_mass_res = bkg_mass_res_;
  }
  
  void setParMC(float MCN_, float MCAlpha_, float MCSigma1S_, float MCM0_, float MCf_, float MCX_)
  {
    MCN = MCN_ ;  MCAlpha = MCAlpha_ ;  MCSigma1S = MCSigma1S_ ;  MCSigma2S = MCSigma1S*(pdgMass.Y2S/pdgMass.Y1S) ;  MCSigma3S = MCSigma1S*(pdgMass.Y3S/pdgMass.Y1S) ; 
    MCM0 = MCM0_ ;   MCf = MCf_ ;  MCX = MCX_ ;
  }

  void setNAlphaSigmaCB2(float n1_1, float n1_2, float n2_1, float n2_2, float n3_1, float n3_2, float a1_1, float a1_2, float a2_1, float a2_2, float a3_1, float a3_2, float s1_1, float s1_2, float s2_1, float s2_2, float s3_1, float s3_2)
  {
    N1s_1 = n1_1 ;  N1s_2 = n1_2 ;  N2s_1 = n2_1 ;  N2s_2 = n2_2 ;  N3s_1 = n3_1 ;  N3s_2 = n3_2 ;  
    Alpha1s_1 = a1_1 ;  Alpha1s_2 = a1_2 ;  Alpha2s_1 = a2_1 ;  Alpha2s_2 = a2_2 ;  Alpha3s_1 = a3_1 ;  Alpha3s_2 = a3_2 ;
    Sigma1s_1 = s1_1 ;  Sigma1s_2 = s1_2 ;  Sigma2s_1 = s2_1 ;  Sigma2s_2 = s2_2 ;  Sigma3s_1 = s3_1 ;  Sigma3s_2 = s3_2 ;
  }

  void setNAlphaSigma(float n1_, float n2_, float n3_, float a1_, float a2_, float a3_, float s1_, float s2_, float s3_) {
    n1s = n1_ ;    n2s = n2_ ;    n3s = n3_ ;
    alpha1s = a1_ ;    alpha2s = a2_ ;    alpha3s = a3_ ;
    sigma1s = s1_ ;    sigma2s = s2_ ;    sigma3s = s3_ ;
  }
  
  void setParBkg2ErrExp(float bkg_mu1_, float bkg_sigma1_, float bkg_lambda1_, float bkg_mu2_, float bkg_sigma2_, float bkg_lambda2_, float rBkg2over1_)
  {
    bkg_mu1 = bkg_mu1_;  bkg_sigma1 = bkg_sigma1_;  bkg_lambda1 = bkg_lambda1_;
    bkg_mu2 = bkg_mu2_;  bkg_sigma2 = bkg_sigma2_;  bkg_lambda2 = bkg_lambda2_; rBkg2over1 = rBkg2over1_;
  }

  void setParBkgErrExpExp(float bkg4_mu_, float bkg4_sigma_, float bkg4_lambda_, float bkg4_lambda2_, float rBkg42over1_)
  {
    bkg4_mu = bkg4_mu_;  bkg4_sigma = bkg4_sigma_;  bkg4_lambda = bkg4_lambda_;
    bkg4_lambda2 = bkg4_lambda2_; rBkg42over1 = rBkg42over1_;
  }
  
  void setParBkgPol3(float k1_, float k2_, float k3_) 
  {
    ch3_k1 = k1_ ;      ch3_k2 = k2_ ;       ch3_k3 = k3_ ;
  }

  void setParBkgPol4(float k1_, float k2_, float k3_, float k4_) 
  {
    ch4_k1 = k1_ ;      ch4_k2 = k2_ ;       ch4_k3 = k3_ ;   ch4_k4 = k4_;
  }

  void setSig1sF21NBkg(float sig1s_, float f21_, float nBkg_)
  {
    nSignal1s = sig1s_;
    nSignal2s = sig1s_ * f21_;
    nBkg = nBkg_;
  }
  
  void setSig1sF321NBkg(float sig1s_, float f21_, float f31_, float nBkg_)
  {
    nSignal1s = sig1s_;
    nSignal2s = sig1s_ * f21_;
    nSignal3s = sig1s_ * f31_;
    nBkg = nBkg_;
  }
  
  
  
  
  void reset() {
    n1s = 0 ;    n2s = 0 ;    n3s = 0 ;
    alpha1s = 0 ;    alpha2s = 0 ;    alpha3s = 0 ;
    sigma1s = 0 ;    sigma2s = 0 ;    sigma3s = 0 ;
    N1s_1 = 0 ;  N1s_2 = 0 ;  N2s_1 = 0 ;  N2s_2 = 0 ;  N3s_1 = 0 ;  N3s_2 =0 ; 
    Alpha1s_1 = 0 ;  Alpha1s_2 = 0 ;  Alpha2s_1 = 0 ;  Alpha2s_2 = 0 ;  Alpha3s_1 = 0 ;  Alpha3s_2 = 0 ;
    Sigma1s_1 = 0 ;  Sigma1s_2 = 0 ;  Sigma2s_1 = 0 ;  Sigma2s_2 = 0 ;  Sigma3s_1 = 0 ;  Sigma3s_2 = 0 ;
    MCN = 0 ;  MCAlpha = 0 ;  MCSigma1S = 0 ;  MCSigma2S = 0 ;  MCSigma3S = 0 ;  MCM0 = 0 ;  MCf = 0 ;  MCX = 0 ;
    bkg_mu = 0; bkg_sigma = 0; bkg_lambda = 0;
    bkg_mu_res = 0; bkg_sigma_res = 0; bkg_lambda_res = 0; bkg_mass_res=0;
    bkg_mu1 = 0; bkg_sigma1 = 0; bkg_lambda1 = 0; bkg_mu2 = 0; bkg_sigma2 = 0; bkg_lambda2=0; rBkg2over1=0; // double ErrFunction
    ch3_k1 = 0; ch3_k2 = 0; ch3_k3=0 ;
    ch4_k1 = 0; ch4_k2 = 0; ch4_k3=0 ; ch4_k4=0;
    nSignal1s=0; nSignal2s=0; nSignal3s=0; nBkg=0;
    bkg4_mu = 0; bkg4_sigma = 0; bkg4_lambda = 0; bkg4_lambda2 = 0; rBkg42over1=0;
    
  } 
  
};


// Upsilon nominal bins
const int nPtBinsUps = 2;   double ptBinUps[nPtBinsUps+1] = {0, 5,     100};
const int  nYBinsUps = 2;   double yBinUps[nYBinsUps+1] =   {0, 1.2,   2.4};
const int nPBinsUps  = 3;   double pBinUps[nPBinsUps+1] =   {0, 0.167, 0.333,  0.5};

int const nPtBin =4 ;
float ptBin[nPtBin+1] =  {3,6.5,12,20,30};
int const nYBin =2 ;
float yBin[nYBin+1] =  {0, 1.6, 2.4};

TString getKineLabel(int collId, float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut, int cLow, int cHigh, float dphiEp2Low, float dphiEp2High) {
  TString kineLabel = Form("%s_pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f",getCollID(collId).Data(), ptLow,ptHigh, yLow, yHigh, (float)muPtCut) ;
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATAPeri ) || ( collId == kAADATACentL3) || (collId == kAAMCUps1S) || ( collId==kAAMCUps2S) || (collId == kAAMCUps3S) || (collId == kPPAADATASIMUL) || (collId == kPPAADATAPeriSIMUL))
    kineLabel = kineLabel+ Form("_centrality%d-%d_dphiEp_%.2fPI_%.2fPI",(int)cLow, (int)cHigh, (float)dphiEp2Low, (float)dphiEp2High ) ;
  return kineLabel;
}

#endif
