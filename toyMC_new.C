#include "rootFitHeaders.h"
#include "commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
#include "PsetCollection.h"
#include "RooRandom.h"

int kChPol3 = 1 ;
int kErrExp = 2 ;
int kChPol4 = 3 ;
int kErrExpExp = 4 ;

using namespace std;
using namespace RooFit;
void toyMC(
	   int collId = kAADATA,
	   float ptLow=0, float ptHigh=5,
	   float yLow=0, float yHigh=2.4,
	   int cLow=0, int cHigh=200,
	   float muPtCut=4.0,
	   int inputOption=kChPol4, //kChPol3,
	   int nGen = 10000,
	   int useCentIntBkgShape = 1,
     int nToys = 1000
	    ) 
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(111);
  gStyle->SetEndErrorSize(0);
  float Val_2S_1S_nom = 0;
  float Val_2S_1S_alt = 0;
  float Dev_2S_1S = 0;
  

  TString fcoll;
  TString finput;
  if(collId == kAADATA) fcoll = "AA";
  else if(collId == kPPDATA) fcoll = "PP";
  if(inputOption == 3) finput = "4th poly";
  else if(inputOption == 4) finput = "Nominal+Exp";
  
  TFile *wf = new TFile(Form("%s_fit_pt%.1f-%.1f_rap%.1f-%.1f_cent%d-%d_Gen%d_input%d_useCentBkg%d_nToys%d.root",fcoll.Data(),ptLow,ptHigh,yLow,yHigh,cLow,cHigh,nGen,inputOption,useCentIntBkgShape,nToys),"recreate");
  
  TH1D *h1 = new TH1D("h1",Form("SR Nominal, %d toys, %d events, cent %d-%d;2S/1S nom;Counts",nToys,nGen,cLow,cHigh),100,0,1);
  TH1D *h2 = new TH1D("h2",Form("SR %s, %d toys, %d events, cent %d-%d;2S/1S nom;Counts",finput.Data(),nToys,nGen,cLow,cHigh),100,0,1);
  TH1D *h3 = new TH1D("h3","Deviation;2S/1S dev;Counts",1000,0,100);

  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  
  for(int i=0;i<nToys;i++){
  
  float massLow = 8. ;
  float massHigh = 14.;
  int   nMassBin  = (massHigh-massLow)*10;
  
  RooWorkspace *ws = new RooWorkspace("ws");
  RooWorkspace *wsinp = new RooWorkspace("wsinp");

  RooRealVar mass("mass","mass", massLow, massHigh);
  ws->import(mass);
  wsinp->import(mass);
  mass.Print();
  

  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  
  PSet3SingleCB InitialSetUpsilons = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut) ; 
  RooRealVar sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.14);
  RooRealVar sigma2s_1("sigma2s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.14);
  RooRealVar sigma1s_2("sigma1s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.14);
  RooRealVar sigma2s_2("sigma2s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.14);

  RooRealVar alpha1s_1("alpha1s_1","tail shift", 5. , 1.0, 9.8);
  RooRealVar alpha2s_1("alpha2s_1","tail shift", 5. , 1.15, 9.2);
  RooRealVar alpha1s_2("alpha1s_2","tail shift", 5. , 1.0, 9.2);
  RooRealVar alpha2s_2("alpha2s_2","tail shift", 2.5, 1.10, 10.);

  RooRealVar n1s_1("n1s_1","power order", 5. , 1.4, 10.);
  RooRealVar n2s_1("n2s_1","power order", 6. , 1.1, 9.5);
  RooRealVar n1s_2("n1s_2","power order", 5. , 1.4, 10.);
  RooRealVar n2s_2("n2s_2","power order", 6. , 1.1, 9.5);
    
  RooRealVar *f1S = new RooRealVar("f1S","1S CB fraction", InitialSetUpsilons.MCf, InitialSetUpsilons.MCf*0.9, InitialSetUpsilons.MCf*1.1);
  f1S->setVal(InitialSetUpsilons.MCf);  f1S->setConstant();
  RooRealVar X1S("X1S","sigma fraction 1S 2nd CB", InitialSetUpsilons.MCX, InitialSetUpsilons.MCX*0.9, InitialSetUpsilons.MCX*1.1);
  
  // Fix the parameters 
  n1s_1.setVal(InitialSetUpsilons.MCN);  n1s_1.setConstant();  
  n1s_2.setVal(InitialSetUpsilons.MCN);  n1s_2.setConstant();  
  n2s_1.setVal(InitialSetUpsilons.MCN);  n2s_1.setConstant();  
  n2s_2.setVal(InitialSetUpsilons.MCN);  n2s_2.setConstant(); 
  alpha1s_1.setVal(InitialSetUpsilons.MCAlpha);  alpha1s_1.setConstant(); 
  alpha1s_2.setVal(InitialSetUpsilons.MCAlpha);  alpha1s_2.setConstant();
  alpha2s_1.setVal(InitialSetUpsilons.MCAlpha);  alpha2s_1.setConstant();
  alpha2s_2.setVal(InitialSetUpsilons.MCAlpha);  alpha2s_2.setConstant();  
  sigma1s_1.setVal(InitialSetUpsilons.MCSigma1S);  sigma1s_1.setConstant();  
  sigma1s_2.setVal(InitialSetUpsilons.MCSigma1S);  sigma1s_2.setConstant();  
  sigma2s_1.setVal(InitialSetUpsilons.MCSigma1S * InitialSetUpsilons.MCX );  sigma2s_1.setConstant();  
  sigma2s_2.setVal(InitialSetUpsilons.MCSigma1S * InitialSetUpsilons.MCX );  sigma2s_2.setConstant();  
  mean1s.setVal(InitialSetUpsilons.bkg_mass_res); mean1s.setConstant();
  
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );
  
  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", mass, mean1s, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape* cb2s_1 = new RooCBShape("cball2s_1", "cystal Ball", mass, mean2s, sigma2s_1, alpha2s_1, n2s_1);
  
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", mass, mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape* cb2s_2 = new RooCBShape("cball2s_2", "cystal Ball", mass, mean2s, sigma2s_2, alpha2s_2, n2s_2);
 
  RooAddPdf*  cb1s = new RooAddPdf();
  RooAddPdf*  cb2s = new RooAddPdf();
  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1S) );
  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(*cb2s_1,*cb2s_2), RooArgList(*f1S) );


  // Input model 
  PSet3SingleCB bkgParm = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut) ; 
  PSet3SingleCB bkgParmCentInt;
  if ( !( (cLow==0) && (cHigh==200) ) && (collId==kAADATA)   ) {
    bkgParmCentInt      = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh,    0,   200, muPtCut) ; 
    cout << " ok done " << endl;
  }

  // if ( inputOption == kErrExp ) 
  RooRealVar err_mu1("#mu1","err_mu1",  bkgParm.bkg_mu1 ) ;
  RooRealVar err_sigma1("#sigma1","err_sigma1", bkgParm.bkg_sigma1);
  RooRealVar m_decay1("#lambda1","m_decay1", bkgParm.bkg_lambda1);
  RooRealVar err_mu2("#mu2","err_mu2",  bkgParm.bkg_mu2 ) ;
  RooRealVar err_sigma2("#sigma2","err_sigma2", bkgParm.bkg_sigma2);
  RooRealVar m_decay2("#lambda2","m_decay2", bkgParm.bkg_lambda2);
  
  
  float the_ch3_k1 = bkgParm.ch3_k1 ;   float the_ch3_k2 = bkgParm.ch3_k2 ;  float the_ch3_k3 = bkgParm.ch3_k3 ;
  float the_ch4_k1 = bkgParm.ch4_k1 ;   float the_ch4_k2 = bkgParm.ch4_k2 ;  float the_ch4_k3 = bkgParm.ch4_k3 ; float the_ch4_k4 = bkgParm.ch4_k4 ;
  float the_bkg4_mu = bkgParm.bkg4_mu ;           float the_bkg4_sigma = bkgParm.bkg4_sigma; 
  float the_bkg4_lambda = bkgParm.bkg4_lambda ;   float the_bkg4_lambda2 = bkgParm.bkg4_lambda2 ;
  if ( !( (cLow==0) && (cHigh==200) ) && (collId==kAADATA) && useCentIntBkgShape  ) {
    the_ch3_k1 = bkgParmCentInt.ch3_k1 ;    the_ch3_k2 = bkgParmCentInt.ch3_k2 ;   the_ch3_k3 = bkgParmCentInt.ch3_k3 ;
    the_ch4_k1 = bkgParmCentInt.ch4_k1 ;    the_ch4_k2 = bkgParmCentInt.ch4_k2 ;   the_ch4_k3 = bkgParmCentInt.ch4_k3 ; the_ch4_k4 = bkgParmCentInt.ch4_k4 ;
    the_bkg4_mu = bkgParmCentInt.bkg4_mu ;            bkgParmCentInt.bkg4_sigma =bkgParmCentInt.bkg4_sigma;
    the_bkg4_lambda = bkgParmCentInt.bkg4_lambda ;    the_bkg4_lambda2 = bkgParmCentInt.bkg4_lambda2 ;
  }
  // if ( inputOption == kChPol3 ) 
  RooRealVar ch3_k1("pol3_k1","pol3_k1", the_ch3_k1 ) ;
  RooRealVar ch3_k2("pol3_k2","pol3_k2", the_ch3_k2 ) ;
  RooRealVar ch3_k3("pol3_k3","pol3_k3", the_ch3_k3 ) ;
  // if ( inputOption == kChPol4 ) 
  RooRealVar ch4_k1("pol4_k1","pol4_k1", the_ch4_k1 , the_ch4_k1*0.3, the_ch4_k1*1.6) ;
  RooRealVar ch4_k2("pol4_k2","pol4_k2", the_ch4_k2 , the_ch4_k2*0.3, the_ch4_k2*1.6) ;
  RooRealVar ch4_k3("pol4_k3","pol4_k3", the_ch4_k3 , the_ch4_k3*0.3, the_ch4_k3*1.6) ;
  RooRealVar ch4_k4("pol4_k4","pol4_k4", the_ch4_k4 , the_ch4_k4*0.3, the_ch4_k4*1.6) ;


  // if (inputOption == kErrExpExp ) 
  RooRealVar err4_mu("err4_mu","err4_mu",  the_bkg4_mu , the_bkg4_mu*0.4,the_bkg4_mu*1.4) ;
  RooRealVar err4_sigma("err4_sigma","err4_sigma", the_bkg4_sigma, the_bkg4_sigma*0.4, the_bkg4_sigma*1.4);
  RooRealVar m4_decay("err4_lambda","m4_decay", the_bkg4_lambda, the_bkg4_lambda*0.4, the_bkg4_lambda*1.4);
  RooRealVar m4_decay2("err4_lambda2","m4_decay2", the_bkg4_lambda2, the_bkg4_lambda2*0.4, the_bkg4_lambda2*1.4);

  RooGenericPdf *bkgErrExp1;
  RooGenericPdf *bkgErrExp2;

  RooGenericPdf *bkg4ErrExp ; // kErrExpExp
  RooGenericPdf *bkg4Exp = new RooGenericPdf("bkg4Exp","bkg4Exp","TMath::Exp(-@0/@1)",RooArgList(mass,m4_decay2));


 
 if ( ptLow == 0)  { 
   bkg4ErrExp = new RooGenericPdf("bkg4ErrExp","bkg4ErrExp","(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5*TMath::Exp(-@0/@3)",RooArgList(mass,err4_mu,err4_sigma,m4_decay));
   bkgErrExp1 = new RooGenericPdf("bkgErrExp1","Background1","(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5*TMath::Exp(-@0/@3)",RooArgList(mass,err_mu1,err_sigma1,m_decay1));
   bkgErrExp2 = new RooGenericPdf("bkgErrExp2","Background2","(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5*TMath::Exp(-@0/@3)",RooArgList(mass,err_mu2,err_sigma2,m_decay2));
 }
  else  { // if ptLow >= 5 
    bkg4ErrExp = new RooGenericPdf("bkg4ErrExp","bkg4ErrExp", "TMath::Exp(-@0/@1)",RooArgList(mass,m4_decay));
    bkgErrExp1 = new RooGenericPdf("bkgErrExp1","Background1","TMath::Exp(-@0/@1)",RooArgList(mass,m_decay1));
    bkgErrExp2 = new RooGenericPdf("bkgErrExp2","Background2","TMath::Exp(-@0/@1)",RooArgList(mass,m_decay2));
  }




  RooRealVar* rBkg2nd = new RooRealVar("rBkg2over1","rBkg2over1", bkgParm.rBkg42over1); // bkgParm.rBkgErr2over1
  RooAddPdf* bkgDblErr = new RooAddPdf("bkgDblErrExp","Bkg Only",RooArgList(*bkgErrExp2, *bkgErrExp1),RooArgList(*rBkg2nd));   // if ( inputOption == kErrExp )
  RooAddPdf* bkgComp4 = new RooAddPdf("bkgComp4","bkgComp4",RooArgList(*bkg4Exp, *bkg4ErrExp),RooArgList(*rBkg2nd));   // if ( inputOption == kErrExp )

  RooChebychev * bkgChPol3 = new RooChebychev("cPolBkg","Background1",mass,RooArgSet(ch3_k1,ch3_k2,ch3_k3));  // if ( inputOption == kChPol3 )
  RooChebychev * bkgChPol4 = new RooChebychev("cPol4Bkg","Background4",mass,RooArgSet(ch4_k1,ch4_k2,ch4_k3,ch4_k4));  // if ( inputOption == kChPol3 )

  float r1S_overTot = bkgParm.nSignal1s / ( bkgParm.nSignal1s + bkgParm.nSignal2s + bkgParm.nBkg ) ; // Numbers obtained from the real data
  float r2S_overTot = bkgParm.nSignal2s / ( bkgParm.nSignal1s + bkgParm.nSignal2s + bkgParm.nBkg ) ; 
  float rBkg_overTot = bkgParm.nBkg / ( bkgParm.nSignal1s + bkgParm.nSignal2s + bkgParm.nBkg ) ; 
  
  RooRealVar *nSig1sInp  = new RooRealVar("nSig1sInp","nSig1sInp", nGen * r1S_overTot,  0,   nGen);
  RooRealVar *nSig2sInp  = new RooRealVar("nSig2sInp","nSig2sInp", nGen * r2S_overTot, 0,   nGen);
  RooRealVar *nBkgInp  = new RooRealVar("nBkgInp","n_bkgInp",      nGen * rBkg_overTot,  0,   nGen);
  

  //----------------------------------------------------------------------------------------
  //Generating function from nominal fit
  
  RooRealVar err_mu_gen("err_mu_gen","err_mu_gen",  bkgParm.bkg_mu_res) ;
  RooRealVar err_sigma_gen("err_sigma_gen","err_sigma_gen", bkgParm.bkg_sigma_res);
  RooRealVar m_decay_gen("err_lambda_gen","m_decay_gen", bkgParm.bkg_lambda_res);

  err_mu_gen.setVal(bkgParm.bkg_mu_res); err_mu_gen.setConstant();
  err_sigma_gen.setVal(bkgParm.bkg_sigma_res); err_sigma_gen.setConstant();
  m_decay_gen.setVal(bkgParm.bkg_lambda_res); m_decay_gen.setConstant();
  
  RooGenericPdf* bkgInp_gen;
  RooGenericPdf *bkgInp_in;
  if ( ptLow == 0)  { 
    bkgInp_in = new RooGenericPdf("bkgInp_gen","Background Gen","(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5*TMath::Exp(-@0/@3)",RooArgList(mass,err_mu_gen,err_sigma_gen,m_decay_gen));
  }
  else {
    bkgInp_in = new RooGenericPdf("bkgInp_gen","Background Gen","TMath::Exp(-@0/@1)",RooArgList(mass,m_decay_gen));
  }
  
  bkgInp_gen = bkgInp_in;
  
  
  RooAddPdf* modelInput_gen; 
  modelInput_gen = new RooAddPdf("modelInput_gen","1S+2S + Bkg",RooArgList(*cb1s, *cb2s, *bkgInp_gen),RooArgList(*nSig1sInp,*nSig2sInp,*nBkgInp));
  
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  
  RooAddPdf* modelInput; 
  RooGenericPdf* bkgInp;
  if ( inputOption == kErrExp ) { 
    bkgInp = (RooGenericPdf*) bkgDblErr;
  }
  else if ( inputOption == kChPol3 ) { 
    bkgInp = (RooGenericPdf*) bkgChPol3;
  }
  else if ( inputOption == kChPol4 ) { 
    bkgInp = (RooGenericPdf*) bkgChPol4;
  }
  else if ( inputOption == kErrExpExp ) { 
    bkgInp = (RooGenericPdf*) bkgComp4;
  }

  modelInput = new RooAddPdf("modelInput","1S+2S + Bkg",RooArgList(*cb1s, *cb2s, *bkgInp),RooArgList(*nSig1sInp,*nSig2sInp,*nBkgInp));
  wsinp->import(*modelInput);

  Val_2S_1S_nom=0;
  Val_2S_1S_alt=0;
  Dev_2S_1S=0;


  RooDataSet *data = modelInput_gen->generate(mass,nGen) ;

  RooPlot* xframe  = ws->var("mass")->frame(nMassBin); // bins
  xframe->SetXTitle("mass (Gev/c^{2})");
  xframe->GetXaxis()->CenterTitle();
  xframe->GetYaxis()->CenterTitle();
  RooPlot* xframe2 = (RooPlot*)xframe->Clone("xframe2");
  
  RooFitResult* fitResInput = wsinp->pdf("modelInput")->fitTo(*data,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE));
  data->plotOn(xframe,Name("dataHist"),MarkerSize(0.7)) ;
  wsinp->pdf("modelInput")->plotOn(xframe, Name("inputModelHist"));
  wsinp->pdf("modelInput")->plotOn(xframe, Components(RooArgSet(*bkgInp)),LineColor(kBlack),LineStyle(kDashed));
  if ( inputOption == kErrExp ) 
    {  
      modelInput->plotOn(xframe,Components(RooArgSet(*bkgDblErr)),LineColor(kRed),LineStyle(kDashed));
      modelInput->plotOn(xframe,Components(RooArgSet(*bkgErrExp1)),LineColor(kBlack),LineStyle(kDashed));
      modelInput->plotOn(xframe,Components(RooArgSet(*bkgErrExp2)),LineColor(kBlack),LineStyle(kDashed));
    }
  else if  ( inputOption == kChPol3 )
    modelInput->plotOn(xframe,Components(RooArgSet(*bkgChPol3)),LineColor(kBlack),LineStyle(kDashed));
  else if  ( inputOption == kChPol4 )
    modelInput->plotOn(xframe,Components(RooArgSet(*bkgChPol4)),LineColor(kBlack),LineStyle(kDashed));
  else if (inputOption == kErrExpExp ) {
    modelInput->plotOn(xframe,Components(RooArgSet(*bkgComp4)),LineColor(kBlack),LineStyle(kDashed));
    modelInput->plotOn(xframe,Components(RooArgSet(*bkg4ErrExp)),LineColor(kBlack),LineStyle(kDashed));
    modelInput->plotOn(xframe,Components(RooArgSet(*bkg4Exp)),LineColor(kBlack),LineStyle(kDashed));
  }


  // New fit 
  
  float the_bkg_mu = bkgParm.bkg_mu ;
  float the_bkg_sigma = bkgParm.bkg_sigma ;
  float the_bkg_lambda = bkgParm.bkg_lambda ;
  if ( !( (cLow==0) && (cHigh==200) ) && (collId==kAADATA) && useCentIntBkgShape  ) {
    the_bkg_mu = bkgParmCentInt.bkg_mu ;
    the_bkg_sigma = bkgParmCentInt.bkg_sigma ;
    the_bkg_lambda = bkgParmCentInt.bkg_lambda ;
  }

  //RooRealVar err_mu("err_mu","err_mu", the_bkg_mu, 0.0, 40);
  RooRealVar err_mu("err_mu","err_mu", the_bkg_mu, the_bkg_mu*0.4, the_bkg_mu*1.4);
  //RooRealVar err_mu("err_mu","err_mu", 1., 0.0, 30);
  //RooRealVar err_sigma("err_sigma","err_sigma", 1.2, 1.1,55);
  //RooRealVar err_sigma("err_sigma","err_sigma", 10.,0,20);
  RooRealVar err_sigma("err_sigma","err_sigma", the_bkg_sigma, the_bkg_sigma*0.4, the_bkg_sigma*1.4);
  //RooRealVar m_decay("m_decay","m_decay", 10., 6.5, 30);
  RooRealVar m_decay("m_decay","m_decay",the_bkg_lambda, the_bkg_lambda*0.4, the_bkg_lambda*1.4);
  
  if( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4) && collId==kPPDATA) 
  {
    err_sigma.setVal(1.055); 
    err_sigma.setConstant();
  }
  if( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4) && collId==kAADATA) 
  {
    err_sigma.setVal(1.103); 
    err_sigma.setConstant();
  }

  RooGenericPdf *bkgFitOut;
  if ( ptLow == 0)  { 
    bkgFitOut = new RooGenericPdf("bkgFitOut","BackgroundOut","(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5*TMath::Exp(-@0/@3)",RooArgList(mass,err_mu,err_sigma,m_decay));
  }
  else {
    bkgFitOut = new RooGenericPdf("bkgFitOut","BackgroundOut","TMath::Exp(-@0/@1)",RooArgList(mass,m_decay));
  }

  

  RooRealVar *nSig1sOut  = new RooRealVar("nSig1sOut","nSig1sOut", r1S_overTot*nGen, 0,  r1S_overTot*2.*nGen);
  RooRealVar *nSig2sOut  = new RooRealVar("nSig2sOut","nSig2sOut", r2S_overTot*nGen, 0, r2S_overTot*2.*nGen);
  RooRealVar *nBkgOut  = new RooRealVar("nBkgOut","n_bkgOut",nGen * rBkg_overTot, 0, nGen);

  RooAddPdf*  cb1sOut = (RooAddPdf*)cb1s->Clone("cb1sOutput");
  RooAddPdf*  cb2sOut = (RooAddPdf*)cb2s->Clone("cb2sOutput");
  RooAddPdf* modelOutput = new RooAddPdf("modelOutput","1S+2S + Bkg",RooArgList(*cb1sOut, *cb2sOut, *bkgFitOut),RooArgList(*nSig1sOut,*nSig2sOut,*nBkgOut));
  ws->import(*modelOutput);
  
  RooFitResult* fitRes = ws->pdf("modelOutput")->fitTo(*data,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE));
  data->plotOn(xframe2,Name("dataHist2"),MarkerSize(0.7)) ;
  ws->pdf("modelOutput")->plotOn(xframe2, Name("outputModelHist"));
  ws->pdf("modelOutput")->plotOn(xframe2, Components(RooArgSet(*bkgFitOut)),LineColor(kBlack),LineStyle(kDashed));
  
  Val_2S_1S_nom = (float)(ws->var("nSig2sOut")->getVal() / ws->var("nSig1sOut")->getVal());
  Val_2S_1S_alt = (float)(wsinp->var("nSig2sInp")->getVal() / wsinp->var("nSig1sInp")->getVal());
  Dev_2S_1S = (Val_2S_1S_alt/Val_2S_1S_nom - 1) * 100;
  h1->Fill(Val_2S_1S_nom);
  h2->Fill(Val_2S_1S_alt);
  h3->Fill(Dev_2S_1S);

  // DRAW! 
  if(i == 0){
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,800,400);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.49, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  pad1->SetBottomMargin(0); // Upper and lower plot are joined

  xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
  drawText(Form("#Upsilon(2S)/#Upsilon(1S) = %.5f",(float)(wsinp->var("nSig2sInp")->getVal() / wsinp->var("nSig1sInp")->getVal())),0.2,0.54,1,16) ;
  
  if (inputOption==kChPol4 ) 
    drawText("4th order poly. Bkg.",0.2,0.62,2,15) ;
  if (inputOption==kErrExpExp ) 
    drawText("Erf*exp + exp Bkg.",0.2,0.62,2,15) ;
  
  if(collId == kAADATA)
    drawText("PbPb",0.4,0.45,1,15);
  if(collId == kPPDATA)
    drawText("pp", 0.4,0.45,1,15);

  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV",ptLow,ptHigh ),0.5,0.60,1,12);
  drawText(Form("%.1f < y^{#mu#mu} < %.1f",yLow,yHigh ), 0.5,0.55,1,12);
  TString perc = "%";
  if(collId == kAADATA)
    drawText(Form("Cent %d-%d%s",cLow/2,cHigh/2,perc.Data()),0.5,0.5,4,12);
  
  TLatex *tex = new TLatex(0.4,0.88,"Toy MC generated");
  tex->SetTextFont(43);
  tex->SetTextSize(15);
  tex->SetNDC();
  //  tex->SetTextAngle(180);
  tex->Draw();

  RooArgList paramListinp = fitResInput->floatParsFinal();
  paramListinp.Print("v");
  RooPlot* legFrameinp = wsinp->var("mass")->frame(Name("Fit Results"), Title("Fit Results"));
  wsinp->pdf("modelInput")->paramOn(legFrameinp,Layout(.6,.9, .5),Parameters(paramListinp));
  legFrameinp->getAttText()->SetTextAlign(11);
  legFrameinp->getAttText()->SetTextSize(0.028);
  TPaveText* hhinp = (TPaveText*)legFrameinp->findObject(Form("%s_paramBox",wsinp->pdf("modelInput")->GetName()));
  hhinp->SetY1(0.35); hhinp->SetY2(0.83);
  hhinp->Draw();
  // PULL

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.49, 0.25);
  c1->cd();
  pad2->Draw();
  pad2->cd();
  RooHist* hpull = xframe->pullHist("dataHist","inputModelHist");
  RooPlot* pullFrame = wsinp->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(2.57);
  pullFrame->GetYaxis()->SetTitleOffset(1.8) ;
  pullFrame->GetYaxis()->SetLabelSize(0.16) ;
  pullFrame->GetYaxis()->SetRange(-10,10) ;
  pullFrame->GetXaxis()->SetTitleOffset(0.7) ;
  pullFrame->GetXaxis()->SetLabelSize(0.1) ;
  pullFrame->GetXaxis()->SetTitleSize(0.13) ;
  pullFrame->Draw() ;


  TPad *pad3 = new TPad("pad3", "pad3", 0.51, 0.25, 0.99, 1);
  pad3->SetTicks(1,1);
  pad3->SetBottomMargin(0); // Upper and lower plot are joined
  c1->cd();
  pad3->Draw(); pad3->cd();

  xframe2->GetYaxis()->SetTitleOffset(1.4) ; xframe2->Draw() ;
  TLatex *tex2 = new TLatex(0.4,0.9,"Fitted by Nominal function");
  tex2->SetTextFont(43);
  tex2->SetTextSize(15);
  tex2->SetTextColor(2);
  tex2->SetNDC();
  tex2->Draw();
  drawText(Form("#Upsilon(2S)/#Upsilon(1S) = %.5f",(float)(ws->var("nSig2sOut")->getVal() / ws->var("nSig1sOut")->getVal())), 0.4,0.85,1,16 );

  // *~*~*~*~*~*~*~* Draw the parameters in the plot  *~*~*~*~*~*~*~* //
  RooArgList paramList = fitRes->floatParsFinal();
  paramList.Print("v");
  RooPlot* legFrame = ws->var("mass")->frame(Name("Fit Results"), Title("Fit Results"));
  ws->pdf("modelOutput")->paramOn(legFrame,Layout(.6,.9, .5),Parameters(paramList));
  legFrame->getAttText()->SetTextAlign(11);
  legFrame->getAttText()->SetTextSize(0.028);
  TPaveText* hh = (TPaveText*)legFrame->findObject(Form("%s_paramBox",ws->pdf("modelOutput")->GetName()));
  hh->SetY1(0.35); hh->SetY2(0.83);
  hh->Draw();

  TPad *pad4 = new TPad("pad4", "pad4", 0.51, 0.05, 0.99, 0.25);
  // pad4->SetBottomMargin(0); // Upper and lower plot are joined
  c1->cd();
  pad4->Draw();
  pad4->cd();
  RooHist* hpullOut = xframe2->pullHist("dataHist2","outputModelHist");
  RooPlot* pullOutFrm = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullOutFrm->addPlotable(hpullOut,"P") ;
  pullOutFrm->SetTitleSize(2.57);
  pullOutFrm->GetYaxis()->SetTitleOffset(1.8) ;
  pullOutFrm->GetYaxis()->SetLabelSize(0.16) ;
  pullOutFrm->GetYaxis()->SetRange(-10,10) ;
  pullOutFrm->GetXaxis()->SetTitleOffset(0.7) ;
  pullOutFrm->GetXaxis()->SetLabelSize(0.1) ;
  pullOutFrm->GetXaxis()->SetTitleSize(0.13) ;
  pullOutFrm->Draw() ;



  // *~*~*~*~*~*~*~* Print the results *~*~*~*~*~*~*~* //

  //cout << "nSig2sInp/nSig1sInp = " << nSig2sInp->getVal() / nSig1sInp->getVal() << endl;
  cout << "input fit ratio = " << wsinp->var("nSig2sInp")->getVal() / wsinp->var("nSig1sInp")->getVal() << endl;
  cout << "output fit ratio = " << ws->var("nSig2sOut")->getVal() / ws->var("nSig1sOut")->getVal() << endl;
  c1->SaveAs(Form( "toyMCFit_collId%d_pt%.0f-%.0fGeV_y%.0f-%.0f_cBin%d-%d_muPtCut%.0fGeV_BkgPDFOpt%d_nGen%d_useCentIntBkgShape%d.png",
		   collId, ptLow, ptHigh, yLow*10, yHigh*10, cLow, cHigh, muPtCut, inputOption, nGen,useCentIntBkgShape) );
  
  float r1 =  wsinp->var("nSig2sInp")->getVal() / wsinp->var("nSig1sInp")->getVal() ; 
  float r2 =  ws->var("nSig2sOut")->getVal() / ws->var("nSig1sOut")->getVal() ; 
  cout << Form( "collId: %d,    pt: %.0f - %.0fGeV,   y: %.1f - %.1f,  cBin: %d - %d", collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh ) << endl;
  cout << "Uncertainty = "  << (r2 - r1 ) / r1 << endl;
  }  
  
  }

  wf->cd();
  h1->Write();
  h2->Write();
  h3->Write();

   
} 

