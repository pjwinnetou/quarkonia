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
using namespace std;
using namespace RooFit;

RooWorkspace* doFitJpsiPsi( float ptLow=6.5, float ptHigh=30, float yLow=1.6, float yHigh=2.4) {
  
  using namespace RooFit;

  int nMassBins = 60;
  float massLow   = 2.6;
  float massHigh = 4.2;
  float pdgJpsiMass = pdgMass.JPsi;
  //  pdgJpsiMass = pdgJpsiMass - 0.01;
  float pdgJpsiWidth = 0.0961;
  float pdgPsi2sMass = pdgMass.Psi2S;
  float pdgPsi2sWidth = 0.286;

  TString kineCut = Form("pt>%f && pt<%f && abs(y)>%f && abs(y)<%f",ptLow, ptHigh, yLow, yHigh);
  //  TString muonCut = "pt1>3 && pt2>3";   kineCut = kineCut + " && " + muonCut ;
  

  //TFile* f1 = new TFile("skimmedFiles/yskimPbPb_20151242343.root");
  TFile* f1 = new TFile("skimmedFiles/yskimPP_20151242342.root");
  TTree* tree = (TTree*) f1->Get("dimu");
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace(Form("workspace_pt%f-%f_y%f-%f",ptLow,ptHigh, yLow, yHigh));


  ws->import(*dataset);
  ws->data("dataset")->Print();
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  RooRealVar* mass = (RooRealVar*)ws->var("mass");
  RooRealVar pt("pt","pt",0,100);
  mass->Print();
  ws->var("mass")->setRange(massLow, massHigh);
  mass->Print();
  
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,480,530);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0);
  pad1->Draw(); pad1->cd();

  RooPlot* myPlot = ws->var("mass")->frame(nMassBins); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));
  RooRealVar sig_mean1s("sig_mean1s","mean of the signal gaussian mass PDF",pdgJpsiMass,pdgJpsiMass-0.05, pdgJpsiMass+0.05);
  RooRealVar sig_sigma1s("sig_sigma1s","width/sigma of the signal gaussian mass PDF",pdgJpsiWidth, 0.0001, 0.1);
  RooRealVar sig_mean2s("sig_mean2s","mean of the signal gaussian mass PDF",pdgPsi2sMass,massLow,massHigh);
  RooRealVar sig_sigma2s("sig_sigma2s","width/sigma of the signal gaussian mass PDF",pdgPsi2sMass, 0.0001, 0.1);
  //  RooGaussian* sigModel0 = new RooGaussian("sigModel0","signal mass PDF (gaussian)",mass,sig_mean,sig_sigma);
  RooRealVar alpha1s("alpha1s","tail shift", 5,1,10);
  RooRealVar n1s("n1s","power order", 5,1,10);
  RooCBShape* cb1s = new RooCBShape("cball1s", "cystal Ball", *mass, sig_mean1s, sig_sigma1s, alpha1s, n1s);
  RooRealVar gSigma1s("gSigma1s","sigma of gaussian", 1.6943e-02, 0, 0.1);
  RooGaussian* gau1s = new RooGaussian("gau1s","Signal component 2",*mass,sig_mean1s,gSigma1s) ;
  RooRealVar alpha2s("alpha2s","tail shift", 5, 1, 10);
  RooRealVar n2s("n2s","power order", 5, 1, 10);
  RooCBShape* cb2s = new RooCBShape("cball2s", "cystal Ball", *mass, sig_mean2s, sig_sigma2s, alpha2s, n2s);
  RooRealVar gSigma2s("gSigma2s","sigma of gaussian", 1.6943e-02, 0, 0.1);
  RooGaussian* gau2s = new RooGaussian("gau2s","Signal component 2",*mass,sig_mean2s,gSigma2s) ;
  
  RooRealVar sigCBfrac1s("sigCBfrac1s","fraction of CB component in signal", 0.1, 0, 1) ; // fraction of cb out of cb+gaus
  RooRealVar sigCBfrac2s("sigCBfrac2s","fraction of CB component in signal", 0.1, 0, 1) ; // fraction of cb out of cb+gaus
  RooAddPdf* sig1s = new RooAddPdf("sig1s","Signal 1s", RooArgList(*cb1s,*gau1s), sigCBfrac1s);  
  RooAddPdf* sig2s = new RooAddPdf("sig2s","Signal 2s", RooArgList(*cb2s,*gau2s), sigCBfrac2s);  
  
  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0","a0",-4.7516e-01,-100,100.) ;
  RooRealVar a1("a1","a1",-1.7154e-02,-100,100)  ;
  RooChebychev* bkg = new RooChebychev("bkg","Background",*mass,RooArgSet(a0,a1)) ;
  RooRealVar *nSig1s= new RooRealVar("nSig1s","number of 1S signals",1000,0,50000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s","number of 2S signals",1000,0,10000);
  RooRealVar *nBkg= new RooRealVar("nBkg","fraction of component 1 in bkg",1000,0,10000);
  //  RooAddPdf* model = new RooAddPdf("model","sig+bkg",RooArgList(*cb, *bkg),RooArgList(*nSig1s,*nBkg));
  RooAddPdf* model = new RooAddPdf("model","sig+bkg",RooArgList(*sig1s, *sig2s, *bkg),RooArgList(*nSig1s,*nSig2s,*nBkg));
  // RooAddPdf* model = new RooAddPdf("model","sig+bkg",RooArgList(*sig1s, *bkg),RooArgList(*nSig1s,*nBkg));

  ws->import(*model);
  //  ws->import(*cb);
  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2);

  //  ws->var("sig_mean1s")->setVal(pdgJpsiMass);
  //  ws->var("sig_mean1s")->setConstant(kTRUE);
  //  ws->var("sig_mean2s")->setVal(pdgPsi2sMass);
  // ws->var("sig_mean2s")->setConstant(kTRUE);

  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE),Extended(kTRUE));
  //  RooFitResult* fitRes2 = ws->pdf("cball")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb1s)),LineColor(kRed),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*gau1s)),LineColor(kMagenta),LineStyle(kDashed));
  //ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kRed),LineStyle(kDashed));
  //ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*gau2s)),LineColor(kMagenta),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*bkg)),LineColor(kGreen),LineStyle(kDashed));
  ws->pdf("model")->paramOn(myPlot2,Layout(0.5,0.9,0.9));
  myPlot2->getAttText()->SetTextSize(0.035);
  ///  ws->data("reducedDS")->statOn(myPlot2,Layout(0.75,0.99,1));
  myPlot2->SetTitle("sig + background");
  myPlot2->SetAxisRange(massLow,massHigh,"X");
  myPlot2->Draw();
  fitRes2->Print("v");
  
    
  // PULL 

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
  pad2->SetBottomMargin(0); // Upper and lower plot are joined
  c1->cd();  
  pad2->Draw(); 
  pad2->cd();
  
  
  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  RooPlot* pullFrame = mass->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->GetYaxis()->SetTitleOffset(1.6) ;
  pullFrame->Draw() ;
  c1->SaveAs(Form("gifs/fit_pt%f-%f_y%f-%f.gif",ptLow,ptHigh, yLow, yHigh));
  
  return ws ; 
} 
 
