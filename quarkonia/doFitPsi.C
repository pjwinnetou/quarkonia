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

using namespace std;
using namespace RooFit;
void doFitPsi( 
		   int collId = kAADATA,  
		   float ptLow=3.5, float ptHigh=6,
		   float yLow=1.6, float yHigh=2.4,
		   int cLow=20, int cHigh=120,
		   float dphiEp2Low=0, float dphiEp2High=0.5,   // In unit of PI!!
		   int accFlag=0, 
		   bool fixParameters=true
		    ) 
{
  using namespace RooFit;
  
  float massLow        = 2.5;  float massHigh = 4.3;
  float massLowForPlot = massLow;    float massHighForPlot = massHigh;

  int   nMassBin  = (int)((massHigh-massLow)/0.02);
  //  int   nMassBin  = (int)((massHigh-massLow)/0.015);  # Andre's bin

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, dphiEp2Low, dphiEp2High, accFlag) ;
  TString kineCut = Form("pt>%.1f && pt<%.1f && abs(y)>%.1f && abs(y)<%.1f",ptLow, ptHigh, yLow, yHigh);
  if (accFlag>0) kineCut = kineCut + Form(" && (accFlag==%d)",(int)accFlag );
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) )
    kineCut = kineCut + Form(" && (cBin>=%d && cBin<%d) && ( abs(abs(dphiEp2/3.141592)-0.5)>%.3f && abs(abs(dphiEp2/3.141592)-0.5)<%.3f )",cLow, cHigh, dphiEp2Low, dphiEp2High);
  
  //  kineCut = kineCut + " && pt1>4 && pt2>4";
  TFile* f1;
  if      ( collId == kPPDATA) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPP_20161201910.root");
  else if ( collId == kAADATA) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPbPb_20161201914.root");
  
  TTree* tree = (TTree*) f1->Get("dimu");
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace(Form("workspace_%s",kineLabel.Data()));
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
  
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,780,500);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.65, 1.0);
  pad1->Draw(); pad1->cd();

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"), Layout(0,1,0.95));
  // crystal ball functions 
  RooRealVar mean1s("mean1s","mean of the signal gaussian mass PDF",pdgMass.JPsi, pdgMass.JPsi-0.1, pdgMass.JPsi+0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Psi2S / pdgMass.JPsi );
  RooFormulaVar mean2s("mean2s","mean1s*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooRealVar sigma1s("sigma1s","width/sigma of the signal gaussian mass PDF",0.05, 0.01, 0.1);
  RooFormulaVar sigma2s("sigma2s","sigma1s*mRatio21", RooArgSet(sigma1s,mRatio21) );
  RooRealVar n1s("n1s","power order", 2,1,5);
  RooRealVar n2s("n2s","power order", 2,1,5);
  RooRealVar alpha1s("alpha1s","tail shift", 4, 1, 10);
  RooRealVar alpha2s("alpha2s","tail shift", 4, 1, 10);
  
  // Gaussian parameters ... 
  RooRealVar gausWidth1s(    "gausWidth1s","width of gaussian funciton of 1s peak", 0.05, 0.01, 0.1);
  RooFormulaVar gausWidth2s("gausWidth2s","gausWidth1s*mRatio21",RooArgSet(gausWidth1s,mRatio21));
 
  // Fix the parameters 
  /*  
      PSet3SingleCB psetUpsilons = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh) ; 
      if (fixParameters) { 
      if ( psetUpsilons.n1s == -1)  {
      cout << endl << endl << endl << "#########################  ERROR!!!! ##################" << endl;
      cout << "No parameter sets for " << kineLabel << endl;
      cout << "fitting macro is stopped!" << endl << endl << endl;
      return;
      }
      if ( psetUpsilons.n1s != 0     )   {    n1s.setVal(psetUpsilons.n1s);  n1s.setConstant();  }
      if ( psetUpsilons.n2s != 0     )   {    n2s.setVal(psetUpsilons.n2s);  n2s.setConstant();  }
      if ( psetUpsilons.n3s != 0     )   {    n3s.setVal(psetUpsilons.n3s);  n3s.setConstant();  }
      if ( psetUpsilons.alpha1s != 0 )   {    alpha1s.setVal(psetUpsilons.alpha1s);  alpha1s.setConstant();  }
      if ( psetUpsilons.alpha2s != 0 )   {    alpha2s.setVal(psetUpsilons.alpha2s);  alpha2s.setConstant();  }
      if ( psetUpsilons.alpha3s != 0 )   {    alpha3s.setVal(psetUpsilons.alpha3s);  alpha3s.setConstant();  }
      if ( psetUpsilons.sigma1s != 0 )   {    sigma1s.setVal(psetUpsilons.sigma1s);  sigma1s.setConstant();  }
      if ( psetUpsilons.sigma2s != 0 )   {    sigma2s.setVal(psetUpsilons.sigma2s);  sigma2s.setConstant();  }
      if ( psetUpsilons.sigma3s != 0 )   {    sigma3s.setVal(psetUpsilons.sigma3s);  sigma3s.setConstant();  }
      }
  */
  
  
  RooCBShape* cb1s = new RooCBShape("cball1s", "cystal Ball", *mass, mean1s, sigma1s, alpha1s, n1s);
  RooGaussian* gauss1s = new RooGaussian("gauss1s","gaussian for 1s",*mass, mean1s, gausWidth1s);
  RooRealVar* fracCB1s = new RooRealVar("fracCB1s", "fracCB1s", 0.5, 0,1);
  RooAddPdf*  pdf1s = new RooAddPdf("pdf1s","Singal of 1S", RooArgList(*cb1s, *gauss1s),RooArgList(*fracCB1s));
  
  RooCBShape* cb2s = new RooCBShape("cball2s", "cystal Ball", *mass, mean2s, sigma2s, alpha2s, n2s);
  RooGaussian* gauss2s = new RooGaussian("gauss2s","gaussian for 2s",*mass, mean2s, gausWidth2s);
  RooRealVar* fracCB2s = new RooRealVar("fracCB2s", "fracCB2s", 0.5, 0,1);
  RooAddPdf*  pdf2s = new RooAddPdf("pdf2s","Singal of 2S", RooArgList(*cb2s, *gauss2s),RooArgList(*fracCB2s));
  
  RooRealVar* frac1s = new RooRealVar("frac1over12","1S over 1S+2S", 0.8,0,1);
  RooAddPdf*  cb12s = new RooAddPdf("sig12S","Signal 1S+2S",RooArgList(*pdf1s,*pdf2s), RooArgList(*frac1s) );
    
  /*  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", *mass, mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape* cb2s_2 = new RooCBShape("cball2s_2", "cystal Ball", *mass, mean2s, sigma2s_2, alpha2s_2, n2s_2);
  RooCBShape* cb3s_2 = new RooCBShape("cball3s_2", "cystal Ball", *mass, mean3s, sigma3s_2, alpha3s_2, n3s_2); */
  
  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("bkg_a0","bkg_a0",-4.7516e-01,-1,1.) ;
  RooRealVar a1("bkg_a1","bkg_a1",-1.7154e-02,-1,1)  ;
  RooRealVar a2("bkg_a2","bkg_a2",-1.7154e-02,-1,1)  ;

  RooChebychev *bkg = new RooChebychev("bkg","Background",*mass,RooArgSet(a0,a1,a2) );
  //  RooRealVar *nSig1s= new RooRealVar("nSig1s"  ,"number of 1S signals"   ,3000,0,500000); 
  
  RooRealVar *nSig12s= new RooRealVar("nSig12s","number of 1S+2S signals",10000,0,1000000);
  RooRealVar *nBkg= new RooRealVar("nBkg","fraction of component 1 in bkg",1000,0,100000);
  RooAddPdf  *model = new RooAddPdf("model","sig+bkg",RooArgList(*cb12s, *bkg),RooArgList(*nSig12s, *nBkg));

  ws->import(*model);
  
  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2);
  
  /* 
     ws->var("fracCB1_1s")->setVal(1);   ws->var("fracCB1_1s")->setConstant(kTRUE);
     ws->var("fracCB1_2s")->setVal(1);   ws->var("fracCB1_2s")->setConstant(kTRUE);  
     ws->var("fracCB1_3s")->setVal(1);   ws->var("fracCB1_3s")->setConstant(kTRUE);
     ws->var("alpha1s_2")->setVal(1);   
     ws->var("alpha1s_2")->setConstant(kTRUE);
     ws->var("n1s_2")->setVal(1); 
     ws->var("n1s_2")->setConstant(kTRUE);
     ws->var("sigma1s_2")->setVal(1);  
     ws->var("sigma1s_2")->setConstant(kTRUE);
     ws->var("alpha2s_2")->setVal(1);   
     ws->var("alpha2s_2")->setConstant(kTRUE);
     ws->var("n2s_2")->setVal(1); 
     ws->var("n2s_2")->setConstant(kTRUE);
     ws->var("sigma2s_2")->setVal(1);
     ws->var("sigma2s_2")->setConstant(kTRUE);
     ws->var("alpha3s_2")->setVal(1); 
     ws->var("alpha3s_2")->setConstant(kTRUE);
     ws->var("n3s_2")->setVal(1); 
     ws->var("n3s_2")->setConstant(kTRUE);
     ws->var("sigma3s_2")->setVal(1);  
     ws->var("sigma3s_2")->setConstant(kTRUE);
  */
  

  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE),Extended(kTRUE));
  //  RooFitResult* fitRes2 = ws->pdf("cball")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*pdf1s)),LineColor(kRed),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*pdf2s)),LineColor(kRed),LineStyle(kDashed));

  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*bkg)), LineColor(kAzure+7) );
  
    

  //  ws->pdf("model")->paramOn(myPlot2,Layout(0.5,0.9,0.9));
  //  ws->data("reducedDS")->statOn(myPlot2,Layout(0.75,0.99,1));
  myPlot2->SetTitle("sig + background");
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->Draw();
  fitRes2->Print("v");
  
    
  // PULL 

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.65, 0.25);
  pad2->SetBottomMargin(0); // Upper and lower plot are joined
  c1->cd();  
  pad2->Draw(); 
  pad2->cd();
  
  
  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  RooPlot* pullFrame = mass->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->GetYaxis()->SetTitleOffset(1.6) ;
  pullFrame->Draw() ;

  TPad *pad3 = new TPad("pad3", "pad3", 0.65, 0.05, 1, 1);
  pad3->SetBottomMargin(0);
  c1->cd();  
  pad3->Draw(); 
  pad3->cd();

  RooPlot* legFrame = mass->frame(Name("Fit Results"), Title("Fit Results"));
  
  ws->pdf("model")->paramOn(legFrame,Layout(0,.95, .97));
  legFrame->getAttText()->SetTextAlign(11);
  legFrame->getAttText()->SetTextSize(0.05);

  TPaveText* hh = (TPaveText*)legFrame->findObject(Form("%s_paramBox",ws->pdf("model")->GetName()));
  hh->SetY1(0.01); hh->SetY2(0.95);
  hh->Draw();
  //legFrame->findObject(Form("%s_paramBox",ws->pdf("model")->GetName()))->Draw();
				      
  c1->SaveAs(Form("fitresults_upsilon_singleCB_%s.gif",kineLabel.Data()));
  return ;

  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);
  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  outh->GetXaxis()->SetBinLabel(2,"Upsilon2S");
  outh->GetXaxis()->SetBinLabel(3,"Upsilon3S");
  outh->GetXaxis()->SetBinLabel(4,"(2S+3S)/(1S+2S+3S)");
  outh->GetXaxis()->SetBinLabel(5,"(2S+3S)/(1S)");

  outh->SetBinContent(1,  ws->var("nSig123s")->getVal() * ws->var("frac1over123")->getVal() );
  outh->SetBinError(1,    ws->var("nSig123s")->getError() * ws->var("frac1over123")->getVal() );
  outh->SetBinContent(2,  ws->var("nSig123s")->getVal() * (1 - ws->var("frac1over123")->getVal()) * ws->var("frac2over23")->getVal() );
  outh->SetBinError(2,    ws->var("nSig123s")->getError() * (1 - ws->var("frac1over123")->getVal()) * ws->var("frac2over23")->getVal() );
  outh->SetBinContent(3,  ws->var("nSig123s")->getVal() * (1 - ws->var("frac1over123")->getVal()) * ( 1- ws->var("frac2over23")->getVal() ) );
  outh->SetBinError(3,    ws->var("nSig123s")->getError() * (1 - ws->var("frac1over123")->getVal()) * ( 1- ws->var("frac2over23")->getVal() ) );
  
  float r1o123 = ws->var("frac1over123")->getVal();
  float e1o123 = ws->var("frac1over123")->getError();
 
  outh->SetBinContent(4,  1 - r1o123 ) ;
  outh->SetBinError  (4,  e1o123 ) ;
  
  outh->SetBinContent(5, (1-r1o123) / r1o123 ) ;
  outh->SetBinError  (5, e1o123 / ( r1o123*r1o123) );

  cout << "1S signal    =  " << outh->GetBinContent(1) << " +/- " << outh->GetBinError(1) << endl;
  cout << "2S signal    =  " << outh->GetBinContent(2) << " +/- " << outh->GetBinError(2) << endl;
  cout << "3S signal    =  " << outh->GetBinContent(3) << " +/- " << outh->GetBinError(3) << endl;
  cout << "(2S+3S)/(1S+2S+3S) = " << outh->GetBinContent(4) << " +/- " << outh->GetBinError(4) << endl;
  cout << "(2S+3S)/(1S)       = " << outh->GetBinContent(5) << " +/- " << outh->GetBinError(5) << endl;
    
  TFile* outf = new TFile(Form("fitresults_upsilon_singleCB_%s.root",kineLabel.Data()),"recreate");
  outh->Write();
  c1->Write();
  outf->Close();
} 

