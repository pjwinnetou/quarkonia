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
void doFitUpsilonOneCB( 
		   int collId = kPPDATA,  
		   float ptLow=5, float ptHigh=6, 
		   float yLow=0, float yHigh=1.2,
		   int cLow=0, int cHigh=200,
		   float dphiEp2Low=0, float dphiEp2High=0.5,   // In unit of PI!!
		   float muPtCut=4.0,
		   bool fixParameters=false
		    ) 
{
  using namespace RooFit;
  
  float massLow        = 8;  float massHigh = 13.5;
  float massLowForPlot = massLow;    float massHighForPlot = massHigh;
  int   nMassBin  = 65;

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f)", (float)muPtCut, (float)muPtCut );
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATACentL3) || (collId==kAADATAPeri) )
    kineCut = kineCut + Form(" && (cBin>=%d && cBin<%d) && ( abs(abs(dphiEp2/3.141592)-0.5)>%.3f && abs(abs(dphiEp2/3.141592)-0.5)<%.3f )",cLow, cHigh, dphiEp2Low, dphiEp2High);
  
  TFile* f1;
  if      ( collId == kPPDATA) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_201622292.root"); // Updated on Feb 22th.  Major change is the addition of highPt muon cut flag
  else if ( collId == kAADATA) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPbPb_L1DoubleMu0PD_Trig-L1DoubleMu0_EP-OppositeHF_201622297.root");   // Updated on Feb 22th.  Major change is the addition of highPt muon cut flag
  else if ( collId == kAADATACentL3) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPbPb_CentralPD_Trig-L3UpsilonCentral_EP-OppositeHF_2016221913.root");
  else if ( collId == kAADATAPeri) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPbPb_PeripheralPD_Trig-L1DoubleMu0Peripheral_EP-OppositeHF_2016222920.root");
  
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
  RooRealVar mean1s("mean1s","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","mean1s*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","mean1s*mRatio31", RooArgSet(mean1s,mRatio31) );
  
  RooRealVar sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.14);
  RooRealVar sigma2s_1("sigma2s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.14);
  RooRealVar sigma3s_1("sigma3s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.14);
  //  RooRealVar sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",0.0769);
  //  RooRealVar sigma2s_1("sigma2s_1","width/sigma of the signal gaussian mass PDF",0.0816);
  //  RooRealVar sigma3s_1("sigma3s_1","width/sigma of the signal gaussian mass PDF",0.0771);
  
  
  RooRealVar n1s_1("n1s_1","power order", 4,2,10);
  RooRealVar n2s_1("n2s_1","power order", 4,2,10);
  RooRealVar n3s_1("n3s_1","power order", 4,2,10);
    
  RooRealVar alpha1s_1("alpha1s_1","tail shift", 4, 1.2, 10);
  RooRealVar alpha2s_1("alpha2s_1","tail shift", 4, 1.2, 10);  
  RooRealVar alpha3s_1("alpha3s_1","tail shift", 4, 1.2, 10);
  
  // Fix the parameters 
  PSet3SingleCB psetUpsilons = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, muPtCut) ; 
  if (fixParameters) { 
    if ( psetUpsilons.n1s == -1)  {
      cout << endl << endl << endl << "#########################  ERROR!!!! ##################" << endl;
      cout << "No parameter sets for " << kineLabel << endl;
      cout << "fitting macro is stopped!" << endl << endl << endl;
      return;
    }
    if ( psetUpsilons.n1s != 0     )   {    n1s_1.setVal(psetUpsilons.n1s);  n1s_1.setConstant();  }
    if ( psetUpsilons.n2s != 0     )   {    n2s_1.setVal(psetUpsilons.n2s);  n2s_1.setConstant();  }
    if ( psetUpsilons.n3s != 0     )   {    n3s_1.setVal(psetUpsilons.n3s);  n3s_1.setConstant();  }
    if ( psetUpsilons.alpha1s != 0 )   {    alpha1s_1.setVal(psetUpsilons.alpha1s);  alpha1s_1.setConstant();  }
    if ( psetUpsilons.alpha2s != 0 )   {    alpha2s_1.setVal(psetUpsilons.alpha2s);  alpha2s_1.setConstant();  }
    if ( psetUpsilons.alpha3s != 0 )   {    alpha3s_1.setVal(psetUpsilons.alpha3s);  alpha3s_1.setConstant();  }
    if ( psetUpsilons.sigma1s != 0 )   {    sigma1s_1.setVal(psetUpsilons.sigma1s);  sigma1s_1.setConstant();  }
    if ( psetUpsilons.sigma2s != 0 )   {    sigma2s_1.setVal(psetUpsilons.sigma2s);  sigma2s_1.setConstant();  }
    if ( psetUpsilons.sigma3s != 0 )   {    sigma3s_1.setVal(psetUpsilons.sigma3s);  sigma3s_1.setConstant();  }
  }
  
  
  
  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", *mass, mean1s, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape* cb2s_1 = new RooCBShape("cball2s_1", "cystal Ball", *mass, mean2s, sigma2s_1, alpha2s_1, n2s_1);
  RooCBShape* cb3s_1 = new RooCBShape("cball3s_1", "cystal Ball", *mass, mean3s, sigma3s_1, alpha3s_1, n3s_1);
  RooRealVar a0("bkg_a0","bkg_a0",0,-10,10.) ;
  RooRealVar a1("bkg_a1","bkg_a1",0,-10,10.) ;
  RooRealVar a2("bkg_a2","bkg_a2",0,-10,10.) ;
  RooChebychev* bkg = new RooChebychev("bkg","Background",*mass,RooArgSet(a0,a1,a2));
  
  /*  RooRealVar* frac2over23 = new RooRealVar("frac2over23","2S over 2S+3S", 0.8,0,1);
  RooAddPdf*  cb23s_1 = new RooAddPdf("cb23s_1","Signal 2s+3s",RooArgList(*cb2s_1,*cb3s_1), RooArgList(*frac2over23) );
  RooRealVar* frac1s = new RooRealVar("frac1over123","1S over 1S+2S+3S", 0.8,0,1);
  RooAddPdf*  cb123s_1 = new RooAddPdf("cb123s_1","Signal 1S+2S+3S",RooArgList(*cb1s_1,*cb23s_1), RooArgList(*frac1s) );
  RooRealVar a0("bkg_a0","bkg_a0",0,-10,10.) ;
  RooRealVar a1("bkg_a1","bkg_a1",0,-10,10.) ;
  RooRealVar a2("bkg_a2","bkg_a2",0,-10,10.) ;
  RooChebychev* bkg = new RooChebychev("bkg","Background",*mass,RooArgSet(a0,a1,a2));
  RooRealVar *nSig123s= new RooRealVar("nSig123s","number of 1S+2S+3S signals",30000,0,100000);
  RooRealVar *nBkg= new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,1000000);   */

  RooRealVar* frac2over1 = new RooRealVar("frac2over1","2S/1S", 0.2,0,1);
  RooAddPdf*  cb12s_1 = new RooAddPdf("cb12s","Signal 1S+2S",RooArgList(*cb2s_1,*cb1s_1), RooArgList(*frac2over1) );
  RooRealVar *nSig12s= new RooRealVar("nSig12s"," 1S+2S signals",30000,0,100000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",300,0,10000);
  RooRealVar *nBkg= new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,1000000);  

  RooAddPdf* model = new RooAddPdf("model","(1S+2S)+3S + Bkg",RooArgList(*cb12s_1, *cb3s_1, *bkg),RooArgList(*nSig12s, *nSig3s, *nBkg));

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
  
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb1s_1)),LineColor(kRed),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s_1)),LineColor(kRed),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s_1)),LineColor(kRed),LineStyle(kDashed));

  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*bkg)), LineColor(kAzure+7) );
  
    

  //  ws->pdf("model")->paramOn(myPlot2,Layout(0.5,0.9,0.9));
  //  ws->data("reducedDS")->statOn(myPlot2,Layout(0.75,0.99,1));
  myPlot2->SetTitle("sig + background");
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->Draw();
  fitRes2->Print("v");
  
  drawText(getCollID(collId),0.55,0.85,1,20);
  drawText(Form("p_{T} : %.2f  -  %.2f GeV",ptLow,ptHigh ),0.55,0.78,2,20);
  drawText(Form("y   : %.2f  -  %.2f ",yLow,yHigh ), 0.55,0.71,2,20);
  drawText(Form("(p_{T}^{#mu} > %.2f GeV)", muPtCut ), 0.55,0.63,1,15);


    
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
  
  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);
  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  outh->GetXaxis()->SetBinLabel(2,"Upsilon2S");
  outh->GetXaxis()->SetBinLabel(3,"Upsilon3S");
  outh->GetXaxis()->SetBinLabel(4,"2S/1S");
  outh->GetXaxis()->SetBinLabel(5,"3S/1S");
  
  float temp12 = ws->var("nSig12s")->getVal();  
  float temp12err = ws->var("nSig12s")->getError();
  float farc2o1 = ws->var("frac2over1")->getVal();  
  float farc2o1err = ws->var("frac2over1")->getError();  
  float temp3 = ws->var("nSig3s")->getVal();  
  float temp3err = ws->var("nSig3s")->getError();  
  
  outh->SetBinContent(1,  temp12 * (1.-farc2o1 ) ) ;
  outh->SetBinError  (1,  sqrt( TMath::Power(temp12err * (1.-farc2o1), 2) + TMath::Power(temp12 * farc2o1err, 2)  )  ) ;
  outh->SetBinContent(2,  temp12 * farc2o1  ) ;
  outh->SetBinError  (2,  sqrt( TMath::Power(temp12err * farc2o1, 2) + TMath::Power(temp12 * farc2o1err, 2)  )  ) ;
  outh->SetBinContent(3,  temp3 ) ;
  outh->SetBinError  (3,  temp3err ) ;
  outh->SetBinContent(4,  farc2o1 ); 
  outh->SetBinError  (4,  farc2o1err ); 
  outh->SetBinContent(5,  temp3 / outh->GetBinContent(1) ) ;
  outh->SetBinError  (5,  temp3 / outh->GetBinContent(1) * sqrt( TMath::Power(outh->GetBinError(1) / outh->GetBinContent(1), 2) + TMath::Power(temp3err/temp3,2)) );
  
  cout << "1S signal    =  " << outh->GetBinContent(1) << " +/- " << outh->GetBinError(1) << endl;
  cout << "2S signal    =  " << outh->GetBinContent(2) << " +/- " << outh->GetBinError(2) << endl;
  cout << "3S signal    =  " << outh->GetBinContent(3) << " +/- " << outh->GetBinError(3) << endl;
  cout << "2S/1S        =  " << outh->GetBinContent(4) << " +/- " << outh->GetBinError(4) << endl;
  cout << "3S/1S        =  " << outh->GetBinContent(5) << " +/- " << outh->GetBinError(5) << endl;

  
  TFile* outf = new TFile(Form("fitresults_upsilon_singleCB_%s.root",kineLabel.Data()),"recreate");
  outh->Write();
  c1->Write();
  outf->Close();


  ///  cout parameters :
  cout << "N1,2,3, Alpha1,2,3, Sigma1,2,3 " << endl;
  cout << "if ( ( ptLow == (float)"<< ptLow <<" ) && (ptHigh == (float)"<<ptHigh<<" ) && (yLow == (float)"<<yLow<<" ) && (yHigh == (float)"<<yHigh<<" ) )" << endl;
  cout << " {ret.setNAlphaSigma( " ;
  cout <<  ws->var("n1s_1")->getVal() << ", " <<  ws->var("n2s_1")->getVal() << ", "<<  ws->var("n2s_1")->getVal() << ", " << endl;
  cout <<  ws->var("alpha1s_1")->getVal() << ", " <<  ws->var("alpha2s_1")->getVal() << ", "<<  ws->var("alpha2s_1")->getVal() << ", " << endl;
  cout <<  ws->var("sigma1s_1")->getVal() << ", " <<  ws->var("sigma2s_1")->getVal() << ", "<<  ws->var("sigma2s_1")->getVal() << " );} " << endl;
  






} 
 
