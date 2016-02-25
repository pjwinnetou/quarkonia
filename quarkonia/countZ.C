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
void countZ(
		   int collId = kPPDATA,  
		   float ptLow=0, float ptHigh=100, 
		   float yLow=0, float yHigh=1.2,
		   int cLow=0, int cHigh=200,
		   float dphiEp2Low=0, float dphiEp2High=0.5,   // In unit of PI!!
		   float muPtCut=20.0
		    ) 
{
  using namespace RooFit;
  
  float massLow        = 60;  float massHigh = 120;
  float massLowForPlot = massLow;    float massHighForPlot = massHigh;
  int   nMassBin  = massHigh - massLow;

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f)", (float)muPtCut, (float)muPtCut );
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATACentL3) || (collId==kAADATAPeri) )
    kineCut = kineCut + Form(" && (cBin>=%d && cBin<%d) && ( abs(abs(dphiEp2/3.141592)-0.5)>%.3f && abs(abs(dphiEp2/3.141592)-0.5)<%.3f )",cLow, cHigh, dphiEp2Low, dphiEp2High);
  kineCut = kineCut + (" && highPtFlag" );
  TFile* f1;
  if      ( collId == kPPDATA) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_201622292.root"); // Updated on Feb 22th.  Major change is the addition of highPt muon cut flag
  else if ( collId == kAADATA) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPbPb_L1DoubleMu0PD_Trig-L1DoubleMu0_EP-OppositeHF_201622297.root");   // Updated on Feb 22th.  Major change is the addition of highPt muon cut flag
  else if ( collId == kAADATACentL3) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPbPb_CentralPD_Trig-L3UpsilonCentral_EP-OppositeHF_2016221913.root");
  else if ( collId == kAADATAPeri) f1 = new TFile("/home/jazzitup/analysis/quarkonia/skimmedFiles/yskimPbPb_PeripheralPD_Trig-L1DoubleMu0Peripheral_EP-OppositeHF_2016222920.root");
  
  TTree* tree = (TTree*) f1->Get("mm");
  TH1D* hmass = new TH1D("hmass",";M_{#mu#mu} (GeV)",nMassBin, massLow, massHigh) ;
  hmass->Sumw2();
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,400,400);
  tree->Draw(Form("mass>>%s",hmass->GetName()),kineCut);
  handsomeTH1(hmass,1);
  cleverRange(hmass, 1.5);
  hmass->Draw();
  
  drawText(getCollID(collId),0.55,0.85,1,20);
  drawText(Form("p_{T} : %.0f - %.0f GeV",ptLow,ptHigh ),0.6,0.78,2,16);
  drawText(Form("y   : %.1f - %.1f ",yLow,yHigh ), 0.6,0.71,2,16);
  drawText(Form("(p_{T}^{#mu} > %.0f GeV)", muPtCut ), 0.65,0.63,1,15);

  
  c1->SaveAs(Form("Zcounts_%s.gif",kineLabel.Data()));
  
  TFile* outf = new TFile(Form("Zcounts_%s.root",kineLabel.Data()),"recreate");
  hmass->Write();
  c1->Write();
  outf->Close();
  cout << "Cut      : " << kineCut << endl;
  cout << "Z counts : " <<  hmass->Integral() << " +/- " << sqrt(hmass->Integral()) << endl;
  
} 
 
