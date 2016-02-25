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
#include "doFitJpsiPsi.C"
using namespace std;
using namespace RooFit;

void extractSignal() {
  
  RooWorkspace* ws[10][10];
  
  for ( int ipt=1; ipt<=nPtBin ; ipt++) {
    for ( int iy=1; iy<=nYBin ; iy++) { 
      float ptLow(ptBin[ipt-1]), ptHigh(ptBin[ipt]);
      float yLow(yBin[iy-1]), yHigh(yBin[iy]);
      ws[ipt][iy] = doFitJpsiPsi(ptLow, ptHigh, yLow, yHigh);
      ws[ipt][iy]->var("nSig1s")->Print();
      changeLine();
    }
  }
  
  TH1D* hpt[nYBin];
  for ( int iy=1; iy<=nYBin ; iy++) {
    hpt[iy] = new TH1D(Form("hpt_iy%d",iy),";p_{T} (GeV); Entries",nPtBin,ptBin);
    handsomeTH1(hpt[iy],iy);
  }
  TH1D* hy[nPtBin];
  for ( int ipt=1; ipt<=nPtBin ; ipt++) {
    hy[ipt] = new TH1D(Form("hy_ipt%d",ipt),"; rapidity ; Entries",nYBin,yBin);
    handsomeTH1(hy[ipt],ipt);
  }  
  
  for ( int ipt=1; ipt<=nPtBin ; ipt++) {
    for ( int iy=1; iy<=nYBin ; iy++) {
      hpt[iy]->SetBinContent(ipt, ws[ipt][iy]->var("nSig1s")->getValV());
      hpt[iy]->SetBinError(ipt, ws[ipt][iy]->var("nSig1s")->getError());
      hy[ipt]->SetBinContent(iy, ws[ipt][iy]->var("nSig1s")->getValV());
      hy[ipt]->SetBinError(iy, ws[ipt][iy]->var("nSig1s")->getError());
    }
  }  
  
  TCanvas* cpt = new TCanvas("cpt","",500,500);
  for ( int iy=1; iy<=nYBin ; iy++) {
    hpt[iy]->SetAxisRange(0.1,10000,"Y");
    hpt[iy]->Draw( (iy==1 ? "" : "same"));
  }

  TCanvas* cy = new TCanvas("cy","",500,500);
  for ( int ipt=1; ipt<=nPtBin ; ipt++) {
    hy[ipt]->SetAxisRange(0.1,10000,"Y");
    hy[ipt]->Draw( (ipt==1 ? "" : "same"));
  }
  
  
}


