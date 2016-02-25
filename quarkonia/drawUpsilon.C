#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
using namespace std;

struct UpsilonYields { double yield1s, yield2s, yield3s, r23o123, r23o1, err1s, err2s, err3s, err23o123, err23o1; };

UpsilonYields getValFromTH1D(TString kineLabel="");

void drawUpsilon() {

  TH1::SetDefaultSumw2();
  const int nPtBins = 2;   double ptBin[nPtBins+1] = {0, 5,     100};
  const int  nYBins = 2;   double yBin[nYBins+1] =   {0, 1.2,   2.4};
  const int nPBins  = 3;   double pBin[nPBins+1] =   {0, 0.167, 0.333,  0.5};
  float muPtCut = 4;
  
  // pp VS PbPb - 100 %
  TH1D* hPt[10];

  for (int ii=0; ii<10; ii++) {
    hPt[ii] = new TH1D(Form("hPt%d",ii),";p_{T} (GeV/c);Signals",nPtBins,0+ii*0.3,10+ii*0.3); 
    hPt[ii]->SetMarkerSize(1.5);
  }
 
  for ( int ipt = 1 ; ipt<=nPtBins ; ipt++) {
    hPt[0]->GetXaxis()->SetBinLabel(ipt,  Form("%.1f - %.1f", ptBin[ipt-1], ptBin[ipt]) );
    hPt[0]->GetXaxis()->SetLabelSize(0.08);
  }
  handsomeTH1(hPt[0],1);
  TH1D* hpt_R23o1_PP[nYBins+1]; 
  TH1D* hpt_R23o1_AA[nYBins+1]; 
  TH1D* hpt_R23o1_AAcent[nYBins+1]; 
  TH1D* hpt_R23o1_AAperi[nYBins+1]; 
  TH1D* hpt_R23o1_AAplane[nYBins+1][nPBins+1];  // In-plane 



  for ( int iy = 1 ; iy<=nYBins ; iy++ )   {  float yLow = yBin[iy-1] ; float yHigh = yBin[iy] ;
    hpt_R23o1_PP[iy] = (TH1D*)hPt[0]->Clone(Form("hR23o1_PP_iy%d",iy));
    hpt_R23o1_AA[iy] = (TH1D*)hPt[0]->Clone(Form("hR23o1_AA_iy%d",iy));
    hpt_R23o1_AAcent[iy] = (TH1D*)hPt[2]->Clone(Form("hR23o1_AAcent_iy%d",iy));
    hpt_R23o1_AAperi[iy] = (TH1D*)hPt[1]->Clone(Form("hR23o1_AAperi_iy%d",iy));
    for ( int ip = 1 ; ip<=nPBins ; ip++ )   {  
      hpt_R23o1_AAplane[iy][ip] = (TH1D*)hPt[ip]->Clone(Form("hR23o1_AAplane_iy%d_ieventplane%d",iy,ip));
	}
  }
  
  for ( int icoll = 0 ; icoll <=2 ; icoll++) {
    for ( int ipt = 1 ; ipt<=nPtBins ; ipt++ ) {  float ptLow = ptBin[ipt-1];     float ptHigh = ptBin[ipt]; 
      for ( int iy = 1 ; iy<=nYBins ; iy++ )   {  float yLow = yBin[iy-1] ; float yHigh = yBin[iy] ;
  	int cLow = 0 ;  int cHigh = 200;	float dphiEp2Low=0 ; float dphiEp2High = 0.5;
	UpsilonYields uu = getValFromTH1D( getKineLabel (icoll, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High))  ;

  	if ( icoll == kPPDATA ) { 
	  hpt_R23o1_PP[iy]->SetBinContent( ipt, uu.r23o1);
	  hpt_R23o1_PP[iy]->SetBinError( ipt, uu.err23o1);
	}
  	if ( icoll == kAADATA) {
	  hpt_R23o1_AA[iy]->SetBinContent( ipt, uu.r23o1);
	  hpt_R23o1_AA[iy]->SetBinError( ipt, uu.err23o1);
  	}
      }
    }
    // peripheral :
    for ( int ipt = 1 ; ipt<=nPtBins ; ipt++ ) {  float ptLow = ptBin[ipt-1];     float ptHigh = ptBin[ipt]; 
      for ( int iy = 1 ; iy<=nYBins ; iy++ )   {  float yLow = yBin[iy-1] ; float yHigh = yBin[iy] ;
  	int cLow = 60 ;  int cHigh = 200;	
	float dphiEp2Low=0 ; float dphiEp2High = 0.5;
	UpsilonYields uu = getValFromTH1D( getKineLabel (icoll, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High))  ;
  	if ( icoll == kAADATA) {
	  hpt_R23o1_AAperi[iy]->SetBinContent( ipt, uu.r23o1);
	  hpt_R23o1_AAperi[iy]->SetBinError( ipt, uu.err23o1);
  	}
      }
    }
    // Central :
    for ( int ipt = 1 ; ipt<=nPtBins ; ipt++ ) {  float ptLow = ptBin[ipt-1];     float ptHigh = ptBin[ipt]; 
      for ( int iy = 1 ; iy<=nYBins ; iy++ )   {  float yLow = yBin[iy-1] ; float yHigh = yBin[iy] ;
  	int cLow = 0 ;  int cHigh = 60;
	float dphiEp2Low=0 ; float dphiEp2High = 0.5;
	UpsilonYields uu = getValFromTH1D( getKineLabel (icoll, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High))  ;
  	if ( icoll == kAADATA) {
	  hpt_R23o1_AAcent[iy]->SetBinContent( ipt, uu.r23o1);
	  hpt_R23o1_AAcent[iy]->SetBinError( ipt, uu.err23o1);
  	}
      }
    }
    // Event Plane dependence : 
    for ( int ipt = 1 ; ipt<=nPtBins ; ipt++ ) {  float ptLow = ptBin[ipt-1];     float ptHigh = ptBin[ipt]; 
      for ( int iy = 1 ; iy<=nYBins ; iy++ )   {  float yLow = yBin[iy-1] ; float yHigh = yBin[iy] ;
	for ( int ip = 1 ; ip<=nPBins ; ip++ )   {  float dphiEp2Low = pBin[ip-1] ; float dphiEp2High = pBin[ip] ;
	  int cLow = 20 ;  int cHigh = 120;
	UpsilonYields uu = getValFromTH1D( getKineLabel (icoll, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High))  ;
	  if ( icoll == kAADATA) {
	    hpt_R23o1_AAplane[iy][ip]->SetBinContent( ipt, uu.r23o1);
	    hpt_R23o1_AAplane[iy][ip]->SetBinError( ipt, uu.err23o1);
	  }
	}
      }
    }
      
    
  }
  
  TCanvas* c1 = new TCanvas("c1","",800,400);
  c1->Divide(2,1);
  c1->cd(1);   // upsilon 1s
  hPt[0]->SetAxisRange(0,0.9,"Y");
  hPt[0]->SetYTitle("#Upsilon (2S+3S) / 1S");
  hPt[0]->DrawCopy();
  handsomeTH1(hpt_R23o1_PP[1],1);
  handsomeTH1(hpt_R23o1_AA[1],2);
  hpt_R23o1_PP[1]->Draw("same");
  hpt_R23o1_AA[1]->Draw("same");
  drawText("Mid-rapidty |y|<1.2", 0.45,0.88,1,20);
  TLegend* leg1 = new TLegend(0.2046176,0.7200982,0.6492568,0.9504435,NULL,"brNDC");
  easyLeg(leg1,"");
  leg1->AddEntry(hpt_R23o1_PP[1], "pp");
  leg1->AddEntry(hpt_R23o1_AA[1], "PbPb 0-100%");
  leg1->Draw();

  c1->cd(2); // upsilon 2s
  hPt[0]->DrawCopy();
  handsomeTH1(hpt_R23o1_PP[2],1);
  handsomeTH1(hpt_R23o1_AA[2],2);
  hpt_R23o1_PP[2]->Draw("same");
  hpt_R23o1_AA[2]->Draw("same");
  drawText("Forward 1.2<|y|<2.4", 0.45,0.88,1,20);

    c1->SaveAs("singleRatio_1.png");  
  // trivial change
  

  TCanvas* c2 = new TCanvas("c2","",800,400);
  c2->Divide(2,1);
  c2->cd(1);   // upsilon 1s
  hPt[0]->DrawCopy();
  handsomeTH1(hpt_R23o1_AAcent[1],2,1,21);
  handsomeTH1(hpt_R23o1_AAperi[1],2,1,25);
  hpt_R23o1_PP[1]->Draw("same");
  hpt_R23o1_AAcent[1]->Draw("same");
  hpt_R23o1_AAperi[1]->Draw("same");
  drawText("Mid-rapidty |y|<1.2", 0.45,0.88,1,20);
  TLegend* leg2 = new TLegend(0.2046176,0.7200982,0.6492568,0.9504435,NULL,"brNDC");
  easyLeg(leg2,"");
  leg2->AddEntry(hpt_R23o1_PP[1], "pp");
  leg2->AddEntry(hpt_R23o1_AAperi[1], "PbPb 30-100%");
  leg2->AddEntry(hpt_R23o1_AAcent[1], "PbPb 0-30%");
  leg2->Draw();

  c2->cd(2);   // upsilon 1s
  hPt[0]->DrawCopy();
  handsomeTH1(hpt_R23o1_AAcent[2],2,1,21);
  handsomeTH1(hpt_R23o1_AAperi[2],2,1,25);
  hpt_R23o1_PP[2]->Draw("same");
  hpt_R23o1_AAcent[2]->Draw("same");
  hpt_R23o1_AAperi[2]->Draw("same");
  drawText("Forward 1.2<|y|<2.4", 0.45,0.88,1,20);
  c2->cd();
    c2->SaveAs("singleRatio_2.png");

  TCanvas* c3 = new TCanvas("c3","",800,400);
  c3->Divide(2,1);
  c3->cd(1);   // upsilon 1s
  hPt[0]->DrawCopy();
  handsomeTH1(hpt_R23o1_AAplane[1][1],4,1,22);   // [rapidity][event plane]
  handsomeTH1(hpt_R23o1_AAplane[1][2],4,1,34);
  handsomeTH1(hpt_R23o1_AAplane[1][3],4,1,23);

  hpt_R23o1_PP[1]->Draw("same");
  hpt_R23o1_AAplane[1][1]->Draw("same");
  hpt_R23o1_AAplane[1][2]->Draw("same");
  hpt_R23o1_AAplane[1][3]->Draw("same");
  drawText("Mid-rapidty |y|<1.2", 0.45,0.88,1,20);
  TLegend* leg3 = new TLegend(0.2046176,0.6800982,0.6492568,0.9504435,NULL,"brNDC");
  easyLeg(leg3,"");
  leg3->AddEntry(hpt_R23o1_PP[1], "pp");
  leg3->AddEntry(hpt_R23o1_AAplane[1][1], "PbPb #Delta#phi  [0,#pi/6]");
  leg3->AddEntry(hpt_R23o1_AAplane[1][2], "               [#pi/6,2#pi/6]");
  leg3->AddEntry(hpt_R23o1_AAplane[1][3], "               [2#pi/6,3#pi/6]");
  leg3->Draw();

  c3->cd(2);
  hPt[0]->DrawCopy();
  handsomeTH1(hpt_R23o1_AAplane[2][1],4,1,22);   // [rapidity][event plane]
  handsomeTH1(hpt_R23o1_AAplane[2][2],4,1,34);
  handsomeTH1(hpt_R23o1_AAplane[2][3],4,1,23);

  hpt_R23o1_PP[2]->Draw("same");
  hpt_R23o1_AAplane[2][1]->Draw("same");
  hpt_R23o1_AAplane[2][2]->Draw("same");
  hpt_R23o1_AAplane[2][3]->Draw("same");
  drawText("Forward |y|<1.2", 0.45,0.88,1,20);

  c3->SaveAs("singleRatio_3.png");
 
}
UpsilonYields getValFromTH1D(TString kineLabel) {
  
  UpsilonYields yields;
  TString fname = Form("fitResults/fixParam1MuPt4_Feb02/fitresults_upsilon_singleCB_%s.root",kineLabel.Data());
  cout << "opening " << fname << " ....." <<endl;
  TFile* f1 = new TFile(fname.Data() );
  if ( f1->IsZombie() )  {
    cout << " no such a file exists..." << endl << endl;
    yields.yield1s = -1 ;     yields.err1s = -1 ; 
    yields.yield2s = -1 ;     yields.err2s = -1 ; 
    yields.yield3s = -1 ;     yields.err3s = -1 ; 
  }
  else { 
    TH1D* h = (TH1D*)f1->Get("fitResults");
    
    yields.yield1s =    h->GetBinContent(1);
    yields.err1s = h->GetBinError (1);
    yields.yield2s =    h->GetBinContent(2);
    yields.err2s = h->GetBinError (2);
    yields.yield3s =    h->GetBinContent(3);
    yields.err3s = h->GetBinError (3);
    yields.r23o123 = h->GetBinContent(4);
    yields.err23o123 = h->GetBinError(4);
    yields.r23o1 = h->GetBinContent(5);
    yields.err23o1 = h->GetBinError(5);

    cout << " 1s = " << yields.yield1s << " +/- " << yields.err1s << endl;
    cout << " 2s = " << yields.yield2s << " +/- " << yields.err2s << endl;
    cout << " 3s = " << yields.yield3s << " +/- " << yields.err3s << endl;
    cout << " (2S+3S)/(1S+2S+3S) = " << yields.r23o123 << " +/- " <<  yields.err23o123 << endl;
    cout << " (2S+3S)/1S = " << yields.r23o1 << " +/- " <<  yields.err23o1 << endl;
    f1->Close();
  }
  cout << "*===*===*===*===*===*===*===*===*===*===*===*===*===*===*"<<endl;
  return yields;
    
}
