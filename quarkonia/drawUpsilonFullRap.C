#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
#include <TGraphErrors.h>
using namespace std;

struct UpsilonYields { double yield1s, yield2s, yield3s, r2o1, r3o1, err1s, err2s, err3s, err2o1, err3o1; };
struct ValErr{ double value, error; };

UpsilonYields getValFromTH1D(TString kineLabel="");
ValErr getZyield(TString kineLabel="");
void drawUpsilonFullRap(float muPtCut = 4, int nCBins=9) {

  TH1::SetDefaultSumw2();
  //  const int nPtBins = 1;   double ptBin[nPtBins+1] = {0,  100};
  float ptLow = 0;   float ptHigh = 100;
  float yLow = 0;   float yHigh = 2.4;
  const int  nCBins9 = 9;
  double cBin9[nCBins9+1]  =  {0, 10, 20, 40, 60, 80, 100, 120,140,200}; 
  double nPart9[nCBins9+1] =  {0, 8.30, 30.59, 53.85, 86.94, 131.4, 189.2, 264.3, 333.4, 384.4};
  const int  nCBins7 = 7;
  double cBin7[nCBins7+1]  =  {0, 10, 20, 40, 60, 80, 100, 200};
  double nPart7[nCBins7+1] =  {0, 21.85, 86.94, 131.4, 189.2, 264.3, 333.4, 384.4};
  

  double* cBin;
  double* nPart;

  if     ( nCBins==7 )  
    { cBin = cBin7;     nPart = nPart7;  }
  else if ( nCBins==9 ) 
    { cBin = cBin9;     nPart = nPart9;  }
  
  
  
  double nPartBoundary[nCBins+1];
  for ( int ii=1; ii<= nCBins ; ii++) {
    nPartBoundary[0] = nPart[0]; // = 0;
    nPartBoundary[ii] = 2* nPart[ii] - nPartBoundary[ii-1];
    cout << "     nPartBoundary_" << ii<< " = " <<     nPartBoundary[ii] << endl;
  }

  for ( int ii=1; ii<= nCBins ; ii++) { // In order to avoid the overlap
    nPartBoundary[ii] = nPartBoundary[ii] + 0;
  }
  
  
  // Ncoll from https://twiki.cern.ch/twiki/pub/CMS/HI2015DailyMeetings/Ncoll_Npart_04Dec2015.pdf
  //  double nColl[nCBins+1] =  { 0, 1626 *(0.1/392.5) , (1005+606)*(0.1/392.5), (348.3+186.2)*(0.1/392.5), (90.69+40.14)*(0.1/392.5), (15.87+5.502+1.642)*(0.1/392.5) };  // 392.5 = average Ncoll in 0-100%
  // Npart 
  
  // pp VS PbPb - 100 %
  TH1D* h1sAA = new TH1D("h1sAA", ";N_{Part}; PbPb/pp #Upsilon(1S) in fit", nCBins, nPartBoundary);
  TH1D* h1sRatio_nPart = (TH1D*)h1sAA->Clone("h1sRatio_nPart");
  //  for ( int ic = 1 ; ic<=nCBins ; ic++) {
  //    h1sAA->GetXaxis()->SetBinLabel( ic,  Form("%d - %d %%", (int)cBin[ic-1]/2, (int)cBin[ic]/2 ) );
  //    h1sAAnPart->GetXaxis()->SetBinLabel( nCBins - ic + 1,  Form("%d - %d %%", (int)cBin[ic-1]/2, (int)cBin[ic]/2 ) );
  //  }
  //  h1sAA->GetXaxis()->SetLabelSize(0.06); 
  //  h1sAAnPart->GetXaxis()->SetLabelSize(0.06); 
  
  TH1D* hR2o1 = (TH1D*)h1sAA->Clone("hR2o1");
  TH1D* hR2o1_nPart = (TH1D*)h1sRatio_nPart->Clone("hR2o1_nPart");
  TH1D* hDoubleR_nPart = (TH1D*)h1sRatio_nPart->Clone("hDoubleR_nPart");
  TH1D* zRatio_nPart  = (TH1D*)h1sRatio_nPart->Clone("zRatio");

  
  UpsilonYields uuPP;
  float ppVal, ppErr ; 
  float r23o1pp = 0;
  ValErr zPP;
  ValErr Y1sPP;
  for ( int icoll = 0 ; icoll <=2 ; icoll++) {
    float dphiEp2Low=0 ; float dphiEp2High = 0.5;
    // pp 
    if ( icoll == kPPDATA ) { 
      int cLow = 0 ; int cHigh = 0 ;
      uuPP = getValFromTH1D( getKineLabel (icoll, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High))  ;
      ppVal = uuPP.r2o1 ;
      ppErr = uuPP.err2o1 ;
      
      Y1sPP.value = uuPP.yield1s ;      Y1sPP.error = uuPP.err1s ;
      zPP = getZyield( getKineLabel (icoll, ptLow, ptHigh, yLow, yHigh, 20, cLow, cHigh, dphiEp2Low, dphiEp2High));
      cout << "zpp = " << zPP.value << " +/- " << zPP.error << endl;
    }
    // pbpb 
    for ( int ic = 1 ; ic<=nCBins ; ic++ )   {  int cLow = cBin[ic-1] ; int cHigh = cBin[ic] ;

      if ( icoll == kAADATA) {


	bool usePeripheralPD = true;
	int tempDataset = icoll; 
	if (usePeripheralPD && (cLow >= 100) ) 
	  tempDataset = kAADATAPeri ;

	UpsilonYields uu = getValFromTH1D( getKineLabel (tempDataset, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High))  ;
	ValErr zAA = getZyield( getKineLabel (tempDataset, ptLow, ptHigh, yLow, yHigh, 20, cLow, cHigh, dphiEp2Low, dphiEp2High));
	
	cout << "z in AA = " << zAA.value << " +/- " << zAA.error << endl;
	
	zRatio_nPart->SetBinContent( nCBins-ic+1, zAA.value / zPP.value) ;
	zRatio_nPart->SetBinError( nCBins-ic+1, zAA.value / zPP.value * sqrt ( zAA.error*zAA.error/zAA.value/zAA.value + zPP.error*zPP.error/zPP.value/zPP.value  ) );
	h1sRatio_nPart->SetBinContent( nCBins-ic+1, uu.yield1s / Y1sPP.value);
	h1sRatio_nPart->SetBinError( nCBins-ic+1, uu.yield1s / Y1sPP.value  * sqrt (uu.err1s*uu.err1s/uu.yield1s/uu.yield1s + Y1sPP.error*Y1sPP.error/Y1sPP.value/Y1sPP.value  ) ) ; 
	
	
	
	float AAVal = uu.r2o1 ;
	float AAErr = uu.err2o1 ; 
	hR2o1_nPart->SetBinContent( nCBins-ic+1, AAVal) ;
	hR2o1_nPart->SetBinError( nCBins-ic+1, AAErr );
	
	hDoubleR_nPart->SetBinContent( nCBins-ic+1, AAVal / ppVal ) ;
	hDoubleR_nPart->SetBinError( nCBins-ic+1, AAErr / ppVal ); //pp error will be in error band


	
	
      }
    }
  }

  TCanvas* c2 = new TCanvas("c2","",400,400);
  handsomeTH1(hR2o1_nPart,2);
  TH1D* htemp = new TH1D("htemp",";N_{part} ;#Upsilon(2S)/#Upsilon(1S)",420,0,420);
  htemp->SetAxisRange(0,0.7,"Y");
  htemp->SetYTitle("#Upsilon(2S) / #Upsilon(1S)");
  htemp->DrawCopy();
  hR2o1_nPart->Draw("same");
  
  drawText(Form("|y| < %.2f",yHigh), 0.3,0.88,2,20);
  drawText(Form("p_{T}^{#mu} > %.1f GeV", muPtCut ), 0.3,0.80,2,18);
  drawText("pp value w/ uncertainty band)", 0.35,0.63,1,18);
  
  TBox *ppBand = new TBox(0, ppVal-ppErr, 400, ppVal+ppErr);
  ppBand->SetLineColor(kYellow+2);
  ppBand->SetFillColor(kYellow+2);
  ppBand->SetFillStyle(3001);
  ppBand->Draw();
  c2->SaveAs(Form("singleRatio_%fGeV.png",muPtCut));
  
  TCanvas* c3 = new TCanvas("c3","",400,400);
  handsomeTH1(hR2o1_nPart,2);
  htemp->SetAxisRange(0,1.4,"Y");
  htemp->SetYTitle("(PbPb/pp) of #Upsilon(2S)/#Upsilon(1S)");
  htemp->DrawCopy();
  hDoubleR_nPart->Draw("same");


  drawText(Form("|y| < %.2f",yHigh), 0.3,0.88,2,20);
  drawText(Form("p_{T}^{#mu} > %.1f GeV", muPtCut ), 0.3,0.80,2,18);
  TBox *ppBandDR = new TBox(0, 1 - ppErr/ppVal, 420, 1 + ppErr/ppVal);
  ppBandDR->SetLineColor(kYellow+2);
  ppBandDR->SetFillColor(kYellow+2);
  ppBandDR->SetFillStyle(3001);
  ppBandDR->Draw();
  /*  TLegend *leg = new TLegend(0.5527638,0.5534759,0.8894472,0.7566845,NULL,"brNDC");
  easyLeg(leg,"");
  leg->AddEntry(hDoubleR_nPart,"Yongsun","p");
  leg->AddEntry(gDoubleRatioChad,"Chad","p");*/
  //  leg->Draw();

  c3->SaveAs(Form("doubleRatio_%fGeV.png",muPtCut));
  TGraphErrors* gR2o1_nPart =  new TGraphErrors(hR2o1_nPart);
  gR2o1_nPart->SetName("gR2o1_nPart");
  TGraphErrors* gDoubleR_nPart =  new TGraphErrors(hDoubleR_nPart);
  gDoubleR_nPart->SetName("gDoubleR_nPart");
  
  
  // Z normalied RAA 
  TCanvas* c4 = new TCanvas("c4","",400,400);
  handsomeTH1(h1sRatio_nPart);
  h1sRatio_nPart->Divide(zRatio_nPart);
  h1sRatio_nPart->SetYTitle("R normalized by Z");
  h1sRatio_nPart->Draw();
  
  TFile* fout = new TFile(Form("upsilon_yongsun_%dCentralityBins_20160217.root",(int)nCBins),"recreate");
  gR2o1_nPart->Write();
  gDoubleR_nPart->Write();
  
  fout->Close();
}


UpsilonYields getValFromTH1D(TString kineLabel) {
  
  UpsilonYields yields;
  TString fname = Form("fitResults/fixParam1MuPt4_fullRap_feb22/fitresults_upsilon_singleCB_%s.root",kineLabel.Data());

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
    yields.r2o1 = h->GetBinContent(4);
    yields.err2o1 = h->GetBinError(4);
    yields.r3o1 = h->GetBinContent(5);
    yields.err3o1 = h->GetBinError(5);

    cout << " 1s = " << yields.yield1s << " +/- " << yields.err1s << endl;
    cout << " 2s = " << yields.yield2s << " +/- " << yields.err2s << endl;
    cout << " 3s = " << yields.yield3s << " +/- " << yields.err3s << endl;
    cout << " 2s/1s = " << yields.r2o1 << " +/- " <<  yields.err2o1 << endl;
    cout << " 3s/1s = " << yields.r3o1 << " +/- " <<  yields.err3o1 << endl;
    f1->Close();
  }
  cout << "*===*===*===*===*===*===*===*===*===*===*===*===*===*===*"<<endl;
  return yields;
    
}

ValErr getZyield(TString kineLabel ) { 

  ValErr yields;
  TString fname = Form("zCountResults/MuPt20_fullRap_feb22/Zcounts_%s.root",kineLabel.Data() );
  
  cout << "opening " << fname << " ....." <<endl;
  TFile* f1 = new TFile(fname.Data() );
  if ( f1->IsZombie() )  {
    cout << " no such a file exists..." << endl << endl;
    yields.value = -1 ;  yields.error = -1;
  }
  else {
    TH1D* h = (TH1D*)f1->Get("hmass");
    yields.value = h->Integral();
    yields.error = sqrt( yields.value ) ;
    
  }
  
  return yields;

}
