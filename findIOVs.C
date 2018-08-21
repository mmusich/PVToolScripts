#include "TFile.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TObjArray.h"
#include "TProcPool.h"
#include "TList.h"
#include "TMath.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <TStopwatch.h>
#include "TArrow.h"
#include "TCanvas.h"
#include "TObjString.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <functional>
#include <iterator>
#include <fstream>
#include <bitset>
#include <sstream>

double square ( const double a )
{
  return a*a;
}

/*--------------------------------------------------------------------*/
TH1F* DrawConstant(TH1F *hist,Int_t iter,Double_t theConst)
/*--------------------------------------------------------------------*/
{ 
  
  Int_t nbins       = hist->GetNbinsX();
  Double_t lowedge  = hist->GetBinLowEdge(1);
  Double_t highedge = hist->GetBinLowEdge(nbins+1);
  
  TH1F *hzero = new TH1F(Form("hconst_%s_%i",hist->GetName(),iter),Form("hconst_%s_%i",hist->GetName(),iter),nbins,lowedge,highedge);
  for (Int_t i=0;i<=hzero->GetNbinsX();i++){
    hzero->SetBinContent(i,theConst);
    hzero->SetBinError(i,0.);
  }
  hzero->SetLineWidth(2);
  hzero->SetLineStyle(9);
  hzero->SetLineColor(kMagenta);
  
  return hzero;
}

/*--------------------------------------------------------------------*/
TH1F* DrawWithErrors(TH1F *hist,Int_t iter,float err)
/*--------------------------------------------------------------------*/
{ 
  
  Int_t nbins       = hist->GetNbinsX();
  Double_t lowedge  = hist->GetBinLowEdge(1);
  Double_t highedge = hist->GetBinLowEdge(nbins+1);
  
  TH1F *hzero = new TH1F(Form("hwithErr_%s_%i",hist->GetName(),iter),Form("hwithErr_%s_%i",hist->GetName(),iter),nbins,lowedge,highedge);
  for (Int_t i=0;i<=hzero->GetNbinsX();i++){
    hzero->SetBinContent(i,hist->GetBinContent(i));
    hzero->SetBinError(i,err);
  }  
  return hzero;
}

/*--------------------------------------------------------------------*/
TH1F* DrawConstantWithErr(TH1F *hist,Int_t iter,Double_t theConst)
/*--------------------------------------------------------------------*/
{ 

  Int_t nbins       = hist->GetNbinsX();
  Double_t lowedge  = hist->GetBinLowEdge(1);
  Double_t highedge = hist->GetBinLowEdge(nbins+1);


  TH1F *hzero = new TH1F(Form("hconstwithErr_%s_%i",hist->GetName(),iter),Form("hconstwithErr_%s_%i",hist->GetName(),iter),nbins,lowedge,highedge);
  for (Int_t i=0;i<=hzero->GetNbinsX();i++){
    hzero->SetBinContent(i,theConst);
    hzero->SetBinError(i,hist->GetBinError(i));
  }
  hzero->SetLineWidth(2);
  hzero->SetLineStyle(9);
  hzero->SetLineColor(kMagenta);
  
  return hzero;
}

/*--------------------------------------------------------------------*/
void makeNiceTrendPlotStyle(TH1 *hist,bool hasexceeded)
/*--------------------------------------------------------------------*/
{ 
  TString myTitle = hist->GetName();
  hist->SetStats(kFALSE);  
  hist->SetLineWidth(2);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.065);
  hist->GetYaxis()->SetTitleSize(0.065);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.05);
  hist->GetXaxis()->SetLabelSize(.07);
  if(myTitle.Contains("chi2Score")){
    hist->SetStats(kTRUE);
    hist->GetXaxis()->SetLabelSize(.05);
  }
  //hist->GetXaxis()->SetNdivisions(505);
  hist->SetMarkerSize(1.);
  hist->SetMarkerStyle(21);
  if(hasexceeded){
    hist->SetLineColor(kRed);
    hist->SetMarkerColor(kRed);
  } else {
    hist->SetLineColor(kBlack);
    hist->SetMarkerColor(kBlack);
  }
}

/*--------------------------------------------------------------------*/
void makeNewXAxis (TH1 *h)
/*--------------------------------------------------------------------*/
{
  
  TString myTitle = h->GetName();
  float axmin = -999;
  float axmax = 999.;
  int ndiv = 510;
  if(myTitle.Contains("eta")){
    axmin = -2.7;
    axmax = 2.7;
    ndiv = 505;
  } else if (myTitle.Contains("phi")){
    axmin = -TMath::Pi();
    axmax = TMath::Pi();
    ndiv = 510;
  } else if (myTitle.Contains("pT")) {
    axmin = 0;
    axmax = 19.99;
    ndiv = 510;
  } else if (myTitle.Contains("ladder")) {
    axmin = 0;
    axmax = 12;
    ndiv  = 510;
  } else if (myTitle.Contains("modZ")){
    axmin = 0;
    axmax = 8;
    ndiv  = 510;
  } else {
    std::cout<<"unrecognized variable"<<std::endl;
  }
  
  // Remove the current axis
  h->GetXaxis()->SetLabelOffset(999);
  h->GetXaxis()->SetTickLength(0);
  
   // Redraw the new axis
  gPad->Update();
  
  TGaxis *newaxis = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),
			       gPad->GetUxmax(),gPad->GetUymin(),
			       axmin,
			       axmax,
			       ndiv,"SDH");
  
  TGaxis *newaxisup =  new TGaxis(gPad->GetUxmin(),gPad->GetUymax(),
                                  gPad->GetUxmax(),gPad->GetUymax(),
                                  axmin,
                                  axmax,                          
                                  ndiv,"-SDH");
    
  newaxis->SetLabelOffset(0.02);
  newaxis->SetLabelFont(42);
  newaxis->SetLabelSize(0.05);
  
  newaxisup->SetLabelOffset(-0.02);
  newaxisup->SetLabelFont(42);
  newaxisup->SetLabelSize(0);
  
  newaxis->Draw();
  newaxisup->Draw();

}

/*--------------------------------------------------------------------*/
void makeNiceLegend(TLegend* &lego)
/*--------------------------------------------------------------------*/
{
  lego->SetFillColor(10);
  lego->SetTextSize(0.042);
  lego->SetTextFont(42);
  lego->SetFillColor(10);
  lego->SetLineColor(10);
  lego->SetShadowColor(10);
}


/*--------------------------------------------------------------------*/
bool exceedsLimit(TH1F* hist,float limit=10.)
/*--------------------------------------------------------------------*/
{
  bool result = ( (hist->GetMaximum() > limit ) || ( hist->GetMinimum() < -limit ) );
  if(result){
    std::cout << hist->GetName() <<  "| maximum: "  << hist->GetMaximum()
	      << " minumum:   "<< hist->GetMinimum() << std::endl;
  }
  return result;
}

/*--------------------------------------------------------------------*/
bool doKSTest(TH1F* hist, double& theScore)
/*--------------------------------------------------------------------*/
{
  TH1F* theZero = DrawConstantWithErr(hist,1,1.);
  TH1F* displaced = (TH1F*)hist->Clone("displaced");
  displaced->Add(theZero);
  theScore   = displaced->KolmogorovTest(theZero);
  //Double_t chi2Score = displaced->Chi2Test(theZero);
  std::cout << hist->GetName() << "KS prob:" << theScore << std::endl;
  if(theScore<0.5) {
    return true;
  } else {
    return false;
  }
}

/*--------------------------------------------------------------------*/
double doChi2Test(TH1F* &hist,TH1F* &hist2)
/*--------------------------------------------------------------------*/
{
  /*
  double theMin = (hist->GetMinimum() <  hist2->GetMinimum()) ? std::abs(hist->GetMinimum()) :  std::abs(hist2->GetMinimum());
  TH1F* theZero = DrawConstant(hist,1,theMin);
  TH1F* displaced1 = (TH1F*)hist->Clone("displaced1");
  TH1F* displaced2 = (TH1F*)hist2->Clone("displaced2");
  if(theMin<0.){
    displaced1->Add(theZero);
    displaced2->Add(theZero);
  }
  */
  //theScore   = displaced->KolmogorovTest(theZero);
  //theScore = hist->Chi2Test(hist2);
  //theScore   = hist->KolmogorovTest(hist2);
  //theScore = displaced1->Chi2Test(displaced2);
  //std::cout << hist->GetName() << "chi2 test:" << theScore << std::endl;

  TH1F* hist_new  = DrawWithErrors(hist,1,10);
  TH1F* hist2_new = DrawWithErrors(hist2,2,10);

  double theScore = hist_new->Chi2Test(hist2_new);

  delete hist_new;
  delete hist2_new;
  return theScore;
  
}

/*--------------------------------------------------------------------*/
double buildMetric(TH1F* &hist,TH1F* &hist2)
/*--------------------------------------------------------------------*/
{
  double theScore   = 0.;
  Int_t nbins       = hist->GetNbinsX();
  Double_t lowedge  = hist->GetBinLowEdge(1);
  Double_t highedge = hist->GetBinLowEdge(nbins+1);
  for (Int_t i=1;i<=hist->GetNbinsX();i++){
    theScore+= std::abs(hist->GetBinContent(i)-hist2->GetBinContent(i))/sqrt(square(std::max(5.,hist->GetBinError(i))) + square(std::max(5.,hist2->GetBinError(i))));
  }

  std::cout << hist->GetName() << "chi2 test:" << theScore << std::endl;

  return theScore;
}


std::vector<int> list_files(const char *dirname=".", const char *ext=".root");

/////////////////////////////////////////////////////////////////////////
//
//  MAIN FUNCTION
//
/////////////////////////////////////////////////////////////////////////

/*--------------------------------------------------------------------*/
void findIOVs(const char* dirname)
/*--------------------------------------------------------------------*/
{

  gStyle->SetOptTitle(0);

  std::vector<int> currentList = list_files(dirname);
  std::sort(currentList.begin(),currentList.end());

  const unsigned int nFiles = currentList.size(); 

  TFile *fins[nFiles];

  TH1F* dxyPhiMeanTrend[nFiles]; 
  TH1F* dzPhiMeanTrend[nFiles];  
  TH1F* dxyEtaMeanTrend[nFiles]; 
  TH1F* dzEtaMeanTrend[nFiles];  
  
  TH1F* dxyPhiMeanDiff[nFiles]; 
  TH1F* dzPhiMeanDiff[nFiles];  
  TH1F* dxyEtaMeanDiff[nFiles]; 
  TH1F* dzEtaMeanDiff[nFiles];  

  TH1F* chi2Score_dxy_phi = new TH1F("chi2Score_dxy_phi","#chi^{2} score d_{xy} vs #phi;#chi^{2} score;n. runs",100,0.,100.);
  TH1F* chi2Score_dxy_eta = new TH1F("chi2Score_dxy_eta","#chi^{2} score d_{xy} vs #eta;#chi^{2} score;n. runs",100,0.,100.);
  TH1F* chi2Score_dz_phi  = new TH1F("chi2Score_dz_phi" ,"#chi^{2} score d_{z} vs #eta ;#chi^{2} score;n. runs",100,0.,100.);
  TH1F* chi2Score_dz_eta  = new TH1F("chi2Score_dz_eta" ,"#chi^{2} score d_{z} vs #eta ;#chi^{2} score;n. runs",100,0.,100.);

  for(unsigned int j=0;j<nFiles;j++){
    std::cout << "count:" << j << std::endl;
    auto run = currentList[j];
    fins[j] = TFile::Open(Form("%s/PVValidation_%s_%i.root",dirname,dirname,run),"R");
    dxyPhiMeanTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_phi");
    dzPhiMeanTrend[j]   = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_phi");
    dxyEtaMeanTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_eta");
    dzEtaMeanTrend[j]   = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_eta");
  }

  TH1F* theConst10      = DrawConstant(dxyPhiMeanTrend[0],1, 10.);
  TH1F* theConstMinus10 = DrawConstant(dxyPhiMeanTrend[0],2,-10.);
  TH1F* theConst20      = DrawConstant(dxyPhiMeanTrend[0],1, 20.);
  TH1F* theConstMinus20 = DrawConstant(dxyPhiMeanTrend[0],2,-20.);

  std::ofstream outfile ("log_IOVs.txt"); 

  TCanvas dummyC("dummyC","dummyC",1200,1200);  
  dummyC.Print("diff.pdf[");

  TCanvas dummyC2("dummyC2","dummyC2",1200,1200);  
  dummyC2.Print("superimpose.pdf[");

  for(unsigned int i=1;i<nFiles;i++){

    bool writeIOV(false);
    bool exceeds_dxy_phi(false);
    bool exceeds_dxy_eta(false);
    bool exceeds_dz_phi(false);
    bool exceeds_dz_eta(false);

    TCanvas *cv = new TCanvas(Form("%s_%i",dirname,i),dirname,1200,1200);
    cv->Divide(2,2);

    TCanvas *csup = new TCanvas(Form("sup_%s_%i",dirname,i),dirname,1200,1200);
    csup->Divide(2,2);

    TLegend *lego = new TLegend(0.62,0.83,0.92,0.93);
    makeNiceLegend(lego);

    TLegend *lego_dxy_phi = new TLegend(0.62,0.80,0.92,0.93);
    makeNiceLegend(lego_dxy_phi);

    TLegend *lego_dxy_eta = new TLegend(0.62,0.80,0.92,0.93);
    makeNiceLegend(lego_dxy_eta);

    TLegend *lego_dz_phi = new TLegend(0.62,0.80,0.92,0.93);
    makeNiceLegend(lego_dz_phi);

    TLegend *lego_dz_eta = new TLegend(0.62,0.80,0.92,0.93);
    makeNiceLegend(lego_dz_eta);

    TPaveText *ptDate =new TPaveText(0.19,0.95,0.69,0.99,"blNDC");
    ptDate->SetFillColor(kYellow);
    //ptDate->SetFillColor(10);
    ptDate->SetBorderSize(1);
    ptDate->SetLineColor(kBlue);
    ptDate->SetLineWidth(1);
    ptDate->SetTextFont(42);
    TText *textDate = ptDate->AddText(Form("#Delta: %i - %i",currentList[i],currentList[i-1]));
    textDate->SetTextSize(0.04);
    textDate->SetTextColor(kBlue);
    textDate->SetTextAlign(22);

    for(Int_t k=0; k<4; k++){
      cv->cd(k+1)->SetBottomMargin(0.14);
      cv->cd(k+1)->SetLeftMargin(0.18);
      cv->cd(k+1)->SetRightMargin(0.01);
      cv->cd(k+1)->SetTopMargin(0.06);

      csup->cd(k+1)->SetBottomMargin(0.14);
      csup->cd(k+1)->SetLeftMargin(0.18);
      csup->cd(k+1)->SetRightMargin(0.01);
      csup->cd(k+1)->SetTopMargin(0.06);
    }

    // dxy vs phi

    dxyPhiMeanDiff[i] = (TH1F*)(dxyPhiMeanTrend[i]->Clone(Form("%s_clone",dxyPhiMeanTrend[i]->GetName())));
    dxyPhiMeanDiff[i]->Add(dxyPhiMeanTrend[i-1],-1); 
    
    if(exceedsLimit(dxyPhiMeanDiff[i])){
      writeIOV=true;
      exceeds_dxy_phi=true;
    }

    dxyPhiMeanDiff[i]->GetYaxis()->SetRangeUser(-30.,30.);
    cv->cd(1);
    
    makeNiceTrendPlotStyle(dxyPhiMeanDiff[i],exceeds_dxy_phi);
    dxyPhiMeanDiff[i]->Draw();
    makeNewXAxis(dxyPhiMeanDiff[i]);

    lego->AddEntry(dxyPhiMeanDiff[i],Form("Run %i",currentList[i]),"L");
    lego->Draw("same");

    theConst10->Draw("PLsame");
    theConstMinus10->Draw("PLsame");
    ptDate->Draw("same");

    //------------------------------------------------------------------------------

    csup->cd(1);
    makeNiceTrendPlotStyle(dxyPhiMeanTrend[i],false);
    makeNiceTrendPlotStyle(dxyPhiMeanTrend[i-1],true);
    dxyPhiMeanTrend[i]->GetYaxis()->SetRangeUser(-30.,30.);
    dxyPhiMeanTrend[i]->Draw();
    dxyPhiMeanTrend[i-1]->Draw("same");
    makeNewXAxis(dxyPhiMeanTrend[i]);
    makeNewXAxis(dxyPhiMeanTrend[i-1]);

    double chi2Score = buildMetric(dxyPhiMeanTrend[i],dxyPhiMeanTrend[i-1]);
    chi2Score_dxy_phi->Fill(chi2Score);

    lego_dxy_phi->SetHeader(Form("#chi^{2} score %.3f",chi2Score),"C");
    lego_dxy_phi->AddEntry(dxyPhiMeanTrend[i],Form("%i",currentList[i]),"L");
    lego_dxy_phi->AddEntry(dxyPhiMeanTrend[i-1],Form("%i",currentList[i-1]),"L");
    lego_dxy_phi->Draw("same");
      
    // dxy vs eta

    dxyEtaMeanDiff[i] = (TH1F*)(dxyEtaMeanTrend[i]->Clone(Form("%s_clone",dxyEtaMeanTrend[i]->GetName())));
    dxyEtaMeanDiff[i]->Add(dxyEtaMeanTrend[i-1],-1); 
    if(exceedsLimit(dxyEtaMeanDiff[i])){
      writeIOV=true;
      exceeds_dxy_eta=true;
    }
    dxyEtaMeanDiff[i]->GetYaxis()->SetRangeUser(-30.,30.);
    cv->cd(2);

    makeNiceTrendPlotStyle(dxyEtaMeanDiff[i],exceeds_dxy_eta);
    dxyEtaMeanDiff[i]->Draw();
    makeNewXAxis(dxyEtaMeanDiff[i]);

    theConst10->Draw("PLsame");
    theConstMinus10->Draw("PLsame");
    ptDate->Draw("same");
    lego->Draw("same");

    //------------------------------------------------------------------------------

    csup->cd(2);
    makeNiceTrendPlotStyle(dxyEtaMeanTrend[i],false);
    makeNiceTrendPlotStyle(dxyEtaMeanTrend[i-1],true);
    dxyEtaMeanTrend[i]->GetYaxis()->SetRangeUser(-30.,30.);
    dxyEtaMeanTrend[i]->Draw();
    dxyEtaMeanTrend[i-1]->Draw("same");
    makeNewXAxis(dxyEtaMeanTrend[i]);
    makeNewXAxis(dxyEtaMeanTrend[i-1]);

    chi2Score = buildMetric(dxyEtaMeanTrend[i],dxyEtaMeanTrend[i-1]);    
    chi2Score_dxy_eta->Fill(chi2Score); 	

    lego_dxy_eta->SetHeader(Form("#chi^{2} score %.3f",chi2Score),"C");
    lego_dxy_eta->AddEntry(dxyPhiMeanTrend[i],Form("%i",currentList[i]),"L");
    lego_dxy_eta->AddEntry(dxyPhiMeanTrend[i-1],Form("%i",currentList[i-1]),"L");
    lego_dxy_eta->Draw("same");

    // dz vs phi

    dzPhiMeanDiff[i] = (TH1F*)(dzPhiMeanTrend[i]->Clone(Form("%s_clone",dzPhiMeanTrend[i]->GetName())));
    dzPhiMeanDiff[i]->Add(dzPhiMeanTrend[i-1],-1); 
    if(exceedsLimit(dzPhiMeanDiff[i])){
      writeIOV=true;
      exceeds_dz_phi=true;
    }
    dzPhiMeanDiff[i]->GetYaxis()->SetRangeUser(-30.,30.);
    cv->cd(3);

    makeNiceTrendPlotStyle(dzPhiMeanDiff[i],exceeds_dz_phi);
    dzPhiMeanDiff[i]->Draw();
    makeNewXAxis(dzPhiMeanDiff[i]);

    theConst10->Draw("PLsame");
    theConstMinus10->Draw("PLsame");
    ptDate->Draw("same");
    lego->Draw("same");
    
    //------------------------------------------------------------------------------

    csup->cd(3);
    makeNiceTrendPlotStyle(dzPhiMeanTrend[i],false);
    makeNiceTrendPlotStyle(dzPhiMeanTrend[i-1],true);
    dzPhiMeanTrend[i]->GetYaxis()->SetRangeUser(-30.,30.);
    dzPhiMeanTrend[i]->Draw();
    dzPhiMeanTrend[i-1]->Draw("same");
    makeNewXAxis(dzPhiMeanTrend[i]);
    makeNewXAxis(dzPhiMeanTrend[i-1]);
    
    chi2Score = buildMetric(dzPhiMeanTrend[i],dzPhiMeanTrend[i-1]);
    chi2Score_dz_phi->Fill(chi2Score);  

    lego_dz_phi->SetHeader(Form("#chi^{2} score %.3f",chi2Score),"C");
    lego_dz_phi->AddEntry(dxyPhiMeanTrend[i],Form("%i",currentList[i]),"L");
    lego_dz_phi->AddEntry(dxyPhiMeanTrend[i-1],Form("%i",currentList[i-1]),"L");
    lego_dz_phi->Draw("same");

    // dz vs eta

    dzEtaMeanDiff[i] = (TH1F*)(dzEtaMeanTrend[i]->Clone(Form("%s_clone",dzEtaMeanTrend[i]->GetName())));
    dzEtaMeanDiff[i]->Add(dzEtaMeanTrend[i-1],-1); 
    if(exceedsLimit(dzEtaMeanDiff[i],20)){ 
      writeIOV=true;
      exceeds_dz_eta=true;
    }

    //double ksScore;
    //doKSTest(dzEtaMeanDiff[i],ksScore);

    dzEtaMeanDiff[i]->GetYaxis()->SetRangeUser(-30.,30.);
    cv->cd(4);

    makeNiceTrendPlotStyle(dzEtaMeanDiff[i],exceeds_dz_eta);
    dzEtaMeanDiff[i]->Draw();
    makeNewXAxis(dzEtaMeanDiff[i]);

    theConst20->Draw("PLsame");
    theConstMinus20->Draw("PLsame");
    ptDate->Draw("same");
    lego->Draw("same");

    //------------------------------------------------------------------------------

    csup->cd(4);
    makeNiceTrendPlotStyle(dzEtaMeanTrend[i],false);
    makeNiceTrendPlotStyle(dzEtaMeanTrend[i-1],true);
    dzEtaMeanTrend[i]->GetYaxis()->SetRangeUser(-30.,30.);
    dzEtaMeanTrend[i]->Draw();
    dzEtaMeanTrend[i-1]->Draw("same");
    makeNewXAxis(dzEtaMeanTrend[i]);
    makeNewXAxis(dzEtaMeanTrend[i-1]);
   
    chi2Score = buildMetric(dzEtaMeanTrend[i],dzEtaMeanTrend[i-1]);
    chi2Score_dz_eta->Fill(chi2Score);  

    lego_dz_eta->SetHeader(Form("#chi^{2} score %.3f",chi2Score),"C");
    lego_dz_eta->AddEntry(dxyEtaMeanTrend[i],Form("%i",currentList[i]),"L");
    lego_dz_eta->AddEntry(dxyEtaMeanTrend[i-1],Form("%i",currentList[i-1]),"L");
    lego_dz_eta->Draw("same");

    if(writeIOV) outfile<< currentList[i] << std::endl;

    cv->Print("diff.pdf");
    csup->Print("superimpose.pdf");

    if (cv) delete cv;
    if(csup) delete csup;
  }

  dummyC.Print("diff.pdf]");
  dummyC2.Print("superimpose.pdf]");

  TCanvas *scores = new TCanvas("scores","scores",1200,1200);
  scores->Divide(2,2);
  for(Int_t k=0; k<4; k++){
    scores->cd(k+1)->SetBottomMargin(0.14);
    scores->cd(k+1)->SetLeftMargin(0.18);
    scores->cd(k+1)->SetRightMargin(0.01);
    scores->cd(k+1)->SetTopMargin(0.06);
  }	

  scores->cd(1);
  makeNiceTrendPlotStyle(chi2Score_dxy_phi,false);
  chi2Score_dxy_phi->Draw(); 
  scores->cd(2);
  makeNiceTrendPlotStyle(chi2Score_dxy_eta,false);
  chi2Score_dxy_eta->Draw(); 
  scores->cd(3);
  makeNiceTrendPlotStyle(chi2Score_dz_phi,false);
  chi2Score_dz_phi->Draw(); 
  scores->cd(4);
  makeNiceTrendPlotStyle(chi2Score_dz_eta,false);
  chi2Score_dz_eta->Draw();  

  scores->SaveAs("scores.png");
  
}


/*--------------------------------------------------------------------*/
std::vector<int> list_files(const char *dirname, const char *ext)
/*--------------------------------------------------------------------*/
{
  std::vector<int> theRunNumbers;
  
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext) && fname.BeginsWith("PVValidation")) {
	//std::cout << fname.Data() << std::endl;
	TObjArray *bits = fname.Tokenize("_");
	TString theRun = bits->At(2)->GetName();
	//std::cout << theRun << std::endl;
	TString formatRun = (theRun.ReplaceAll(".root","")).ReplaceAll("_","");
	//std::cout << dirname << " "<< formatRun.Atoi() << std::endl;
	theRunNumbers.push_back(formatRun.Atoi());
      }
    }
  }
  return theRunNumbers;
}
