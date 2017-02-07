#include "TFile.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TObjArray.h"
#include "TList.h"
#include "TMath.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>

// forward declarations
void arrangeOutCanvas(TCanvas *canv,
		      TH1F* m_11Trend[100],
		      TH1F* m_12Trend[100],
		      TH1F* m_21Trend[100],
		      TH1F* m_22Trend[100],
		      Int_t nFiles, 
		      TString LegLabels[10],
		      unsigned int theRun);

TH1F* DrawConstant(TH1F *hist,Int_t nbins,Double_t lowedge,Double_t highedge,Int_t iter,Double_t theConst);
std::vector<int> list_files(const char *dirname=".", const char *ext=".root");
void MakeNiceTrendPlotStyle(TH1 *hist,Int_t color);
void cmsPrel(TPad* pad);
void makeNewXAxis (TH1 *h);
void setStyle();

// main function

void RunAndPlotPVValidation(TString namesandlabels){
  
  setStyle();

  TList *DirList   = new TList();
  TList *LabelList = new TList();
  
  TObjArray *nameandlabelpairs = namesandlabels.Tokenize(",");
  for (Int_t i = 0; i < nameandlabelpairs->GetEntries(); ++i) {
    TObjArray *aFileLegPair = TString(nameandlabelpairs->At(i)->GetName()).Tokenize("=");
    
    if(aFileLegPair->GetEntries() == 2) {
      DirList->Add(aFileLegPair->At(0)); 
      LabelList->Add(aFileLegPair->At(1));
    }
    else {
      std::cout << "Please give file name and legend entry in the following form:\n" 
		<< " filename1=legendentry1,filename2=legendentry2\n"; 
    }    
  }

  const Int_t nDirs_ = DirList->GetSize();
  TString LegLabels[10];  
  const char* dirs[10];

  std::vector<int> intersection;
 
  for(Int_t j=0; j < nDirs_; j++) {
    
    // Retrieve labels
    TObjString* legend = (TObjString*)LabelList->At(j);
    TObjString* dir    = (TObjString*)DirList->At(j);
    LegLabels[j] = legend->String();
    dirs[j] = (dir->String()).Data();
    cout<<"RunAndPlotPVValidation(): label["<<j<<"]"<<LegLabels[j]<<endl;
    
    std::vector<int> currentList = list_files(dirs[j]);
    std::vector<int> tempSwap;
    
    std::sort(currentList.begin(),currentList.end());

    if(j==0){
      intersection = currentList;
    }

    std::sort(intersection.begin(),intersection.end());

    std::set_intersection(currentList.begin(),currentList.end(),
			  intersection.begin(),intersection.end(),
			  std::back_inserter(tempSwap));
    
    intersection.clear();
    intersection = tempSwap;
    tempSwap.clear();
  }
  

  // loop over the runs in the intersection
  //for(unsigned int n=0; n<intersection.size();n++){
  // in case of debug, use only 2 
  for(unsigned int n=0; n<2;n++){
    std::cout << n << " "<<intersection.at(n) << std::endl;
    
    TFile *fins[nDirs_];

    TH1F* dxyPhiMeanTrend[nDirs_];  
    TH1F* dxyPhiWidthTrend[nDirs_]; 
    TH1F* dzPhiMeanTrend[nDirs_];   
    TH1F* dzPhiWidthTrend[nDirs_];  
    
    TH1F* dxyEtaMeanTrend[nDirs_];  
    TH1F* dxyEtaWidthTrend[nDirs_]; 
    TH1F* dzEtaMeanTrend[nDirs_];   
    TH1F* dzEtaWidthTrend[nDirs_];  
    
    TH1F* dxyNormPhiWidthTrend[nDirs_]; 
    TH1F* dxyNormEtaWidthTrend[nDirs_]; 
    TH1F* dzNormPhiWidthTrend[nDirs_]; 
    TH1F* dzNormEtaWidthTrend[nDirs_]; 

    TH1F* dxyNormPtWidthTrend[nDirs_];     
    TH1F* dzNormPtWidthTrend[nDirs_];      
    TH1F* dxyPtWidthTrend[nDirs_];
    TH1F* dzPtWidthTrend[nDirs_]; 

    bool areAllFilesOK = true;
    Int_t lastOpen = 0;

    // loop over the objects
    for(Int_t j=0; j < nDirs_; j++) {
      fins[j] = TFile::Open(Form("%s/PVValidation_%s_%i.root",dirs[j],dirs[j],intersection[n]));
      std::cout<< Form("%s/PVValidation_%s_%i.root",dirs[j],dirs[j],intersection[n]) 
	       << " has size: "<<fins[j]->GetSize() << " b ";
      
      // sanity check
      TH1F* h_tracks = (TH1F*)fins[j]->Get("PVValidation/EventFeatures/h_nTracks");
      Double_t numEvents = h_tracks->GetEntries();
      if(numEvents<1000){
	areAllFilesOK = false;
	lastOpen=j;
	break;
      }

      dxyPhiMeanTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_phi");
      dxyPhiWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_phi");
      dzPhiMeanTrend[j]   = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_phi");
      dzPhiWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_phi");
                              
      dxyEtaMeanTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_eta");
      dxyEtaWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_eta");
      dzEtaMeanTrend[j]   = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_eta");
      dzEtaWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_eta");
      
      dxyNormPhiWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dxy_phi");
      dxyNormEtaWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dxy_eta");
      dzNormPhiWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dz_phi");
      dzNormEtaWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dz_eta");

      dxyNormPtWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dxy_pTCentral");
      dzNormPtWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dz_pTCentral");
      dxyPtWidthTrend[j]     = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_pTCentral");
      dzPtWidthTrend[j]      = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_pTCentral");

      MakeNiceTrendPlotStyle(dxyPhiMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dxyPhiWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzPhiMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dzPhiWidthTrend[j],j);
      
      MakeNiceTrendPlotStyle(dxyEtaMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dxyEtaWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzEtaMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dzEtaWidthTrend[j],j);  

      MakeNiceTrendPlotStyle(dxyNormPhiWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dxyNormEtaWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzNormPhiWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzNormEtaWidthTrend[j],j);

      MakeNiceTrendPlotStyle(dxyNormPtWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzNormPtWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dxyPtWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzPtWidthTrend[j],j);

    }

    if(!areAllFilesOK){
       
      // do all the necessary deletions
      std::cout<<"\n====> not all files are OK"<<std::endl;

      for(int i=0;i<lastOpen;i++){
	fins[i]->Close();
      }
      continue;
    }

    std::cout<<"I am still here"<<std::endl;

    // Bias plots

    TCanvas *BiasesCanvas = new TCanvas(Form("Biases_%i",intersection.at(n)),"Biases",1200,1200);
    arrangeOutCanvas(BiasesCanvas,dxyPhiMeanTrend,dzPhiMeanTrend,dxyEtaMeanTrend,dzEtaMeanTrend,nDirs_,LegLabels,intersection.at(n));

    BiasesCanvas->SaveAs(Form("Biases_%i.pdf",intersection.at(n)));
    BiasesCanvas->SaveAs(Form("Biases_%i.png",intersection.at(n)));
   
    // Resolution plots

    TCanvas *ResolutionsCanvas = new TCanvas(Form("Resolutions_%i",intersection.at(n)),"Resolutions",1200,1200);
    arrangeOutCanvas(ResolutionsCanvas,dxyPhiWidthTrend,dzPhiWidthTrend,dxyEtaWidthTrend,dzEtaWidthTrend,nDirs_,LegLabels,intersection.at(n));

    ResolutionsCanvas->SaveAs(Form("Resolutions_%i.pdf",intersection.at(n)));
    ResolutionsCanvas->SaveAs(Form("Resolutions_%i.png",intersection.at(n)));
    
     // Pull plots

    TCanvas *PullsCanvas = new TCanvas(Form("Pulls_%i",intersection.at(n)),"Pulls",1200,1200);
    arrangeOutCanvas(PullsCanvas,dxyNormPhiWidthTrend,dzNormPhiWidthTrend,dxyNormEtaWidthTrend,dzNormEtaWidthTrend,nDirs_,LegLabels,intersection.at(n));

    PullsCanvas->SaveAs(Form("Pulls_%i.pdf",intersection.at(n)));
    PullsCanvas->SaveAs(Form("Pulls_%i.png",intersection.at(n)));

    // pT plots

    TCanvas *ResolutionsVsPt = new TCanvas(Form("ResolutionsVsPT_%i",intersection.at(n)),"ResolutionsVsPt",1200,1200);
    arrangeOutCanvas(ResolutionsVsPt,dxyPtWidthTrend,dzPtWidthTrend,dxyNormPtWidthTrend,dzNormPtWidthTrend,nDirs_,LegLabels,intersection.at(n));

    ResolutionsVsPt->SaveAs(Form("ResolutionsVsPt_%i.pdf",intersection.at(n)));
    ResolutionsVsPt->SaveAs(Form("ResolutionsVsPt_%i.png",intersection.at(n)));

    // do all the necessary deletions

    for(int i=0;i<nDirs_;i++){

      delete dxyPhiMeanTrend[i];
      delete dzPhiMeanTrend[i];
      delete dxyEtaMeanTrend[i];
      delete dzEtaMeanTrend[i];
      
      delete dxyPhiWidthTrend[i];
      delete dzPhiWidthTrend[i];
      delete dxyEtaWidthTrend[i];
      delete dzEtaWidthTrend[i];

      delete dxyNormPhiWidthTrend[i];
      delete dxyNormEtaWidthTrend[i];
      delete dzNormPhiWidthTrend[i];
      delete dzNormEtaWidthTrend[i];

      delete dxyNormPtWidthTrend[i]; 
      delete dzNormPtWidthTrend[i];  
      delete dxyPtWidthTrend[i];
      delete dzPtWidthTrend[i]; 

      fins[i]->Close();
    }
    
    std::cout<<std::endl;
  }   
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

//*************************************************************
void arrangeOutCanvas(TCanvas *canv, TH1F* m_11Trend[100],TH1F* m_12Trend[100],TH1F* m_21Trend[100],TH1F* m_22Trend[100],Int_t nDirs, TString LegLabels[10],unsigned int theRun){
//*************************************************************

  TLegend *lego = new TLegend(0.19,0.80,0.79,0.93);
  //lego-> SetNColumns(2);
  lego->SetFillColor(10);
  lego->SetTextSize(0.042);
  lego->SetTextFont(42);
  lego->SetFillColor(10);
  lego->SetLineColor(10);
  lego->SetShadowColor(10);
  
  TPaveText *ptDate =new TPaveText(0.19,0.95,0.45,0.99,"blNDC");
  ptDate->SetFillColor(kYellow);
  //ptDate->SetFillColor(10);
  ptDate->SetBorderSize(1);
  ptDate->SetLineColor(kBlue);
  ptDate->SetLineWidth(1);
  ptDate->SetTextFont(42);
  TText *textDate = ptDate->AddText(Form("Run: %i",theRun));
  textDate->SetTextSize(0.04);
  textDate->SetTextColor(kBlue);
  textDate->SetTextAlign(22);

  canv->SetFillColor(10);  
  canv->Divide(2,2);
 
  TH1F *dBiasTrend[4][nDirs]; 
  
  for(Int_t i=0;i<nDirs;i++){
    dBiasTrend[0][i] = m_11Trend[i];
    dBiasTrend[1][i] = m_12Trend[i];
    dBiasTrend[2][i] = m_21Trend[i];
    dBiasTrend[3][i] = m_22Trend[i];
  }

  Double_t absmin[4]={999.,999.,999.,999.};
  Double_t absmax[4]={-999.,-999.-999.,-999.};

  for(Int_t k=0; k<4; k++){

    canv->cd(k+1)->SetBottomMargin(0.14);
    canv->cd(k+1)->SetLeftMargin(0.18);
    canv->cd(k+1)->SetRightMargin(0.01);
    canv->cd(k+1)->SetTopMargin(0.06);
    
    canv->cd(k+1);
    
    for(Int_t i=0; i<nDirs; i++){
      if(dBiasTrend[k][i]->GetMaximum()>absmax[k]) absmax[k] = dBiasTrend[k][i]->GetMaximum();
      if(dBiasTrend[k][i]->GetMinimum()<absmin[k]) absmin[k] = dBiasTrend[k][i]->GetMinimum();
    }

    Double_t safeDelta=(absmax[k]-absmin[k])/8.;
    Double_t theExtreme=std::max(absmax[k],TMath::Abs(absmin[k]));

    for(Int_t i=0; i<nDirs; i++){
      if(i==0){

	TString theTitle = dBiasTrend[k][i]->GetName();
	Int_t nbins =  dBiasTrend[k][i]->GetNbinsX();
	Double_t lowedge  = dBiasTrend[k][i]->GetBinLowEdge(1);
	Double_t highedge = dBiasTrend[k][i]->GetBinLowEdge(nbins+1);
	
	if(theTitle.Contains("norm")){
	  //dBiasTrend[k][i]->GetYaxis()->SetRangeUser(std::min(-0.48,absmin[k]-safeDelta/2.),std::max(0.48,absmax[k]+safeDelta/2.));
	  dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,1.8);
	} else {
	  if(!theTitle.Contains("width")){
	    //dBiasTrend[k][i]->GetYaxis()->SetRangeUser(-theExtreme-(safeDelta/2.),theExtreme+(safeDelta/2.));
	    dBiasTrend[k][i]->GetYaxis()->SetRangeUser(-40.,40.);
	  } else {
	    // dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,theExtreme+(safeDelta/2.));
	    dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,300.);
	  }
	} 
	dBiasTrend[k][i]->Draw("Le1");
	makeNewXAxis(dBiasTrend[k][i]);
      
	Double_t theC = -1.;
	
	if(theTitle.Contains("width")){ 
	  if(theTitle.Contains("norm") ){
	    theC = 1.;
	  } else {
	    theC = -1.;
	  }
	} else {
	  theC = 0.;
	}
	
	TH1F* theConst = DrawConstant(dBiasTrend[k][i],nbins,lowedge,highedge,1,theC);
	theConst->Draw("PLsame");

      } else { 
	dBiasTrend[k][i]->Draw("Le1sames");
	makeNewXAxis(dBiasTrend[k][i]);
      }
      TPad *current_pad = static_cast<TPad*>(canv->GetPad(k+1));
      cmsPrel(current_pad);
      ptDate->Draw("same");

      if(k==0){
	lego->AddEntry(dBiasTrend[k][i],LegLabels[i]);
      } 
    }  
  
    lego->Draw();
  } 
}

/*--------------------------------------------------------------------*/
void  MakeNiceTrendPlotStyle(TH1 *hist,Int_t color)
/*--------------------------------------------------------------------*/
{ 

  Int_t markers[8] = {kFullSquare,kFullCircle,kOpenSquare,kOpenCircle,kFullTriangleDown,kFullTriangleUp,kOpenTriangleDown,kOpenTriangleUp};
  Int_t colors[8]  = {kBlack,kBlue,kRed,kGreen+2,kOrange,kMagenta,kCyan,kViolet};
  
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
  //hist->GetXaxis()->SetNdivisions(505);
  if(color!=8){
    hist->SetMarkerSize(1.5);
  } else {
    hist->SetLineWidth(3);
    hist->SetMarkerSize(0.0);    
  }
  hist->SetMarkerStyle(markers[color]);
  hist->SetLineColor(colors[color]);
  hist->SetMarkerColor(colors[color]);
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
    axmin = -2.5;
    axmax = 2.5;
    ndiv = 505;
  } else if (myTitle.Contains("phi")){
    axmin = -TMath::Pi();
    axmax = TMath::Pi();
    ndiv = 510;
  } else if (myTitle.Contains("pT")) {
    axmin = 0;
    axmax = 19.99;
    ndiv = 510;
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
void setStyle(){
/*--------------------------------------------------------------------*/
  TH1::StatOverflows(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat("e");
  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadBottomMargin(0.15);
  //gStyle->SetPadLeftMargin(0.17);
  //gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBorderMode(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.05);///---> gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(510);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

/*--------------------------------------------------------------------*/
void cmsPrel(TPad* pad) {
/*--------------------------------------------------------------------*/
  
  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  float relPosX = 0.015;
  float relPosY = 0.045;
  float lumiTextOffset = 0.8;

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42); //22

  float posX_ =  1-r - relPosX*(1-l-r);
  float posY_ =  1-t + 0.05; /// - relPosY*(1-t-b);

  latex->SetTextAlign(33);
  latex->DrawLatex(posX_,posY_,"CMS Preliminary (13 TeV)");
  
}

/*--------------------------------------------------------------------*/
TH1F* DrawConstant(TH1F *hist,Int_t nbins,Double_t lowedge,Double_t highedge,Int_t iter,Double_t theConst)
/*--------------------------------------------------------------------*/
{ 

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
