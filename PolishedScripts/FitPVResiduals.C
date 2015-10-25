#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include "TFile.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TCut.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TString.h"
#include "TMath.h"
#include <TDatime.h>
#include <TSpectrum.h>
#include <TSystem.h>

void arrangeCanvas(TCanvas *canv,TH1F* meanplots[100],TH1F* widthplots[100],Int_t nFiles,TString LegLabels[10],TString theLeg);
void arrangeFitCanvas(TCanvas *canv,TH1F* meanplots[100],Int_t nFiles, TString LegLabels[10],TString theLeg);

std::pair<Double_t,Double_t> getMedian(TH1F *histo);
std::pair<Double_t,Double_t> getMAD(TH1F *histo);

std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > fitResiduals(TH1 *hist);

Double_t DoubleSidedCB(double* x, double* par);
std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > fitResidualsCB(TH1 *hist);

Double_t tp0Fit( Double_t *x, Double_t *par5 );
std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > fitStudentTResiduals(TH1 *hist);

void FillTrendPlot(TH1F* trendPlot, TH1F* residualsPlot[100], TString fitPar_, TString var_,Int_t nbins);
void MakeNiceTrendPlotStyle(TH1 *hist,Int_t color);
void MakeNiceTF1Style(TF1 *f1,Int_t color);

void FitPVResiduals(TString namesandlabels,TString theLeg);

void setStyle();

// ancillary fitting functions
Double_t fULine(Double_t *x, Double_t *par);
Double_t fDLine(Double_t *x, Double_t *par);
void FitULine(TH1 *hist);
void FitDLine(TH1 *hist);

// global variables

ofstream outfile("FittedDeltaZ.txt");
Int_t my_colors[10]={kBlack,kRed,kBlue,kMagenta,kBlack,kRed,kBlue,kGreen};

const Int_t nBins_ = 24;
Float_t _boundMin   = -0.5;
Float_t _boundSx    = (nBins_/4.)-0.5;
Float_t _boundDx    = 3*(nBins_/4.)-0.5;
Float_t _boundMax   = nBins_-0.5;

//*************************************************************
void FitPVResiduals(TString namesandlabels,TString theLeg){
//*************************************************************
 
  setStyle();

  TList *FileList  = new TList();
  TList *LabelList = new TList();
  
  TObjArray *nameandlabelpairs = namesandlabels.Tokenize(",");
  for (Int_t i = 0; i < nameandlabelpairs->GetEntries(); ++i) {
    TObjArray *aFileLegPair = TString(nameandlabelpairs->At(i)->GetName()).Tokenize("=");
    
    if(aFileLegPair->GetEntries() == 2) {
      FileList->Add( TFile::Open(aFileLegPair->At(0)->GetName())  );  // 2
      LabelList->Add( aFileLegPair->At(1) );
    }
    else {
      std::cout << "Please give file name and legend entry in the following form:\n" 
		<< " filename1=legendentry1,filename2=legendentry2\n";
      
    }    
  }

  const Int_t nFiles_ = FileList->GetSize();
  TString LegLabels[10];  
  TFile *fins[nFiles_]; 

  for(Int_t j=0; j < nFiles_; j++) {
    
    // Retrieve files
    fins[j] = (TFile*)FileList->At(j);    
 
    // Retrieve labels
    TObjString* legend = (TObjString*)LabelList->At(j);
    LegLabels[j] = legend->String();
    LegLabels[j].ReplaceAll("_"," ");
    cout<<"FitPVResiduals::FitPVResiduals(): label["<<j<<"]"<<LegLabels[j]<<endl;
    
  }
  
  TString append_ = "run"+theLeg+"_"+LegLabels[0];
  if(nFiles_>1){
    for(Int_t i=1;i<nFiles_;i++){
      append_+=("_vs_"+LegLabels[i]);
    }
  }

  // absolute residuals
  TH1F* dxyPhiResiduals[nFiles_][nBins_];
  TH1F* dxyEtaResiduals[nFiles_][nBins_];
  				        
  TH1F* dzPhiResiduals[nFiles_][nBins_];
  TH1F* dzEtaResiduals[nFiles_][nBins_];

  // normalized residuals
  TH1F* dxyNormPhiResiduals[nFiles_][nBins_];
  TH1F* dxyNormEtaResiduals[nFiles_][nBins_];
  				        
  TH1F* dzNormPhiResiduals[nFiles_][nBins_];
  TH1F* dzNormEtaResiduals[nFiles_][nBins_];
  
  for(Int_t i=0;i<nFiles_;i++){
    for(Int_t j=0;j<nBins_;j++){
      
      // absolute residuals

      //std::cout<<Form("PVValidation/Abs_Transv_Phi_Residuals/histo_dxy_phi_plot%i",j)<<std::endl;   
      dxyPhiResiduals[i][j] = (TH1F*)fins[i]->Get(Form("PVValidation/Abs_Transv_Phi_Residuals/histo_dxy_phi_plot%i",j));
      //std::cout<<Form("PVValidation/Abs_Transv_Ets_Residuals/histo_dxy_eta_plot%i",j)<<std::endl;
      dxyEtaResiduals[i][j] = (TH1F*)fins[i]->Get(Form("PVValidation/Abs_Transv_Eta_Residuals/histo_dxy_eta_plot%i",j));
      //std::cout<<Form("PVValidation/Abs_Long_Phi_Residuals/histo_dz_phi_plot%i",j)<<std::endl;
      dzPhiResiduals[i][j]  = (TH1F*)fins[i]->Get(Form("PVValidation/Abs_Long_Phi_Residuals/histo_dz_phi_plot%i",j));
      //std::cout<<Form("PVValidation/Abs_Long_Eta_Residuals/histo_dz_eta_plot%i",j)<<std::endl;
      dzEtaResiduals[i][j]  = (TH1F*)fins[i]->Get(Form("PVValidation/Abs_Long_Eta_Residuals/histo_dz_eta_plot%i",j));
      
      // normalized residuals
      
      //std::cout<<Form("PVValidation/Norm_Transv_Phi_Residuals/histo_dxy_phi_plot%i",j)<<std::endl;   
      dxyNormPhiResiduals[i][j] = (TH1F*)fins[i]->Get(Form("PVValidation/Norm_Transv_Phi_Residuals/histo_norm_dxy_phi_plot%i",j));
      //std::cout<<Form("PVValidation/Norm_Transv_Ets_Residuals/histo_dxy_eta_plot%i",j)<<std::endl;
      dxyNormEtaResiduals[i][j] = (TH1F*)fins[i]->Get(Form("PVValidation/Norm_Transv_Eta_Residuals/histo_norm_dxy_eta_plot%i",j));
      //std::cout<<Form("PVValidation/Norm_Long_Phi_Residuals/histo_dz_phi_plot%i",j)<<std::endl;
      dzNormPhiResiduals[i][j]  = (TH1F*)fins[i]->Get(Form("PVValidation/Norm_Long_Phi_Residuals/histo_norm_dz_phi_plot%i",j));
      //std::cout<<Form("PVValidation/Norm_Long_Eta_Residuals/histo_dz_eta_plot%i",j)<<std::endl;
      dzNormEtaResiduals[i][j]  = (TH1F*)fins[i]->Get(Form("PVValidation/Norm_Long_Eta_Residuals/histo_norm_dz_eta_plot%i",j));

    }
  }
 
  Double_t highedge=nBins_-0.5;
  Double_t lowedge=-0.5;
  
  // absoulute

  TH1F* dxyPhiMeanTrend[nFiles_];  
  TH1F* dxyPhiWidthTrend[nFiles_]; 
  TH1F* dzPhiMeanTrend[nFiles_];   
  TH1F* dzPhiWidthTrend[nFiles_];  
  		       	 
  TH1F* dxyEtaMeanTrend[nFiles_];  
  TH1F* dxyEtaWidthTrend[nFiles_]; 
  TH1F* dzEtaMeanTrend[nFiles_];   
  TH1F* dzEtaWidthTrend[nFiles_];  

  // normalized

  TH1F* dxyNormPhiMeanTrend[nFiles_];  
  TH1F* dxyNormPhiWidthTrend[nFiles_]; 
  TH1F* dzNormPhiMeanTrend[nFiles_];   
  TH1F* dzNormPhiWidthTrend[nFiles_];  
  		       	 
  TH1F* dxyNormEtaMeanTrend[nFiles_];  
  TH1F* dxyNormEtaWidthTrend[nFiles_]; 
  TH1F* dzNormEtaMeanTrend[nFiles_];   
  TH1F* dzNormEtaWidthTrend[nFiles_];  

  for(Int_t i=0;i<nFiles_;i++){

    dxyPhiMeanTrend[i]  = new TH1F(Form("means_dxy_phi_%i",i),"#LT d_{xy} #GT vs #phi sector;#varphi (sector) [degrees];#LT d_{xy} #GT [#mum]",nBins_,lowedge,highedge); 
    dxyPhiWidthTrend[i] = new TH1F(Form("widths_dxy_phi_%i",i),"#sigma(d_{xy}) vs #phi sector;#varphi (sector) [degrees];#sigma(d_{xy}) [#mum]",nBins_,lowedge,highedge);
    dzPhiMeanTrend[i]   = new TH1F(Form("means_dz_phi_%i",i),"#LT d_{z} #GT vs #phi sector;#varphi (sector) [degrees];#LT d_{z} #GT [#mum]",nBins_,lowedge,highedge); 
    dzPhiWidthTrend[i]  = new TH1F(Form("widths_dz_phi_%i",i),"#sigma(d_{z}) vs #phi sector;#varphi (sector) [degrees];#sigma(d_{z}) [#mum]",nBins_,lowedge,highedge);
    
    dxyEtaMeanTrend[i]  = new TH1F(Form("means_dxy_eta_%i",i),"#LT d_{xy} #GT vs #eta sector;#eta (sector);#LT d_{xy} #GT [#mum]",nBins_,lowedge,highedge);
    dxyEtaWidthTrend[i] = new TH1F(Form("widths_dxy_eta_%i",i),"#sigma(d_{xy}) vs #eta sector;#eta (sector);#sigma(d_{xy}) [#mum]",nBins_,lowedge,highedge);
    dzEtaMeanTrend[i]   = new TH1F(Form("means_dz_eta_%i",i),"#LT d_{z} #GT vs #eta sector;#eta (sector);#LT d_{z} #GT [#mum]",nBins_,lowedge,highedge); 
    dzEtaWidthTrend[i]  = new TH1F(Form("widths_dz_eta_%i",i),"#sigma(d_{xy}) vs #eta sector;#eta (sector);#sigma(d_{z}) [#mum]",nBins_,lowedge,highedge);

    dxyNormPhiMeanTrend[i] = new TH1F(Form("means_dxyNorm_phi_%i",i),"#LT d_{xy}/#sigma_{d_{xy}} #GT vs #phi sector;#varphi (sector) [degrees];#LT d_{xy}/#sigma_{d_{xy}} #GT [#mum]",nBins_,lowedge,highedge); 
    dxyNormPhiWidthTrend[i]= new TH1F(Form("widths_dxyNorm_phi_%i",i),"#sigma(d_{xy}/#sigma_{d_{xy}}) vs #phi sector;#varphi (sector) [degrees];#sigma(d_{xy}/#sigma_{d_{xy}}) [#mum]",nBins_,lowedge,highedge);
    dzNormPhiMeanTrend[i]  = new TH1F(Form("means_dzNorm_phi_%i",i),"#LT d_{z}/#sigma_{d_{z}} #GT vs #phi sector;#varphi (sector) [degrees];#LT d_{z}/#sigma_{d_{z}} #GT [#mum]",nBins_,lowedge,highedge); 
    dzNormPhiWidthTrend[i] = new TH1F(Form("widths_dzNorm_phi_%i",i),"#sigma(d_{z}/#sigma_{d_{z}}) vs #phi sector;#varphi (sector) [degrees];#sigma(d_{z}/#sigma_{d_{z}}) [#mum]",nBins_,lowedge,highedge);
    
    dxyNormEtaMeanTrend[i] = new TH1F(Form("means_dxyNorm_eta_%i",i),"#LT d_{xy}/#sigma_{d_{xy}} #GT vs #eta sector;#eta (sector);#LT d_{xy}/#sigma_{d_{xy}} #GT [#mum]",nBins_,lowedge,highedge);
    dxyNormEtaWidthTrend[i]= new TH1F(Form("widths_dxyNorm_eta_%i",i),"#sigma(d_{xy}/#sigma_{d_{xy}}) vs #eta sector;#eta (sector);#sigma(d_{xy}/#sigma_{d_{xy}}) [#mum]",nBins_,lowedge,highedge);
    dzNormEtaMeanTrend[i]  = new TH1F(Form("means_dzNorm_eta_%i",i),"#LT d_{z}/#sigma_{d_{z}} #GT vs #eta sector;#eta (sector);#LT d_{z}/#sigma_{d_{z}} #GT [#mum]",nBins_,lowedge,highedge); 
    dzNormEtaWidthTrend[i] = new TH1F(Form("widths_dzNorm_eta_%i",i),"#sigma(d_{z}/#sigma_{d_{z}}) vs #eta sector;#eta (sector);#sigma(d_{z}/#sigma_{d_{z}}) [#mum]",nBins_,lowedge,highedge);

    // dxyPhiMeanBiasTrend  = new TH1F("norm_means_dxy_phi","#LT d_{xy}/#sigma_{d_{xy}} #GT vs #phi sector;#varphi (sector) [degrees];#LT d_{xy}/#sigma_{d_{xy}} #GT",nBins_,lowedge,highedge);
    // dxyPhiWidthBiasTrend = new TH1F("norm_widths_dxy_phi","width(d_{xy}/#sigma_{d_{xy}}) vs #phi sector;#varphi (sector) [degrees]; width(d_{xy}/#sigma_{d_{xy}})",nBins_,lowedge,highedge);
    // dzPhiMeanBiasTrend   = new TH1F("norm_means_dz_phi","#LT d_{z}/#sigma_{d_{z}} #GT vs #phi sector;#varphi (sector) [degrees];#LT d_{z}/#sigma_{d_{z}} #GT",nBins_,lowedge,highedge); 
    // dzPhiWidthBiasTrend  = new TH1F("norm_widths_dz_phi","width(d_{z}/#sigma_{d_{z}}) vs #phi sector;#varphi (sector) [degrees];width(d_{z}/#sigma_{d_{z}})",nBins_,lowedge,highedge);
    			             
    // dxyEtaMeanBiasTrend  = new TH1F("norm_means_dxy_eta","#LT d_{xy}/#sigma_{d_{xy}} #GT vs #eta sector;#eta (sector);#LT d_{xy}/#sigma_{d_{z}} #GT",nBins_,lowedge,highedge);
    // dxyEtaWidthBiasTrend = new TH1F("norm_widths_dxy_eta","width(d_{xy}/#sigma_{d_{xy}}) vs #eta sector;#eta (sector);width(d_{xy}/#sigma_{d_{z}})",nBins_,lowedge,highedge);
    // dzEtaMeanBiasTrend   = new TH1F("norm_means_dz_eta","#LT d_{z}/#sigma_{d_{z}} #GT vs #eta sector;#eta (sector);#LT d_{z}/#sigma_{d_{z}} #GT",nBins_,lowedge,highedge);  
    // dzEtaWidthBiasTrend  = new TH1F("norm_widths_dz_eta","width(d_{z}/#sigma_{d_{z}}) vs #eta sector;#eta (sector);width(d_{z}/#sigma_{d_{z}})",nBins_,lowedge,highedge);    
    
    // absolute

    FillTrendPlot(dxyPhiMeanTrend[i] ,dxyPhiResiduals[i],"mean","phi",nBins_);  
    FillTrendPlot(dxyPhiWidthTrend[i],dxyPhiResiduals[i],"width","phi",nBins_);
    FillTrendPlot(dzPhiMeanTrend[i]  ,dzPhiResiduals[i] ,"mean","phi",nBins_);   
    FillTrendPlot(dzPhiWidthTrend[i] ,dzPhiResiduals[i] ,"width","phi",nBins_);  
    
    FillTrendPlot(dxyEtaMeanTrend[i] ,dxyEtaResiduals[i],"mean","eta",nBins_); 
    FillTrendPlot(dxyEtaWidthTrend[i],dxyEtaResiduals[i],"width","eta",nBins_);
    FillTrendPlot(dzEtaMeanTrend[i]  ,dzEtaResiduals[i] ,"mean","eta",nBins_); 
    FillTrendPlot(dzEtaWidthTrend[i] ,dzEtaResiduals[i] ,"width","eta",nBins_);

    MakeNiceTrendPlotStyle(dxyPhiMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dxyPhiWidthTrend[i],i);
    MakeNiceTrendPlotStyle(dzPhiMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dzPhiWidthTrend[i],i);
  
    MakeNiceTrendPlotStyle(dxyEtaMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dxyEtaWidthTrend[i],i);
    MakeNiceTrendPlotStyle(dzEtaMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dzEtaWidthTrend[i],i);
    
    // normalized

    FillTrendPlot(dxyNormPhiMeanTrend[i] ,dxyNormPhiResiduals[i],"mean","phi",nBins_);  
    FillTrendPlot(dxyNormPhiWidthTrend[i],dxyNormPhiResiduals[i],"width","phi",nBins_);
    FillTrendPlot(dzNormPhiMeanTrend[i]  ,dzNormPhiResiduals[i] ,"mean","phi",nBins_);   
    FillTrendPlot(dzNormPhiWidthTrend[i] ,dzNormPhiResiduals[i] ,"width","phi",nBins_);  
    
    FillTrendPlot(dxyNormEtaMeanTrend[i] ,dxyNormEtaResiduals[i],"mean","eta",nBins_); 
    FillTrendPlot(dxyNormEtaWidthTrend[i],dxyNormEtaResiduals[i],"width","eta",nBins_);
    FillTrendPlot(dzNormEtaMeanTrend[i]  ,dzNormEtaResiduals[i] ,"mean","eta",nBins_); 
    FillTrendPlot(dzNormEtaWidthTrend[i] ,dzNormEtaResiduals[i] ,"width","eta",nBins_);

    MakeNiceTrendPlotStyle(dxyNormPhiMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dxyNormPhiWidthTrend[i],i);
    MakeNiceTrendPlotStyle(dzNormPhiMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dzNormPhiWidthTrend[i],i);
  
    MakeNiceTrendPlotStyle(dxyNormEtaMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dxyNormEtaWidthTrend[i],i);
    MakeNiceTrendPlotStyle(dzNormEtaMeanTrend[i],i);
    MakeNiceTrendPlotStyle(dzNormEtaWidthTrend[i],i);

  }
  
  // absolute

  TCanvas *dxyPhiTrend = new TCanvas("dxyPhiTrend","dxyPhiTrend",1200,600);
  arrangeCanvas(dxyPhiTrend,dxyPhiMeanTrend,dxyPhiWidthTrend,nFiles_,LegLabels,theLeg);

  dxyPhiTrend->SaveAs("dxyPhiTrend_"+append_+".pdf");
  dxyPhiTrend->SaveAs("dxyPhiTrend_"+append_+".png");

  TCanvas *dzPhiTrend = new TCanvas("dzPhiTrend","dzPhiTrend",1200,600);
  arrangeCanvas(dzPhiTrend,dzPhiMeanTrend,dzPhiWidthTrend,nFiles_,LegLabels,theLeg);

  dzPhiTrend->SaveAs("dzPhiTrend_"+append_+".pdf");
  dzPhiTrend->SaveAs("dzPhiTrend_"+append_+".png");

  TCanvas *dxyEtaTrend = new TCanvas("dxyEtaTrend","dxyEtaTrend",1200,600);
  arrangeCanvas(dxyEtaTrend,dxyEtaMeanTrend,dxyEtaWidthTrend,nFiles_,LegLabels,theLeg);

  dxyEtaTrend->SaveAs("dxyEtaTrend_"+append_+".pdf");
  dxyEtaTrend->SaveAs("dxyEtaTrend_"+append_+".png");

  TCanvas *dzEtaTrend = new TCanvas("dzEtaTrend","dzEtaTrend",1200,600);
  arrangeCanvas(dzEtaTrend,dzEtaMeanTrend,dzEtaWidthTrend,nFiles_,LegLabels,theLeg);

  dzEtaTrend->SaveAs("dzEtaTrend_"+append_+".pdf");
  dzEtaTrend->SaveAs("dzEtaTrend_"+append_+".png");

  // fit dz vs phi
  TCanvas *dzPhiTrendFit = new TCanvas("dzPhiTrendFit","dzPhiTrendFit",1200,600);
  arrangeFitCanvas(dzPhiTrendFit,dzPhiMeanTrend,nFiles_,LegLabels,theLeg);

  dzPhiTrendFit->SaveAs("dzPhiTrendFit_"+append_+".pdf");
  dzPhiTrendFit->SaveAs("dzPhiTrendFit_"+append_+".png");

  // normalized

  TCanvas *dxyNormPhiTrend = new TCanvas("dxyNormPhiTrend","dxyNormPhiTrend",1200,600);
  arrangeCanvas(dxyNormPhiTrend,dxyNormPhiMeanTrend,dxyNormPhiWidthTrend,nFiles_,LegLabels,theLeg);

  dxyNormPhiTrend->SaveAs("dxyNormPhiTrend_"+append_+".pdf");
  dxyNormPhiTrend->SaveAs("dxyNormPhiTrend_"+append_+".png");

  TCanvas *dzNormPhiTrend = new TCanvas("dzNormPhiTrend","dzNormPhiTrend",1200,600);
  arrangeCanvas(dzNormPhiTrend,dzNormPhiMeanTrend,dzNormPhiWidthTrend,nFiles_,LegLabels,theLeg);

  dzNormPhiTrend->SaveAs("dzNormPhiTrend_"+append_+".pdf");
  dzNormPhiTrend->SaveAs("dzNormPhiTrend_"+append_+".png");

  TCanvas *dxyNormEtaTrend = new TCanvas("dxyNormEtaTrend","dxyNormEtaTrend",1200,600);
  arrangeCanvas(dxyNormEtaTrend,dxyNormEtaMeanTrend,dxyNormEtaWidthTrend,nFiles_,LegLabels,theLeg);

  dxyNormEtaTrend->SaveAs("dxyNormEtaTrend_"+append_+".pdf");
  dxyNormEtaTrend->SaveAs("dxyNormEtaTrend_"+append_+".png");

  TCanvas *dzNormEtaTrend = new TCanvas("dzNormEtaTrend","dzNormEtaTrend",1200,600);
  arrangeCanvas(dzNormEtaTrend,dzNormEtaMeanTrend,dzNormEtaWidthTrend,nFiles_,LegLabels,theLeg);

  dzNormEtaTrend->SaveAs("dzNormEtaTrend_"+append_+".pdf");
  dzNormEtaTrend->SaveAs("dzNormEtaTrend_"+append_+".png");

}

//*************************************************************
void arrangeCanvas(TCanvas *canv,TH1F* meanplots[100],TH1F* widthplots[100],Int_t nFiles, TString LegLabels[10],TString theLeg){
//*************************************************************

  TLegend *lego = new TLegend(0.18,0.75,0.58,0.92);
  lego->SetFillColor(10);
  lego->SetTextSize(0.05);
  lego->SetTextFont(42);
  lego->SetFillColor(10);
  lego->SetLineColor(10);
  lego->SetShadowColor(10);

  TPaveText *pt = new TPaveText(0.148,0.95,0.89,0.97,"NDC");
  pt->SetFillColor(10);
  pt->SetTextColor(1);
  pt->SetTextFont(42);
  pt->SetTextAlign(11);
  TText *text1 = pt->AddText("CMS preliminary 2015    p-p data, #sqrt{s}=8 TeV");
  text1->SetTextSize(0.04);
 
  TPaveText *pt2 = new TPaveText(0.73,0.77,0.88,0.92,"NDC");
  pt2->SetFillColor(10);
  pt2->SetTextColor(1);
  //pt->SetTextSize(0.05);
  pt2->SetTextFont(42);
  pt2->SetTextAlign(11);
  TText *text2 = pt2->AddText("run: "+theLeg);
  text2->SetTextSize(0.04); 

  canv->SetFillColor(10);  
  canv->Divide(2,1);
 
  canv->cd(1)->SetBottomMargin(0.12);
  canv->cd(1)->SetLeftMargin(0.17);
  canv->cd(1)->SetRightMargin(0.02);
  canv->cd(1)->SetTopMargin(0.06);  

  canv->cd(2)->SetBottomMargin(0.12);
  canv->cd(2)->SetLeftMargin(0.17);
  canv->cd(2)->SetRightMargin(0.02);
  canv->cd(2)->SetTopMargin(0.06);  
  
  canv->cd(1);
  Double_t absmin(999.);
  Double_t absmax(-999.);

  for(Int_t i=0; i<nFiles; i++){
    if(meanplots[i]->GetMaximum()>absmax) absmax = meanplots[i]->GetMaximum();
    if(meanplots[i]->GetMinimum()<absmin) absmin = meanplots[i]->GetMinimum();
  }

  Double_t safeDelta=(absmax-absmin)/2.;

  for(Int_t i=0; i<nFiles; i++){
    if(i==0){
      meanplots[i]->GetYaxis()->SetRangeUser(absmin-safeDelta/2.,absmax+safeDelta);
      meanplots[i]->Draw("e1");
    }
    else meanplots[i]->Draw("e1sames");
    lego->AddEntry(meanplots[i],LegLabels[i]); 
  }  
  
  lego->Draw();
  pt->Draw("same");
  pt2->Draw("same");

  canv->cd(2);
  Double_t absmax2(-999.);

  for(Int_t i=0; i<nFiles; i++){
    if(widthplots[i]->GetMaximum()>absmax2) absmax2 = widthplots[i]->GetMaximum();
  }

  Double_t safeDelta2=absmax2/3.;

  for(Int_t i=0; i<nFiles; i++){
     if(i==0) widthplots[i]->Draw("e1");
     else widthplots[i]->Draw("e1sames");
     widthplots[i]->SetMinimum(0.5);
     widthplots[i]->SetMaximum(absmax2+safeDelta2);
  }
  
  lego->Draw();
  pt->Draw("same");
  pt2->Draw("same");

}

//*************************************************************
void arrangeFitCanvas(TCanvas *canv,TH1F* meanplots[100],Int_t nFiles, TString LegLabels[10],TString theLeg)
//*************************************************************
{
  canv->SetBottomMargin(0.12);
  canv->SetLeftMargin(0.1);
  canv->SetRightMargin(0.02);
  canv->SetTopMargin(0.06);  

  TLegend *lego = new TLegend(0.12,0.73,0.52,0.90);
  lego->SetFillColor(10);
  lego->SetTextSize(0.05);
  lego->SetTextFont(42);
  lego->SetFillColor(10);
  lego->SetLineColor(10);
  lego->SetShadowColor(10);

  TPaveText *pt = new TPaveText(0.148,0.95,0.89,0.97,"NDC");
  pt->SetFillColor(10);
  pt->SetTextColor(1);
  pt->SetTextFont(42);
  pt->SetTextAlign(11);
  TText *text1 = pt->AddText("CMS preliminary 2015    p-p data, #sqrt{s}=8 TeV");
  text1->SetTextSize(0.04);

  TPaveText *pt2 = new TPaveText(0.83,0.75,0.94,0.90,"NDC");
  pt2->SetFillColor(10);
  pt2->SetTextColor(1);
  //pt->SetTextSize(0.05);
  pt2->SetTextFont(42);
  pt2->SetTextAlign(11);
  TText *text2 = pt2->AddText("run: "+theLeg);
  text2->SetTextSize(0.04); 


  TF1 *fleft[nFiles]; 
  TF1 *fright[nFiles];
  TF1 *fall[nFiles];  

  TF1 *FitDzUp[nFiles];
  TF1 *FitDzDown[nFiles];

  for(Int_t j=0;j<nFiles;j++){
    
    std::cout<<"my_colors["<<j<<"]="<<my_colors[j]<<std::endl;

    Double_t deltaZ(0);
    Double_t sigmadeltaZ(-1);

    TCanvas *theNewCanvas2 = new TCanvas("NewCanvas2","Fitting Canvas 2",800,600);
    theNewCanvas2->Divide(2,1);

    TH1F *hnewUp   = (TH1F*)meanplots[j]->Clone("hnewUp");
    TH1F *hnewDown = (TH1F*)meanplots[j]->Clone("hnewDown");
    
    fleft[j]  = new TF1(Form("fleft_%i",j),fULine,_boundMin,_boundSx,1);
    fright[j] = new TF1(Form("fright_%i",j),fULine,_boundDx,_boundMax,1);
    fall[j]   = new TF1(Form("fall_%i",j),fDLine,_boundSx,_boundDx,1);
    
    FitULine(hnewUp);  
    FitDzUp[j]   = (TF1*)hnewUp->GetListOfFunctions()->FindObject("lineUp"); 
    if(FitDzUp[j]){
      fleft[j]->SetParameters(FitDzUp[j]->GetParameters());
      fleft[j]->SetParErrors(FitDzUp[j]->GetParErrors());
      hnewUp->GetListOfFunctions()->Add(fleft[j]);
      fright[j]->SetParameters(FitDzUp[j]->GetParameters());
      fright[j]->SetParErrors(FitDzUp[j]->GetParErrors());
      hnewUp->GetListOfFunctions()->Add(fright[j]);
      FitDzUp[j]->Delete();

      theNewCanvas2->cd(1);
      MakeNiceTF1Style(fright[j],my_colors[j]);
      MakeNiceTF1Style(fleft[j],my_colors[j]);
      fright[j]->Draw("same");
      fleft[j]->Draw("same");
    }
    
    FitDLine(hnewDown);  
    FitDzDown[j] = (TF1*)hnewDown->GetListOfFunctions()->FindObject("lineDown");    
    
    if(FitDzDown[j]){
      fall[j]->SetParameters(FitDzDown[j]->GetParameters());
      fall[j]->SetParErrors(FitDzDown[j]->GetParErrors());
      hnewDown->GetListOfFunctions()->Add(fall[j]);
      FitDzDown[j]->Delete();
      theNewCanvas2->cd(2);
      MakeNiceTF1Style(fall[j],my_colors[j]);
      fall[j]->Draw("same");
      canv->cd();
      hnewUp->GetYaxis()->SetTitleOffset(0.8);
      if(j==0)
	hnewUp->Draw();
      else 
	hnewUp->Draw("same");
      fright[j]->Draw("sames");
      fleft[j]->Draw("same");
      fall[j]->Draw("same");
    }
    
    if(j==nFiles-1){
      theNewCanvas2->Close();
    }
    
    deltaZ=(fright[j]->GetParameter(0) - fall[j]->GetParameter(0))/2;
    sigmadeltaZ=0.5*TMath::Sqrt(fright[j]->GetParError(0)*fright[j]->GetParError(0) + fall[j]->GetParError(0)*fall[j]->GetParError(0));
    TString COUT = Form(" #Delta z = %.f #pm %.f #mum",deltaZ,sigmadeltaZ);
    
    lego->AddEntry(meanplots[j],LegLabels[j]+COUT); 

    if(j==nFiles-1){ 
      outfile <<deltaZ<<"|"<<sigmadeltaZ<<endl;
    }
    
    delete theNewCanvas2;

  }
 
  lego->Draw("same");
  pt->Draw("same");
  pt2->Draw("same");

}

//*************************************************************
std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > fitStudentTResiduals(TH1 *hist)
//*************************************************************
{

  hist->SetMarkerStyle(21);
  hist->SetMarkerSize(0.8);
  hist->SetStats(1);
 
  double dx = hist->GetBinWidth(1);
  double nmax = hist->GetBinContent(hist->GetMaximumBin());
  double xmax = hist->GetBinCenter(hist->GetMaximumBin());
  double nn = 7*nmax;
  
  int nb = hist->GetNbinsX();
  double n1 = hist->GetBinContent(1);
  double n9 = hist->GetBinContent(nb);
  double bg = 0.5*(n1+n9);
  
  double x1 = hist->GetBinCenter(1);
  double x9 = hist->GetBinCenter(nb);
  
  // create a TF1 with the range from x1 to x9 and 5 parameters
  
  TF1 *tp0Fcn = new TF1("tmp", tp0Fit, x1, x9, 5 );
  
  tp0Fcn->SetParName( 0, "mean" );
  tp0Fcn->SetParName( 1, "sigma" );
  tp0Fcn->SetParName( 2, "nu" );
  tp0Fcn->SetParName( 3, "area" );
  tp0Fcn->SetParName( 4, "BG" );
  
  tp0Fcn->SetNpx(500);
  tp0Fcn->SetLineWidth(2);
  //tp0Fcn->SetLineColor(kMagenta);
  //tp0Fcn->SetLineColor(kGreen);
  tp0Fcn->SetLineColor(kRed);

  // set start values for some parameters:
    
  tp0Fcn->SetParameter( 0, xmax ); // peak position
  tp0Fcn->SetParameter( 1, 4*dx ); // width
  tp0Fcn->SetParameter( 2, 2.2 ); // nu
  tp0Fcn->SetParameter( 3, nn ); // N
  tp0Fcn->SetParameter( 4, bg );
  
  hist->Fit( "tmp", "R", "ep" );
  // h->Fit("tmp","V+","ep");
  
  hist->Draw("histepsame");  // data again on top
  
  float res_mean  = tp0Fcn->GetParameter(0);
  float res_width = tp0Fcn->GetParameter(1);
  
  float res_mean_err  = tp0Fcn->GetParError(0);
  float res_width_err = tp0Fcn->GetParError(1);

  std::pair<Double_t,Double_t> resultM;
  std::pair<Double_t,Double_t> resultW;

  resultM = std::make_pair(res_mean,res_mean_err);
  resultW = std::make_pair(res_width,res_width_err);

  std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > result;
  
  result = std::make_pair(resultM,resultW);
  return result;


}

//*************************************************************
Double_t tp0Fit( Double_t *x, Double_t *par5 ) 
//*************************************************************
{
  static int nn = 0;
  nn++;
  static double dx = 0.1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) dx = x[0] - b1;
  //
  //--  Mean and width:
  //
  double xm = par5[0];
  double t = ( x[0] - xm ) / par5[1];
  double tt = t*t;
  //
  //--  exponent:
  //
  double rn = par5[2];
  double xn = 0.5 * ( rn + 1.0 );
  //
  //--  Normalization needs Gamma function:
  //
  double pk = 0.0;

  if( rn > 0.0 ) {

    double pi = 3.14159265358979323846;
    double aa = dx / par5[1] / sqrt(rn*pi) * TMath::Gamma(xn) / TMath::Gamma(0.5*rn);

    pk = par5[3] * aa * exp( -xn * log( 1.0 + tt/rn ) );
  }

  return pk + par5[4];
}

//*************************************************************
std::pair<Double_t,Double_t> getMedian(TH1F *histo)
//*************************************************************
{
  Double_t median = 999;
  int nbins = histo->GetNbinsX();

  //extract median from histogram
  double *x = new double[nbins];
  double *y = new double[nbins];
  for (int j = 0; j < nbins; j++) {
    x[j] = histo->GetBinCenter(j+1);
    y[j] = histo->GetBinContent(j+1);
  }
  median = TMath::Median(nbins, x, y);
  
  delete[] x; x = 0;
  delete [] y; y = 0;  

  std::pair<Double_t,Double_t> result;
  result = std::make_pair(median,median/TMath::Sqrt(histo->GetEntries()));

  return result;

}

//*************************************************************
std::pair<Double_t,Double_t> getMAD(TH1F *histo)
//*************************************************************
{

  int nbins = histo->GetNbinsX();
  Double_t median = getMedian(histo).first;
  Double_t x_lastBin = histo->GetBinLowEdge(nbins+1);
  const char *HistoName =histo->GetName();
  TString Finalname = Form("resMed%s",HistoName);
  TH1F *newHisto = new TH1F(Finalname,Finalname,nbins,0.,x_lastBin);
  Double_t *residuals = new Double_t[nbins];
  Double_t *weights = new Double_t[nbins];

  for (int j = 0; j < nbins; j++) {
    residuals[j] = TMath::Abs(median - histo->GetBinCenter(j+1));
    weights[j]=histo->GetBinContent(j+1);
    newHisto->Fill(residuals[j],weights[j]);
  }
  
  Double_t theMAD = (getMedian(newHisto).first)*1.4826;
  newHisto->Delete("");
  
  std::pair<Double_t,Double_t> result;
  result = std::make_pair(theMAD,theMAD/histo->GetEntries());

  return result;

}

//*************************************************************
std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > fitResiduals(TH1 *hist)
//*************************************************************
{
  //float fitResult(9999);
  //if (hist->GetEntries() < 20) return ;
  
  float mean  = hist->GetMean();
  float sigma = hist->GetRMS();
  
  if(TMath::IsNaN(mean) || TMath::IsNaN(sigma)){  
    mean=0;
    sigma= - hist->GetXaxis()->GetBinLowEdge(1) + hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX()+1) ;
  }

  TF1 func("tmp", "gaus", mean - 2.*sigma, mean + 2.*sigma); 
  if (0 == hist->Fit(&func,"QNR")) { // N: do not blow up file by storing fit!
    mean  = func.GetParameter(1);
    sigma = func.GetParameter(2);
    // second fit: three sigma of first fit around mean of first fit
    func.SetRange(mean - 2*sigma, mean + 2*sigma);
      // I: integral gives more correct results if binning is too wide
      // L: Likelihood can treat empty bins correctly (if hist not weighted...)
    if (0 == hist->Fit(&func, "Q0LR")) {
      if (hist->GetFunction(func.GetName())) { // Take care that it is later on drawn:
	hist->GetFunction(func.GetName())->ResetBit(TF1::kNotDraw);
      }
    }
  }

  float res_mean  = func.GetParameter(1);
  float res_width = func.GetParameter(2);
  
  float res_mean_err  = func.GetParError(1);
  float res_width_err = func.GetParError(2);

  std::pair<Double_t,Double_t> resultM;
  std::pair<Double_t,Double_t> resultW;

  resultM = std::make_pair(res_mean,res_mean_err);
  resultW = std::make_pair(res_width,res_width_err);

  std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > result;
  
  result = std::make_pair(resultM,resultW);
  return result;
}

//*************************************************************
Double_t DoubleSidedCB(double* x, double* par){
//*************************************************************

  double m      = x[0];
  double m0     = par[0]; 
  double sigma  = par[1];
  double alphaL = par[2];
  double alphaR = par[3];
  double nL     = par[4];
  double nR     = par[5];
  double N      = par[6];

  Double_t arg = m - m0;
  
  if (arg < 0.0) {
    Double_t t = (m-m0)/sigma; //t < 0
    Double_t absAlpha = fabs((Double_t)alphaL); //slightly redundant since alpha > 0 anyway, but never mind
    if (t >= -absAlpha) { //-absAlpha <= t < 0
      return N*exp(-0.5*t*t);
    } else {
      Double_t a = TMath::Power(nL/absAlpha,nL)*exp(-0.5*absAlpha*absAlpha);
      Double_t b = nL/absAlpha - absAlpha;
      return N*(a/TMath::Power(b - t, nL)); //b - t
    }
  } else {
    Double_t t = (m-m0)/sigma; //t > 0
    Double_t absAlpha = fabs((Double_t)alphaR);
    if (t <= absAlpha) { //0 <= t <= absAlpha
      return N*exp(-0.5*t*t);
    } else {
      Double_t a = TMath::Power(nR/absAlpha,nR)*exp(-0.5*absAlpha*absAlpha);
      Double_t b = nR/absAlpha - absAlpha;   
      return N*(a/TMath::Power(b + t, nR)); //b + t
    }
  }
}

//*************************************************************
std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > fitResidualsCB(TH1 *hist)
//*************************************************************
{
  
  //hist->Rebin(2);

  float mean  = hist->GetMean();
  float sigma = hist->GetRMS();
  //int   nbinsX   = hist->GetNbinsX();
  float nentries = hist->GetEntries();
  float meanerr  = sigma/TMath::Sqrt(nentries);
  float sigmaerr = TMath::Sqrt(sigma*sigma*TMath::Sqrt(2/nentries));

  float lowBound  = hist->GetXaxis()->GetBinLowEdge(1);
  float highBound = hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX()+1);

  if(TMath::IsNaN(mean) || TMath::IsNaN(sigma)){  
    mean=0;
    sigma= - lowBound + highBound;
  }
  
  TF1 func("tmp", "gaus", mean - 1.*sigma, mean + 1.*sigma); 
  if (0 == hist->Fit(&func,"QNR")) { // N: do not blow up file by storing fit!
    mean  = func.GetParameter(1);
    sigma = func.GetParameter(2);
  }
  
  // first round
  TF1 *doubleCB = new TF1("myDoubleCB",DoubleSidedCB,lowBound,highBound,7);
  doubleCB->SetParameters(mean,sigma,1.5,1.5,2.5,2.5,100);
  doubleCB->SetParLimits(0,mean-meanerr,mean+meanerr);
  doubleCB->SetParLimits(1,0.,sigma+2*sigmaerr);
  doubleCB->SetParLimits(2,0.,30.);
  doubleCB->SetParLimits(3,0.,30.);
  doubleCB->SetParLimits(4,0.,50.);
  doubleCB->SetParLimits(5,0.,50.);
  doubleCB->SetParLimits(6,0.,100*nentries);

  doubleCB->SetParNames("#mu","#sigma","#alpha_{L}","#alpha_{R}","n_{L}","n_{R}","N");
  doubleCB->SetLineColor(kRed);
  doubleCB->SetNpx(1000);
  // doubleCB->SetRange(0.8*lowBound,0.8*highBound);

  hist->Fit(doubleCB,"QM");

  // second round

  float p0 = doubleCB->GetParameter(0);
  float p1 = doubleCB->GetParameter(1);
  float p2 = doubleCB->GetParameter(2);
  float p3 = doubleCB->GetParameter(3);
  float p4 = doubleCB->GetParameter(4);
  float p5 = doubleCB->GetParameter(5);
  float p6 = doubleCB->GetParameter(6);
  
  float p0err = doubleCB->GetParError(0);
  float p1err = doubleCB->GetParError(1);
  float p2err = doubleCB->GetParError(2);
  float p3err = doubleCB->GetParError(3);
  float p4err = doubleCB->GetParError(4);
  float p5err = doubleCB->GetParError(5);
  float p6err = doubleCB->GetParError(6);

  if( (doubleCB->GetChisquare()/doubleCB->GetNDF()) >5){

    std::cout<<"------------------------"<<std::endl;
    std::cout<<"chi2 1st:"<<doubleCB->GetChisquare()<<std::endl;

    //std::cout<<"p0: "<<p0<<"+/-"<<p0err<<std::endl;
    //std::cout<<"p1: "<<p1<<"+/-"<<p1err<<std::endl;
    //std::cout<<"p2: "<<p2<<"+/-"<<p2err<<std::endl;
    //std::cout<<"p3: "<<p3<<"+/-"<<p3err<<std::endl;
    //std::cout<<"p4: "<<p4<<"+/-"<<p4err<<std::endl;
    //std::cout<<"p5: "<<p5<<"+/-"<<p5err<<std::endl;
    //std::cout<<"p6: "<<p6<<"+/-"<<p6err<<std::endl;

    doubleCB->SetParameters(p0,p1,3,3,6,6,p6);
    doubleCB->SetParLimits(0,p0-2*p0err,p0+2*p0err);
    doubleCB->SetParLimits(1,p1-2*p1err,p0+2*p1err);
    doubleCB->SetParLimits(2,p2-2*p2err,p0+2*p2err);
    doubleCB->SetParLimits(3,p3-2*p3err,p0+2*p3err);
    doubleCB->SetParLimits(4,p4-2*p4err,p0+2*p4err);
    doubleCB->SetParLimits(5,p5-2*p5err,p0+2*p5err);
    doubleCB->SetParLimits(6,p6-2*p6err,p0+2*p6err);

    hist->Fit(doubleCB,"MQ");

    //gMinuit->Command("SCAn 1");
    //TGraph *gr = (TGraph*)gMinuit->GetPlot();
    //gr->SetMarkerStyle(21);
    //gr->Draw("alp"); 

    std::cout<<"chi2 2nd:"<<doubleCB->GetChisquare()<<std::endl;
    
  }

  float res_mean  = doubleCB->GetParameter(0);
  float res_width = doubleCB->GetParameter(1);
  
  float res_mean_err  = doubleCB->GetParError(0);
  float res_width_err = doubleCB->GetParError(1);

  std::pair<Double_t,Double_t> resultM;
  std::pair<Double_t,Double_t> resultW;

  resultM = std::make_pair(res_mean,res_mean_err);
  resultW = std::make_pair(res_width,res_width_err);

  std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > result;
  
  result = std::make_pair(resultM,resultW);
  return result;

}

//*************************************************************
void FillTrendPlot(TH1F* trendPlot, TH1F* residualsPlot[100], TString fitPar_, TString var_,Int_t nBins)
//*************************************************************
{

  std::cout<<"trendPlot name: "<<trendPlot->GetName()<<std::endl;

  float phiInterval = (360.)/nBins;
  float etaInterval = 5./nBins;
 
  for ( int i=0; i<nBins; ++i ) {
    
    int binn = i+1;

    char phipositionString[129];
    float phiposition = (-180+i*phiInterval)+(phiInterval/2);
    sprintf(phipositionString,"%.f",phiposition);
    
    char etapositionString[129];
    float etaposition = (-2.5+i*etaInterval)+(etaInterval/2);
    sprintf(etapositionString,"%.1f",etaposition);

    std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t> > myFit = std::make_pair(std::make_pair(0.,0.),std::make_pair(0.,0.));

    if ( ((TString)trendPlot->GetName()).Contains("Norm") ) {
      //std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > myFit = fitResiduals(residualsPlot[i]);
      myFit = fitResiduals(residualsPlot[i]);
    } else {
      //std::pair<std::pair<Double_t,Double_t>, std::pair<Double_t,Double_t>  > myFit = fitStudentTResiduals(residualsPlot[i]);
      //myFit = fitStudentTResiduals(residualsPlot[i]);
      myFit = fitResiduals(residualsPlot[i]); 
    }

    if(fitPar_=="mean"){
      float mean_      = myFit.first.first;
      float meanErr_   = myFit.first.second;
      trendPlot->SetBinContent(i+1,mean_);
      trendPlot->SetBinError(i+1,meanErr_);
    } else if (fitPar_=="width"){
      float width_     = myFit.second.first;
      float widthErr_  = myFit.second.second;
      trendPlot->SetBinContent(i+1,width_);
      trendPlot->SetBinError(i+1,widthErr_);
    } else if (fitPar_=="median"){
      float median_    = getMedian(residualsPlot[i]).first;
      float medianErr_ = getMedian(residualsPlot[i]).second;
      trendPlot->SetBinContent(i+1,median_);
      trendPlot->SetBinError(i+1,medianErr_);
    } else if (fitPar_=="mad"){
      float mad_       = getMAD(residualsPlot[i]).first; 
      float madErr_    = getMAD(residualsPlot[i]).second;
      trendPlot->SetBinContent(i+1,mad_);
      trendPlot->SetBinError(i+1,madErr_);
    } else {
      std::cout<<"PrimaryVertexValidation::FillTrendPlot() "<<fitPar_<<" unknown estimator!"<<std::endl;
    }

    if(var_=="phi"){
      if( ((binn%2==0)&&(binn<=12)) || ((binn%2==1)&&(binn>12)) ) trendPlot->GetXaxis()->SetBinLabel(binn,phipositionString); 
    } else if(var_=="eta"){
      if( ((binn%2==0)&&(binn<=12)) || ((binn%2==1)&&(binn>12)) ) trendPlot->GetXaxis()->SetBinLabel(binn,etapositionString); 
    } else {
      //std::cout<<"PrimaryVertexValidation::FillTrendPlot() "<<var_<<" unknown track parameter!"<<std::endl;
    }
  }

  if(fitPar_=="mean" || fitPar_=="median"){

    TString res;
    if(TString(residualsPlot[0]->GetName()).Contains("dxy")) res="dxy";
    else if(TString(residualsPlot[0]->GetName()).Contains("dz")) res="dz";
    
    TCanvas *fitOutput = new TCanvas(Form("fitOutput_%s_%s_%s",res.Data(),var_.Data(),trendPlot->GetName()),Form("fitOutput_%s_%s",res.Data(),var_.Data()),1200,1200);
    fitOutput->Divide(6,4);
    
    TCanvas *fitPulls = new TCanvas(Form("fitPulls_%s_%s_%s",res.Data(),var_.Data(),trendPlot->GetName()),Form("fitPulls_%s_%s",res.Data(),var_.Data()),1200,1200);
    fitPulls->Divide(6,4);

    TH1F* residualsPull[nBins];

    for(Int_t i=0;i<nBins;i++){
      
      TF1 *tmp1 = (TF1*)residualsPlot[i]->GetListOfFunctions()->FindObject("tmp");
      fitOutput->cd(i+1)->SetLogy();
      fitOutput->cd(i+1)->SetBottomMargin(0.15);
      //residualsPlot[i]->Sumw2();
      residualsPlot[i]->SetMarkerStyle(20);
      residualsPlot[i]->SetMarkerSize(1.);
      residualsPlot[i]->SetStats(0);
      //residualsPlot[i]->GetXaxis()->SetRangeUser(-3*(tmp1->GetParameter(1)),3*(tmp1->GetParameter(1)));
      residualsPlot[i]->Draw("e1");
      residualsPlot[i]->GetYaxis()->UnZoom();

      //std::cout<<"*********************"<<std::endl;
      //std::cout<<"fitOutput->cd("<<i+1<<")"<<std::endl;
      //std::cout<<"residualsPlot["<<i<<"]->GetTitle() = "<<residualsPlot[i]->GetTitle()<<std::endl;
      
      // -- for chi2 ----
      TPaveText *pt = new TPaveText(0.15,0.73,0.35,0.83,"NDC");
      pt->SetFillColor(10);
      pt->SetTextColor(1);
      pt->SetTextSize(0.07);
      pt->SetTextFont(42);
      pt->SetTextAlign(11);

      //TF1 *tmp1 = (TF1*)residualsPlot[i]->GetListOfFunctions()->FindObject("tmp");
      TString COUT = Form("#chi^{2}/ndf=%.1f",tmp1->GetChisquare()/tmp1->GetNDF());
      
      TText *text1 = pt->AddText(COUT);
      text1->SetTextFont(72);
      text1->SetTextColor(kBlue);
      pt->Draw("same");
    
      // -- for bins --
     
      TPaveText *title = new TPaveText(0.1,0.93,0.8,0.95,"NDC");
      title->SetFillColor(10);
      title->SetTextColor(1);
      title->SetTextSize(0.07);
      title->SetTextFont(42);
      title->SetTextAlign(11);
 
      //TText *text2 = title->AddText(residualsPlot[i]->GetTitle());
      //text2->SetTextFont(72);
      //text2->SetTextColor(kBlue);

      title->Draw("same");

      fitPulls->cd(i+1);
      fitPulls->cd(i+1)->SetBottomMargin(0.15);
      
      residualsPull[i]=(TH1F*)residualsPlot[i]->Clone(Form("pull_%s",residualsPlot[i]->GetName()));
      for(Int_t nbin=1;nbin< residualsPull[i]->GetNbinsX(); nbin++){
      	if(residualsPlot[i]->GetBinContent(nbin)!=0){ 
	  residualsPull[i]->SetBinContent(nbin,(residualsPlot[i]->GetBinContent(nbin) - tmp1->Eval(residualsPlot[i]->GetBinCenter(nbin)))/residualsPlot[i]->GetBinContent(nbin));
	  residualsPull[i]->SetBinError(nbin,0.1);
	}
      }

      TF1* toDel = (TF1*)residualsPull[i]->FindObject("tmp");
      if(toDel) residualsPull[i]->GetListOfFunctions()->Remove(toDel); 
      residualsPull[i]->SetMarkerStyle(20);
      residualsPull[i]->SetMarkerSize(1.);
      residualsPull[i]->SetStats(0);
      // residualsPull[i]->SetOptTitle(1);
      residualsPull[i]->GetXaxis()->SetLabelFont(42);
      residualsPull[i]->GetYaxis()->SetLabelFont(42);
      residualsPull[i]->GetYaxis()->SetLabelSize(.07);
      residualsPull[i]->GetXaxis()->SetLabelSize(.07);
      residualsPull[i]->GetYaxis()->SetTitleSize(.07);
      residualsPull[i]->GetXaxis()->SetTitleSize(.07);
      residualsPull[i]->GetXaxis()->SetTitleOffset(0.9);
      residualsPull[i]->GetYaxis()->SetTitleOffset(1.2);
      residualsPull[i]->GetXaxis()->SetTitleFont(42);
      residualsPull[i]->GetYaxis()->SetTitleFont(42);
      
      residualsPull[i]->Draw("e1");
      residualsPull[i]->GetYaxis()->UnZoom();
    }
    
    
    TString tpName =trendPlot->GetName();

    TString FitNameToSame  = Form("fitOutput_%s",(tpName.ReplaceAll("means_","").Data()));
    fitOutput->SaveAs(FitNameToSame+".pdf");
    //fitOutput->SaveAs(FitNameToSame+".png");
    TString PullNameToSave = Form("fitPulls_%s",(tpName.ReplaceAll("means_","").Data()));
    fitPulls->SaveAs(PullNameToSave+".pdf");
    //fitPulls->SaveAs(PullNameToSave+".png");
    
    //fitOutput->SaveAs(Form("fitOutput_%s_%s_%s.pdf",res.Data(),var_.Data(),trendPlot->GetName()));
    //fitOutput->SaveAs(Form("fitOutput_%s.pdf",(((TString)trendPlot->GetName()).ReplaceAll("means_","")).Data()));
    //fitPulls->SaveAs(Form("fitPulls_%s.pdf",(((TString)trendPlot->GetName()).ReplaceAll("means_","")).Data()));
    //fitOutput->SaveAs(Form("fitOutput_%s.png",(((TString)trendPlot->GetName()).ReplaceAll("means_","")).Data()));

  }
}

/*--------------------------------------------------------------------*/
void MakeNiceTrendPlotStyle(TH1 *hist,Int_t color)
/*--------------------------------------------------------------------*/
{ 
  Int_t markers[4] = {kFullCircle,kFullSquare,kOpenCircle,kOpenSquare};
  Int_t colors[10]={kBlack,kRed,kBlue,kMagenta,kBlack,kRed,kBlue,kGreen};
  hist->SetStats(kFALSE);  
  hist->SetLineWidth(2);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(.05);
  hist->GetXaxis()->SetLabelSize(.05);
  hist->SetMarkerSize(1.4);
  hist->SetMarkerStyle(markers[color]);
  hist->SetLineColor(colors[color]);
  hist->SetMarkerColor(colors[color]);
}

/*--------------------------------------------------------------------*/
void setStyle()
/*--------------------------------------------------------------------*/
{
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
}

/*--------------------------------------------------------------------*/
Double_t fDLine(Double_t *x, Double_t *par)
/*--------------------------------------------------------------------*/
{
  if (x[0] < _boundSx && x[0] > _boundDx) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0];
}

/*--------------------------------------------------------------------*/
Double_t fULine(Double_t *x, Double_t *par)
/*--------------------------------------------------------------------*/
{
  if (x[0] >= _boundSx && x[0] <= _boundDx) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0];
}

/*--------------------------------------------------------------------*/
void FitULine(TH1 *hist)
/*--------------------------------------------------------------------*/
{ 
  // define fitting function
  TF1 func1("lineUp",fULine,_boundMin,_boundMax,1);
  //TF1 func1("lineUp","pol0",-0.5,11.5);
  
  if (0 == hist->Fit(&func1,"QR")) {
    if (hist->GetFunction(func1.GetName())) { // Take care that it is later on drawn:
      hist->GetFunction(func1.GetName())->ResetBit(TF1::kNotDraw);
    }
    cout<<"AugoFitResiduals() fit Up done!"<<endl;
  }
  
}

/*--------------------------------------------------------------------*/
void FitDLine(TH1 *hist)
/*--------------------------------------------------------------------*/
{
  // define fitting function
  // TF1 func1("lineDown",fDLine,-0.5,11.5,1);
  
  TF1 func2("lineDown","pol0",_boundSx,_boundDx);
  func2.SetRange(_boundSx,_boundDx);
  
  if (0 == hist->Fit(&func2,"QR")) {
    if (hist->GetFunction(func2.GetName())) { // Take care that it is later on drawn:
      hist->GetFunction(func2.GetName())->ResetBit(TF1::kNotDraw);
    }
    cout<<"FitPVResiduals() fit Down done!"<<endl;
  } 
}

/*--------------------------------------------------------------------*/
void MakeNiceTF1Style(TF1 *f1,Int_t color)
/*--------------------------------------------------------------------*/
{
  f1->SetLineColor(color);
  f1->SetLineWidth(3);
  f1->SetLineStyle(2);
}
