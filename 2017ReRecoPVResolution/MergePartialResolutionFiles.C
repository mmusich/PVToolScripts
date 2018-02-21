#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TStyle.h"
#include "TClass.h"
#include <iostream>
#include "TLegend.h"
#include <string>
#include "TPaveText.h"

#include <fstream>      // std::ofstream
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <regex>

#include "Alignment/OfflineValidation/interface/PVValidationHelpers.h"

namespace statmode{
  using fitParams = std::pair<std::pair<double,double>, std::pair<double,double> >;
}

// forward declarations

void fillTrendPlotByIndex(TH1F* trendPlot,std::map<std::string,TH1F*>& h, std::regex toMatch, PVValHelper::estimator fitPar_);
statmode::fitParams fitResiduals(TH1 *hist,bool singleTime=false);

void MergePartialResolutionFiles(TString FileName){
  
  TFile *fin = TFile::Open(FileName,"READ");
  
  // max sumPt
  const int max_sum_pt  = 30;
  std::array<float, max_sum_pt+1>  mypT_bins_ = PVValHelper::makeLogBins<float,max_sum_pt>(1.,1e3); 

  // max track
  const int max_n_tracks = 60;
  std::array<float,max_n_tracks+1> myNTrack_bins_;
  for(float i=0; i<=max_n_tracks; i++){
    myNTrack_bins_[i]=-0.5+i*2;
  }
  
  // summary plots

  TH1F*  p_resolX_vsSumPt_    = new TH1F("p_resolX_vsSumPt"  , "x-resolution vs #Sigma p_{T};#sum p_{T} [GeV]; x vertex resolution [#mum]", mypT_bins_.size()-1 , mypT_bins_.data() );  			  
  TH1F*  p_resolY_vsSumPt_    = new TH1F("p_resolY_vsSumPt"  , "y-resolution vs #Sigma p_{T};#sum p_{T} [GeV]; y vertex resolution [#mum]", mypT_bins_.size()-1 , mypT_bins_.data() );  			  
  TH1F*  p_resolZ_vsSumPt_    = new TH1F("p_resolZ_vsSumPt"  , "z-resolution vs #Sigma p_{T};#sum p_{T} [GeV]; z vertex resolution [#mum]", mypT_bins_.size()-1 , mypT_bins_.data() );  			  
  			                                                                                                                                                                                     
  TH1F*  p_resolX_vsNtracks_  = new TH1F("p_resolX_vsNtracks"  , "x-resolution vs n_{tracks};n_{tracks} in vertex; x vertex resolution [#mum]", myNTrack_bins_.size()-1 , myNTrack_bins_.data() );
  TH1F*  p_resolY_vsNtracks_  = new TH1F("p_resolY_vsNtracks"  , "y-resolution vs n_{tracks};n_{tracks} in vertex; y vertex resolution [#mum]", myNTrack_bins_.size()-1 , myNTrack_bins_.data() );
  TH1F*  p_resolZ_vsNtracks_  = new TH1F("p_resolZ_vsNtracks"  , "z-resolution vs n_{tracks};n_{tracks} in vertex; z vertex resolution [#mum]", myNTrack_bins_.size()-1 , myNTrack_bins_.data() );
  
  TH1F*  p_pullX_vsSumPt_     = new TH1F("p_pullX_vsSumPt"  , "x-pull vs #Sigma p_{T};#sum p_{T} [GeV]; x vertex pull", mypT_bins_.size()-1 , mypT_bins_.data() ); 		     
  TH1F*  p_pullY_vsSumPt_     = new TH1F("p_pullY_vsSumPt"  , "y-pull vs #Sigma p_{T};#sum p_{T} [GeV]; y vertex pull", mypT_bins_.size()-1 , mypT_bins_.data() ); 		     
  TH1F*  p_pullZ_vsSumPt_     = new TH1F("p_pullZ_vsSumPt"  , "z-pull vs #Sigma p_{T};#sum p_{T} [GeV]; z vertex pull", mypT_bins_.size()-1 , mypT_bins_.data() ); 		     
  			                                                                                                                                                                       
  TH1F*  p_pullX_vsNtracks_   = new TH1F("p_pullX_vsNtracks"  , "x-pull vs n_{tracks};n_{tracks} in vertex; x vertex pull", myNTrack_bins_.size()-1 , myNTrack_bins_.data() );
  TH1F*  p_pullY_vsNtracks_   = new TH1F("p_pullY_vsNtracks"  , "y-pull vs n_{tracks};n_{tracks} in vertex; y vertex pull", myNTrack_bins_.size()-1 , myNTrack_bins_.data() );
  TH1F*  p_pullZ_vsNtracks_   = new TH1F("p_pullZ_vsNtracks"  , "z-pull vs n_{tracks};n_{tracks} in vertex; z vertex pull", myNTrack_bins_.size()-1 , myNTrack_bins_.data() );

  // define the maps
  std::map<std::string, TH1F*> hpulls_;
  std::map<std::string, TH1F*> hdiffs_;
  
  for (int i = 0; i < max_n_tracks; i ++){
    hpulls_[Form("pullX_%dTrks",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/xPullNtracks/histo_pullX_Ntracks_plot%i",i));
    hpulls_[Form("pullY_%dTrks",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/yPullNtracks/histo_pullY_Ntracks_plot%i",i));
    hpulls_[Form("pullZ_%dTrks",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/zPullNtracks/histo_pullZ_Ntracks_plot%i",i));
    hdiffs_[Form("diffX_%dTrks",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/xResolNtracks/histo_resolX_Ntracks_plot%i",i));
    hdiffs_[Form("diffY_%dTrks",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/yResolNtracks/histo_resolY_Ntracks_plot%i",i));
    hdiffs_[Form("diffZ_%dTrks",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/zResolNtracks/histo_resolZ_Ntracks_plot%i",i));
  }


  for (int i = 0; i < max_sum_pt; i ++){
    hpulls_[Form("pullX_%dsumPt",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/xPullSumPt/histo_pullX_sumPt_plot%i",i));
    hpulls_[Form("pullY_%dsumPt",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/yPullSumPt/histo_pullY_sumPt_plot%i",i));
    hpulls_[Form("pullZ_%dsumPt",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/zPullSumPt/histo_pullZ_sumPt_plot%i",i));
    hdiffs_[Form("diffX_%dsumPt",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/xResolSumPt/histo_resolX_sumPt_plot%i",i));
    hdiffs_[Form("diffY_%dsumPt",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/yResolSumPt/histo_resolY_sumPt_plot%i",i));
    hdiffs_[Form("diffZ_%dsumPt",i)] = (TH1F*)fin->Get(Form("PrimaryVertexResolution/zResolSumPt/histo_resolZ_sumPt_plot%i",i));
  }

  // diffs

  fillTrendPlotByIndex(p_resolX_vsSumPt_  , hdiffs_,std::regex("diffX_(.*)sumPt"),PVValHelper::WIDTH);
  fillTrendPlotByIndex(p_resolY_vsSumPt_  , hdiffs_,std::regex("diffY_(.*)sumPt"),PVValHelper::WIDTH); 
  fillTrendPlotByIndex(p_resolZ_vsSumPt_  , hdiffs_,std::regex("diffZ_(.*)sumPt"),PVValHelper::WIDTH);
  
  fillTrendPlotByIndex(p_resolX_vsNtracks_, hdiffs_,std::regex("diffX_(.*)Trks"),PVValHelper::WIDTH);
  fillTrendPlotByIndex(p_resolY_vsNtracks_, hdiffs_,std::regex("diffY_(.*)Trks"),PVValHelper::WIDTH); 
  fillTrendPlotByIndex(p_resolZ_vsNtracks_, hdiffs_,std::regex("diffZ_(.*)Trks"),PVValHelper::WIDTH);

  // pulls

  fillTrendPlotByIndex(p_pullX_vsSumPt_   , hpulls_,std::regex("pullX_(.*)sumPt"),PVValHelper::WIDTH);
  fillTrendPlotByIndex(p_pullY_vsSumPt_   , hpulls_,std::regex("pullY_(.*)sumPt"),PVValHelper::WIDTH); 
  fillTrendPlotByIndex(p_pullZ_vsSumPt_   , hpulls_,std::regex("pullZ_(.*)sumPt"),PVValHelper::WIDTH);
  
  fillTrendPlotByIndex(p_pullX_vsNtracks_ , hpulls_,std::regex("pullX_(.*)Trks"),PVValHelper::WIDTH);
  fillTrendPlotByIndex(p_pullY_vsNtracks_ , hpulls_,std::regex("pullY_(.*)Trks"),PVValHelper::WIDTH); 
  fillTrendPlotByIndex(p_pullZ_vsNtracks_ , hpulls_,std::regex("pullZ_(.*)Trks"),PVValHelper::WIDTH);  

  // write out to file

  TFile *MyFile = new TFile("merged"+FileName,"RECREATE");
  if ( MyFile->IsOpen() ) cout << "File opened successfully" << endl;
  MyFile->cd();
  gDirectory->mkdir("PrimaryVertexResolution");
  MyFile->cd("PrimaryVertexResolution");


  p_resolX_vsSumPt_->Write();   
  p_resolY_vsSumPt_->Write();   
  p_resolZ_vsSumPt_->Write();   
  		     
  p_resolX_vsNtracks_->Write(); 
  p_resolY_vsNtracks_->Write(); 
  p_resolZ_vsNtracks_->Write(); 
                      
  p_pullX_vsSumPt_->Write();    
  p_pullY_vsSumPt_->Write();    
  p_pullZ_vsSumPt_->Write();    
  		     
  p_pullX_vsNtracks_->Write();  
  p_pullY_vsNtracks_->Write();  
  p_pullZ_vsNtracks_->Write();  

  MyFile->Close();
  fin->Close();
}



/*--------------------------------------------------------------------*/
void fillTrendPlotByIndex(TH1F* trendPlot,std::map<std::string,TH1F*>& h, std::regex toMatch, PVValHelper::estimator fitPar_)
/*--------------------------------------------------------------------*/
{  

  for(const auto &iterator: h) {
    
    statmode::fitParams myFit = fitResiduals(iterator.second);

    int bin=-1;
    std::string result;
    try {
      std::smatch match;
      if (std::regex_search(iterator.first, match,toMatch) && match.size() > 1) {
	result = match.str(1);
	bin = std::stoi(result)+1;
      } else {
	result = std::string("");
	continue;
      } 
    } catch (std::regex_error& e) {
      // Syntax error in the regular expression
    }
    
    switch(fitPar_)
      {
      case PVValHelper::MEAN: 
	{   
	  float mean_      = myFit.first.first;
	  float meanErr_   = myFit.first.second;
	  trendPlot->SetBinContent(bin,mean_);
	  trendPlot->SetBinError(bin,meanErr_);
	  break;
	}
      case PVValHelper::WIDTH:
	{
	  float width_     = myFit.second.first;
	  float widthErr_  = myFit.second.second;
	  trendPlot->SetBinContent(bin,width_);
	  trendPlot->SetBinError(bin,widthErr_);
	  break;
	}
      case PVValHelper::MEDIAN:
	{
	  float median_    = PVValHelper::getMedian(iterator.second).value();
	  float medianErr_ = PVValHelper::getMedian(iterator.second).error();
	  trendPlot->SetBinContent(bin,median_);
	  trendPlot->SetBinError(bin,medianErr_);
	  break;
	}
      case PVValHelper::MAD:
	{
	  float mad_       = PVValHelper::getMAD(iterator.second).value(); 
	  float madErr_    = PVValHelper::getMAD(iterator.second).error();
	  trendPlot->SetBinContent(bin,mad_);
	  trendPlot->SetBinError(bin,madErr_);
	  break;
	}
      default:
	edm::LogWarning("PrimaryVertexResolution")<<"fillTrendPlotByIndex() "<<fitPar_<<" unknown estimator!"<<std::endl;
	break;
      }
  }
}

/*--------------------------------------------------------------------*/
statmode::fitParams fitResiduals(TH1 *hist,bool singleTime)
/*--------------------------------------------------------------------*/
{
  if (hist->GetEntries() < 10){ 
    // std::cout<<"hist name: "<<hist->GetName() << std::endl;
    return std::make_pair(std::make_pair(0.,0.),std::make_pair(0.,0.));
  }
  
  float maxHist = hist->GetXaxis()->GetXmax();
  float minHist = hist->GetXaxis()->GetXmin();
  float mean  = hist->GetMean();
  float sigma = hist->GetRMS();
  
  if(TMath::IsNaN(mean) || TMath::IsNaN(sigma)){  
    mean=0;
    //sigma= - hist->GetXaxis()->GetBinLowEdge(1) + hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX()+1);
    sigma = - minHist + maxHist;
    edm::LogWarning("PrimaryVertexResolution")<< "FitPVResiduals::fitResiduals(): histogram" << hist->GetName()  << " mean or sigma are NaN!!"<< std::endl;
  }

  TF1 func("tmp", "gaus", mean - 2.*sigma, mean + 2.*sigma); 
  if (0 == hist->Fit(&func,"QNR")) { // N: do not blow up file by storing fit!
    mean  = func.GetParameter(1);
    sigma = func.GetParameter(2);

    if(!singleTime){
      // second fit: three sigma of first fit around mean of first fit
      func.SetRange(std::max(mean - 3*sigma,minHist),std::min(mean + 3*sigma,maxHist));
      // I: integral gives more correct results if binning is too wide
      // L: Likelihood can treat empty bins correctly (if hist not weighted...)
      if (0 == hist->Fit(&func, "Q0LR")) {
	if (hist->GetFunction(func.GetName())) { // Take care that it is later on drawn:
	  hist->GetFunction(func.GetName())->ResetBit(TF1::kNotDraw);
	}
      }
    }
  }

  return std::make_pair(std::make_pair(func.GetParameter(1),func.GetParError(1)),std::make_pair(func.GetParameter(2),func.GetParError(2)));

}
