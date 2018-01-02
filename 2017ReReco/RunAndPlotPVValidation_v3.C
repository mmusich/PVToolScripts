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
#include "TSystem.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include <TStopwatch.h>
#include "TArrow.h"
#include "TCanvas.h"
#include "TObjString.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <iterator>
#include <fstream>
#include <sstream>

#define DEBUG false

namespace pv{
  enum view {
    dxyphi,
    dzphi,
    dxyeta,
    dzeta,
    pT,
    generic
  };

  // brief method to find first value that doesn not compare lett
  int closest(std::vector<int> const& vec, int value) {
    auto const it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end()) { return -1; }
    return *it;
  }
}

// auxilliary struct to store
// histogram features
struct unrolledHisto {
  double m_y_min;
  double m_y_max;
  unsigned int m_n_bins;
  std::vector<double> m_bincontents;

  unrolledHisto(const double& y_min,const double& y_max, const unsigned int& n_bins, const std::vector<double>& bincontents){
    m_y_min  = y_min;
    m_y_max  = y_max;
    m_n_bins = n_bins;
    m_bincontents = bincontents;
  } //look, a constructor
  
  double get_y_min(){
    return m_y_min;
  }

  double get_y_max(){
    return m_y_max;
  }
  
  unsigned int get_n_bins(){
    return m_n_bins;
  }
  
  std::vector<double> get_bin_contents(){
    return m_bincontents;
  }

  double get_integral(){
    double ret(0.);
    for(const auto &binc: m_bincontents){
      ret+=binc;
    }
    return ret;
  }

};


// all marker and style types
Int_t markers[8] = {kFullSquare,kFullCircle,kFullTriangleDown,kOpenSquare,kOpenCircle,kFullTriangleUp,kOpenTriangleDown,kOpenTriangleUp};
Int_t colors[8]  = {kBlack,kBlue,kRed,kGreen+2,kOrange,kMagenta,kCyan,kViolet};

// forward declarations
void RunAndPlotPVValidation_v3(TString namesandlabels="",bool lumi_axis_format=false,bool time_axis_format=false,bool useRMS=true);
void arrangeOutCanvas(TCanvas *canv,
		      TH1F* m_11Trend[100],
		      TH1F* m_12Trend[100],
		      TH1F* m_21Trend[100],
		      TH1F* m_22Trend[100],
		      Int_t nFiles, 
		      TString LegLabels[10],
		      unsigned int theRun);

std::pair<std::pair<Double_t,Double_t>, Double_t> getBiases(TH1F* hist,bool useRMS_);
unrolledHisto getUnrolledHisto(TH1F* hist);

TH1F* DrawConstant(TH1F *hist,Int_t iter,Double_t theConst);
TH1F* DrawConstantGraph(TGraph *graph,Int_t iter,Double_t theConst);
std::vector<int> list_files(const char *dirname=".", const char *ext=".root");
TH1F* checkTH1AndReturn(TFile *f,TString address);
void MakeNiceTrendPlotStyle(TH1 *hist,Int_t color);
void cmsPrel(TPad* pad);
void makeNewXAxis (TH1 *h);
void beautify(TGraph *g);
void beautify(TH1 *h);
void adjustmargins(TCanvas *canv);
void adjustmargins(TVirtualPad*canv);
void setStyle();
pv::view checkTheView(const TString &toCheck);
template<typename T> void timify(T *mgr);
Double_t getMaximumFromArray(TObjArray *array);
void superImposeIOVBoundaries(TCanvas *c,bool lumi_axis_format,bool time_axis_format,const std::map<int,double> &lumiMapByRun,const std::map<int,TDatime>& timeMap);

typedef std::map<TString, std::vector<double> > alignmentTrend; 

// utility function to split strings
/*--------------------------------------------------------------------*/
std::vector<std::string> split(const std::string& s,char delimiter)
/*--------------------------------------------------------------------*/
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)){
    tokens.push_back(token);
  }
  return tokens;
}

///////////////////////////////////
//
//  Main function
//
///////////////////////////////////

void RunAndPlotPVValidation_v3(TString namesandlabels,bool lumi_axis_format,bool time_axis_format,bool useRMS){
  
  TStopwatch timer; 	 
  timer.Start();

  // consistency check, we cannot do plot vs lumi if time_axis
  if(lumi_axis_format && time_axis_format){
    std::cout<<"##########################################################################################"<<std::endl;
    std::cout<<"msg-i: RunAndPlotPVValidation_v3(): you're requesting both summary vs lumi and vs time, "<< std::endl;
    std::cout<<"       this combination is inconsistent --> exiting!"<< std::endl;
    return;
  }

  // preload the dates from file
  std::map<int,TDatime> times;
  
  if(time_axis_format){
    
    std::ifstream infile("times.txt");

    if(!infile){
      std::cout<<"missing input file :("<<std::endl;
      std::cout<<" -- exiting" << std::endl;
      return;
    }

    std::string line;
    while (std::getline(infile, line)){

      std::istringstream iss(line);
      std::string a,b,c;
      if (!(iss >> a >> b >> c)) { break; } // error
      
      //std::cout<<a<<"  "<<b<<"   "<<c<<"   "<<std::endl;
      
      int run  = std::stoi(a);
      auto tokens_b = split(b,'-');
      int year  = std::stoi(tokens_b[0]);
      int month = std::stoi(tokens_b[1]);
      int day   = std::stoi(tokens_b[2]);
      
      auto tokens_c  = split(c,'.');
      auto tokens_c1 = split(tokens_c[0],':');

      int hour   = std::stoi(tokens_c1[0]);
      int minute = std::stoi(tokens_c1[2]);
      int second = std::stoi(tokens_c1[2]);

      //std::cout<<run<<" "<<year<<" "<<month<<" "<<day<<" "<<hour<<" "<<minute<<" "<<second<<" "<<std::endl;
	
      TDatime da(year,month,day,hour,minute,second);
      times[run]=da;
    }    
  } // if time axis in the plots

  std::ofstream outfile ("lumiByRun.txt"); 
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
  std::vector<double> runs;
  std::vector<double> lumiByRun;
  std::map<int,double> lumiMapByRun;
  std::vector<double> x_ticks;

  std::vector<double> runtimes;
  if(time_axis_format){
    for(const auto &element : times){
      runtimes.push_back((element.second).Convert());
    }
  }

  for(Int_t j=0; j < nDirs_; j++) {
    
    // Retrieve labels
    TObjString* legend = (TObjString*)LabelList->At(j);
    TObjString* dir    = (TObjString*)DirList->At(j);
    LegLabels[j] = legend->String();
    dirs[j] = (dir->String()).Data();
    cout<<"RunAndPlotPVValidation_v3(): label["<<j<<"]"<<LegLabels[j]<<endl;
    
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
  
  // debug only
  for(UInt_t index=0;index<intersection.size();index++){
    std::cout<<index<<" "<<intersection[index]<<std::endl;
  }


  // book the vectors of values
  alignmentTrend dxyPhiMeans_;
  alignmentTrend dxyPhiHi_;
  alignmentTrend dxyPhiLo_;
 
  alignmentTrend dxyEtaMeans_;
  alignmentTrend dxyEtaHi_;
  alignmentTrend dxyEtaLo_;
  
  alignmentTrend dzPhiMeans_;
  alignmentTrend dzPhiHi_;
  alignmentTrend dzPhiLo_;
  
  alignmentTrend dzEtaMeans_;
  alignmentTrend dzEtaHi_;
  alignmentTrend dzEtaLo_;

  // unrolled histos
  
  std::map<TString,std::vector<unrolledHisto> > dxyVect;
  std::map<TString,std::vector<unrolledHisto> > dzVect;

  double lumiSoFar=0.0;

  // loop over the runs in the intersection
  unsigned int last = (DEBUG==true) ? 50 : intersection.size();

  for(unsigned int n=0; n<last;n++){
  //in case of debug, use only 50
  //for(unsigned int n=0; n<50;n++){

    //if(intersection.at(n)!=283946) 
    //  continue;

    std::cout << n << " "<<intersection.at(n) << std::endl;
    
    TFile *fins[nDirs_];

    TH1F* dxyPhiMeanTrend[nDirs_];  
    TH1F* dxyPhiWidthTrend[nDirs_]; 
    TH1F* dzPhiMeanTrend[nDirs_];   
    TH1F* dzPhiWidthTrend[nDirs_];  

    TH1F* dxyLadderMeanTrend[nDirs_];    
    TH1F* dxyLadderWidthTrend[nDirs_]; 
    TH1F* dzLadderWidthTrend[nDirs_];  
    TH1F* dzLadderMeanTrend[nDirs_];  

    TH1F* dxyModZMeanTrend[nDirs_];  
    TH1F* dxyModZWidthTrend[nDirs_]; 
    TH1F* dzModZMeanTrend[nDirs_];   
    TH1F* dzModZWidthTrend[nDirs_];  

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
        
    TH1F* dxyIntegralTrend[nDirs_];
    TH1F* dzIntegralTrend[nDirs_];

    bool areAllFilesOK = true;
    Int_t lastOpen = 0;
 
    // loop over the objects
    for(Int_t j=0; j < nDirs_; j++) {

      //fins[j] = TFile::Open(Form("%s/PVValidation_%s_%i.root",dirs[j],dirs[j],intersection[n]));

      size_t position = std::string(dirs[j]).find("/");   
      string stem = std::string(dirs[j]).substr(position+1);     // get from position to the end
      
      fins[j] = new TFile(Form("%s/PVValidation_%s_%i.root",dirs[j],stem.c_str(),intersection[n]));
      if(fins[j]->IsZombie()){
	std::cout<< Form("%s/PVValidation_%s_%i.root",dirs[j],stem.c_str(),intersection[n]) << " is a Zombie! cannot combine" << std::endl;
	areAllFilesOK = false;
	lastOpen=j;
	break;
      }

      std::cout<< Form("%s/PVValidation_%s_%i.root",dirs[j],stem.c_str(),intersection[n]) 
	       << " has size: "<<fins[j]->GetSize() << " b ";
      
      // sanity check
      TH1F* h_tracks = (TH1F*)fins[j]->Get("PVValidation/EventFeatures/h_nTracks");
      if(j==0){
	TH1F* h_lumi   = (TH1F*)fins[j]->Get("PVValidation/EventFeatures/h_lumiFromConfig");
	double lumi = h_lumi->GetBinContent(1);
	lumiSoFar+=lumi;
	//std::cout<<"lumi: "<<lumi
	//		 <<" ,lumi so far: "<<lumiSoFar<<std::endl;

	outfile<<"run "<<intersection[n]<<" lumi: "<<lumi
	       <<" ,lumi so far: "<<lumiSoFar<<std::endl;

      }

      Double_t numEvents = h_tracks->GetEntries();
      if(numEvents<2500){
	std::cout<<"excluding " << intersection[n] << "because it has less than 2.5k events" << std::endl;
	areAllFilesOK = false;
	lastOpen=j;
	break;
      }

      dxyPhiMeanTrend[j]     = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_phi");
      //dxyPhiMeanTrend[j]     = checkTH1AndReturn(fins[j],"PVValidation/MeanTrends/means_dxy_phi");
      dxyPhiWidthTrend[j]    = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_phi");
      dzPhiMeanTrend[j]      = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_phi");
      dzPhiWidthTrend[j]     = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_phi");

      dxyLadderMeanTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_ladder");
      dxyLadderWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_ladder");
      dzLadderMeanTrend[j]   = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_ladder");
      dzLadderWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_ladder");
                              
      dxyEtaMeanTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_eta");
      dxyEtaWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_eta");
      dzEtaMeanTrend[j]   = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_eta");
      dzEtaWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_eta");
      
      dxyModZMeanTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dxy_modZ");
      dxyModZWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_modZ");
      dzModZMeanTrend[j]   = (TH1F*)fins[j]->Get("PVValidation/MeanTrends/means_dz_modZ");
      dzModZWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_modZ");
      
      dxyNormPhiWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dxy_phi");
      dxyNormEtaWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dxy_eta");
      dzNormPhiWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dz_phi");
      dzNormEtaWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dz_eta");

      dxyNormPtWidthTrend[j] = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dxy_pTCentral");
      dzNormPtWidthTrend[j]  = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/norm_widths_dz_pTCentral");
      dxyPtWidthTrend[j]     = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dxy_pTCentral");
      dzPtWidthTrend[j]      = (TH1F*)fins[j]->Get("PVValidation/WidthTrends/widths_dz_pTCentral");

      dxyIntegralTrend[j]    = (TH1F*)fins[j]->Get("PVValidation/ProbeTrackFeatures/h_probedxyRefitV");
      dzIntegralTrend[j]     = (TH1F*)fins[j]->Get("PVValidation/ProbeTrackFeatures/h_probedzRefitV");

      // fill the vectors of biases
      
      auto dxyPhiBiases = getBiases(dxyPhiMeanTrend[j],useRMS);
      
      //std::cout<<"\n" <<j<<" "<< LegLabels[j] << " dxy(phi) mean: "<< dxyPhiBiases.second
      //	       <<" dxy(phi) max: "<< dxyPhiBiases.first.first
      //       <<" dxy(phi) min: "<< dxyPhiBiases.first.second
      //       << std::endl;

      dxyPhiMeans_[LegLabels[j]].push_back(dxyPhiBiases.second);
      dxyPhiLo_[LegLabels[j]].push_back(dxyPhiBiases.first.first);
      dxyPhiHi_[LegLabels[j]].push_back(dxyPhiBiases.first.second);

      auto dxyEtaBiases = getBiases(dxyEtaMeanTrend[j],useRMS);
      dxyEtaMeans_[LegLabels[j]].push_back(dxyEtaBiases.second);
      dxyEtaLo_[LegLabels[j]].push_back(dxyEtaBiases.first.first);
      dxyEtaHi_[LegLabels[j]].push_back(dxyEtaBiases.first.second);

      auto dzPhiBiases = getBiases(dzPhiMeanTrend[j],useRMS);
      dzPhiMeans_[LegLabels[j]].push_back(dzPhiBiases.second);
      dzPhiLo_[LegLabels[j]].push_back(dzPhiBiases.first.first);
      dzPhiHi_[LegLabels[j]].push_back(dzPhiBiases.first.second);

      auto dzEtaBiases = getBiases(dzEtaMeanTrend[j],useRMS);
      dzEtaMeans_[LegLabels[j]].push_back(dzEtaBiases.second);
      dzEtaLo_[LegLabels[j]].push_back(dzEtaBiases.first.first);
      dzEtaHi_[LegLabels[j]].push_back(dzEtaBiases.first.second);

      // unrolled histograms
      dxyVect[LegLabels[j]].push_back(getUnrolledHisto(dxyIntegralTrend[j]));
      dzVect[LegLabels[j]].push_back(getUnrolledHisto(dzIntegralTrend[j]));
      
      std::cout<<std::endl;
      std::cout<<" n. bins: "<< dxyVect[LegLabels[j]].back().get_n_bins() 
	       <<" y-min:   "<< dxyVect[LegLabels[j]].back().get_y_min()
	       <<" y-max:   "<< dxyVect[LegLabels[j]].back().get_y_max() << std::endl;

      // beautify the histograms
      MakeNiceTrendPlotStyle(dxyPhiMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dxyPhiWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzPhiMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dzPhiWidthTrend[j],j);

      MakeNiceTrendPlotStyle(dxyLadderMeanTrend[j],j); 
      MakeNiceTrendPlotStyle(dxyLadderWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzLadderMeanTrend[j],j); 
      MakeNiceTrendPlotStyle(dzLadderWidthTrend[j],j);

      MakeNiceTrendPlotStyle(dxyEtaMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dxyEtaWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzEtaMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dzEtaWidthTrend[j],j);  

      MakeNiceTrendPlotStyle(dxyModZMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dxyModZWidthTrend[j],j);
      MakeNiceTrendPlotStyle(dzModZMeanTrend[j],j);
      MakeNiceTrendPlotStyle(dzModZWidthTrend[j],j);  

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
    } else {
      runs.push_back(intersection.at(n));
      // push back the vector of lumi (in fb at this point)
      lumiByRun.push_back(lumiSoFar/1000.);
      lumiMapByRun[intersection.at(n)]=lumiSoFar/1000.;
    }

    //std::cout<<"I am still here"<<std::endl;

    // Bias plots

    TCanvas *BiasesCanvas = new TCanvas(Form("Biases_%i",intersection.at(n)),"Biases",1200,1200);
    arrangeOutCanvas(BiasesCanvas,dxyPhiMeanTrend,dzPhiMeanTrend,dxyEtaMeanTrend,dzEtaMeanTrend,nDirs_,LegLabels,intersection.at(n));

    BiasesCanvas->SaveAs(Form("Biases_%i.pdf",intersection.at(n)));
    BiasesCanvas->SaveAs(Form("Biases_%i.png",intersection.at(n)));
   
    // Bias vs L1 modules position

    TCanvas *BiasesL1Canvas = new TCanvas(Form("BiasesL1_%i",intersection.at(n)),"BiasesL1",1200,1200);
    arrangeOutCanvas(BiasesL1Canvas,dxyLadderMeanTrend,dzLadderMeanTrend,dxyModZMeanTrend,dzModZMeanTrend,nDirs_,LegLabels,intersection.at(n));

    BiasesL1Canvas->SaveAs(Form("BiasesL1_%i.pdf",intersection.at(n)));
    BiasesL1Canvas->SaveAs(Form("BiasesL1_%i.png",intersection.at(n)));

    // Resolution plots

    TCanvas *ResolutionsCanvas = new TCanvas(Form("Resolutions_%i",intersection.at(n)),"Resolutions",1200,1200);
    arrangeOutCanvas(ResolutionsCanvas,dxyPhiWidthTrend,dzPhiWidthTrend,dxyEtaWidthTrend,dzEtaWidthTrend,nDirs_,LegLabels,intersection.at(n));

    ResolutionsCanvas->SaveAs(Form("Resolutions_%i.pdf",intersection.at(n)));
    ResolutionsCanvas->SaveAs(Form("Resolutions_%i.png",intersection.at(n)));
    
    // Resolution plots vs L1 modules position

    TCanvas *ResolutionsL1Canvas = new TCanvas(Form("ResolutionsL1_%i",intersection.at(n)),"Resolutions",1200,1200);
    arrangeOutCanvas(ResolutionsL1Canvas,dxyLadderWidthTrend,dzLadderWidthTrend,dxyModZWidthTrend,dzModZWidthTrend,nDirs_,LegLabels,intersection.at(n));

    ResolutionsL1Canvas->SaveAs(Form("ResolutionsL1_%i.pdf",intersection.at(n)));
    ResolutionsL1Canvas->SaveAs(Form("ResolutionsL1_%i.png",intersection.at(n)));

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
    
    delete BiasesCanvas;
    delete BiasesL1Canvas;
    delete ResolutionsCanvas; 
    delete ResolutionsL1Canvas;
    delete PullsCanvas;
    delete ResolutionsVsPt;

    std::cout<<std::endl;
  }

  // do the trend-plotting!

  TCanvas *c_dxy_phi_vs_run = new TCanvas("c_dxy_phi_vs_run","dxy(#phi) bias vs run number",1600,800);
  TCanvas *c_dxy_eta_vs_run = new TCanvas("c_dxy_eta_vs_run","dxy(#eta) bias vs run number",1600,800);
  TCanvas *c_dz_phi_vs_run  = new TCanvas("c_dz_phi_vs_run" ,"dz(#phi) bias vs run number" ,1600,800);
  TCanvas *c_dz_eta_vs_run  = new TCanvas("c_dz_eta_vs_run" ,"dz(#eta) bias vs run number" ,1600,800);

  TCanvas *c_RMS_dxy_phi_vs_run = new TCanvas("c_RMS_dxy_phi_vs_run","dxy(#phi) bias vs run number",1600,800);
  TCanvas *c_RMS_dxy_eta_vs_run = new TCanvas("c_RMS_dxy_eta_vs_run","dxy(#eta) bias vs run number",1600,800);
  TCanvas *c_RMS_dz_phi_vs_run  = new TCanvas("c_RMS_dz_phi_vs_run","dxy(#phi) bias vs run number",1600,800);
  TCanvas *c_RMS_dz_eta_vs_run  = new TCanvas("c_RMS_dz_eta_vs_run","dxy(#eta) bias vs run number",1600,800);

  TCanvas *c_mean_dxy_phi_vs_run = new TCanvas("c_mean_dxy_phi_vs_run","dxy(#phi) bias vs run number",1600,800);
  TCanvas *c_mean_dxy_eta_vs_run = new TCanvas("c_mean_dxy_eta_vs_run","dxy(#eta) bias vs run number",1600,800);
  TCanvas *c_mean_dz_phi_vs_run  = new TCanvas("c_mean_dz_phi_vs_run","dxy(#phi) bias vs run number",1600,800);
  TCanvas *c_mean_dz_eta_vs_run  = new TCanvas("c_mean_dz_eta_vs_run","dxy(#eta) bias vs run number",1600,800);

  TCanvas *Scatter_dxy_vs_run = new TCanvas("Scatter_dxy_vs_run","dxy bias vs run number",1600,800);
  Scatter_dxy_vs_run->Divide(1,nDirs_);
  TCanvas *Scatter_dz_vs_run  = new TCanvas("Scatter_dz_vs_run","dxy bias vs run number",1600,800);
  Scatter_dz_vs_run->Divide(1,nDirs_);    

  // bias on the mean

  TGraph *g_dxy_phi_vs_run[nDirs_];
  TGraph *gprime_dxy_phi_vs_run[nDirs_];  
  TGraph *g_dxy_phi_hi_vs_run[nDirs_];
  TGraph *g_dxy_phi_lo_vs_run[nDirs_];
  
  TGraph *g_dxy_eta_vs_run[nDirs_];
  TGraph *gprime_dxy_eta_vs_run[nDirs_];
  TGraph *g_dxy_eta_hi_vs_run[nDirs_];
  TGraph *g_dxy_eta_lo_vs_run[nDirs_];

  TGraph *g_dz_phi_vs_run[nDirs_];
  TGraph *gprime_dz_phi_vs_run[nDirs_];
  TGraph *g_dz_phi_hi_vs_run[nDirs_];
  TGraph *g_dz_phi_lo_vs_run[nDirs_];

  TGraph *g_dz_eta_vs_run[nDirs_];
  TGraph *gprime_dz_eta_vs_run[nDirs_];
  TGraph *g_dz_eta_hi_vs_run[nDirs_];
  TGraph *g_dz_eta_lo_vs_run[nDirs_];

  // resolutions 

  TH1F *h_RMS_dxy_phi_vs_run[nDirs_];
  TH1F *h_RMS_dxy_eta_vs_run[nDirs_];
  TH1F *h_RMS_dz_phi_vs_run[nDirs_];
  TH1F *h_RMS_dz_eta_vs_run[nDirs_];   

  // scatters of integrated bias

  TH2F *h2_scatter_dxy_vs_run[nDirs_];
  TH2F *h2_scatter_dz_vs_run[nDirs_];

  // decide the type

  TString theType="";
  TString theTypeLabel="";
  if(lumi_axis_format){
    theType="luminosity";
    theTypeLabel="processed luminosity (1/fb)";
    x_ticks = lumiByRun;
  } else {
    if(!time_axis_format){
      theType="run number";
      theTypeLabel="run number";
      x_ticks = runs;
    } else {
      theType="date";
      theTypeLabel="date";
      for(const auto &run : runs){
	x_ticks.push_back(times[run].Convert());
      }
    }
  }
  
  TLegend *my_lego = new TLegend(0.12,0.80,0.25,0.93);
  //my_lego-> SetNColumns(2);
  my_lego->SetFillColor(10);
  my_lego->SetTextSize(0.042);
  my_lego->SetTextFont(42);
  my_lego->SetFillColor(10);
  my_lego->SetLineColor(10);
  my_lego->SetShadowColor(10);

  // arrays for storing RMS histograms
  TObjArray *arr_dxy_phi = new TObjArray();
  TObjArray *arr_dz_phi  = new TObjArray();
  TObjArray *arr_dxy_eta = new TObjArray();
  TObjArray *arr_dz_eta  = new TObjArray();

  arr_dxy_phi->Expand(nDirs_);
  arr_dz_phi->Expand(nDirs_);
  arr_dxy_eta->Expand(nDirs_);
  arr_dz_eta->Expand(nDirs_);

  for(Int_t j=0; j < nDirs_; j++) {

    // check on the sanity
    std::cout<<"x_ticks.size()= "<<x_ticks.size()<<"d xyPhiMeans_[LegLabels["<<j<<"]].size()="<<dxyPhiMeans_[LegLabels[j]].size()<<std::endl;

    // *************************************
    // dxy vs phi
    // *************************************
    
    g_dxy_phi_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dxyPhiMeans_[LegLabels[j]])[0]));
    g_dxy_phi_hi_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dxyPhiHi_[LegLabels[j]])[0]));
    g_dxy_phi_lo_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dxyPhiLo_[LegLabels[j]])[0]));

    adjustmargins(c_dxy_phi_vs_run);
    c_dxy_phi_vs_run->cd();
    g_dxy_phi_vs_run[j]->SetMarkerStyle(markers[j]);
    g_dxy_phi_vs_run[j]->SetMarkerColor(colors[j]);
    g_dxy_phi_vs_run[j]->SetLineColor(colors[j]);
    g_dxy_phi_hi_vs_run[j]->SetLineColor(colors[j]);
    g_dxy_phi_lo_vs_run[j]->SetLineColor(colors[j]);

    g_dxy_phi_vs_run[j]->SetName(Form("g_bias_dxy_phi_%s",LegLabels[j].Data()));
    g_dxy_phi_vs_run[j]->SetTitle(Form("Bias of d_{xy}(#varphi) vs %s",theType.Data()));
    g_dxy_phi_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_dxy_phi_vs_run[j]->GetYaxis()->SetTitle("#LT d_{xy}(#phi) #GT [#mum]");
    g_dxy_phi_vs_run[j]->GetYaxis()->SetRangeUser(-50,50);
    beautify(g_dxy_phi_vs_run[j]);
 
    my_lego->AddEntry(g_dxy_phi_vs_run[j],LegLabels[j],"PL");

    if(j==0){
      g_dxy_phi_vs_run[j]->Draw("AP");
    } else {
      g_dxy_phi_vs_run[j]->Draw("Psame");
    }
    g_dxy_phi_hi_vs_run[j]->Draw("Lsame");
    g_dxy_phi_lo_vs_run[j]->Draw("Lsame");

    if(time_axis_format){
      timify(g_dxy_phi_vs_run[j]);
      timify(g_dxy_phi_hi_vs_run[j]);
      timify(g_dxy_phi_lo_vs_run[j]);
    }

    if(j==nDirs_-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(g_dxy_phi_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    TPad *current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_dxy_phi_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // mean only
    adjustmargins(c_mean_dxy_phi_vs_run);
    c_mean_dxy_phi_vs_run->cd();
    gprime_dxy_phi_vs_run[j] =  (TGraph*)g_dxy_phi_vs_run[j]->Clone();
    if(j==0){
      gprime_dxy_phi_vs_run[j]->GetYaxis()->SetRangeUser(-10.,10.);
      gprime_dxy_phi_vs_run[j]->Draw("APL");
    } else {
      gprime_dxy_phi_vs_run[j]->Draw("PLsame");
    }
    
    if(j==nDirs_-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(gprime_dxy_phi_vs_run[j],2,0.);
      theZero->Draw("E1same");

      current_pad = static_cast<TPad*>(gPad);
      cmsPrel(current_pad);

      superImposeIOVBoundaries(c_mean_dxy_phi_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);
    }

    // scatter or RMS TH1

    h_RMS_dxy_phi_vs_run[j] = new TH1F(Form("h_RMS_dxy_phi_%s",LegLabels[j].Data()),Form("scatter of d_{xy}(#varphi) vs %s;%s;maximum scatter of d_{xy}(#phi) [#mum]",theType.Data(),theTypeLabel.Data()),x_ticks.size()-1,&(x_ticks[0]));
    h_RMS_dxy_phi_vs_run[j]->SetStats(kFALSE);

    int bincounter=0;
    for(const auto &tick : x_ticks ){
      bincounter++;
      h_RMS_dxy_phi_vs_run[j]->SetBinContent(bincounter,std::abs(dxyPhiHi_[LegLabels[j]][bincounter-1]-dxyPhiLo_[LegLabels[j]][bincounter-1]));
      h_RMS_dxy_phi_vs_run[j]->SetBinError(bincounter,0.01);
    }

    h_RMS_dxy_phi_vs_run[j]->SetLineColor(colors[j]);
    h_RMS_dxy_phi_vs_run[j]->SetLineWidth(2);
    h_RMS_dxy_phi_vs_run[j]->SetMarkerStyle(markers[j]);
    h_RMS_dxy_phi_vs_run[j]->SetMarkerColor(colors[j]);
    adjustmargins(c_RMS_dxy_phi_vs_run);
    c_RMS_dxy_phi_vs_run->cd();
    beautify(h_RMS_dxy_phi_vs_run[j]);

    if(time_axis_format){
      timify(h_RMS_dxy_phi_vs_run[j]);
    }

    arr_dxy_phi->Add(h_RMS_dxy_phi_vs_run[j]);

    // at the last file re-loop
    if(j==nDirs_-1){

      auto theMax = getMaximumFromArray(arr_dxy_phi);

      for(Int_t k=0; k< nDirs_; k++){
	h_RMS_dxy_phi_vs_run[k]->GetYaxis()->SetRangeUser(-theMax*0.45,theMax*1.3);
	if(k==0){
	  h_RMS_dxy_phi_vs_run[k]->Draw("L");
	} else {
	  h_RMS_dxy_phi_vs_run[k]->Draw("Lsame");
	}
      }
      
      my_lego->Draw("same");
      TH1F* theConst = DrawConstant(h_RMS_dxy_phi_vs_run[j],1,0.);
      theConst->Draw("same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_RMS_dxy_phi_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // dxy vs eta
    // *************************************
    g_dxy_eta_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dxyEtaMeans_[LegLabels[j]])[0]));
    g_dxy_eta_hi_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dxyEtaHi_[LegLabels[j]])[0]));
    g_dxy_eta_lo_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dxyEtaLo_[LegLabels[j]])[0]));

    adjustmargins(c_dxy_eta_vs_run);
    c_dxy_eta_vs_run->cd();
    g_dxy_eta_vs_run[j]->SetMarkerStyle(markers[j]);
    g_dxy_eta_vs_run[j]->SetMarkerColor(colors[j]);
    g_dxy_eta_vs_run[j]->SetLineColor(colors[j]);
    g_dxy_eta_hi_vs_run[j]->SetLineColor(colors[j]);
    g_dxy_eta_lo_vs_run[j]->SetLineColor(colors[j]);
    
    g_dxy_eta_vs_run[j]->SetName(Form("g_bias_dxy_eta_%s",LegLabels[j].Data()));
    g_dxy_eta_vs_run[j]->SetTitle(Form("Bias of d_{xy}(#eta) vs %s",theType.Data()));
    g_dxy_eta_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_dxy_eta_vs_run[j]->GetYaxis()->SetTitle("#LT d_{xy}(#eta) #GT [#mum]");
    beautify(g_dxy_eta_vs_run[j]);

    g_dxy_eta_vs_run[j]->GetYaxis()->SetRangeUser(-20,20);
    if(j==0){
      g_dxy_eta_vs_run[j]->Draw("AP");
    } else {
      g_dxy_eta_vs_run[j]->Draw("Psame");
    }
    g_dxy_eta_hi_vs_run[j]->Draw("Lsame");
    g_dxy_eta_lo_vs_run[j]->Draw("Lsame");

    if(time_axis_format){
      timify(g_dxy_eta_vs_run[j]);
      timify(g_dxy_eta_hi_vs_run[j]);
      timify(g_dxy_eta_lo_vs_run[j]);
    }

    if(j==nDirs_-1){
      my_lego->Draw("same");
      TH1F* theZero = DrawConstantGraph(g_dxy_eta_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }
	
    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_dxy_eta_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // mean only
    adjustmargins(c_mean_dxy_eta_vs_run);
    c_mean_dxy_eta_vs_run->cd();
    gprime_dxy_eta_vs_run[j] =  (TGraph*)g_dxy_eta_vs_run[j]->Clone();
    if(j==0){
      gprime_dxy_eta_vs_run[j]->GetYaxis()->SetRangeUser(-10.,10.);
      gprime_dxy_eta_vs_run[j]->Draw("APL");
    } else {
      gprime_dxy_eta_vs_run[j]->Draw("PLsame");
    }
    
    if(j==nDirs_-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(gprime_dxy_eta_vs_run[j],2,0.);
      theZero->Draw("E1same");

      current_pad = static_cast<TPad*>(gPad);
      cmsPrel(current_pad);

      superImposeIOVBoundaries(c_mean_dxy_eta_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);
    }

    // scatter or RMS TH1

    h_RMS_dxy_eta_vs_run[j] = new TH1F(Form("h_RMS_dxy_eta_%s",LegLabels[j].Data()),Form("scatter of d_{xy}(#eta) vs %s;%s;maximum scatter of d_{xy}(#eta) [#mum]",theType.Data(),theTypeLabel.Data()),x_ticks.size()-1,&(x_ticks[0]));
    h_RMS_dxy_eta_vs_run[j]->SetStats(kFALSE);

    bincounter=0;
    for(const auto &tick : x_ticks ){
      bincounter++;
      h_RMS_dxy_eta_vs_run[j]->SetBinContent(bincounter,std::abs(dxyEtaHi_[LegLabels[j]][bincounter-1]-dxyEtaLo_[LegLabels[j]][bincounter-1]));
      h_RMS_dxy_eta_vs_run[j]->SetBinError(bincounter,0.01);
    }

    h_RMS_dxy_eta_vs_run[j]->SetLineColor(colors[j]);
    h_RMS_dxy_eta_vs_run[j]->SetLineWidth(2);
    h_RMS_dxy_eta_vs_run[j]->SetMarkerStyle(markers[j]);
    h_RMS_dxy_eta_vs_run[j]->SetMarkerColor(colors[j]);
    adjustmargins(c_RMS_dxy_eta_vs_run);
    c_RMS_dxy_eta_vs_run->cd();
    beautify(h_RMS_dxy_eta_vs_run[j]);

    if(time_axis_format){
      timify(h_RMS_dxy_eta_vs_run[j]);
    }

    arr_dxy_eta->Add(h_RMS_dxy_eta_vs_run[j]);

    // at the last file re-loop
    if(j==nDirs_-1){

      auto theMax = getMaximumFromArray(arr_dxy_eta);

      for(Int_t k=0; k< nDirs_; k++){
	h_RMS_dxy_eta_vs_run[k]->GetYaxis()->SetRangeUser(-theMax*0.45,theMax*1.30);
	if(k==0){
	  h_RMS_dxy_eta_vs_run[k]->Draw("L");
	} else {
	  h_RMS_dxy_eta_vs_run[k]->Draw("Lsame");
	}
      }
      my_lego->Draw("same");
      TH1F* theConst = DrawConstant(h_RMS_dxy_eta_vs_run[j],1,0.);
      theConst->Draw("same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_RMS_dxy_eta_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // dz vs phi
    // *************************************
    g_dz_phi_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dzPhiMeans_[LegLabels[j]])[0]));
    g_dz_phi_hi_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dzPhiHi_[LegLabels[j]])[0]));
    g_dz_phi_lo_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dzPhiLo_[LegLabels[j]])[0]));

    adjustmargins(c_dz_phi_vs_run);
    c_dz_phi_vs_run->cd();
    g_dz_phi_vs_run[j]->SetMarkerStyle(markers[j]);
    g_dz_phi_vs_run[j]->SetMarkerColor(colors[j]);
    g_dz_phi_vs_run[j]->SetLineColor(colors[j]);
    g_dz_phi_hi_vs_run[j]->SetLineColor(colors[j]);
    g_dz_phi_lo_vs_run[j]->SetLineColor(colors[j]);
    beautify(g_dz_phi_vs_run[j]);

    g_dz_phi_vs_run[j]->SetName(Form("g_bias_dz_phi_%s",LegLabels[j].Data()));
    g_dz_phi_vs_run[j]->SetTitle(Form("Bias of d_{z}(#varphi) vs %s",theType.Data()));
    g_dz_phi_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_dz_phi_vs_run[j]->GetYaxis()->SetTitle("#LT d_{z}(#phi) #GT [#mum]");
    
    g_dz_phi_vs_run[j]->GetYaxis()->SetRangeUser(-20,20);
    if(j==0){
      g_dz_phi_vs_run[j]->Draw("AP");
    } else {
      g_dz_phi_vs_run[j]->Draw("Psame");
    }
    g_dz_phi_hi_vs_run[j]->Draw("Lsame");
    g_dz_phi_lo_vs_run[j]->Draw("Lsame");

    if(time_axis_format){
      timify(g_dz_phi_vs_run[j]);
      timify(g_dz_phi_hi_vs_run[j]);
      timify(g_dz_phi_lo_vs_run[j]);
    }

    if(j==nDirs_-1){
      my_lego->Draw("same");
      TH1F* theZero = DrawConstantGraph(g_dz_phi_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);
    
    superImposeIOVBoundaries(c_dz_phi_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // mean only
    adjustmargins(c_mean_dz_phi_vs_run);
    c_mean_dz_phi_vs_run->cd();
    gprime_dz_phi_vs_run[j] =  (TGraph*)g_dz_phi_vs_run[j]->Clone();
    if(j==0){
      gprime_dz_phi_vs_run[j]->GetYaxis()->SetRangeUser(-10.,10.);
      gprime_dz_phi_vs_run[j]->Draw("APL");
    } else {
      gprime_dz_phi_vs_run[j]->Draw("PLsame");
    }
    
    if(j==nDirs_-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(gprime_dz_phi_vs_run[j],2,0.);
      theZero->Draw("E1same");

      current_pad = static_cast<TPad*>(gPad);
      cmsPrel(current_pad);

      superImposeIOVBoundaries(c_mean_dz_phi_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);
    }

    // scatter or RMS TH1

    h_RMS_dz_phi_vs_run[j] = new TH1F(Form("h_RMS_dz_phi_%s",LegLabels[j].Data()),Form("scatter of d_{xy}(#varphi) vs %s;%s;maximum scatter of d_{z}(#phi) [#mum]",theType.Data(),theTypeLabel.Data()),x_ticks.size()-1,&(x_ticks[0]));
    h_RMS_dz_phi_vs_run[j]->SetStats(kFALSE);

    bincounter=0;
    for(const auto &tick : x_ticks ){
      bincounter++;
      h_RMS_dz_phi_vs_run[j]->SetBinContent(bincounter,std::abs(dzPhiHi_[LegLabels[j]][bincounter-1]-dzPhiLo_[LegLabels[j]][bincounter-1]));
      h_RMS_dz_phi_vs_run[j]->SetBinError(bincounter,0.01);
    }

    h_RMS_dz_phi_vs_run[j]->SetLineColor(colors[j]);
    h_RMS_dz_phi_vs_run[j]->SetLineWidth(2);
    h_RMS_dz_phi_vs_run[j]->SetMarkerStyle(markers[j]);
    h_RMS_dz_phi_vs_run[j]->SetMarkerColor(colors[j]);
    adjustmargins(c_RMS_dz_phi_vs_run);
    c_RMS_dz_phi_vs_run->cd();
    beautify(h_RMS_dz_phi_vs_run[j]);

    if(time_axis_format){
      timify(h_RMS_dz_phi_vs_run[j]);
    }

    arr_dz_phi->Add(h_RMS_dz_phi_vs_run[j]);

    // at the last file re-loop
    if(j==nDirs_-1){

      auto theMax = getMaximumFromArray(arr_dz_phi);
      
      for(Int_t k=0; k< nDirs_; k++){
	h_RMS_dz_phi_vs_run[k]->GetYaxis()->SetRangeUser(-theMax*0.45,theMax*1.30);
	if(k==0){
	  h_RMS_dz_phi_vs_run[k]->Draw("L");
	} else {
	  h_RMS_dz_phi_vs_run[k]->Draw("Lsame");
	}
      }
      my_lego->Draw("same");
      TH1F* theConst = DrawConstant(h_RMS_dz_phi_vs_run[j],1,0.);
      theConst->Draw("same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_RMS_dz_phi_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // dz vs eta
    // *************************************
    g_dz_eta_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dzEtaMeans_[LegLabels[j]])[0]));
    g_dz_eta_hi_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dzEtaHi_[LegLabels[j]])[0]));
    g_dz_eta_lo_vs_run[j] = new TGraph(x_ticks.size(),&(x_ticks[0]),&((dzEtaLo_[LegLabels[j]])[0]));

    adjustmargins(c_dz_eta_vs_run);
    c_dz_eta_vs_run->cd();
    g_dz_eta_vs_run[j]->SetMarkerStyle(markers[j]);
    g_dz_eta_vs_run[j]->SetMarkerColor(colors[j]);
    g_dz_eta_vs_run[j]->SetLineColor(colors[j]);
    g_dz_eta_hi_vs_run[j]->SetLineColor(colors[j]);
    g_dz_eta_lo_vs_run[j]->SetLineColor(colors[j]);
    beautify(g_dz_eta_vs_run[j]);

    g_dz_eta_vs_run[j]->SetName(Form("g_bias_dz_eta_%s",LegLabels[j].Data()));
    g_dz_eta_vs_run[j]->SetTitle(Form("Bias of d_{z}(#eta) vs %s",theType.Data()));
    g_dz_eta_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_dz_eta_vs_run[j]->GetYaxis()->SetTitle("#LT d_{z}(#eta) #GT [#mum]");

    g_dz_eta_vs_run[j]->GetYaxis()->SetRangeUser(-100,100);
    if(j==0){
      g_dz_eta_vs_run[j]->Draw("AP");
    } else {
      g_dz_eta_vs_run[j]->Draw("Psame");
    }
    g_dz_eta_hi_vs_run[j]->Draw("Lsame");
    g_dz_eta_lo_vs_run[j]->Draw("Lsame");

    if(time_axis_format){
      timify(g_dz_eta_vs_run[j]);
      timify(g_dz_eta_hi_vs_run[j]);
      timify(g_dz_eta_lo_vs_run[j]);
    }

    if(j==nDirs_-1){ 
      my_lego->Draw("same");
      TH1F* theZero = DrawConstantGraph(g_dz_eta_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }
	
    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_dz_eta_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // mean only
    adjustmargins(c_mean_dz_eta_vs_run);
    c_mean_dz_eta_vs_run->cd();
    gprime_dz_eta_vs_run[j] =  (TGraph*)g_dz_eta_vs_run[j]->Clone();
    if(j==0){
      gprime_dz_eta_vs_run[j]->GetYaxis()->SetRangeUser(-20.,20.);
      gprime_dz_eta_vs_run[j]->Draw("APL");
    } else {
      gprime_dz_eta_vs_run[j]->Draw("PLsame");
    }
    
    if(j==nDirs_-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(gprime_dz_eta_vs_run[j],2,0.);
      theZero->Draw("E1same");

      current_pad = static_cast<TPad*>(gPad);
      cmsPrel(current_pad);

      superImposeIOVBoundaries(c_mean_dz_eta_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);
    }

    // scatter or RMS TH1
    h_RMS_dz_eta_vs_run[j] = new TH1F(Form("h_RMS_dz_eta_%s",LegLabels[j].Data()),Form("scatter of d_{xy}(#eta) vs %s;%s;maximum scatter of d_{z}(#eta) [#mum]",theType.Data(),theTypeLabel.Data()),x_ticks.size()-1,&(x_ticks[0]));
    h_RMS_dz_eta_vs_run[j]->SetStats(kFALSE);

    bincounter=0;
    for(const auto &tick : x_ticks ){
      bincounter++;
      h_RMS_dz_eta_vs_run[j]->SetBinContent(bincounter,std::abs(dzEtaHi_[LegLabels[j]][bincounter-1]-dzEtaLo_[LegLabels[j]][bincounter-1]));
      h_RMS_dz_eta_vs_run[j]->SetBinError(bincounter,0.01);
    }

    h_RMS_dz_eta_vs_run[j]->SetLineColor(colors[j]);
    h_RMS_dz_eta_vs_run[j]->SetLineWidth(2);
    h_RMS_dz_eta_vs_run[j]->SetMarkerStyle(markers[j]);
    h_RMS_dz_eta_vs_run[j]->SetMarkerColor(colors[j]);
    adjustmargins(c_RMS_dz_eta_vs_run);
    c_RMS_dz_eta_vs_run->cd();
    beautify(h_RMS_dz_eta_vs_run[j]);

    if(time_axis_format){
      timify(h_RMS_dz_eta_vs_run[j]);
    }

    arr_dz_eta->Add(h_RMS_dz_eta_vs_run[j]);

    // at the last file re-loop
    if(j==nDirs_-1){

      auto theMax = getMaximumFromArray(arr_dz_eta);

      for(Int_t k=0; k< nDirs_; k++){
	h_RMS_dz_eta_vs_run[k]->GetYaxis()->SetRangeUser(-theMax*0.45,theMax*1.30);
	if(k==0){
	  h_RMS_dz_eta_vs_run[k]->Draw("L");
	} else {
	  h_RMS_dz_eta_vs_run[k]->Draw("Lsame");
	}
      }
      my_lego->Draw("same");
      TH1F* theConst = DrawConstant(h_RMS_dz_eta_vs_run[j],1,0.);
      theConst->Draw("same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_RMS_dz_eta_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);
    
    // *************************************
    // Integrated bias dxy scatter plots
    // *************************************

    h2_scatter_dxy_vs_run[j] = new TH2F(Form("h2_scatter_dxy_%s",LegLabels[j].Data()),Form("scatter of d_{xy} vs %s;%s;d_{xy} [cm]",theType.Data(),theTypeLabel.Data()),x_ticks.size()-1,&(x_ticks[0]),dxyVect[LegLabels[j]][0].get_n_bins(),dxyVect[LegLabels[j]][0].get_y_min(),dxyVect[LegLabels[j]][0].get_y_max());
    h2_scatter_dxy_vs_run[j]->SetStats(kFALSE);

    for(unsigned int runindex=0;runindex<x_ticks.size();runindex++){
      for(unsigned int binindex=0;binindex<dxyVect[LegLabels[j]][runindex].get_n_bins();binindex++){
	h2_scatter_dxy_vs_run[j]->SetBinContent(runindex+1,binindex+1,dxyVect[LegLabels[j]][runindex].get_bin_contents().at(binindex)/dxyVect[LegLabels[j]][runindex].get_integral());
      }
    }
  
    //Scatter_dxy_vs_run->cd();
    h2_scatter_dxy_vs_run[j]->SetFillColorAlpha(colors[j],0.3);
    h2_scatter_dxy_vs_run[j]->SetMarkerColor(colors[j]);
    h2_scatter_dxy_vs_run[j]->SetLineColor(colors[j]);
    h2_scatter_dxy_vs_run[j]->SetMarkerStyle(markers[j]);

    auto h_dxypfx_tmp = (TProfile*)(((TH2F*)h2_scatter_dxy_vs_run[j])->ProfileX(Form("_apfx_%i",j),1,-1,"o"));
    h_dxypfx_tmp->SetName(TString(h2_scatter_dxy_vs_run[j]->GetName())+"_pfx");
    h_dxypfx_tmp->SetStats(kFALSE);
    h_dxypfx_tmp->SetMarkerColor(colors[j]);
    h_dxypfx_tmp->SetLineColor(colors[j]);
    h_dxypfx_tmp->SetLineWidth(2); 
    h_dxypfx_tmp->SetMarkerSize(1); 
    h_dxypfx_tmp->SetMarkerStyle(markers[j]);

    beautify(h2_scatter_dxy_vs_run[j]);
    beautify(h_dxypfx_tmp);

    Scatter_dxy_vs_run->cd(j+1);
    adjustmargins(Scatter_dxy_vs_run->cd(j+1));
    //h_dxypfx_tmp->GetYaxis()->SetRangeUser(-0.01,0.01);
    //h2_scatter_dxy_vs_run[j]->GetYaxis()->SetRangeUser(-0.5,0.5);
    h2_scatter_dxy_vs_run[j]->Draw("colz");  
    h_dxypfx_tmp->Draw("same");

    // *************************************
    // Integrated bias dz scatter plots
    // *************************************

    h2_scatter_dz_vs_run[j] = new TH2F(Form("h2_scatter_dz_%s",LegLabels[j].Data()),Form("scatter of d_{z} vs %s;%s;d_{z} [cm]",theType.Data(),theTypeLabel.Data()),x_ticks.size()-1,&(x_ticks[0]),dzVect[LegLabels[j]][0].get_n_bins(),dzVect[LegLabels[j]][0].get_y_min(),dzVect[LegLabels[j]][0].get_y_max());
    h2_scatter_dz_vs_run[j]->SetStats(kFALSE);

    for(unsigned int runindex=0;runindex<x_ticks.size();runindex++){
      for(unsigned int binindex=0;binindex<dzVect[LegLabels[j]][runindex].get_n_bins();binindex++){
	h2_scatter_dz_vs_run[j]->SetBinContent(runindex+1,binindex+1,dzVect[LegLabels[j]][runindex].get_bin_contents().at(binindex)/dzVect[LegLabels[j]][runindex].get_integral());
      }
    }
    
    //Scatter_dz_vs_run->cd();
    h2_scatter_dz_vs_run[j]->SetFillColorAlpha(colors[j],0.3);
    h2_scatter_dz_vs_run[j]->SetMarkerColor(colors[j]);
    h2_scatter_dz_vs_run[j]->SetLineColor(colors[j]);
    h2_scatter_dz_vs_run[j]->SetMarkerStyle(markers[j]);

    auto h_dzpfx_tmp = (TProfile*)(((TH2F*)h2_scatter_dz_vs_run[j])->ProfileX(Form("_apfx_%i",j),1,-1,"o"));
    h_dzpfx_tmp->SetName(TString(h2_scatter_dz_vs_run[j]->GetName())+"_pfx");
    h_dzpfx_tmp->SetStats(kFALSE);
    h_dzpfx_tmp->SetMarkerColor(colors[j]);
    h_dzpfx_tmp->SetLineColor(colors[j]);
    h_dzpfx_tmp->SetLineWidth(2); 
    h_dzpfx_tmp->SetMarkerSize(1); 
    h_dzpfx_tmp->SetMarkerStyle(markers[j]);

    beautify(h2_scatter_dz_vs_run[j]);
    beautify(h_dzpfx_tmp);

    Scatter_dz_vs_run->cd(j+1);
    adjustmargins(Scatter_dz_vs_run->cd(j+1));
    //h_dzpfx_tmp->GetYaxis()->SetRangeUser(-0.01,0.01);
    //h2_scatter_dz_vs_run[j]->GetYaxis()->SetRangeUser(-0.5,0.5); 
    h2_scatter_dz_vs_run[j]->Draw("colz");  
    h_dzpfx_tmp->Draw("same");

 
  }
  
  // delete the array for the maxima
  delete arr_dxy_phi;
  delete arr_dz_phi;
  delete arr_dxy_eta;
  delete arr_dz_eta;

  TString append;
  if(lumi_axis_format){ 
    append = "lumi";
  } else{ 
    if(time_axis_format){
      append = "date";
    } else {
      append = "run";
    }
  }   

  c_dxy_phi_vs_run->SaveAs("dxy_phi_vs_"+append+".pdf");
  c_dxy_phi_vs_run->SaveAs("dxy_phi_vs_"+append+".png");
  c_dxy_phi_vs_run->SaveAs("dxy_phi_vs_"+append+".root");

  c_dxy_eta_vs_run->SaveAs("dxy_eta_vs_"+append+".pdf");
  c_dxy_eta_vs_run->SaveAs("dxy_eta_vs_"+append+".png");
  c_dxy_eta_vs_run->SaveAs("dxy_eta_vs_"+append+".root");

  c_dz_phi_vs_run->SaveAs("dz_phi_vs_"+append+".pdf");
  c_dz_phi_vs_run->SaveAs("dz_phi_vs_"+append+".png");
  c_dz_phi_vs_run->SaveAs("dz_phi_vs_"+append+".root");

  c_dz_eta_vs_run->SaveAs("dz_eta_vs_"+append+".pdf");
  c_dz_eta_vs_run->SaveAs("dz_eta_vs_"+append+".png");
  c_dz_eta_vs_run->SaveAs("dz_eta_vs_"+append+".root");

  // mean

  c_mean_dxy_phi_vs_run->SaveAs("mean_dxy_phi_vs_"+append+".pdf");
  c_mean_dxy_phi_vs_run->SaveAs("mean_dxy_phi_vs_"+append+".png");

  c_mean_dxy_eta_vs_run->SaveAs("mean_dxy_eta_vs_"+append+".pdf");
  c_mean_dxy_eta_vs_run->SaveAs("mean_dxy_eta_vs_"+append+".png");

  c_mean_dz_phi_vs_run->SaveAs("mean_dz_phi_vs_"+append+".pdf");
  c_mean_dz_phi_vs_run->SaveAs("mean_dz_phi_vs_"+append+".png");

  c_mean_dz_eta_vs_run->SaveAs("mean_dz_eta_vs_"+append+".pdf");
  c_mean_dz_eta_vs_run->SaveAs("mean_dz_eta_vs_"+append+".png");

  // RMS

  c_RMS_dxy_phi_vs_run->SaveAs("RMS_dxy_phi_vs_"+append+".pdf");
  c_RMS_dxy_phi_vs_run->SaveAs("RMS_dxy_phi_vs_"+append+".png");

  c_RMS_dxy_eta_vs_run->SaveAs("RMS_dxy_eta_vs_"+append+".pdf");
  c_RMS_dxy_eta_vs_run->SaveAs("RMS_dxy_eta_vs_"+append+".png");

  c_RMS_dz_phi_vs_run->SaveAs("RMS_dz_phi_vs_"+append+".pdf");
  c_RMS_dz_phi_vs_run->SaveAs("RMS_dz_phi_vs_"+append+".png");

  c_RMS_dz_eta_vs_run->SaveAs("RMS_dz_eta_vs_"+append+".pdf");
  c_RMS_dz_eta_vs_run->SaveAs("RMS_dz_eta_vs_"+append+".png");

  // scatter

  Scatter_dxy_vs_run->SaveAs("Scatter_dxy_vs_"+append+".pdf");
  Scatter_dxy_vs_run->SaveAs("Scatter_dxy_vs_"+append+".png");

  Scatter_dz_vs_run->SaveAs("Scatter_dz_vs_"+append+".pdf");
  Scatter_dz_vs_run->SaveAs("Scatter_dz_vs_"+append+".png");


  // do all the deletes

  for(int iDir=0;iDir<nDirs_;iDir++){

   delete g_dxy_phi_vs_run[iDir];    
   delete gprime_dxy_phi_vs_run[iDir];    
   delete g_dxy_phi_hi_vs_run[iDir]; 
   delete g_dxy_phi_lo_vs_run[iDir]; 
                               
   delete g_dxy_eta_vs_run[iDir];    
   delete gprime_dxy_eta_vs_run[iDir];    
   delete g_dxy_eta_hi_vs_run[iDir]; 
   delete g_dxy_eta_lo_vs_run[iDir]; 
                               
   delete g_dz_phi_vs_run[iDir];     
   delete gprime_dz_phi_vs_run[iDir];     
   delete g_dz_phi_hi_vs_run[iDir];  
   delete g_dz_phi_lo_vs_run[iDir];  
                               
   delete g_dz_eta_vs_run[iDir];     
   delete gprime_dz_eta_vs_run[iDir];     
   delete g_dz_eta_hi_vs_run[iDir];  
   delete g_dz_eta_lo_vs_run[iDir];  
                                                             
   delete h_RMS_dxy_phi_vs_run[iDir];  
   delete h_RMS_dxy_eta_vs_run[iDir];  
   delete h_RMS_dz_phi_vs_run[iDir];   
   delete h_RMS_dz_eta_vs_run[iDir];   

  }

  // mv the run-by-run plots into the folders

  gSystem->mkdir("Biases");
  TString processline = ".! mv Bias*.p* ./Biases/";
  std::cout<<"Executing: \n"
	   <<processline<< "\n"<<std::endl;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
  processline.Clear();

  gSystem->mkdir("ResolutionsVsPt");
  processline = ".! mv ResolutionsVsPt*.p* ./ResolutionsVsPt/";
  std::cout<<"Executing: \n"
	   <<processline<< "\n"<<std::endl;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
  processline.Clear();

  gSystem->mkdir("Resolutions");
  processline = ".! mv Resolutions*.p* ./Resolutions/";
  std::cout<<"Executing: \n"
	   <<processline<< "\n"<<std::endl;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
  processline.Clear();

  gSystem->mkdir("Pulls");
  processline = ".! mv Pulls*.p* ./Pulls/";
  std::cout<<"Executing: \n"
	   <<processline<< "\n" <<std::endl;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
  processline.Clear();

  timer.Stop(); 	 
  timer.Print();

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
	
	if(theTitle.Contains("norm")){
	  //dBiasTrend[k][i]->GetYaxis()->SetRangeUser(std::min(-0.48,absmin[k]-safeDelta/2.),std::max(0.48,absmax[k]+safeDelta/2.));
	  dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,1.8);
	} else {
	  if(!theTitle.Contains("width")){
	    //dBiasTrend[k][i]->GetYaxis()->SetRangeUser(-theExtreme-(safeDelta/2.),theExtreme+(safeDelta/2.));
	    dBiasTrend[k][i]->GetYaxis()->SetRangeUser(-40.,40.);
	  } else {
	    // dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,theExtreme+(safeDelta/2.));
	    // if(theTitle.Contains("eta")) {
	    //   dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,500.);
	    // } else {
	    //   dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,200.);
	    // }
	    auto my_view = checkTheView(theTitle);

	    //std::cout<<" ----------------------------------> " << theTitle << " view: " << my_view <<  std::endl;

	    switch(my_view){
	    case pv::dxyphi  : 
	      dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,200.); 
	      break;
	    case pv::dzphi   : 
	      dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,400.);
	      break;
	    case pv::dxyeta  : 
	      dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,300.);
	      break;
	    case pv::dzeta   : 
	      dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,1.e3);
	      break;
	    case pv::pT : 
	      dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,200.);
	      break;
	    case pv::generic : 
	      dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,350.);
	      break;
	    default : 
	      dBiasTrend[k][i]->GetYaxis()->SetRangeUser(0.,300.);
	    }
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
	
	TH1F* theConst = DrawConstant(dBiasTrend[k][i],1,theC);
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
void setStyle(){
/*--------------------------------------------------------------------*/

  TGaxis::SetMaxDigits(6);
  
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

  //gStyle->SetPalette(kInvertedDarkBodyRadiator);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  /*
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  */

  Double_t stops[NRGBs] = {0.00, 0.01, 0.05, 0.09, 0.1};
  Double_t red[NRGBs]   = {1.00, 0.84, 0.61, 0.34, 0.00};
  Double_t green[NRGBs] = {1.00, 0.84, 0.61, 0.34, 0.00};
  Double_t blue[NRGBs]  = {1.00, 0.84, 0.61, 0.34, 0.00};

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

  float posX_ =  1-r - relPosX*(1-l-r);
  //  float posXCMS_ = 1-r -15*relPosX*(1-l-r);
  float posY_ =  1-t + 0.05; /// - relPosY*(1-t-b);

  latex->SetTextAlign(33);
  //latex->SetTextFont(61);
  //latex->DrawLatex(posXCMS_,posY_,"CMS");
  latex->SetTextFont(42); //22
  latex->DrawLatex(posX_,posY_,"CMS Internal (13 TeV)");
  //latex->DrawLatex(posX_,posY_,"CMS Preliminary (13 TeV)");
  //latex->DrawLatex(posX_,posY_,"CMS 2017 Work in progress (13 TeV)");
  
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
TH1F* DrawConstantGraph(TGraph *graph,Int_t iter,Double_t theConst)
/*--------------------------------------------------------------------*/
{ 
 
  Double_t xmin = graph->GetXaxis()->GetXmin(); //TMath::MinElement(graph->GetN(),graph->GetX());
  Double_t xmax = graph->GetXaxis()->GetXmax(); //TMath::MaxElement(graph->GetN(),graph->GetX()); 

  //std::cout<<xmin<<" : "<<xmax<<std::endl;

  TH1F *hzero = new TH1F(Form("hconst_%s_%i",graph->GetName(),iter),Form("hconst_%s_%i",graph->GetName(),iter),graph->GetN(),xmin,xmax);
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
unrolledHisto getUnrolledHisto(TH1F* hist)
/*--------------------------------------------------------------------*/
{
  
  /*
    Double_t y_min = hist->GetBinLowEdge(1);
    Double_t y_max = hist->GetBinLowEdge(hist->GetNbinsX()+1);
  */    

  Double_t y_min = -0.1;
  Double_t y_max = 0.1;

  std::vector<Double_t> contents;
  for (int j = 0; j < hist->GetNbinsX(); j++) {
    if(std::abs(hist->GetXaxis()->GetBinCenter(j))<=0.1) contents.push_back(hist->GetBinContent(j+1));
  }
  
  auto ret = unrolledHisto(y_min,y_max,contents.size(),contents);
  return ret;

}

/*--------------------------------------------------------------------*/
std::pair<std::pair<Double_t,Double_t>, Double_t> getBiases(TH1F* hist,bool useRMS_)
/*--------------------------------------------------------------------*/
{
  Double_t mean=0;
  Double_t rms=0;

  int nbins = hist->GetNbinsX();

  //extract median from histogram
  double *y   = new double[nbins];
  double *err = new double[nbins];
  for (int j = 0; j < nbins; j++) {
    y[j]   = hist->GetBinContent(j+1);
    err[j] = hist->GetBinError(j+1);
  }
  mean = TMath::Mean(nbins,y,err);
  rms =  TMath::RMS(nbins,y,err);

  Double_t max=hist->GetMaximum();
  Double_t min=hist->GetMinimum();

  /*
  for(Int_t i=1;i<=hist->GetNbinsX();i++){
    mean+=hist->GetBinContent(i);
  }

  mean = mean/hist->GetNbinsX();
  
  //std::pair<Double_t,Double_t> resultBounds = std::make_pair(min,max);
  */
  
  /*
    // in case one would like to use a pol0 fit
    hist->Fit("pol0","Q0+");
    TF1* f = (TF1*)hist->FindObject("pol0");
    f->SetLineColor(hist->GetLineColor());
    f->SetLineStyle(hist->GetLineStyle());
    mean = f->GetParameter(0);
  */

  std::pair<std::pair<Double_t,Double_t>, Double_t> result;
  std::pair<Double_t,Double_t> resultBounds;

  resultBounds = useRMS_ ? std::make_pair(mean-rms,mean+rms) :  std::make_pair(min,max)  ;

  result = make_pair(resultBounds,mean);
  
  return result;
   
}

/*--------------------------------------------------------------------*/
void beautify(TGraph *g){
/*--------------------------------------------------------------------*/
  g->GetXaxis()->SetLabelFont(42);
  g->GetYaxis()->SetLabelFont(42);
  g->GetYaxis()->SetLabelSize(.055);
  g->GetXaxis()->SetLabelSize(.055);
  g->GetYaxis()->SetTitleSize(.055);
  g->GetXaxis()->SetTitleSize(.055);
  g->GetXaxis()->SetTitleOffset(1.1);
  g->GetYaxis()->SetTitleOffset(0.9);
  g->GetXaxis()->SetTitleFont(42);
  g->GetYaxis()->SetTitleFont(42);
  g->GetXaxis()->CenterTitle(true);
  g->GetYaxis()->CenterTitle(true);
  g->GetXaxis()->SetNdivisions(505);
}

/*--------------------------------------------------------------------*/
void beautify(TH1 *h){
/*--------------------------------------------------------------------*/
  h->SetMinimum(0.);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(.055);
  h->GetXaxis()->SetLabelSize(.055);
  h->GetYaxis()->SetTitleSize(.055);
  h->GetXaxis()->SetTitleSize(.055);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->GetXaxis()->SetNdivisions(505);
}

/*--------------------------------------------------------------------*/
void adjustmargins(TCanvas *canv){
/*--------------------------------------------------------------------*/
  canv->cd()->SetBottomMargin(0.14);
  canv->cd()->SetLeftMargin(0.11);
  canv->cd()->SetRightMargin(0.01);
  canv->cd()->SetTopMargin(0.06);
}

/*--------------------------------------------------------------------*/
void adjustmargins(TVirtualPad *canv){
/*--------------------------------------------------------------------*/
  canv->SetBottomMargin(0.14);
  canv->SetLeftMargin(0.08);
  canv->SetRightMargin(0.08);
  canv->SetTopMargin(0.06);
}

/*--------------------------------------------------------------------*/
TH1F* checkTH1AndReturn(TFile *f,TString address){
/*--------------------------------------------------------------------*/
  TH1F* h = NULL; 
  if(f->GetListOfKeys()->Contains(address)){
    h = (TH1F*)f->Get(address);
  } 
  return h;
}

/*--------------------------------------------------------------------*/
pv::view checkTheView(const TString &toCheck){
/*--------------------------------------------------------------------*/
  if (toCheck.Contains("dxy")){
    if (toCheck.Contains("phi") || toCheck.Contains("ladder")){
      return pv::dxyphi;
    } else if (toCheck.Contains("eta") || toCheck.Contains("modZ")){
      return pv::dxyeta;
    } else {
      return pv::pT;
    }
  } else if (toCheck.Contains("dz")){
    if (toCheck.Contains("phi") || toCheck.Contains("ladder")){
      return pv::dzphi;
    } else if (toCheck.Contains("eta") || toCheck.Contains("modZ")){
      return pv::dzeta;
    } else {
      return pv::pT;
    }
  } else {
    return pv::generic;	
  } 
}

/*--------------------------------------------------------------------*/
template<typename T>
void timify(T *mgr)
/*--------------------------------------------------------------------*/
{
  mgr->GetXaxis()->SetTimeDisplay(1);
  mgr->GetXaxis()->SetNdivisions(510);
  mgr->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  mgr->GetXaxis()->SetTimeOffset(0,"gmt");
  mgr->GetXaxis()->SetLabelSize(.035);
}

/*--------------------------------------------------------------------*/
Double_t getMaximumFromArray(TObjArray *array)
/*--------------------------------------------------------------------*/
{

  Double_t theMaximum = (static_cast<TH1*>(array->At(0)))->GetMaximum();
  for(Int_t i = 0; i< array->GetSize(); i++){
    if( (static_cast<TH1*>(array->At(i)))->GetMaximum() > theMaximum){
      theMaximum = (static_cast<TH1*>(array->At(i)))->GetMaximum();
      //cout<<"i= "<<i<<" theMaximum="<<theMaximum<<endl;
    }
  }
  return theMaximum;
}

/*--------------------------------------------------------------------*/
void superImposeIOVBoundaries(TCanvas *c,bool lumi_axis_format,bool time_axis_format,const std::map<int,double> &lumiMapByRun,const std::map<int,TDatime>& timeMap)
/*--------------------------------------------------------------------*/
{
 
  /*
    1   296641       2017-11-07 12:31:49  bf4bdd66393bee0623b730e11f3db799674c4daf  Alignments   
    2   297179       2017-11-07 12:31:49  a5858bccd2e2c031a12aa6cbe38e8adaeca19331  Alignments   
    3   297281       2017-11-07 12:31:49  b423ed0f118c16d99cf92d12a1bd33bb1902f533  Alignments   
    4   298653       2017-11-07 12:31:49  3e598e771471289c59d755ef4ebb993e6c4b7890  Alignments   
    5   299277       2017-11-07 12:31:49  d8631cc034f3971131f7e7551610acf790192cd2  Alignments   
    6   299443       2017-11-07 12:31:49  58961d0e9322c4d7064b31e73d9edfe9a861e1f8  Alignments   
    7   300389       2017-11-07 12:31:49  7c0a2dcfa527e8e336ff6b8b198827a5fcb9c287  Alignments   
    8   301046       2017-11-07 12:31:49  4bf6bfbd7f8a202e8a6bba5e8730b27b00e5142b  Alignments   
    9   302131       2017-11-07 12:31:49  387e5f629f808464f5c69ab8bb7319576d81e768  Alignments   
    10  303790       2017-11-07 12:31:49  21b79cc2a163df4fbc6c059eed8c96b86f8a07ff  Alignments   
    11  304911       2017-11-07 12:31:49  44c7e1e991484d1490ac5ad6b734b6622f1e9d85  Alignments   
    12  305898       2017-11-28 16:03:38  64aae34a7273b374738ed4b9b20c332abf503ab4  Alignments   
  */
 
  // get the vector of runs in the lumiMap
  std::vector<int> vruns;
  for(auto const& imap: lumiMapByRun){
    vruns.push_back(imap.first);
    //std::cout<<" run:" << imap.first << " lumi: "<< imap.second << std::endl;
  }

  //std::vector<vint> truns;
  for(auto const& imap: timeMap){
    std::cout<<" run:" << imap.first << " time: "<< imap.second.Convert() << std::endl;
  }

  static const int nIOVs=12; //     1      2      3      4      5      6      7      8      9     10     11     12      
  int IOVboundaries[nIOVs]  = {296641,297179,297281,298653,299277,299443,300389,301046,302131,303790,304911,305898};
  TArrow* IOV_lines[nIOVs];
  c->cd();
  c->Update();

  TArrow* a_lines[nIOVs];
  for(Int_t IOV=0;IOV<nIOVs;IOV++){

    // check we are not in the RMS histogram to avoid first line
    if(IOVboundaries[IOV]<vruns.front() && ((TString)c->GetName()).Contains("RMS")) continue;
    int closestrun = pv::closest(vruns,IOVboundaries[IOV]); 

    if(lumi_axis_format){

      if(closestrun<0) break;
      //std::cout<< "natural boundary: " << IOVboundaries[IOV] << " closest:" << closestrun << std::endl;

      a_lines[IOV] = new TArrow(lumiMapByRun.at(closestrun),(c->GetUymin()),lumiMapByRun.at(closestrun),0.65*c->GetUymax(),0.5,"|>");
    } else if(time_axis_format){
      
      if(closestrun<0) break;
      std::cout<< "natural boundary: " << IOVboundaries[IOV] << " closest:" << closestrun << std::endl;
      a_lines[IOV] = new TArrow(timeMap.at(closestrun).Convert(),(c->GetUymin()),timeMap.at(closestrun).Convert(),0.65*c->GetUymax(),0.5,"|>");
    } else {
      a_lines[IOV] = new TArrow(IOVboundaries[IOV],(c->GetUymin()),IOVboundaries[IOV],0.65*c->GetUymax(),0.5,"|>"); //(c->GetUymin()+0.2*(c->GetUymax()-c->GetUymin()) ),0.5,"|>");
    }
    a_lines[IOV]->SetLineColor(kRed);
    a_lines[IOV]->SetLineStyle(9);
    a_lines[IOV]->SetLineWidth(1);
    a_lines[IOV]->Draw("same");
  }

  TPaveText* runnumbers[nIOVs];
  
  for(Int_t IOV=0;IOV<nIOVs;IOV++){

    if(IOVboundaries[IOV]<vruns.front() && ((TString)c->GetName()).Contains("RMS")) continue;
    int closestrun = pv::closest(vruns,IOVboundaries[IOV]);
    
    Int_t ix1;
    Int_t ix2;
    Int_t iw = gPad->GetWw();
    Int_t ih = gPad->GetWh();
    Double_t x1p,y1p,x2p,y2p;
    gPad->GetPadPar(x1p,y1p,x2p,y2p);
    ix1 = (Int_t)(iw*x1p);
    ix2 = (Int_t)(iw*x2p);
    Double_t wndc  = TMath::Min(1.,(Double_t)iw/(Double_t)ih);
    Double_t rw    = wndc/(Double_t)iw;
    Double_t x1ndc = (Double_t)ix1*rw;
    Double_t x2ndc = (Double_t)ix2*rw;
    Double_t rx1,ry1,rx2,ry2;
    gPad->GetRange(rx1,ry1,rx2,ry2);
    Double_t rx = (x2ndc-x1ndc)/(rx2-rx1);
    Double_t _sx;
    if(lumi_axis_format){
      if(closestrun<0) break; 
      _sx = rx*(lumiMapByRun.at(closestrun)-rx1)+x1ndc; //-0.05;
    } else if(time_axis_format){
      if(closestrun<0) break; 
      _sx = rx*(timeMap.at(closestrun).Convert()-rx1)+x1ndc; //-0.05;
    } else {
      _sx = rx*(IOVboundaries[IOV]-rx1)+x1ndc; //-0.05
    }
    Double_t _dx = _sx+0.05;
    
    Int_t index;
    if(IOV<5) 
      index=IOV;
    else{
      index=IOV-5;
    }
    
    runnumbers[IOV] = new TPaveText(_sx+0.001,0.14+(0.03*index),_dx,(0.17+0.03*index),"blNDC");
    //runnumbers[IOV]->SetTextAlign(11);
    TText *textRun = runnumbers[IOV]->AddText(Form("%i",int(IOVboundaries[IOV])));
    textRun->SetTextSize(0.028);
    textRun->SetTextColor(kRed);
    runnumbers[IOV]->SetFillColor(10);
    runnumbers[IOV]->SetLineColor(kRed);
    runnumbers[IOV]->SetBorderSize(1);
    runnumbers[IOV]->SetLineWidth(2);
    runnumbers[IOV]->SetTextColor(kRed);
    runnumbers[IOV]->SetTextFont(42);
    runnumbers[IOV]->Draw("same");
  }


}
