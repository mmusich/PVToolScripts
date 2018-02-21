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
#include "TKey.h"
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
#include <functional>
#include <iterator>
#include <fstream>
#include <sstream>

const size_t nWorkers = 10;

typedef std::map<TString, std::vector<double> > resolutionTrend; 

namespace pv {
  
  // brief method to find first value that doesn not compare lett
  int closest(std::vector<int> const& vec, int value) {
    auto const it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end()) { return -1; }
    return *it;
  }

  const Int_t markers[8] = {kFullSquare,kFullCircle,kFullTriangleDown,kOpenSquare,kOpenCircle,kFullTriangleUp,kOpenTriangleDown,kOpenTriangleUp};
  const Int_t colors[8]  = {kBlack,kBlue,kRed,kGreen+2,kOrange,kMagenta,kCyan,kViolet};

  struct resolutions {
    
    // contructor
    resolutions(double resol50,double resol100,double resol200,double resol400){
      m_resol50  = resol50; 
      m_resol100 = resol100;
      m_resol200 = resol200;
      m_resol400 = resol400;
    }
    
    // empty constructor
    resolutions(){
      init();
    }    

    void init(){
      m_resol50  = -999.;
      m_resol100 = -999.;
      m_resol200 = -999.;
      m_resol400 = -999.;
    }

    double getResol50(){ return m_resol50;}
    double getResol100(){ return m_resol100;}
    double getResol200(){ return m_resol200;}
    double getResol400(){ return m_resol400;}    

  private:
    double m_resol50;
    double m_resol100;
    double m_resol200;
    double m_resol400;
  };

}

// auxilliary struct to be returned by the functor
struct outTrends {

  int m_index;
  double m_lumiSoFar;
  std::vector<double> m_runs;
  std::vector<double> m_lumiByRun;
  std::map<int,double> m_lumiMapByRun; 
  
  resolutionTrend m_xResolTrend200;
  resolutionTrend m_yResolTrend200;
  resolutionTrend m_zResolTrend200;

  resolutionTrend m_xResolTrend400;
  resolutionTrend m_yResolTrend400;
  resolutionTrend m_zResolTrend400;


  void init(){
    m_index=-1;
    m_lumiSoFar=0.;
    m_runs.clear();
    m_lumiByRun.clear();
    m_lumiMapByRun.clear();
    
    m_xResolTrend200.clear();
    m_yResolTrend200.clear();
    m_zResolTrend200.clear();

    m_xResolTrend400.clear();
    m_yResolTrend400.clear();
    m_zResolTrend400.clear();

  }
};

///////////////////////////////////
//
//  Forward declarations
//
///////////////////////////////////

void MultiRunPVResolution_withRef(TString namesandlabels="",bool lumi_axis_format=false,bool time_axis_format=false,bool useRMS=true,TString MCRef="");
outTrends doStuff(size_t iter,std::vector<int> intersection,const Int_t nDirs_,const Int_t nMCDirs_,const char* dirs[10], TString LegLabels[10], bool useRMS);

void arrangeOutCanvas(TCanvas *canv,
		      TH1F* m_11Trend[100],
		      TH1F* m_12Trend[100],
		      TH1F* m_13Trend[100],
		      Int_t nFiles, 
		      TString LegLabels[10],
		      unsigned int theRun);

std::vector<int> list_files(const char *dirname=".", const char *ext=".root");
void setStyle();
pv::resolutions getResolutions(TH1F* hist);
TH1F* DrawConstant(TH1F *hist,Int_t iter,Double_t theConst);
TH1F* DrawConstantWithErr(TH1F *hist,Int_t iter,Double_t theConst);
TH1F* DrawConstantGraph(TGraph *graph,Int_t iter,Double_t theConst);
TH1F* checkTH1AndReturn(TFile *f,TString address);
void MakeNiceTrendPlotStyle(TH1 *hist,Int_t color);
void cmsPrel(TPad* pad,size_t ipads=1);
void makeNewXAxis (TH1 *h);
void beautify(TGraph *g);
void beautify(TH1 *h);
void adjustmargins(TCanvas *canv);
void adjustmargins(TVirtualPad*canv);
template<typename T> void timify(T *mgr);
Double_t getMaximumFromArray(TObjArray *array);
void superImposeIOVBoundaries(TCanvas *c,bool lumi_axis_format,bool time_axis_format,const std::map<int,double> &lumiMapByRun,const std::map<int,TDatime>& timeMap);

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
void MultiRunPVResolution_withRef(TString namesandlabels,bool lumi_axis_format,bool time_axis_format,bool useRMS,TString MCRef){
  
  TStopwatch timer; 	 
  timer.Start();

  using namespace std::placeholders;  // for _1, _2, _3...
  gROOT->ProcessLine("gErrorIgnoreLevel = kError;"); 	

  ROOT::EnableThreadSafety();
  TH1::AddDirectory(kFALSE);

  // consistency check, we cannot do plot vs lumi if time_axis
  if(lumi_axis_format && time_axis_format){
    std::cout<<"##########################################################################################"<<std::endl;
    std::cout<<"msg-i: MultiRunPVResolution(): you're requesting both summary vs lumi and vs time, "<< std::endl;
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

  //std::ofstream outfile ("lumiByRun.txt"); 
  std::ofstream outfile ("log.txt"); 
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

  TList *MCDirList   = new TList();
  TList *MCLabelList = new TList();

  if(MCRef.Sizeof()!=0){  
    TObjArray *nameandlabelpairs = MCRef.Tokenize(",");
    for (Int_t i = 0; i < nameandlabelpairs->GetEntries(); ++i) {
      TObjArray *aFileLegPair = TString(nameandlabelpairs->At(i)->GetName()).Tokenize("=");
      
      if(aFileLegPair->GetEntries() == 2) {
	MCDirList->Add(aFileLegPair->At(0)); 
	MCLabelList->Add(aFileLegPair->At(1));
      }
      else {
	std::cout << "Please give the MC reference file names and legend entries in the following form:\n" 
		  << " filename1=legendentry1,filename2=legendentry2\n"; 
      }       
    }
  }

  const Int_t nDirs_   = DirList->GetSize();
  const Int_t nMCDirs_ = MCDirList->GetSize();
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

  // Filling data Info
  for(Int_t j=0; j < nDirs_; j++) {
    
    // Retrieve labels
    TObjString* legend = (TObjString*)LabelList->At(j);
    TObjString* dir    = (TObjString*)DirList->At(j);
    LegLabels[j] = legend->String();
    dirs[j] = (dir->String()).Data();
    cout<<"MultiRunPVResolution(): label["<<j<<"]"<<LegLabels[j]<<endl;
    
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

  // Filling MC info
  for(Int_t jMC=0; jMC< nMCDirs_;jMC++){
    TObjString* legend = (TObjString*)MCLabelList->At(jMC);
    TObjString* dir    = (TObjString*)MCDirList->At(jMC);
    LegLabels[nDirs_+jMC] = legend->String();
    dirs[nDirs_+jMC] = (dir->String()).Data();
  }

  // debug only
  for(UInt_t index=0;index<intersection.size();index++){
    std::cout<<index<<" "<<intersection[index]<<std::endl;
  }

  // book the vector of values
  resolutionTrend xResolTrend200_;
  resolutionTrend yResolTrend200_;
  resolutionTrend zResolTrend200_;

  resolutionTrend xResolTrend400_;
  resolutionTrend yResolTrend400_;
  resolutionTrend zResolTrend400_;

  double lumiSoFar=0.0;
  
  std::cout<<" pre do-stuff: " << runs.size() << std::endl;
  
  //we should use std::bind to create a functor and then pass it to the procPool
  auto f_doStuff = std::bind(doStuff,_1,intersection,nDirs_,nMCDirs_,dirs,LegLabels,useRMS);
  
  TProcPool procPool(nWorkers);
  std::vector<size_t> range(nWorkers);
  std::iota(range.begin(),range.end(),0);
  //procPool.Map([&f_doStuff](size_t a) { f_doStuff(a); },{1,2,3});
  auto extracts = procPool.Map(f_doStuff,range);
  
  // sort the extracts according to the global index
  std::sort(extracts.begin(), extracts.end(), 
	    [](const outTrends & a, const outTrends & b) -> bool
	    { 
	      return a.m_index < b.m_index; 
	    });
  
  // re-assemble everything together
  for (auto extractedTrend : extracts){
    std::cout << "lumiSoFar: " << lumiSoFar <<"/fb" << std::endl;

    runs.insert(std::end(runs),std::begin(extractedTrend.m_runs), std::end(extractedTrend.m_runs));

    // luminosity needs a different treatment
    // we need to re-sum the luminosity so far

    for (const auto &run : extractedTrend.m_runs){

      std::cout<< run << " " << lumiSoFar+extractedTrend.m_lumiMapByRun[run] << std::endl;
      lumiByRun.push_back(lumiSoFar+extractedTrend.m_lumiMapByRun[run]);
      lumiMapByRun[run]=(lumiSoFar+extractedTrend.m_lumiMapByRun[run]);

    }

    lumiSoFar += (extractedTrend.m_lumiSoFar/1000.);
 
    for(const auto &label : LegLabels){

       // at sumPt = 200 GeV
      xResolTrend200_[label].insert(std::end(xResolTrend200_[label]), std::begin(extractedTrend.m_xResolTrend200[label]), std::end(extractedTrend.m_xResolTrend200[label]));
      yResolTrend200_[label].insert(std::end(yResolTrend200_[label]), std::begin(extractedTrend.m_yResolTrend200[label]), std::end(extractedTrend.m_yResolTrend200[label]));
      zResolTrend200_[label].insert(std::end(zResolTrend200_[label]), std::begin(extractedTrend.m_zResolTrend200[label]), std::end(extractedTrend.m_zResolTrend200[label]));

      // at sumPt = 400 GeV
      xResolTrend400_[label].insert(std::end(xResolTrend400_[label]), std::begin(extractedTrend.m_xResolTrend400[label]), std::end(extractedTrend.m_xResolTrend400[label]));
      yResolTrend400_[label].insert(std::end(yResolTrend400_[label]), std::begin(extractedTrend.m_yResolTrend400[label]), std::end(extractedTrend.m_yResolTrend400[label]));
      zResolTrend400_[label].insert(std::end(zResolTrend400_[label]), std::begin(extractedTrend.m_zResolTrend400[label]), std::end(extractedTrend.m_zResolTrend400[label]));

    }
  }

  TCanvas *c_xResol200_vs_run = new TCanvas("c_xResol200_vs_run","x vertex resolution vs run number",2000,800);
  TCanvas *c_yResol200_vs_run = new TCanvas("c_yResol200_vs_run","y vertex resolution vs run number",2000,800);
  TCanvas *c_zResol200_vs_run = new TCanvas("c_zResol200_vs_run" ,"z vertex resolution vs run number",2000,800);

  TCanvas *c_xResol400_vs_run = new TCanvas("c_xResol400_vs_run","x vertex resolution vs run number",2000,800);
  TCanvas *c_yResol400_vs_run = new TCanvas("c_yResol400_vs_run","y vertex resolution vs run number",2000,800);
  TCanvas *c_zResol400_vs_run = new TCanvas("c_zResol400_vs_run" ,"z vertex resolution vs run number",2000,800);

  TGraph *g_xResol200_vs_run[(nDirs_+nMCDirs_)];
  TGraph *g_yResol200_vs_run[(nDirs_+nMCDirs_)];
  TGraph *g_zResol200_vs_run[(nDirs_+nMCDirs_)];

  TGraph *g_xResol400_vs_run[(nDirs_+nMCDirs_)];
  TGraph *g_yResol400_vs_run[(nDirs_+nMCDirs_)];
  TGraph *g_zResol400_vs_run[(nDirs_+nMCDirs_)];

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
      theTypeLabel="UTC date";
      for(const auto &run : runs){
	x_ticks.push_back(times[run].Convert());
      }
    }
  }
  
  TLegend *my_lego = new TLegend(0.75,0.83,0.95,0.93);
  //my_lego-> SetNColumns(2);
  my_lego->SetFillColor(10);
  my_lego->SetTextSize(0.042);
  my_lego->SetTextFont(42);
  my_lego->SetFillColor(10);
  my_lego->SetLineColor(10);
  my_lego->SetShadowColor(10);

  TLatex t1;
  t1.SetTextAlign(21);
  t1.SetTextSize(0.05);
  t1.SetTextFont(42);

  for(Int_t j=0; j < (nDirs_+nMCDirs_); j++) {

    // check on the sanity
    std::cout<<"x_ticks.size()= "<<x_ticks.size()<<"d xyPhiMeans_[LegLabels["<<j<<"]].size()="<<xResolTrend200_[LegLabels[j]].size()<<std::endl;

    // *************************************
    // x resolution @ sumpT=200
    // *************************************
    
    g_xResol200_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((xResolTrend200_[LegLabels[j]])[0]));

    adjustmargins(c_xResol200_vs_run);
    c_xResol200_vs_run->cd();
    g_xResol200_vs_run[j]->SetMarkerStyle(pv::markers[j]);
    g_xResol200_vs_run[j]->SetMarkerColor(pv::colors[j]);
    g_xResol200_vs_run[j]->SetLineColor(pv::colors[j]);
    if(j>=nDirs_)  g_xResol200_vs_run[j]->SetLineWidth(3);

    g_xResol200_vs_run[j]->SetName(Form("g_xResol200_%s",LegLabels[j].Data()));
    g_xResol200_vs_run[j]->SetTitle(Form("Primary Vertex x-Resolution vs %s",theType.Data()));
    g_xResol200_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_xResol200_vs_run[j]->GetYaxis()->SetTitle("Primary Vertex x-Resolution [#mum]");
    g_xResol200_vs_run[j]->GetYaxis()->SetRangeUser(-10.,50.);
    beautify(g_xResol200_vs_run[j]);
 
    j<nDirs_ ? my_lego->AddEntry(g_xResol200_vs_run[j],LegLabels[j],"PL") : my_lego->AddEntry(g_xResol200_vs_run[j],LegLabels[j],"L") ;

    if(j==0){
      g_xResol200_vs_run[j]->Draw("AP");
      t1.DrawLatexNDC(0.20,0.85,"#sump_{T} = 200 GeV");
    } else {
      j<nDirs_ ? g_xResol200_vs_run[j]->Draw("Psame") : g_xResol200_vs_run[j]->Draw("Lsame") ;
    }

    if(time_axis_format){
      timify(g_xResol200_vs_run[j]);
    }

    if(j==(nDirs_+nMCDirs_)-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(g_xResol200_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    TPad *current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_xResol200_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // y resolution @ sumpT=200
    // *************************************
    
    g_yResol200_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((yResolTrend200_[LegLabels[j]])[0]));

    adjustmargins(c_yResol200_vs_run);
    c_yResol200_vs_run->cd();
    g_yResol200_vs_run[j]->SetMarkerStyle(pv::markers[j]);
    g_yResol200_vs_run[j]->SetMarkerColor(pv::colors[j]);
    g_yResol200_vs_run[j]->SetLineColor(pv::colors[j]);
    if(j>=nDirs_)  g_yResol200_vs_run[j]->SetLineWidth(3);

    g_yResol200_vs_run[j]->SetName(Form("g_yResol200_%s",LegLabels[j].Data()));
    g_yResol200_vs_run[j]->SetTitle(Form("Primary Vertex y-Resolution vs %s",theType.Data()));
    g_yResol200_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_yResol200_vs_run[j]->GetYaxis()->SetTitle("Primary Vertex y-Resolution [#mum]");
    g_yResol200_vs_run[j]->GetYaxis()->SetRangeUser(-10.,50.);
    beautify(g_yResol200_vs_run[j]);
 
    if(j==0){
      g_yResol200_vs_run[j]->Draw("AP");
      t1.DrawLatexNDC(0.20,0.85,"#sump_{T} = 200 GeV");
    } else {
      j<nDirs_ ? g_yResol200_vs_run[j]->Draw("Psame") : g_yResol200_vs_run[j]->Draw("Lsame");
    }

    if(time_axis_format){
      timify(g_yResol200_vs_run[j]);
    }

    if(j==(nDirs_+nMCDirs_)-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(g_yResol200_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_yResol200_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // z resolution @ sumpT=200
    // *************************************
    
    g_zResol200_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((zResolTrend200_[LegLabels[j]])[0]));

    adjustmargins(c_zResol200_vs_run);
    c_zResol200_vs_run->cd();
    g_zResol200_vs_run[j]->SetMarkerStyle(pv::markers[j]);
    g_zResol200_vs_run[j]->SetMarkerColor(pv::colors[j]);
    g_zResol200_vs_run[j]->SetLineColor(pv::colors[j]);
    if(j>=nDirs_)  g_zResol200_vs_run[j]->SetLineWidth(3);

    g_zResol200_vs_run[j]->SetName(Form("g_zResol200_%s",LegLabels[j].Data()));
    g_zResol200_vs_run[j]->SetTitle(Form("Primary Vertex z-Resolution vs %s",theType.Data()));
    g_zResol200_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_zResol200_vs_run[j]->GetYaxis()->SetTitle("Primary Vertex z-Resolution [#mum]");
    g_zResol200_vs_run[j]->GetYaxis()->SetRangeUser(-10.,50.);
    beautify(g_zResol200_vs_run[j]);
 
    if(j==0){
      g_zResol200_vs_run[j]->Draw("AP");
      t1.DrawLatexNDC(0.20,0.85,"#sump_{T} = 200 GeV");
    } else {
      j<nDirs_ ?  g_zResol200_vs_run[j]->Draw("Psame") : g_zResol200_vs_run[j]->Draw("Lsame");
    }

    if(time_axis_format){
      timify(g_zResol200_vs_run[j]);
    }

    if(j==(nDirs_+nMCDirs_)-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(g_zResol200_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_zResol200_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // x resolution @ sumpT=400
    // *************************************
    
    g_xResol400_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((xResolTrend400_[LegLabels[j]])[0]));

    adjustmargins(c_xResol400_vs_run);
    c_xResol400_vs_run->cd();
    g_xResol400_vs_run[j]->SetMarkerStyle(pv::markers[j]);
    g_xResol400_vs_run[j]->SetMarkerColor(pv::colors[j]);
    g_xResol400_vs_run[j]->SetLineColor(pv::colors[j]);
    if(j>=nDirs_)  g_xResol400_vs_run[j]->SetLineWidth(3);

    g_xResol400_vs_run[j]->SetName(Form("g_xResol400_%s",LegLabels[j].Data()));
    g_xResol400_vs_run[j]->SetTitle(Form("Primary Vertex x-Resolution vs %s",theType.Data()));
    g_xResol400_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_xResol400_vs_run[j]->GetYaxis()->SetTitle("Primary Vertex x-Resolution [#mum]");
    g_xResol400_vs_run[j]->GetYaxis()->SetRangeUser(-10.,50.);
    beautify(g_xResol400_vs_run[j]);
 
    if(j==0){
      g_xResol400_vs_run[j]->Draw("AP");
      t1.DrawLatexNDC(0.20,0.85,"#sump_{T} = 400 GeV");
    } else {
      j<nDirs_ ? g_xResol400_vs_run[j]->Draw("Psame") : g_xResol400_vs_run[j]->Draw("Lsame");
    }

    if(time_axis_format){
      timify(g_xResol400_vs_run[j]);
    }

    if(j==(nDirs_+nMCDirs_)-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(g_xResol400_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_xResol400_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // y resolution @ sumpT=400
    // *************************************
    
    g_yResol400_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((yResolTrend400_[LegLabels[j]])[0]));

    adjustmargins(c_yResol400_vs_run);
    c_yResol400_vs_run->cd();
    g_yResol400_vs_run[j]->SetMarkerStyle(pv::markers[j]);
    g_yResol400_vs_run[j]->SetMarkerColor(pv::colors[j]);
    g_yResol400_vs_run[j]->SetLineColor(pv::colors[j]);
    if(j>=nDirs_)  g_yResol400_vs_run[j]->SetLineWidth(3);

    g_yResol400_vs_run[j]->SetName(Form("g_yResol400_%s",LegLabels[j].Data()));
    g_yResol400_vs_run[j]->SetTitle(Form("Primary Vertex y-Resolution vs %s",theType.Data()));
    g_yResol400_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_yResol400_vs_run[j]->GetYaxis()->SetTitle("Primary Vertex y-Resolution [#mum]");
    g_yResol400_vs_run[j]->GetYaxis()->SetRangeUser(-10.,50.);
    beautify(g_yResol400_vs_run[j]);
 
    if(j==0){
      g_yResol400_vs_run[j]->Draw("AP");
      t1.DrawLatexNDC(0.20,0.85,"#sump_{T} = 400 GeV");
    } else {
      j<nDirs_ ? g_yResol400_vs_run[j]->Draw("Psame") :  g_yResol400_vs_run[j]->Draw("Lsame") ;
    }

    if(time_axis_format){
      timify(g_yResol400_vs_run[j]);
    }

    if(j==(nDirs_+nMCDirs_)-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(g_yResol400_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_yResol400_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

    // *************************************
    // z resolution @ sumpT=400
    // *************************************
    
    g_zResol400_vs_run[j]    = new TGraph(x_ticks.size(),&(x_ticks[0]),&((zResolTrend400_[LegLabels[j]])[0]));

    adjustmargins(c_zResol400_vs_run);
    c_zResol400_vs_run->cd();
    g_zResol400_vs_run[j]->SetMarkerStyle(pv::markers[j]);
    g_zResol400_vs_run[j]->SetMarkerColor(pv::colors[j]);
    g_zResol400_vs_run[j]->SetLineColor(pv::colors[j]);
    if(j>=nDirs_)  g_zResol400_vs_run[j]->SetLineWidth(3);

    g_zResol400_vs_run[j]->SetName(Form("g_zResol400_%s",LegLabels[j].Data()));
    g_zResol400_vs_run[j]->SetTitle(Form("Primary Vertex z-Resolution vs %s",theType.Data()));
    g_zResol400_vs_run[j]->GetXaxis()->SetTitle(theTypeLabel.Data());
    g_zResol400_vs_run[j]->GetYaxis()->SetTitle("Primary Vertex z-Resolution [#mum]");
    g_zResol400_vs_run[j]->GetYaxis()->SetRangeUser(-10.,50.);
    beautify(g_zResol400_vs_run[j]);
 
    if(j==0){
      g_zResol400_vs_run[j]->Draw("AP");
      t1.DrawLatexNDC(0.20,0.85,"#sump_{T} = 400 GeV");
    } else {
      j<nDirs_ ? g_zResol400_vs_run[j]->Draw("Psame") : g_zResol400_vs_run[j]->Draw("Lsame") ;
    }

    if(time_axis_format){
      timify(g_zResol400_vs_run[j]);
    }

    if(j==(nDirs_+nMCDirs_)-1){
      my_lego->Draw("same");
    }

    if(j==0){
      TH1F* theZero = DrawConstantGraph(g_zResol400_vs_run[j],1,0.);
      theZero->Draw("E1same");
    }

    current_pad = static_cast<TPad*>(gPad);
    cmsPrel(current_pad);

    superImposeIOVBoundaries(c_zResol400_vs_run,lumi_axis_format,time_axis_format,lumiMapByRun,times);

  }

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

  c_xResol200_vs_run->SaveAs("xResol200_vs_"+append+".png");
  c_yResol200_vs_run->SaveAs("yResol200_vs_"+append+".png");
  c_zResol200_vs_run->SaveAs("zResol200_vs_"+append+".png");

  c_xResol200_vs_run->SaveAs("xResol200_vs_"+append+".pdf");
  c_yResol200_vs_run->SaveAs("yResol200_vs_"+append+".pdf");
  c_zResol200_vs_run->SaveAs("zResol200_vs_"+append+".pdf");

  c_xResol400_vs_run->SaveAs("xResol400_vs_"+append+".png");
  c_yResol400_vs_run->SaveAs("yResol400_vs_"+append+".png");
  c_zResol400_vs_run->SaveAs("zResol400_vs_"+append+".png");

  c_xResol400_vs_run->SaveAs("xResol400_vs_"+append+".pdf");
  c_yResol400_vs_run->SaveAs("yResol400_vs_"+append+".pdf");
  c_zResol400_vs_run->SaveAs("zResol400_vs_"+append+".pdf");

  // do all the deletes

  for(int iDir=0;iDir<(nDirs_+nMCDirs_);iDir++){

   delete g_xResol200_vs_run[iDir]; 
   delete g_yResol200_vs_run[iDir];    
   delete g_zResol200_vs_run[iDir];    

   delete g_xResol400_vs_run[iDir]; 
   delete g_yResol400_vs_run[iDir];    
   delete g_zResol400_vs_run[iDir];    

  }

  // mv the run-by-run plots into the folders

  gSystem->mkdir("VertexResolutionsVsPt");
  TString processline = ".! mv VertexResolutionsVsPt*.p* ./VertexResolutionsVsPt/";
  std::cout<<"Executing: \n"
	   <<processline<< "\n"<<std::endl;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
  processline.Clear();

  timer.Stop(); 	 
  timer.Print();

}

///////////////////////////////////////////////////////////////
//
// END OF MAIN FUNCTION
//
///////////////////////////////////////////////////////////////


//*************************************************************
void arrangeOutCanvas(TCanvas *canv, TH1F* m_11Trend[100],TH1F* m_12Trend[100],TH1F* m_13Trend[100],Int_t nDirs, TString LegLabels[10],unsigned int theRun){
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
  canv->Divide(3,1);
 
  TH1F *dResolTrend[3][nDirs]; 
  
  for(Int_t i=0;i<nDirs;i++){
    dResolTrend[0][i] = m_11Trend[i];
    dResolTrend[1][i] = m_12Trend[i];
    dResolTrend[2][i] = m_13Trend[i];
  }

  Double_t absmin[3]={999.,999.,999.};
  Double_t absmax[3]={-999.,-999.-999.};

  for(Int_t k=0; k<3; k++){

    canv->cd(k+1)->SetBottomMargin(0.14);
    canv->cd(k+1)->SetLeftMargin(0.18);
    canv->cd(k+1)->SetRightMargin(0.01);
    canv->cd(k+1)->SetTopMargin(0.06);
    canv->cd(k+1);
    
    for(Int_t i=0; i<nDirs; i++){
      if(dResolTrend[k][i]->GetMaximum()>absmax[k]) absmax[k] = dResolTrend[k][i]->GetMaximum();
      if(dResolTrend[k][i]->GetMinimum()<absmin[k]) absmin[k] = dResolTrend[k][i]->GetMinimum();
    }

    Double_t safeDelta=(absmax[k]-absmin[k])/8.;
    Double_t theExtreme=std::max(absmax[k],TMath::Abs(absmin[k]));

    for(Int_t i=0; i<nDirs; i++){
      if(i==0){

	TString theTitle = dResolTrend[k][i]->GetName();
	dResolTrend[k][i]->GetYaxis()->SetRangeUser(0.,theExtreme+(safeDelta/2.));
       
	dResolTrend[k][i]->Draw("Le1");
	makeNewXAxis(dResolTrend[k][i]);
      
	Double_t theC = 10.;
		
	TH1F* theConst = DrawConstant(dResolTrend[k][i],1,theC);
	theConst->Draw("PLsame");

      } else { 
	dResolTrend[k][i]->Draw("Le1sames");
	makeNewXAxis(dResolTrend[k][i]);
      }
      TPad *current_pad = static_cast<TPad*>(canv->GetPad(k+1));
      cmsPrel(current_pad,2);
      ptDate->Draw("same");

      if(k==0){
	lego->AddEntry(dResolTrend[k][i],LegLabels[i]);
      } 
    }  
  
    lego->Draw();
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
      if (!file->IsDirectory() && fname.EndsWith(ext) && fname.BeginsWith("pvresolution")) {
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
  hist->SetMarkerStyle(pv::markers[color]);
  hist->SetLineColor(pv::colors[color]);
  hist->SetMarkerColor(pv::colors[color]);
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
  } else if (myTitle.Contains("SumPt")){
    axmin = 0.;
    axmax = 1000.;
    ndiv = 505;
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
void cmsPrel(TPad* pad,size_t ipads) {
/*--------------------------------------------------------------------*/
  
  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  float relPosX = 0.009;
  float relPosY = 0.045;
  float lumiTextOffset = 0.8;

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);

  float posX_    = 1-(r/ipads);
  float posY_    = 1-t + 0.05; /// - relPosY*(1-t-b);
  float factor   = 1./0.82;

  latex->SetTextAlign(33);
  latex->SetTextSize(0.045);
  latex->SetTextFont(42); //22
  latex->DrawLatex(posX_,posY_,"Internal (13 TeV)");

  UInt_t w;
  UInt_t h;
  latex->GetTextExtent(w,h,"Internal (13 TeV)");
  float size = w/(W/ipads);
  //std::cout<<w<<" "<<" "<<W<<" "<<size<<std::endl;
  float posXCMS_ = posX_- size*(1+0.025*ipads);

  latex->SetTextAlign(33);
  latex->SetTextFont(61);
  latex->SetTextSize(0.045*factor);
  latex->DrawLatex(posXCMS_,posY_+0.004,"CMS");

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
TH1F* DrawConstantWithErr(TH1F *hist,Int_t iter,Double_t theConst)
/*--------------------------------------------------------------------*/
{ 

  Int_t nbins       = hist->GetNbinsX();
  Double_t lowedge  = hist->GetBinLowEdge(1);
  Double_t highedge = hist->GetBinLowEdge(nbins+1);


  TH1F *hzero = new TH1F(Form("hconst_%s_%i",hist->GetName(),iter),Form("hconst_%s_%i",hist->GetName(),iter),nbins,lowedge,highedge);
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
pv::resolutions getResolutions(TH1F* hist)
/*--------------------------------------------------------------------*/
{
  int nbins = hist->GetNbinsX();

  int bin50  = hist->FindBin(50.);
  int bin100 = hist->FindBin(100.);
  int bin200 = hist->FindBin(200.);
  int bin400 = hist->FindBin(400.);

  float resol50  = hist->GetBinContent(bin50)  !=0. ? hist->GetBinContent(bin50)  : -9999.;
  float resol100 = hist->GetBinContent(bin100) !=0. ? hist->GetBinContent(bin100) : -9999.;
  float resol200 = hist->GetBinContent(bin200) !=0. ? hist->GetBinContent(bin200) : -9999.;
  float resol400 = hist->GetBinContent(bin400) !=0. ? hist->GetBinContent(bin400) : -9999.;
  
  pv::resolutions result(resol50,resol100,resol200,resol400);
  
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
  g->GetYaxis()->SetTitleOffset(0.6);
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
  h->GetYaxis()->SetTitleOffset(0.6);
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
  canv->cd()->SetLeftMargin(0.07);
  canv->cd()->SetRightMargin(0.03);
  canv->cd()->SetTopMargin(0.06);
}

/*--------------------------------------------------------------------*/
void adjustmargins(TVirtualPad *canv){
/*--------------------------------------------------------------------*/
  canv->SetBottomMargin(0.12);
  canv->SetLeftMargin(0.07);
  canv->SetRightMargin(0.01);
  canv->SetTopMargin(0.02);
}

/*--------------------------------------------------------------------*/
TH1F* checkTH1AndReturn(TFile *f,TString address){
/*--------------------------------------------------------------------*/
  TH1F* h(nullptr);
  if(f->GetListOfKeys()->Contains(address)){
    h = (TH1F*)f->Get(address);
  } 
  return h;
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


struct increase
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

/*--------------------------------------------------------------------*/
Double_t getMaximumFromArray(TObjArray *array)
/*--------------------------------------------------------------------*/
{

  Double_t theMaximum = -999.; //(static_cast<TH1*>(array->At(0)))->GetMaximum();

  for(Int_t i = 0; i< array->GetSize(); i++){
    
    double theMaxForThisHist;
    auto hist = static_cast<TH1*>(array->At(i));
    std::vector<double> maxima;
    for (int j=0;j<hist->GetNbinsX();j++) maxima.push_back(hist->GetBinContent(j));
    std::sort(std::begin(maxima), std::end(maxima));//,increase());
    double rms_maxima  = TMath::RMS(hist->GetNbinsX(),&(maxima[0]));
    double mean_maxima = TMath::Mean(hist->GetNbinsX(),&(maxima[0]));

    const Int_t nq = 100;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    Double_t yq[nq];  // array to contain the quantiles
    for (Int_t i=0;i<nq;i++) xq[i] = 0.9+(Float_t(i+1)/nq)*0.10;
    TMath::Quantiles(maxima.size(),nq,&(maxima[0]),yq, xq);
    
    //for(int q=0;q<nq;q++){
    //  std::cout<<q<<" "<<xq[q]<<" "<<yq[q]<<std::endl;
    //}

    //for (const auto &element : maxima){
    //  if(element<1.5*mean_maxima){
    //	theMaxForThisHist=element;
    //	break;
    //  }
    //} 

    theMaxForThisHist=yq[80];

    std::cout<<"rms_maxima["<<i<<"]"<<rms_maxima<<" mean maxima["<<mean_maxima<<"] purged maximum:"<<theMaxForThisHist<<std::endl;

    if(theMaxForThisHist>theMaximum) theMaximum=theMaxForThisHist;

    /*
    if( (static_cast<TH1*>(array->At(i)))->GetMaximum() > theMaximum){
      theMaximum = (static_cast<TH1*>(array->At(i)))->GetMaximum();
      //cout<<"i= "<<i<<" theMaximum="<<theMaximum<<endl;
    }
    */

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

  static const int nIOVs=13; //     1      2      3      4      5      6      7      8      9     10     11     12     13 
  int IOVboundaries[nIOVs]  = {294034,296641,297179,297281,298653,299277,299443,300389,301046,302131,303790,304911,305898};
  int benchmarkruns[nIOVs]  = {296173,297057,297219,297503,299061,299368,300157,300560,301472,302472,304292,305108,305898};
  TArrow* IOV_lines[nIOVs];
  c->cd();
  c->Update();

  TArrow* a_lines[nIOVs];
  TArrow* b_lines[nIOVs];
  for(Int_t IOV=0;IOV<nIOVs;IOV++){

    // check we are not in the RMS histogram to avoid first line
    if(IOVboundaries[IOV]<vruns.front() && ((TString)c->GetName()).Contains("RMS")) continue;
    int closestrun = pv::closest(vruns,IOVboundaries[IOV]); 
    int closestbenchmark = pv::closest(vruns,benchmarkruns[IOV]);

    if(lumi_axis_format){

      if(closestrun<0) continue;
      //std::cout<< "natural boundary: " << IOVboundaries[IOV] << " closest:" << closestrun << std::endl;

      a_lines[IOV] = new TArrow(lumiMapByRun.at(closestrun),(c->GetUymin()),lumiMapByRun.at(closestrun),0.65*c->GetUymax(),0.5,"|>");

      if(closestbenchmark<0) continue;
      b_lines[IOV] = new TArrow(lumiMapByRun.at(closestbenchmark),(c->GetUymin()),lumiMapByRun.at(closestbenchmark),0.65*c->GetUymax(),0.5,"|>");

    } else if(time_axis_format){
      
      if(closestrun<0) continue;
      std::cout<< "natural boundary: " << IOVboundaries[IOV] << " closest:" << closestrun << std::endl;
      a_lines[IOV] = new TArrow(timeMap.at(closestrun).Convert(),(c->GetUymin()),timeMap.at(closestrun).Convert(),0.65*c->GetUymax(),0.5,"|>");

      if(closestbenchmark<0) continue;
      b_lines[IOV] = new TArrow(timeMap.at(closestbenchmark).Convert(),(c->GetUymin()),timeMap.at(closestbenchmark).Convert(),0.65*c->GetUymax(),0.5,"|>");

    } else {
      a_lines[IOV] = new TArrow(IOVboundaries[IOV],(c->GetUymin()),IOVboundaries[IOV],0.65*c->GetUymax(),0.5,"|>"); //(c->GetUymin()+0.2*(c->GetUymax()-c->GetUymin()) ),0.5,"|>");
      b_lines[IOV] = new TArrow(benchmarkruns[IOV],(c->GetUymin()),benchmarkruns[IOV],0.65*c->GetUymax(),0.5,"|>"); //(c->GetUymin()+0.2*(c->GetUymax()-c->GetUymin()) ),0.5,"|>");
      
    }
    a_lines[IOV]->SetLineColor(kRed);
    a_lines[IOV]->SetLineStyle(9);
    a_lines[IOV]->SetLineWidth(1);
    a_lines[IOV]->Draw("same");

    b_lines[IOV]->SetLineColor(kGray);
    b_lines[IOV]->SetLineStyle(1);
    b_lines[IOV]->SetLineWidth(2);
    b_lines[IOV]->Draw("same");

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
    
    Int_t index = IOV%5;
    // if(IOV<5) 
    //   index=IOV;
    // else{
    //   index=IOV-5;
    // }
    
    runnumbers[IOV] = new TPaveText(_sx+0.001,0.14+(0.03*index),_dx,(0.17+0.03*index),"blNDC");
    //runnumbers[IOV]->SetTextAlign(11);
    TText *textRun = runnumbers[IOV]->AddText(Form("%i",int(IOVboundaries[IOV])));
    textRun->SetTextSize(0.028);
    textRun->SetTextColor(kRed);
    runnumbers[IOV]->SetFillColor(10);
    runnumbers[IOV]->SetLineColor(kRed);
    runnumbers[IOV]->SetBorderSize(1);
    runnumbers[IOV]->SetLineWidth(1);
    runnumbers[IOV]->SetTextColor(kRed);
    runnumbers[IOV]->SetTextFont(42);
    runnumbers[IOV]->Draw("same");
  }
}

/*--------------------------------------------------------------------*/
outTrends doStuff(size_t iter,std::vector<int> intersection,const Int_t nDirs_,const Int_t nMCDirs_,const char* dirs[10], TString LegLabels[10],bool useRMS)
/*--------------------------------------------------------------------*/
{
  outTrends ret;
  ret.init();

  unsigned int pitch = std::ceil(intersection.size()/nWorkers);
  unsigned int first = iter*pitch;
  unsigned int last  = std::min((iter+1)*pitch-1,intersection.size());

  std::cout<< "pitch: " << pitch<< " first: "<< first << " last: "<< last<< std::endl;

  ret.m_index=iter;

  for(unsigned int n=first; n<last;n++){
    
    std::cout << n << " "<<intersection.at(n) << std::endl;
    
    TFile *fins[nDirs_+nMCDirs_];

    TH1F* xPVResolVsPT[nDirs_+nMCDirs_];  
    TH1F* yPVResolVsPT[nDirs_+nMCDirs_]; 
    TH1F* zPVResolVsPT[nDirs_+nMCDirs_];   
    
    bool areAllFilesOK = true;
    Int_t lastOpen = 0;
 
    // loop over the objects
    for(Int_t j=0; j < nDirs_; j++) {

      //fins[j] = TFile::Open(Form("%s/pvresolution_%s_%i.root",dirs[j],dirs[j],intersection[n]));

      size_t position = std::string(dirs[j]).find("/");   
      string stem = std::string(dirs[j]).substr(position+1);     // get from position to the end
      
      fins[j] = new TFile(Form("%s/pvresolution_%s_%i.root",dirs[j],stem.c_str(),intersection[n]));
      if(fins[j]->IsZombie()){
	std::cout<< Form("%s/pvresolution_%s_%i.root",dirs[j],stem.c_str(),intersection[n]) << " is a Zombie! cannot combine" << std::endl;
	areAllFilesOK = false;
	lastOpen=j;
	break;
      }

      std::cout<< Form("%s/pvresolution_%s_%i.root",dirs[j],stem.c_str(),intersection[n]) 
	       << " has size: "<<fins[j]->GetSize() << " b ";
      
      // sanity check
      TH1F* h_tracks = (TH1F*)fins[j]->Get("PrimaryVertexResolution/h_nVertices");
      if(j==0){
	//auto h_lumi = checkTH1AndReturn(fins[j],"PrimaryVertexResolution/EventFeatures/h_lumiFromConfig");
	double lumi(0.0);
	if(fins[j]->GetDirectory("PrimaryVertexResolution/EventFeatures/")!= 0){
	  TH1F* h_lumi   = (TH1F*)fins[j]->Get("PrimaryVertexResolution/EventFeatures/h_lumiFromConfig");
	  std::cout<<"lumi: "<<lumi << std::endl;
	  lumi = h_lumi->GetBinContent(1);
	} else {
	  lumi=0.;
	}
	ret.m_lumiSoFar+=lumi;
	std::cout<<"lumi: "<<lumi
		 <<" ,lumi so far: "<<ret.m_lumiSoFar<<std::endl;

	// outfile<<"run "<<intersection[n]<<" lumi: "<<lumi
	//        <<" ,lumi so far: "<<ret.m_lumiSoFar<<std::endl;

      }

      Double_t numEvents = h_tracks->GetEntries();
      if(numEvents<2500){
	std::cout<<"excluding " << intersection[n] << "because it has less than 2.5k events" << std::endl;
	areAllFilesOK = false;
	lastOpen=j;
	break;
      }

      xPVResolVsPT[j] = (TH1F*)fins[j]->Get("PrimaryVertexResolution/p_resolX_vsSumPt");
      yPVResolVsPT[j] = (TH1F*)fins[j]->Get("PrimaryVertexResolution/p_resolY_vsSumPt");
      zPVResolVsPT[j] = (TH1F*)fins[j]->Get("PrimaryVertexResolution/p_resolZ_vsSumPt");
    
      // fill the vectors of resolutions
      
      auto xResol = getResolutions(xPVResolVsPT[j]);
      ret.m_xResolTrend200[LegLabels[j]].push_back(xResol.getResol200());
      ret.m_xResolTrend400[LegLabels[j]].push_back(xResol.getResol400());

      auto yResol = getResolutions(yPVResolVsPT[j]);
      ret.m_yResolTrend200[LegLabels[j]].push_back(yResol.getResol200());
      ret.m_yResolTrend400[LegLabels[j]].push_back(yResol.getResol400());

      auto zResol = getResolutions(zPVResolVsPT[j]);
      ret.m_zResolTrend200[LegLabels[j]].push_back(zResol.getResol200());
      ret.m_zResolTrend400[LegLabels[j]].push_back(zResol.getResol400());

      // beautify the histograms
      MakeNiceTrendPlotStyle(xPVResolVsPT[j],j);
      MakeNiceTrendPlotStyle(yPVResolVsPT[j],j);
      MakeNiceTrendPlotStyle(zPVResolVsPT[j],j);

    }

    // now for the MC
    for(Int_t jMC=nDirs_; jMC<(nDirs_+nMCDirs_); jMC++){
      size_t position = std::string(dirs[jMC]).find("/");   
      string stem = std::string(dirs[jMC]).substr(position+1);     // get from position to the end
      
      fins[jMC] = new TFile(Form("%s/pvresolution_%s_1.root",dirs[jMC],stem.c_str()));
      if(fins[jMC]->IsZombie()){
	std::cout<< Form("%s/pvresolution_%s_1.root",dirs[jMC],stem.c_str()) << " is a Zombie! cannot combine" << std::endl;
	areAllFilesOK = false;
	lastOpen=jMC;
	break;
      }

      std::cout<< Form("%s/pvresolution_%s_1.root",dirs[jMC],stem.c_str())
	       << " has size: "<<fins[jMC]->GetSize() << " b ";

      xPVResolVsPT[jMC] = (TH1F*)fins[jMC]->Get("PrimaryVertexResolution/p_resolX_vsSumPt");
      yPVResolVsPT[jMC] = (TH1F*)fins[jMC]->Get("PrimaryVertexResolution/p_resolY_vsSumPt");
      zPVResolVsPT[jMC] = (TH1F*)fins[jMC]->Get("PrimaryVertexResolution/p_resolZ_vsSumPt");
    
      // fill the vectors of resolutions
      
      auto xResol = getResolutions(xPVResolVsPT[jMC]);
      ret.m_xResolTrend200[LegLabels[jMC]].push_back(xResol.getResol200());
      ret.m_xResolTrend400[LegLabels[jMC]].push_back(xResol.getResol400());

      auto yResol = getResolutions(yPVResolVsPT[jMC]);
      ret.m_yResolTrend200[LegLabels[jMC]].push_back(yResol.getResol200());
      ret.m_yResolTrend400[LegLabels[jMC]].push_back(yResol.getResol400());

      auto zResol = getResolutions(zPVResolVsPT[jMC]);
      ret.m_zResolTrend200[LegLabels[jMC]].push_back(zResol.getResol200());
      ret.m_zResolTrend400[LegLabels[jMC]].push_back(zResol.getResol400());

      // beautify the histograms
      MakeNiceTrendPlotStyle(xPVResolVsPT[jMC],jMC);
      MakeNiceTrendPlotStyle(yPVResolVsPT[jMC],jMC);
      MakeNiceTrendPlotStyle(zPVResolVsPT[jMC],jMC);

    }

    if(!areAllFilesOK){
       
      // do all the necessary deletions
      std::cout<<"\n====> not all files are OK"<<std::endl;

      for(int i=0;i<lastOpen;i++){
	fins[i]->Close();
      }
      continue;
    } else {
      ret.m_runs.push_back(intersection.at(n));
      // push back the vector of lumi (in fb at this point)
      ret.m_lumiByRun.push_back(ret.m_lumiSoFar/1000.);
      ret.m_lumiMapByRun[intersection.at(n)]=ret.m_lumiSoFar/1000.;
    }

    std::cout<<"I am still here - runs.size(): "<<ret.m_runs.size()<<std::endl;

    // resolutions pT plots

    TCanvas *VertexResolutionsVsPt = new TCanvas(Form("VertexResolutionsVsPT_%i",intersection.at(n)),"VertexResolutionsVsPt",1800,800);
    arrangeOutCanvas(VertexResolutionsVsPt,xPVResolVsPT,yPVResolVsPT,zPVResolVsPT,(nDirs_+nMCDirs_),LegLabels,intersection.at(n));

    VertexResolutionsVsPt->SaveAs(Form("VertexResolutionsVsPt_%i.pdf",intersection.at(n)));
    VertexResolutionsVsPt->SaveAs(Form("VertexResolutionsVsPt_%i.png",intersection.at(n)));

    // do all the necessary deletions

    for(int i=0;i<nDirs_+nMCDirs_;i++){

      delete xPVResolVsPT[i];
      delete yPVResolVsPT[i];
      delete zPVResolVsPT[i];
    
      fins[i]->Close();
    }
    
    delete VertexResolutionsVsPt;

    std::cout<<std::endl;
  }
  
  return ret;

}
