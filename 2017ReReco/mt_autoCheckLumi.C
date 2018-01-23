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
#include <TStopwatch.h>
#include "TArrow.h"
#include "TCanvas.h"
#include "TObjString.h"
#include <iostream>
#include <iomanip>      // std::setprecision
#include <vector>
#include <algorithm>
#include <map>
#include <functional>
#include <iterator>
#include <fstream>
#include <sstream>

const size_t nWorkers=10;

//***********************************
// auxilliary struct to store information
//***********************************

namespace check {
  
  struct info{
    
    info(){
      m_events=0;
      m_tracks=0;
      m_lumi=0.;
    }
    
    info(int events,int tracks,double lumi){
      m_events=events;
      m_tracks=tracks;
      m_lumi=lumi;
    }

    int getEvents(){return m_events;}
    int getTracks(){return m_events;}
    double getLumi(){return m_lumi;}
    double getTkEfficiency(){return double(m_tracks/m_events);};

  private:
    int m_events;
    int m_tracks;
    double m_lumi;
  }; 

  using infoPerRun = std::map<int,check::info>;

}

/*--------------------------------------------------------------------*/
check::info checkLumi(TString name)
/*--------------------------------------------------------------------*/
{
  
  gErrorIgnoreLevel=kError;
  TFile *fin = TFile::Open(name,"READ");

  TH1F* h_tracks = (TH1F*)fin->Get("PVValidation/EventFeatures/h_nTracks");
  TH1F* h_eta = (TH1F*)fin->Get("PVValidation/DA/eta_all");
  TH1F* h_lumi = (TH1F*)fin->Get("PVValidation/EventFeatures/h_lumiFromConfig");
  
  double numEvents = h_tracks->GetEntries();
  double numTracks = h_eta->GetEntries();
  double trkEff    = double(numTracks/numEvents);
  double lumi      = h_lumi->GetBinContent(1);
  
  cout.precision(6);
  cout<<"number of events: "<<numEvents<<endl;
  cout<<"number of tracks: "<<numTracks<<endl;
  cout<<"probe trks/event: "<<trkEff<<endl;
  cout<<"lumi: "<<lumi << "/pb" << endl;

  check::info myret = check::info(numEvents,numTracks,lumi);
  fin->Close();
  return myret;
 
}

/*--------------------------------------------------------------------*/
check::infoPerRun checkLumiMultiRun(size_t iter,std::vector<int> runs)
/*--------------------------------------------------------------------*/
{
  const char* dirname="PromptGT";
  unsigned int pitch = std::ceil(runs.size()/nWorkers);
  unsigned int first = iter*pitch;
  unsigned int last  = std::min((iter+1)*pitch-1,runs.size());

  check::infoPerRun ret;
  size_t position = std::string(dirname).find("/");   
  string stem = std::string(dirname).substr(position+1);     // get from position to the end

  for (unsigned int n=first; n<last;n++){
    TString filename = Form("%s/PVValidation_%s_%i.root",dirname,stem.c_str(),n);
    check::info myInfo = checkLumi(filename);
    ret[n]=myInfo;
  }
  
  return ret;

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

//***********************************
//  MAIN METHOD
//***********************************
/*--------------------------------------------------------------------*/
void mt_autoCheckLumi()
/*--------------------------------------------------------------------*/
{
  const char* dirname="PromptGT";
  std::vector<int> runs = list_files(dirname,".root");
  std::sort(runs.begin(),runs.end());
  double lumiSoFar=0.;

  size_t position = std::string(dirname).find("/");   
  string stem = std::string(dirname).substr(position+1);     // get from position to the end
  
  std::ofstream outfile ("lumiInfoByRun.txt"); 

  for (const auto &run : runs){
    TString filename = Form("%s/PVValidation_%s_%i.root",dirname,stem.c_str(),run);
    check::info myInfo = checkLumi(filename);
    lumiSoFar+=myInfo.getLumi();
    std::cout<<"-------------------------------------------------"<<std::endl;
    std::cout<<"run: "<<run<<" lumi: " <<  std::setprecision(5) << myInfo.getLumi() << "/pb" <<" lumiSoFar: "<< std::setprecision(8) << lumiSoFar<< "/pb" << std::endl;

    outfile<<"-------------------------------------------------"<<std::endl;
    outfile<<"run: "<<run<<" lumi: " <<  std::setprecision(5) <<  myInfo.getLumi() << "/pb" <<" lumiSoFar: "<<  std::setprecision(8) << lumiSoFar<< "/pb" << std::endl;

  }

  outfile.close();
  
  /*
  auto f_doCheck = std::bind(checkLumiMultiRun,std::placeholders::_1,runs);

  TProcPool procPool(nWorkers);
  std::vector<size_t> range(nWorkers);
  std::iota(range.begin(),range.end(),0);
  auto extracts = procPool.Map(f_doCheck,range);
  */

}
