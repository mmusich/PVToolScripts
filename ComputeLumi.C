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
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <iterator>

std::vector<int> list_files(const char *dirname=".", const char *ext=".root");

void ComputeLumi(TString dir){
  
  std::vector<int> currentList = list_files(dir.Data());
  const int nfiles = currentList.size();

  TFile *fins[nfiles];
  double lumiSoFar = 0.0;

  for(unsigned int n=0; n<currentList.size();n++){
    fins[n] = TFile::Open(Form("%s/PVValidation_%s_%i.root",dir.Data(),dir.Data(),currentList[n]));
    TH1F* h_lumi = (TH1F*)fins[n]->Get("PVValidation/EventFeatures/h_lumiFromConfig");
    double lumi  = h_lumi->GetBinContent(1);
    lumiSoFar+=lumi;
    std::cout << "Run:"<< currentList[n] << " | lumi: " << lumi << " | total lumi: "<< lumiSoFar << std::endl;
    fins[n]->Close();
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

  // sort the run numbers
  std::sort(theRunNumbers.begin(),theRunNumbers.end());

  return theRunNumbers;
}
