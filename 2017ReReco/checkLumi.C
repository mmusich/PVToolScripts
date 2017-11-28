#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include <iostream>

Double_t checkLumi(TString name){
  
  gErrorIgnoreLevel=kError;
  TFile *fin = TFile::Open(name,"READ");

  TH1F* h_tracks = (TH1F*)fin->Get("PVValidation/EventFeatures/h_nTracks");
  TH1F* h_eta = (TH1F*)fin->Get("PVValidation/DA/eta_all");
  TH1F* h_lumi = (TH1F*)fin->Get("PVValidation/EventFeatures/h_lumiFromConfig");
  
  Double_t numEvents = h_tracks->GetEntries();
  Double_t numTracks = h_eta->GetEntries();
  Double_t trkEff    = Double_t(numTracks/numEvents);
  Double_t lumi      = h_lumi->GetBinContent(1);
  
  cout.precision(6);
  cout<<"number of events: "<<numEvents<<endl;
  cout<<"number of tracks: "<<numTracks<<endl;
  cout<<"probe trks/event: "<<trkEff<<endl;
  cout<<"lumi: "<<lumi << endl;

  fin->Close();
  return lumi;
 
}
