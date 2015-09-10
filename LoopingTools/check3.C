Long64_t check(TString name){
  
  gErrorIgnoreLevel=kError;
  TFile *fin = TFile::Open(name);

  TH1F* h_tracks = (TH1F*)fin->Get("PVValidation/EventFeatures/h_nTracks");
  TH1F* h_eta = (TH1F*)fin->Get("PVValidation/DA/eta_all");
  
  Double_t numEvents = h_tracks->GetEntries();
  Double_t numTracks = h_eta->GetEntries();
  Double_t trkEff = Double_t(numTracks/numEvents);
  
  cout.precision(6);
  cout<<"number of events: "<<numEvents<<endl;
  cout<<"number of tracks: "<<numTracks<<endl;
  cout<<"probe trks/event: "<<trkEff<<endl;

  return numEvents;
}
