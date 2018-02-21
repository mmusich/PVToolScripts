Long64_t checkFile(TString name){
  
  gErrorIgnoreLevel=kError;
  TFile *fin = TFile::Open(name,"READ");
  TH1F* h_nVert = (TH1F*)fin->Get("PrimaryVertexResolution/h_nVertices");
  TH1F* h_diffX = (TH1F*)fin->Get("PrimaryVertexResolution/h_diffX");
  
  Double_t numEvents = h_nVert->GetEntries();
  Double_t numVerts  = h_diffX->GetEntries();
  Double_t vtxEff    = Double_t(numVerts/numEvents);
  
  cout.precision(6);
  cout<<"number of events: "<<numEvents<<endl;
  cout<<"number of vertices: "<<numVerts<<endl;
  cout<<"probe vtx/event: "<<vtxEff<<endl;

  fin->Close();
  delete fin;

  return numEvents;
}
