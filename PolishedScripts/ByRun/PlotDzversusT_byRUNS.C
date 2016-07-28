#include <cmath>
#include "TCanvas.h"
#include "TROOT.h"
#include "TkAlStyle.cc"
#include "TAxis.h"
#include "TArrow.h"
#include "TTreeIndex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
//#include "tdrstyle.C"

void makeNiceTGraph(TGraph *graph);

void PlotDzversusT_byRUNS(string tag="PromptGT"){

  gROOT->ProcessLine(".L TkAlStyle.cc+");
  TkAlStyle::set(INPROGRESS);	// set publication status

  //gROOT->ProcessLine(".L ./setTDRStyle.C++");
  //setTDRStyle();

  //string tag="ORIGINAL";//"TBDkb";
  //ORIGINAL_DzDisparity_byRUNS.root
  //TBDkb_newIOV_DzDisparity_byRUNS_fromJuly1.root

  string inFileName=tag+".root";
  TFile *f = new TFile(inFileName.c_str());
  TTree* t = (TTree*) f->Get("PVValidationTree");

  const int DummyNumEvents = int(t->GetEntries());
  const int NumEvents=DummyNumEvents;
  cout<<"NumEvents = "<<NumEvents <<endl;
  
  double mean_dxy_phi,width_dxy_phi,mean_dz_phi,width_dz_phi,mean_dxy_eta,width_dxy_eta,mean_dz_eta,width_dz_eta;
  double mean_n_dxy_phi,width_n_dxy_phi,mean_n_dz_phi,width_n_dz_phi,mean_n_dxy_eta,width_n_dxy_eta,mean_n_dz_eta,width_n_dz_eta;
  
  double dz_fit,dz_fitError;

  double min_dxy_phi,max_dxy_phi,min_dz_phi,max_dz_phi,min_dxy_eta,max_dxy_eta,min_dz_eta,max_dz_eta;
  double min_n_dxy_phi,max_n_dxy_phi,min_n_dz_phi,max_n_dz_phi,min_n_dxy_eta,max_n_dxy_eta,min_n_dz_eta,max_n_dz_eta;
  
  int run;
  double nevents;
  double ntracks;

  int x[NumEvents];
  double y_dxy_phi[NumEvents];
  double y_fake[NumEvents];
  double y_dxy_phi_max[NumEvents];
  double y_dxy_phi_min[NumEvents];
  double y_dxy_phi_Error[NumEvents];
  double x_dxy_phi[NumEvents];
  double x_dxy_phi_Error[NumEvents];
  
  double y_dz_phi[NumEvents];
  double y_dz_phi_max[NumEvents];
  double y_dz_phi_min[NumEvents];
  double y_dz_phi_Error[NumEvents];
  double x_dz_phi[NumEvents];
  double x_dz_phi_Error[NumEvents];

  double y_dxy_eta[NumEvents];
  double y_dxy_eta_max[NumEvents];
  double y_dxy_eta_min[NumEvents];
  double y_dxy_eta_Error[NumEvents];
  double x_dxy_eta[NumEvents];
  double x_dxy_eta_Error[NumEvents];

  double y_dz_eta[NumEvents];
  double y_dz_eta_max[NumEvents];
  double y_dz_eta_min[NumEvents];  
  double y_dz_eta_Error[NumEvents];
  double x_dz_eta[NumEvents];
  double x_dz_eta_Error[NumEvents];

  double y_dzfit[NumEvents];
  double y_dzfit_Error[NumEvents];

  double noDatax[10000];
  double noDatay[10000];

  t->SetBranchAddress("run"             , &run         );
  t->SetBranchAddress("ntracks"         , &ntracks     );
  t->SetBranchAddress("dz_fit"          , &dz_fit      );
  t->SetBranchAddress("dz_fit_error"    , &dz_fitError );

  t->SetBranchAddress("nevents"         , &nevents     );
  t->SetBranchAddress("mean_dxy_phi"    , &mean_dxy_phi);
  t->SetBranchAddress("width_dxy_phi"   , &width_dxy_phi);
  t->SetBranchAddress("mean_dz_phi"     , &mean_dz_phi);
  t->SetBranchAddress("width_dz_phi"    , &width_dz_phi);
  t->SetBranchAddress("mean_dxy_eta"    , &mean_dxy_eta);
  t->SetBranchAddress("width_dxy_eta"   , &width_dxy_eta);  
  t->SetBranchAddress("mean_dz_eta"     , &mean_dz_eta);
  t->SetBranchAddress("width_dz_eta"    , &width_dz_eta);
  t->SetBranchAddress("mean_n_dxy_phi"  , &mean_n_dxy_phi);
  t->SetBranchAddress("width_n_dxy_phi" , &width_n_dxy_phi);
  t->SetBranchAddress("mean_n_dz_phi"   , &mean_n_dz_phi);
  t->SetBranchAddress("width_n_dz_phi"  , &width_n_dz_phi);
  t->SetBranchAddress("mean_n_dxy_eta"  , &mean_n_dxy_eta);
  t->SetBranchAddress("width_n_dxy_eta" , &width_n_dxy_eta);  
  t->SetBranchAddress("mean_n_dz_eta"   , &mean_n_dz_eta);
  t->SetBranchAddress("width_n_dz_eta"  , &width_n_dz_eta);  

  t->SetBranchAddress("min_dxy_phi"    , &min_dxy_phi);
  t->SetBranchAddress("max_dxy_phi"    , &max_dxy_phi);
  t->SetBranchAddress("min_dz_phi"     , &min_dz_phi);
  t->SetBranchAddress("max_dz_phi"     , &max_dz_phi);
  t->SetBranchAddress("min_dxy_eta"    , &min_dxy_eta);
  t->SetBranchAddress("max_dxy_eta"    , &max_dxy_eta);  
  t->SetBranchAddress("min_dz_eta"     , &min_dz_eta);
  t->SetBranchAddress("max_dz_eta"     , &max_dz_eta);
  t->SetBranchAddress("min_n_dxy_phi"  , &min_n_dxy_phi);
  t->SetBranchAddress("max_n_dxy_phi"  , &max_n_dxy_phi);
  t->SetBranchAddress("min_n_dz_phi"   , &min_n_dz_phi);
  t->SetBranchAddress("max_n_dz_phi"   , &max_n_dz_phi);
  t->SetBranchAddress("min_n_dxy_eta"  , &min_n_dxy_eta);
  t->SetBranchAddress("max_n_dxy_eta"  , &max_n_dxy_eta);  
  t->SetBranchAddress("min_n_dz_eta"   , &min_n_dz_eta);
  t->SetBranchAddress("max_n_dz_eta"   , &max_n_dz_eta); 

  //-----------------------------------------------

  int numEventsGood=0;
  int numEventsBad=0;
  
  t->BuildIndex("run");
  TTreeIndex *index = (TTreeIndex*)t->GetTreeIndex();
  for( int i = index->GetN() - 1; i >=0 ; --i ) {
    Long64_t local = t->LoadTree( index->GetIndex()[i] );

  //  for(int ievt=0; ievt<NumEvents; ievt++){
       
    t->GetEntry(local);

    //cout<<run<<endl;
    //if(run<256000)
    //  continue;
   
    if(nevents<1000)
      continue;
    
    //if(sqrt(innerFitErr*innerFitErr+outerFitErr*outerFitErr)>30) continue;

    if(fabs(min_dxy_phi) > 1000. || 
       fabs(max_dxy_phi) > 1000. || 
       fabs(min_dz_phi)  > 500. || 
       fabs(max_dz_phi)  > 500. ||
       fabs(min_dxy_eta) > 2000. || 
       fabs(max_dxy_eta) > 2000. ||
       fabs(min_dz_eta)  > 2000. || 
       fabs(max_dz_eta)  > 2000.
       ){
      noDatay[numEventsBad] = mean_dxy_phi;
      noDatax[numEventsBad] = i;
      numEventsBad++;

    } else {

      y_fake[numEventsGood]    = 0.;
      y_dxy_phi[numEventsGood] = mean_dxy_phi;
      y_dz_phi[numEventsGood]  = mean_dz_phi;
      y_dxy_eta[numEventsGood] = mean_dxy_eta;
      y_dz_eta[numEventsGood]  = mean_dz_eta;
      
      y_dxy_phi_min[numEventsGood] = min_dxy_phi;
      y_dz_phi_min[numEventsGood]  = min_dz_phi;
      y_dxy_eta_min[numEventsGood] = min_dxy_eta;
      y_dz_eta_min[numEventsGood]  = min_dz_eta;
      
      y_dxy_phi_max[numEventsGood] = max_dxy_phi;
      y_dz_phi_max[numEventsGood]  = max_dz_phi;
      y_dxy_eta_max[numEventsGood] = max_dxy_eta;
      y_dz_eta_max[numEventsGood]  = max_dz_eta;
      
      y_dzfit[numEventsGood]       = dz_fit;
      y_dzfit_Error[numEventsGood]  = dz_fitError;
 
      if(fabs(y_dxy_phi[numEventsGood])>0.01){
	
	y_dxy_phi_Error[numEventsGood] = 0.;//width_dxy_phi/2.;
	y_dz_phi_Error[numEventsGood]  = 0.;//width_dz_phi/2.;
	y_dxy_eta_Error[numEventsGood] = 0.;//width_dxy_eta/2.;
	y_dz_eta_Error[numEventsGood]  = 0.;// width_dz_eta/2.;
	
      } else {
	
	y_dxy_phi_Error[numEventsGood] = 0.;
	y_dz_phi_Error[numEventsGood]  = 0.;
	y_dxy_eta_Error[numEventsGood] = 0.;
	y_dz_eta_Error[numEventsGood]  = 0.;

      }

      x_dxy_phi[numEventsGood]  = run; // numEventsGood
      x_dz_phi[numEventsGood]   = run; // numEventsGood
      x_dxy_eta[numEventsGood]  = run; // numEventsGood
      x_dz_eta[numEventsGood]   = run; // numEventsGood

      x_dxy_phi_Error[numEventsGood] = .5;
      x_dz_phi_Error[numEventsGood]  = .5;
      x_dxy_eta_Error[numEventsGood] = .5;
      x_dz_eta_Error[numEventsGood]  = .5;

      x[numEventsGood] = run;

      cout<<" x_dxy_phi[" << numEventsGood << "]" << x_dxy_phi[numEventsGood] << " y_dxy_phi["<<numEventsGood<<"] " <<  y_dxy_phi[numEventsGood] << endl; 
      numEventsGood++;
      
    }

  }

  //////***************************************************///////

  Double_t max_scale = 30;

  TGraphErrors *dummyGraph_dxy_phi = new TGraphErrors(numEventsGood,x_dxy_phi,y_fake,0,0);
  TGraphErrors *dummyGraph_dz_phi  = new TGraphErrors(numEventsGood,x_dxy_phi,y_fake,0,0);
  TGraphErrors *dummyGraph_dxy_eta = new TGraphErrors(numEventsGood,x_dxy_phi,y_fake,0,0);
  TGraphErrors *dummyGraph_dz_eta  = new TGraphErrors(numEventsGood,x_dxy_phi,y_fake,0,0);
  TGraphErrors *dummyGraph_dz_fit  = new TGraphErrors(numEventsGood,x_dxy_phi,y_fake,0,0);

  TCanvas *c_dxy_phi = new TCanvas("c_dxy_phi","", 1000,500);
  c_dxy_phi->SetFillColor(0);
  c_dxy_phi->SetBottomMargin(0.14);
  c_dxy_phi->SetLeftMargin(0.07);
  c_dxy_phi->SetRightMargin(0.05);
  c_dxy_phi->SetTopMargin(0.08);  
  c_dxy_phi->SetGrid();

  TGraphErrors *graph_dxy_phi = new TGraphErrors(numEventsGood,x_dxy_phi ,y_dxy_phi ,x_dxy_phi_Error ,y_dxy_phi_Error);

  TGraphErrors *graph_dxy_phi_min   = new TGraphErrors(numEventsGood,x_dxy_phi,y_dxy_phi_min,x_dxy_phi_Error ,y_dxy_phi_Error);
  TGraphErrors *graph_dxy_phi_max   = new TGraphErrors(numEventsGood,x_dxy_phi,y_dxy_phi_max,x_dxy_phi_Error ,y_dxy_phi_Error);
  TGraph *graph_dxy_phi_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_dxy_phi_shade->SetPoint(i,x_dxy_phi[i],y_dxy_phi_max[i]);
    graph_dxy_phi_shade->SetPoint(numEventsGood+i,x_dxy_phi[numEventsGood-i-1],y_dxy_phi_min[numEventsGood-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_dxy_phi_shade->SetFillStyle(3013);
  graph_dxy_phi_shade->SetFillColor(16);
  graph_dxy_phi_min->SetLineColor(kGray);
  graph_dxy_phi_max->SetLineColor(kGray);
  graph_dxy_phi_min->SetLineWidth(2);
  graph_dxy_phi_max->SetLineWidth(2);

  makeNiceTGraph(graph_dxy_phi);

  graph_dxy_phi->GetXaxis()->SetTitle("Run Number");
  graph_dxy_phi->GetYaxis()->SetTitle("#LT d_{xy} #GT vs #phi [#mum]");
  graph_dxy_phi->SetTitle("d_{xy} vs #phi average bias vs. time");
  graph_dxy_phi->SetMarkerStyle(20);
  graph_dxy_phi->SetMarkerSize(0.5);
  graph_dxy_phi->SetMarkerColor(kBlack);
  graph_dxy_phi->SetLineColor(kBlack);

  //////***************************************************///////

  TCanvas *c_dz_phi = new TCanvas("c_dz_phi","", 1000,500);
  c_dz_phi->SetFillColor(0);
  c_dz_phi->SetBottomMargin(0.14);
  c_dz_phi->SetLeftMargin(0.07);
  c_dz_phi->SetRightMargin(0.05);
  c_dz_phi->SetTopMargin(0.08);  
  c_dz_phi->SetGrid();

  TGraphErrors *graph_dz_phi  = new TGraphErrors(numEventsGood,x_dz_phi  ,y_dz_phi  ,x_dz_phi_Error  ,y_dz_phi_Error);

  TGraph *graph_dz_phi_min   = new TGraph(numEventsGood,x_dz_phi,y_dz_phi_min);
  TGraph *graph_dz_phi_max   = new TGraph(numEventsGood,x_dz_phi,y_dz_phi_max);
  TGraph *graph_dz_phi_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_dz_phi_shade->SetPoint(i,x_dz_phi[i],y_dz_phi_max[i]);
    graph_dz_phi_shade->SetPoint(numEventsGood+i,x_dz_phi[numEventsGood-i-1],y_dz_phi_min[numEventsGood-i-1]);
  }
  graph_dz_phi_shade->SetFillStyle(3013);
  graph_dz_phi_shade->SetFillColor(16);
  graph_dz_phi_min->SetLineColor(kGray);
  graph_dz_phi_max->SetLineColor(kGray);
  graph_dz_phi_min->SetLineWidth(2);
  graph_dz_phi_max->SetLineWidth(2);
  
  makeNiceTGraph(graph_dz_phi);

  graph_dz_phi->GetXaxis()->SetTitle("Run Number");
  graph_dz_phi->GetYaxis()->SetTitle("#LT d_{z} #GT vs #phi [#mum]");
  graph_dz_phi->SetTitle("d_{z} vs #phi average bias vs. time");
  graph_dz_phi->SetMarkerStyle(20);
  graph_dz_phi->SetMarkerSize(0.5);
  graph_dz_phi->SetMarkerColor(kBlack);
  graph_dz_phi->SetLineColor(kBlack);

  //////***************************************************///////

  TCanvas *c_dxy_eta = new TCanvas("c_dxy_eta","", 1000,500);
  c_dxy_eta->SetFillColor(0);
  c_dxy_eta->SetBottomMargin(0.14);
  c_dxy_eta->SetLeftMargin(0.07);
  c_dxy_eta->SetRightMargin(0.05);
  c_dxy_eta->SetTopMargin(0.08);  
  c_dxy_eta->SetGrid();
  
  TGraphErrors *graph_dxy_eta = new TGraphErrors(numEventsGood,x_dxy_eta ,y_dxy_eta ,x_dxy_eta_Error ,y_dxy_eta_Error);

  TGraph *graph_dxy_eta_min   = new TGraph(numEventsGood,x_dxy_eta,y_dxy_eta_min);
  TGraph *graph_dxy_eta_max   = new TGraph(numEventsGood,x_dxy_eta,y_dxy_eta_max);
  TGraph *graph_dxy_eta_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_dxy_eta_shade->SetPoint(i,x_dxy_eta[i],y_dxy_eta_max[i]);
    graph_dxy_eta_shade->SetPoint(numEventsGood+i,x_dxy_eta[numEventsGood-i-1],y_dxy_eta_min[numEventsGood-i-1]);
  }
  graph_dxy_eta_shade->SetFillStyle(3013);
  graph_dxy_eta_shade->SetFillColor(16);
  graph_dxy_eta_min->SetLineColor(kGray);
  graph_dxy_eta_max->SetLineColor(kGray);
  graph_dxy_eta_min->SetLineWidth(2);
  graph_dxy_eta_max->SetLineWidth(2);

  makeNiceTGraph(graph_dxy_eta);

  graph_dxy_eta->GetXaxis()->SetTitle("Run Number");
  graph_dxy_eta->GetYaxis()->SetTitle("#LT d_{xy} #GT vs #eta [#mum]");
  graph_dxy_eta->SetTitle("d_{xy} vs #eta average bias vs. time");

  graph_dxy_eta->SetMarkerStyle(20);
  graph_dxy_eta->SetMarkerSize(0.5);
  graph_dxy_eta->SetMarkerColor(kBlack);
  graph_dxy_eta->SetLineColor(kBlack);

  //////***************************************************///////

  TCanvas *c_dz_eta = new TCanvas("c_dz_eta","", 1000,500);
  c_dz_eta->SetFillColor(0);
  c_dz_eta->SetBottomMargin(0.14);
  c_dz_eta->SetLeftMargin(0.07);
  c_dz_eta->SetRightMargin(0.05);
  c_dz_eta->SetTopMargin(0.08);  
  c_dz_eta->SetGrid();

  TGraphErrors *graph_dz_eta  = new TGraphErrors(numEventsGood,x_dz_eta  ,y_dz_eta  ,x_dz_eta_Error  ,y_dz_eta_Error); 

  TGraph *graph_dz_eta_min   = new TGraph(numEventsGood,x_dz_eta,y_dz_eta_min);
  TGraph *graph_dz_eta_max   = new TGraph(numEventsGood,x_dz_eta,y_dz_eta_max);
  TGraph *graph_dz_eta_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_dz_eta_shade->SetPoint(i,x_dz_eta[i],y_dz_eta_max[i]);
    graph_dz_eta_shade->SetPoint(numEventsGood+i,x_dz_eta[numEventsGood-i-1],y_dz_eta_min[numEventsGood-i-1]);
  }
  graph_dz_eta_shade->SetFillStyle(3013);
  graph_dz_eta_shade->SetFillColor(16);
  graph_dz_eta_min->SetLineColor(kGray);
  graph_dz_eta_max->SetLineColor(kGray);
  graph_dz_eta_min->SetLineWidth(2);
  graph_dz_eta_max->SetLineWidth(2);

  makeNiceTGraph(graph_dz_eta);

  graph_dz_eta->GetXaxis()->SetTitle("Run Number");
  graph_dz_eta->GetYaxis()->SetTitle("#LT d_{z} #GT vs #eta [#mum]");
  graph_dz_eta->SetTitle("d_{z} vs #eta average bias vs. time");

  graph_dz_eta->SetMarkerStyle(20);
  graph_dz_eta->SetMarkerSize(0.5);
  graph_dz_eta->SetMarkerColor(kBlack);
  graph_dz_eta->SetLineColor(kBlack);
  
  TCanvas *c_dz_fit = new TCanvas("c_dz_fit","", 1000,500);
  c_dz_fit->SetFillColor(0);
  c_dz_fit->SetBottomMargin(0.14);
  c_dz_fit->SetLeftMargin(0.07);
  c_dz_fit->SetRightMargin(0.05);
  c_dz_fit->SetTopMargin(0.08);  
  c_dz_fit->SetGrid();

  TGraphErrors *graph_dz_fit = new TGraphErrors(numEventsGood,x_dxy_phi ,y_dzfit ,x_dxy_phi_Error ,y_dzfit_Error);
  graph_dz_fit->GetXaxis()->SetLabelSize(0.05);
  graph_dz_fit->GetYaxis()->SetLabelSize(0.05);
  graph_dz_fit->GetXaxis()->SetNoExponent();
  graph_dz_fit->GetXaxis()->SetNdivisions(510);
  graph_dz_fit->GetXaxis()->SetTitle("Run Number");
  graph_dz_fit->GetYaxis()->SetTitle("#LT d_{xy} #GT vs #phi [#mum]");
  graph_dz_fit->SetTitle("d_{xy} vs #phi average bias vs. time");
  graph_dz_fit->GetYaxis()->CenterTitle(true);
  graph_dz_fit->GetXaxis()->CenterTitle(true);
  graph_dz_fit->GetYaxis()->SetNdivisions(505);
  graph_dz_fit->GetYaxis()->SetTitleOffset(0.80);
  graph_dz_fit->SetMarkerStyle(20);
  graph_dz_fit->SetMarkerSize(0.5);
  graph_dz_fit->SetMarkerColor(kBlack);
  graph_dz_fit->SetLineColor(kBlack);
  
  TGraph *noDataGraph = new TGraph(NumEvents-numEventsGood,noDatax,noDatay);

  /*
    from Hcal gains tag
    250988 3.8T 3e6bded5a30d31acf6735b37f8834364d396e636
    252039 0T  	c529f437df2fb09fcf6e6f6f64a832389372f848
    253999 3.8T	3e6bded5a30d31acf6735b37f8834364d396e636
    254289 0T  	c529f437df2fb09fcf6e6f6f64a832389372f848
    254656 3.8T	3e6bded5a30d31acf6735b37f8834364d396e636
    255988 0T  	c529f437df2fb09fcf6e6f6f64a832389372f848
    256547 3.8T	3e6bded5a30d31acf6735b37f8834364d396e636
  */

  TBox *magnetCyle252039_x = new TBox(252039,-max_scale,253998,max_scale);
  TBox *magnetCyle254289_x = new TBox(254289,-max_scale,254655,max_scale);
  TBox *magnetCyle255988_x = new TBox(255988,-max_scale,256546,max_scale);
  TBox *magnetCyle257021_x = new TBox(257021,-max_scale,257375,max_scale);
  TBox *magnetCyle259896_x = new TBox(259896,-max_scale,260335,max_scale);
  TBox *magnetCyle260432_x = new TBox(260432,-max_scale,260528,max_scale);  
  TBox *magnetCyle260744_x = new TBox(260744,-max_scale,261084,max_scale);
  TBox *magnetCyle261122_x = new TBox(261122,-max_scale,261914,max_scale);
  TBox *magnetCyle262348_x = new TBox(262348,-max_scale,262442,max_scale);
  TBox *magnetCyle263096_x = new TBox(263096,-max_scale,263230,max_scale);

  magnetCyle252039_x->SetFillColorAlpha(42,0.5);
  magnetCyle254289_x->SetFillColorAlpha(42,0.5);
  magnetCyle255988_x->SetFillColorAlpha(42,0.5);
  magnetCyle257021_x->SetFillColorAlpha(42,0.5);
  magnetCyle259896_x->SetFillColorAlpha(42,0.5);
  magnetCyle260432_x->SetFillColorAlpha(42,0.5);
  magnetCyle260744_x->SetFillColorAlpha(42,0.5);
  magnetCyle261122_x->SetFillColorAlpha(42,0.5);
  magnetCyle262348_x->SetFillColorAlpha(42,0.5);
  magnetCyle263096_x->SetFillColorAlpha(42,0.5);

  TBox *magnetCyle252039_L = new TBox(252039,-max_scale*4,253998,max_scale*4);
  TBox *magnetCyle254289_L = new TBox(254289,-max_scale*4,254655,max_scale*4);
  TBox *magnetCyle255988_L = new TBox(255988,-max_scale*4,256546,max_scale*4);
  TBox *magnetCyle257021_L = new TBox(257021,-max_scale*4,257375,max_scale*4);
  TBox *magnetCyle259896_L = new TBox(259896,-max_scale*4,260335,max_scale*4);
  TBox *magnetCyle260432_L = new TBox(260432,-max_scale*4,260528,max_scale*4);  
  TBox *magnetCyle260744_L = new TBox(260744,-max_scale*4,261084,max_scale*4);
  TBox *magnetCyle261122_L = new TBox(261122,-max_scale*4,261914,max_scale*4);
  TBox *magnetCyle262348_L = new TBox(262348,-max_scale*4,262442,max_scale*4);
  TBox *magnetCyle263096_L = new TBox(263096,-max_scale*4,263230,max_scale*4);

  magnetCyle252039_L->SetFillColorAlpha(42,0.5);
  magnetCyle254289_L->SetFillColorAlpha(42,0.5);
  magnetCyle255988_L->SetFillColorAlpha(42,0.5);
  magnetCyle257021_L->SetFillColorAlpha(42,0.5);
  magnetCyle259896_L->SetFillColorAlpha(42,0.5);
  magnetCyle260432_L->SetFillColorAlpha(42,0.5);
  magnetCyle260744_L->SetFillColorAlpha(42,0.5);
  magnetCyle261122_L->SetFillColorAlpha(42,0.5);
  magnetCyle262348_L->SetFillColorAlpha(42,0.5);
  magnetCyle263096_L->SetFillColorAlpha(42,0.5);

  /*  '256355', # 0T
      '256490', # ramping up
      '256547', # 3.8T
      '257021', # ramping down
      '257042', # 0T
      '257356', # ramping up
      '257376', # 3.8T
      '259896', # ramping down
      '259913', # 0T
      '260325', # ramping up
      '260335', # 3.8T
      '260432', # ramping down
      '260456', # 0T
      '260507', # ramping up
      '260528', # 3.8T
  */

  //** alignment tag updates **//
  /*
    251100 
    251607 
    256715 
    256746 
    257826 
    262922 
   */
  
  TLine *l1 = new TLine(251100,-max_scale,251100,max_scale);
  TLine *l2 = new TLine(251607,-max_scale,251607,max_scale);
  TLine *l3 = new TLine(256715,-max_scale,256715,max_scale);
  TLine *l4 = new TLine(256746,-max_scale,256746,max_scale);
  TLine *l5 = new TLine(257826,-max_scale,257826,max_scale);
  TLine *l6 = new TLine(262922,-max_scale,262922,max_scale);

  TLine *L1 = new TLine(251100,-max_scale*4,251100,max_scale*4);
  TLine *L2 = new TLine(251607,-max_scale*4,251607,max_scale*4);
  TLine *L3 = new TLine(256715,-max_scale*4,256715,max_scale*4);
  TLine *L4 = new TLine(256746,-max_scale*4,256746,max_scale*4);
  TLine *L5 = new TLine(257826,-max_scale*4,257826,max_scale*4);
  TLine *L6 = new TLine(262922,-max_scale*4,262922,max_scale*4);

  l1->SetLineColor(4);
  l1->SetLineWidth(1);
  l1->SetLineStyle(9);  
  l2->SetLineColor(4);
  l2->SetLineWidth(1);
  l2->SetLineStyle(9);  
  l3->SetLineColor(4);
  l3->SetLineWidth(1);
  l3->SetLineStyle(9);
  l4->SetLineColor(4);
  l4->SetLineWidth(1);
  l4->SetLineStyle(9);
  l5->SetLineColor(4);
  l5->SetLineWidth(1);
  l5->SetLineStyle(9);
  l6->SetLineColor(4);
  l6->SetLineWidth(1);
  l6->SetLineStyle(9);

  L1->SetLineColor(4);
  L1->SetLineWidth(1);
  L1->SetLineStyle(9);  
  L2->SetLineColor(4);
  L2->SetLineWidth(1);
  L2->SetLineStyle(9);  
  L3->SetLineColor(4);
  L3->SetLineWidth(1);
  L3->SetLineStyle(9);
  L4->SetLineColor(4);
  L4->SetLineWidth(1);
  L4->SetLineStyle(9);
  L5->SetLineColor(4);
  L5->SetLineWidth(1);
  L5->SetLineStyle(9);
  L6->SetLineColor(4);
  L6->SetLineWidth(1);
  L6->SetLineStyle(9);

  TString boundaries[27] = {"275657","275658","275766","275772","275778","275828","275829","275837","275841","275847","275887","275921","276242","276243","276244","276282","276283","276315","276317","276327","276352","276361","276363","276384","276453","276454","276653"};

  TArrow *lines[27];

  for(Int_t i=0;i<27;i++){
    lines[i]=new TArrow(boundaries[i].Atoi(),-30.,boundaries[i].Atoi(),30,0.5,"|>");
    lines[i]->SetLineColor(kBlue);
    lines[i]->SetLineStyle(9);
    lines[i]->SetLineWidth(2);
  }

  dummyGraph_dxy_phi->SetMarkerColorAlpha(0,0.);
  dummyGraph_dz_phi->SetMarkerColorAlpha(0,0.);
  dummyGraph_dxy_eta->SetMarkerColorAlpha(0,0.);
  dummyGraph_dz_eta->SetMarkerColorAlpha(0,0.);
  
  /* 
     dummyGraph_xPos->GetXaxis()->SetTitle("Run");
     dummyGraph_xPos->GetYaxis()->SetTitle("#Deltax [#mum]");
     dummyGraph_xPos->GetYaxis()->SetRangeUser(-xMax,xMax);
     dummyGraph_xPos->GetXaxis()->SetLimits(usedRuns[0]-1,usedRuns[runRange-1]+100);
     dummyGraph_xPos->GetXaxis()->SetNdivisions(4);
     dummyGraph_xPos->GetXaxis()->SetLabelSize(0.06);
  */


  //******************************************//

  c_dxy_phi->cd();
  makeNiceTGraph(dummyGraph_dxy_phi);

  dummyGraph_dxy_phi->GetXaxis()->SetTitle("Run Number");
  dummyGraph_dxy_phi->GetYaxis()->SetTitle("#LT d_{z} #GT vs #phi [#mum]");
  dummyGraph_dxy_phi->SetTitle("d_{xy} vs #phi average bias vs. time");
  dummyGraph_dxy_phi->GetYaxis()->SetRangeUser(-max_scale,max_scale);
  dummyGraph_dxy_phi->Draw("ap");
  
  //graph_dxy_phi->GetYaxis()->SetRangeUser(-max_scale,max_scale);

  magnetCyle252039_x->Draw("same");
  magnetCyle254289_x->Draw("same");
  magnetCyle255988_x->Draw("same");
  magnetCyle257021_x->Draw("same");
  magnetCyle259896_x->Draw("same");
  magnetCyle260432_x->Draw("same");

  graph_dxy_phi->Draw("lpsame");
  graph_dxy_phi_shade->Draw("fsame");
  graph_dxy_phi_min->SetMarkerSize(0.2);
  graph_dxy_phi_max->SetMarkerSize(0.2);  

  graph_dxy_phi_min->SetLineColor(kGray);
  graph_dxy_phi_max->SetLineColor(kGray);
  graph_dxy_phi_min->Draw("lsame"); 
  graph_dxy_phi_max->Draw("lsame");
 
  l1->Draw("SAME");
  l2->Draw("SAME");
  l3->Draw("SAME");
  l4->Draw("SAME");
  l5->Draw("SAME");
  l6->Draw("SAME");

  for(Int_t i=0;i<27;i++){
    lines[i]->Draw("same");
  }

  c_dxy_phi->SaveAs(Form("dxy_phi_%s.png",tag.c_str()));
  c_dxy_phi->SaveAs(Form("dxy_phi_%s.pdf",tag.c_str()));

  //******************************************//

  c_dz_phi->cd();

  makeNiceTGraph(dummyGraph_dz_phi);
  dummyGraph_dz_phi->GetXaxis()->SetTitle("Run Number");
  dummyGraph_dz_phi->GetYaxis()->SetTitle("#LT d_{z} #GT vs #phi [#mum]");
  dummyGraph_dz_phi->SetTitle("d_{z} vs #phi average bias vs. time");
  dummyGraph_dz_phi->GetYaxis()->SetRangeUser(-max_scale,max_scale);
  dummyGraph_dz_phi->Draw("ap");

  magnetCyle252039_x->Draw("same");
  magnetCyle254289_x->Draw("same");
  magnetCyle255988_x->Draw("same");
  magnetCyle257021_x->Draw("same");
  magnetCyle259896_x->Draw("same");
  magnetCyle260432_x->Draw("same");

  //graph_dz_phi->GetYaxis()->SetRangeUser(-max_scale,max_scale);
  graph_dz_phi->Draw("lpsame");
  graph_dz_phi_shade->Draw("fsame");

  graph_dz_phi_min->SetLineColor(kGray);
  graph_dz_phi_max->SetLineColor(kGray);

  graph_dz_phi_min->Draw("lsame"); 
  graph_dz_phi_max->Draw("lsame");
 
  l1->Draw("SAME");
  l2->Draw("SAME");
  l3->Draw("SAME");
  l4->Draw("SAME");
  l5->Draw("SAME");
  l6->Draw("SAME");

  for(Int_t i=0;i<27;i++){
    lines[i]->Draw("same");
  }

  c_dz_phi->SaveAs(Form("dz_phi_%s.png",tag.c_str()));
  c_dz_phi->SaveAs(Form("dz_phi_%s.pdf",tag.c_str()));

  //******************************************//

  c_dxy_eta->cd();

  makeNiceTGraph(dummyGraph_dxy_eta);
  //graph_dxy_eta->GetYaxis()->SetRangeUser(-max_scale,max_scale);
  dummyGraph_dxy_eta->GetXaxis()->SetTitle("Run Number");
  dummyGraph_dxy_eta->GetYaxis()->SetTitle("#LT d_{xy} #GT vs #eta [#mum]");
  dummyGraph_dxy_eta->SetTitle("d_{xy} vs #eta average bias vs. time");
  dummyGraph_dxy_eta->GetYaxis()->SetRangeUser(-max_scale,max_scale);
  dummyGraph_dxy_eta->Draw("ap");

  magnetCyle252039_x->Draw("same");
  magnetCyle254289_x->Draw("same");
  magnetCyle255988_x->Draw("same");
  magnetCyle257021_x->Draw("same");
  magnetCyle259896_x->Draw("same");
  magnetCyle260432_x->Draw("same");

  graph_dxy_eta->Draw("lpsame");
  graph_dxy_eta_shade->Draw("fsame");
  
  graph_dxy_eta_min->SetLineColor(kGray);
  graph_dxy_eta_max->SetLineColor(kGray);
  
  graph_dxy_eta_min->Draw("lsame"); 
  graph_dxy_eta_max->Draw("lsame");

  l1->Draw("SAME");
  l2->Draw("SAME");
  l3->Draw("SAME");
  l4->Draw("SAME");
  l5->Draw("SAME");
  l6->Draw("SAME");

  for(Int_t i=0;i<27;i++){
    lines[i]->Draw("same");
  }

  c_dxy_eta->SaveAs(Form("dxy_eta_%s.png",tag.c_str()));
  c_dxy_eta->SaveAs(Form("dxy_eta_%s.pdf",tag.c_str()));

   //******************************************//

  c_dz_fit->cd();

  makeNiceTGraph(dummyGraph_dz_fit);
  //graph_dz_fit->GetYaxis()->SetRangeUser(-max_scale,max_scale);
  dummyGraph_dz_fit->GetXaxis()->SetTitle("Run Number");
  dummyGraph_dz_fit->GetYaxis()->SetTitle("#Delta d_{z} [#mum]");
  dummyGraph_dz_fit->SetTitle("#Delta d_{z} bias vs. time");
  dummyGraph_dz_fit->GetYaxis()->SetRangeUser(-max_scale,max_scale);
  dummyGraph_dz_fit->Draw("ap");

  magnetCyle252039_x->Draw("same");
  magnetCyle254289_x->Draw("same");
  magnetCyle255988_x->Draw("same");
  magnetCyle257021_x->Draw("same");
  magnetCyle259896_x->Draw("same");
  magnetCyle260432_x->Draw("same");

  graph_dz_fit->Draw("lpsame");

  l1->Draw("SAME");
  l2->Draw("SAME");
  l3->Draw("SAME");
  l4->Draw("SAME");
  l5->Draw("SAME");
  l6->Draw("SAME");

  for(Int_t i=0;i<27;i++){
    lines[i]->Draw("same");
  }
  
  c_dz_fit->SaveAs(Form("dz_fit_%s.png",tag.c_str()));
  c_dz_fit->SaveAs(Form("dz_fit_%s.pdf",tag.c_str()));


  //******************************************//

  c_dz_eta->cd();
  makeNiceTGraph(dummyGraph_dz_eta);
  dummyGraph_dz_eta->GetXaxis()->SetTitle("Run Number");
  dummyGraph_dz_eta->GetYaxis()->SetTitle("#LT d_{z} #GT vs #eta [#mum]");
  dummyGraph_dz_eta->SetTitle("d_{z} vs #eta average bias vs. time");
  dummyGraph_dz_eta->GetYaxis()->SetRangeUser(-max_scale*4,max_scale*4);
  dummyGraph_dz_eta->Draw("ap");

  magnetCyle252039_L->Draw("same");
  magnetCyle254289_L->Draw("same");
  magnetCyle255988_L->Draw("same");
  magnetCyle257021_L->Draw("same");
  magnetCyle259896_L->Draw("same");
  magnetCyle260432_L->Draw("same");

  graph_dz_eta->Draw("lpsame");
  graph_dz_eta_shade->Draw("fsame");
  graph_dz_eta_min->SetLineColor(kGray);
  graph_dz_eta_max->SetLineColor(kGray);
  graph_dz_eta_min->Draw("lsame"); 
  graph_dz_eta_max->Draw("lsame");

  //graph_dz_eta->GetYaxis()->SetRangeUser(-max_scale*4,max_scale*4);
  
  L1->Draw("SAME");
  L2->Draw("SAME");
  L3->Draw("SAME");
  L4->Draw("SAME");
  L5->Draw("SAME");
  L6->Draw("SAME");
 
  for(Int_t i=0;i<27;i++){
    lines[i]->Draw("same");
  }

  c_dz_eta->SaveAs(Form("dz_eta_%s.pdf",tag.c_str()));
  c_dz_eta->SaveAs(Form("dz_eta_%s.png",tag.c_str()));

  /*
  int lineBegin=-20;
  int lineEnd=40;

  time.Set(2011,3,25,12,0,0);
  
  time.Set(2011,4,19,12,0,0);
  TLine *l2 = new TLine(time.Convert()-offset,lineBegin,time.Convert()-offset,lineEnd);
  l2->SetLineColor(4);
  l2->SetLineWidth(2);
  time.Set(2011,5,4,12,0,0);
  TLine *l3 = new TLine(time.Convert()-offset,lineBegin,time.Convert()-offset,lineEnd);
  l3->SetLineColor(4);
  l3->SetLineWidth(2);
  time.Set(2011,5,20,12,0,0);
  TLine *l4 = new TLine(time.Convert()-offset,lineBegin,time.Convert()-offset,lineEnd);
  l4->SetLineColor(4);
  l4->SetLineWidth(2);
  time.Set(2011,5,28,12,0,0);
  TLine *l5 = new TLine(time.Convert()-offset,lineBegin,time.Convert()-offset,lineEnd);
  l5->SetLineColor(4);
  l5->SetLineWidth(2);
  time.Set(2011,6,2,12,0,0);
  TLine *l6 = new TLine(time.Convert()-offset,lineBegin,time.Convert()-offset,lineEnd);
  l6->SetLineColor(4);
  l6->SetLineWidth(2);
  time.Set(2011,6,15,12,0,0);
  TLine *l7 = new TLine(time.Convert()-offset,lineBegin,time.Convert()-offset,lineEnd);
  l7->SetLineColor(4);
  l7->SetLineWidth(2);
  */

  // noDataGraph->GetXaxis()->SetTimeDisplay(1);
  // noDataGraph->SetMarkerColor(kBlack);

  //graph_dz_phi->GetYaxis()->SetRangeUser(-100,100);
  //graph_dz_phi->Draw("same");

  //graph_dxy_eta->GetYaxis()->SetRangeUser(-100,100);
  //graph_dxy_eta->Draw("same");

  // graph_dz_eta->GetYaxis()->SetRangeUser(-100,100);
  //graph_dz_eta->Draw("same");

  // noDataGraph->Draw("SAME*");
  /*
  l1->Draw("SAME");
  l2->Draw("SAME");
  l3->Draw("SAME");
  l4->Draw("SAME");
  l5->Draw("SAME");
  l6->Draw("SAME");
  l7->Draw("SAME");
  */
  //---------SAVE DATA---------------
  //TFile* dumpFile = new TFile("dzVsTime.root","UPDATE");
  //string name="dzVsTime_"+tag;
  //graph->Write(name.c_str());
  /*
  l1->Write("l1");
  l2->Write("l2");
  l3->Write("l3");
  l4->Write("l4");
  l5->Write("l5");
  l6->Write("l6");
  l7->Write("l7");
  */


  //  noDataGraph->Write("noData");
  //dumpFile->Close();
}

void makeNiceTGraph(TGraph *graph){

  graph->GetXaxis()->SetLabelSize(0.05);
  graph->GetYaxis()->SetLabelSize(0.05);
  graph->GetXaxis()->SetTitleSize(0.05);
  graph->GetYaxis()->SetTitleSize(0.05);
  graph->GetXaxis()->SetNoExponent();
  graph->GetXaxis()->SetNdivisions(510);
  graph->GetYaxis()->CenterTitle(true);
  graph->GetXaxis()->CenterTitle(true);
  graph->GetYaxis()->SetNdivisions(505);
  graph->GetYaxis()->SetTitleOffset(0.6);

}
