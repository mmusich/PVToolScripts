#include <cmath>
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TTreeIndex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <algorithm>

//#include "tdrstyle.C"

void makeNiceTGraph(TGraph *graph);

void replaceAll( string &s, const string &search, const string &replace ) {
  for( size_t pos = 0; ; pos += replace.length() ) {
    // Locate the substring to replace
    pos = s.find( search, pos );
    if( pos == string::npos ) break;
    // Replace by erasing and inserting
    s.erase( pos, search.length() );
    s.insert( pos, replace );
  }
}

void PlotSIPversusT_byLUMI(string inFileName="fit_pTCut_25_PromptGT.root"){

  //gROOT->ProcessLine(".L TkAlStyle.cc+");
  //TkAlStyle::set(INPROGRESS);	// set publication status

  //gROOT->ProcessLine(".L ./setTDRStyle.C++");
  //setTDRStyle();

  //string tag="ORIGINAL";//"TBDkb";
  //ORIGINAL_DzDisparity_byRUNS.root
  //TBDkb_newIOV_DzDisparity_byRUNS_fromJuly1.root


  TFile *f = new TFile(inFileName.c_str());
  TTree* t = (TTree*) f->Get("SIPValidationTree");

  replaceAll(inFileName,".root","");

  string outFileName=inFileName+"_ByLumi.root";
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");

  const int DummyNumEvents = int(t->GetEntries());
  const int numEvents=62; //DummyNumEvents;
  cout<<"NumEvents = "<<numEvents <<endl;
  
  double mean,meanError;
  double rms,rmsError;
  double lanWidth,lanMPV,area,GWidth;
  double err_lanWidth,err_lanMPV,err_area,err_GWidth;
  double fitChisqNdof;

  int run;

  double x[numEvents];
  double x_err[numEvents];
  double y_fake[numEvents];
  double y_mean[numEvents];
  double y_rms[numEvents];
  double y_lanWidth[numEvents];
  double y_lanMPV[numEvents];
  double y_area[numEvents];
  double y_fitChisq[numEvents];
  double y_GWidth[numEvents];
  double y_err_mean[numEvents];
  double y_err_rms[numEvents];
  double y_err_lanWidth[numEvents];
  double y_err_lanMPV[numEvents];
  double y_err_area[numEvents];
  double y_err_GWidth[numEvents];

  t->SetBranchAddress("run"             , &run         );
  t->SetBranchAddress("mean"            , &mean        );
  t->SetBranchAddress("meanError"       , &meanError   );
  t->SetBranchAddress("rms"             , &rms        );
  t->SetBranchAddress("rmsError"        , &rmsError   );
  t->SetBranchAddress("lanWidth"        , &lanWidth    );                      
  t->SetBranchAddress("lanMPV"          , &lanMPV      );                      
  t->SetBranchAddress("area"            , &area        );                         
  t->SetBranchAddress("GWidth"          , &GWidth      );                                                                                                               
  t->SetBranchAddress("err_lanWidth"    , &err_lanWidth);                                                                         
  t->SetBranchAddress("err_lanMPV"      , &err_lanMPV  );                                                                         
  t->SetBranchAddress("err_area"        , &err_area    );                                                                         
  t->SetBranchAddress("err_GWidth"      , &err_GWidth  );                                                                         
  t->SetBranchAddress("fitChisqNdof"    , &fitChisqNdof);  

  //-----------------------------------------------
 
  t->BuildIndex("run");
  
  double lumis[63] = {0.367706,1.006291,1.525403,2.197349,2.697196,3.238763,3.839784,4.491583,5.087056,5.559197,6.274189,6.726759,7.438598,7.664904,8.366887,8.859215,9.561405,9.700322,10.199944,11.02933,11.512876,12.00973,12.58104,13.37714,14.237404,14.966782,15.318177,16.05554,16.789415,17.17412,17.895814,18.548021,19.021169,19.788166,20.194597,20.987246,21.426946,21.895406,22.518841,23.220566,23.769254,24.447126,24.95788,25.820072,26.35347,27.026041,27.54226,27.985273,28.84846,29.686073,30.370196,30.925188,31.462419,31.775375,32.550908,33.159204,34.018259,34.410174,34.973403,35.668558,36.431429,36.768285,36.772695};

  TTreeIndex *index = (TTreeIndex*)t->GetTreeIndex();
  for( int i = 0; i <= index->GetN(); i++ ) {
    
    Long64_t local = t->LoadTree( index->GetIndex()[i] );
   
    t->GetEntry(local);

    if(fitChisqNdof<0.)
      continue;
    
    y_fake   [i]      = 0.;
    y_lanWidth[i]     = lanWidth;
    y_lanMPV[i]       = lanMPV;
    y_area[i]         = area;
    y_GWidth[i]       = GWidth;
    y_fitChisq[i]     = fitChisqNdof;
    y_mean[i]         = mean;
    y_rms[i]          = rms;
    
    y_err_lanWidth[i] = err_lanWidth;
    y_err_lanMPV[i]   = err_lanMPV;
    y_err_area[i]     = err_area;
    y_err_GWidth[i]   = err_GWidth;
    y_err_mean[i]     = meanError;
    y_err_rms[i]      = rmsError;
    
    //x[i] = run;
    //x[i] = (i+1)*0.5;
    x[i] = lumis[i];
    x_err[i] = 0.;

    cout<<" x[" << i << "]" << x[i] << " y_lanMPV["<<i<<"] " <<  y_lanMPV[i] << " +/- y_err_lanMPV["<<i<<"] " <<   y_err_lanMPV[i] << endl; 
    
  }

  //////***************************************************///////

  Double_t max_scale = 0.9;
  Double_t min_scale = 0.4;

  TGraphErrors *dummyGraph = new TGraphErrors(numEvents,x,y_fake,0,0);
  // -- mean

  TCanvas *cMean = new TCanvas("cMean","", 1000,500);
  cMean->SetFillColor(0);
  cMean->SetBottomMargin(0.14);
  cMean->SetLeftMargin(0.07);
  cMean->SetRightMargin(0.05);
  cMean->SetTopMargin(0.08);  
  cMean->SetGrid();

  TGraphErrors *graph_mean = new TGraphErrors(numEvents,x,y_mean,x_err,y_err_mean);
  TGraph *graph_mean_shade = new TGraph(2*numEvents);
  for (Int_t i=0;i<numEvents;i++) {
    graph_mean_shade->SetPoint(i,x[i],y_mean[i]+1*y_err_mean[i]);
    graph_mean_shade->SetPoint(numEvents+i,x[numEvents-i-1],y_mean[numEvents-i-1]-1*y_err_mean[numEvents-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_mean_shade->SetFillStyle(3013);
  graph_mean_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_mean);

  graph_mean->SetName("graph_mean");
  graph_mean->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  graph_mean->GetYaxis()->SetTitle("mean of SIP_{3D}");
  graph_mean->SetTitle("mean of SIP_{3D} vs. time");
  graph_mean->SetMarkerStyle(20);
  graph_mean->SetMarkerSize(1.2);
  graph_mean->SetMarkerColor(kRed);
  graph_mean->SetLineColor(kRed);

  cMean->cd();
  dummyGraph->GetYaxis()->SetRangeUser(min_scale,max_scale);
  graph_mean->Draw("alpsame");
  graph_mean_shade->Draw("fsame");

  fout->cd();
  graph_mean->Write();

  // -- rms 

  TCanvas *cRms = new TCanvas("cRms","", 1000,500);
  cRms->SetFillColor(0);
  cRms->SetBottomMargin(0.14);
  cRms->SetLeftMargin(0.07);
  cRms->SetRightMargin(0.05);
  cRms->SetTopMargin(0.08);  
  cRms->SetGrid();

  TGraphErrors *graph_rms = new TGraphErrors(numEvents,x,y_rms,x_err,y_err_rms);
  TGraph *graph_rms_shade = new TGraph(2*numEvents);
  for (Int_t i=0;i<numEvents;i++) {
    graph_rms_shade->SetPoint(i,x[i],y_rms[i]+1*y_err_rms[i]);
    graph_rms_shade->SetPoint(numEvents+i,x[numEvents-i-1],y_rms[numEvents-i-1]-1*y_err_rms[numEvents-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_rms_shade->SetFillStyle(3013);
  graph_rms_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_rms);

  graph_rms->SetName("graph_rms");
  graph_rms->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  graph_rms->GetYaxis()->SetTitle("rms of SIP_{3D}");
  graph_rms->SetTitle("rms of SIP_{3D} vs. time");
  graph_rms->SetMarkerStyle(20);
  graph_rms->SetMarkerSize(1.2);
  graph_rms->SetMarkerColor(kRed);
  graph_rms->SetLineColor(kRed);

  cRms->cd();
  dummyGraph->GetYaxis()->SetRangeUser(min_scale,max_scale);
  graph_rms->Draw("alpsame");
  graph_rms_shade->Draw("fsame");

  fout->cd();
  graph_rms->Write();

  // -- lan MPV

  TCanvas *c = new TCanvas("c","", 1000,500);
  c->SetFillColor(0);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.07);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.08);  
  c->SetGrid();

  TGraphErrors *graph_lanMPV = new TGraphErrors(numEvents,x,y_lanMPV,x_err,y_err_lanMPV);
  TGraph *graph_lanMPV_shade = new TGraph(2*numEvents);
  for (Int_t i=0;i<numEvents;i++) {
    graph_lanMPV_shade->SetPoint(i,x[i],y_lanMPV[i]+1*y_err_lanMPV[i]);
    graph_lanMPV_shade->SetPoint(numEvents+i,x[numEvents-i-1],y_lanMPV[numEvents-i-1]-1*y_err_lanMPV[numEvents-i-1]);
  }
  graph_lanMPV_shade->SetFillStyle(3013);
  graph_lanMPV_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_lanMPV);

  graph_lanMPV->SetName("graph_lanMPV");
  graph_lanMPV->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  graph_lanMPV->GetYaxis()->SetTitle("MPV of SIP_{3D}");
  graph_lanMPV->SetTitle("MPV(SIP_{3D}) vs. time");
  graph_lanMPV->SetMarkerStyle(20);
  graph_lanMPV->SetMarkerSize(1.2);
  graph_lanMPV->SetMarkerColor(kRed);
  graph_lanMPV->SetLineColor(kRed);

  c->cd();
  dummyGraph->GetYaxis()->SetRangeUser(min_scale,max_scale);
  graph_lanMPV->Draw("alpsame");
  graph_lanMPV_shade->Draw("fsame");

  fout->cd();
  graph_lanMPV->Write();

  //******************************

  TCanvas *c2 = new TCanvas("c2","", 1000,500);
  c2->SetFillColor(0);
  c2->SetBottomMargin(0.14);
  c2->SetLeftMargin(0.07);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.08);  
  c2->SetGrid();

  TGraphErrors *graph_GWidth = new TGraphErrors(numEvents,x,y_GWidth,x_err,y_err_GWidth);
  TGraph *graph_GWidth_shade = new TGraph(2*numEvents);
  for (Int_t i=0;i<numEvents;i++) {
    graph_GWidth_shade->SetPoint(i,x[i],y_GWidth[i]+1*y_err_GWidth[i]);
    graph_GWidth_shade->SetPoint(numEvents+i,x[numEvents-i-1],y_GWidth[numEvents-i-1]-1*y_err_GWidth[numEvents-i-1]);
  }
  graph_GWidth_shade->SetFillStyle(3013);
  graph_GWidth_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_GWidth);

  graph_GWidth->SetName("graph_GWidth");
  graph_GWidth->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  graph_GWidth->GetYaxis()->SetTitle("Gaussian width of SIP_{3D}");
  graph_GWidth->SetTitle("#sigma_{G}(SIP_{3D}) vs. time");
  graph_GWidth->SetMarkerStyle(20);
  graph_GWidth->SetMarkerSize(1.2);
  graph_GWidth->SetMarkerColor(kRed);
  graph_GWidth->SetLineColor(kRed);

  c2->cd();
  dummyGraph->GetYaxis()->SetRangeUser(min_scale,max_scale);
  graph_GWidth->Draw("alpsame");
  graph_GWidth_shade->Draw("fsame");

  fout->cd();
  graph_GWidth->Write();

  //******************************

  TCanvas *c3 = new TCanvas("c3","", 1000,500);
  c3->SetFillColor(0);
  c3->SetBottomMargin(0.14);
  c3->SetLeftMargin(0.07);
  c3->SetRightMargin(0.05);
  c3->SetTopMargin(0.08);  
  c3->SetGrid();

  TGraphErrors *graph_area = new TGraphErrors(numEvents,x,y_area,x_err,y_err_area);
  TGraph *graph_area_shade = new TGraph(2*numEvents);
  for (Int_t i=0;i<numEvents;i++) {
    graph_area_shade->SetPoint(i,x[i],y_area[i]+1*y_err_area[i]);
    graph_area_shade->SetPoint(numEvents+i,x[numEvents-i-1],y_area[numEvents-i-1]-1*y_err_area[numEvents-i-1]);
  }
  graph_area_shade->SetFillStyle(3013);
  graph_area_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_area);

  graph_area->SetName("graph_area");
  graph_area->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  graph_area->GetYaxis()->SetTitle("Area parameter of SIP_{3D}");
  graph_area->SetTitle("A(SIP_{3D}) vs. time");
  graph_area->SetMarkerStyle(20);
  graph_area->SetMarkerSize(1.2);
  graph_area->SetMarkerColor(kRed);
  graph_area->SetLineColor(kRed);

  c3->cd();
  // dummyGraph->GetYaxis()->SetRangeUser(min_scale,max_scale);
  graph_area->Draw("alpsame");
  graph_area_shade->Draw("fsame");

  fout->cd();
  graph_area->Write();

  //******************************

  TCanvas *c4 = new TCanvas("c4","", 1000,500);
  c4->SetFillColor(0);
  c4->SetBottomMargin(0.14);
  c4->SetLeftMargin(0.07);
  c4->SetRightMargin(0.05);
  c4->SetTopMargin(0.08);  
  c4->SetGrid();

  TGraphErrors *graph_lanWidth = new TGraphErrors(numEvents,x,y_lanWidth,x_err,y_err_lanWidth);
  TGraph *graph_lanWidth_shade = new TGraph(2*numEvents);
  for (Int_t i=0;i<numEvents;i++) {
    graph_lanWidth_shade->SetPoint(i,x[i],y_lanWidth[i]+1*y_err_lanWidth[i]);
    graph_lanWidth_shade->SetPoint(numEvents+i,x[numEvents-i-1],y_lanWidth[numEvents-i-1]-1*y_err_lanWidth[numEvents-i-1]);
  }
  graph_lanWidth_shade->SetFillStyle(3013);
  graph_lanWidth_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_lanWidth);

  graph_lanWidth->SetName("graph_lanWidth");
  graph_lanWidth->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  graph_lanWidth->GetYaxis()->SetTitle("Landau width of SIP_{3D}");
  graph_lanWidth->SetTitle("#sigma_{L}(SIP_{3D}) vs. time");
  graph_lanWidth->SetMarkerStyle(20);
  graph_lanWidth->SetMarkerSize(1.2);
  graph_lanWidth->SetMarkerColor(kRed);
  graph_lanWidth->SetLineColor(kRed);

  c4->cd();
  // dummyGraph->GetYaxis()->SetRangeUser(min_scale,max_scale);
  graph_lanWidth->Draw("alpsame");
  graph_lanWidth_shade->Draw("fsame");

  fout->cd();
  graph_lanWidth->Write();

   //******************************

  TCanvas *c5 = new TCanvas("c5","", 1000,500);
  c5->SetFillColor(0);
  c5->SetBottomMargin(0.14);
  c5->SetLeftMargin(0.07);
  c5->SetRightMargin(0.05);
  c5->SetTopMargin(0.08);  
  c5->SetGrid();

  TGraphErrors *graph_fitChisq = new TGraphErrors(numEvents,x,y_fitChisq,x_err,y_fake);
  makeNiceTGraph(graph_fitChisq);

  graph_fitChisq->SetName("graph_fitChisq");
  graph_fitChisq->GetXaxis()->SetTitle("Luminosity (fb^{-1})");
  graph_fitChisq->GetYaxis()->SetTitle("Fit #chi^{2}/ndf Landau fit");
  graph_fitChisq->SetTitle("#chi^{2}/ndf vs. time");
  graph_fitChisq->SetMarkerStyle(20);
  graph_fitChisq->SetMarkerSize(1.2);
  graph_fitChisq->SetMarkerColor(kRed);
  graph_fitChisq->SetLineColor(kRed);

  c5->cd();
  // dummyGraph->GetYaxis()->SetRangeUser(min_scale,max_scale);
  graph_fitChisq->Draw("alpsame");
 
  fout->cd();
  graph_fitChisq->Write();

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

