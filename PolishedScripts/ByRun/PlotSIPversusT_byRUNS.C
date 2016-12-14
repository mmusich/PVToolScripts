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
//#include "tdrstyle.C"

void makeNiceTGraph(TGraph *graph);

void PlotSIPversusT_byRUNS(string tag="PromptGT"){

  //gROOT->ProcessLine(".L TkAlStyle.cc+");
  //TkAlStyle::set(INPROGRESS);	// set publication status

  //gROOT->ProcessLine(".L ./setTDRStyle.C++");
  //setTDRStyle();

  //string tag="ORIGINAL";//"TBDkb";
  //ORIGINAL_DzDisparity_byRUNS.root
  //TBDkb_newIOV_DzDisparity_byRUNS_fromJuly1.root

  string inFileName=tag+".root";
  TFile *f = new TFile(inFileName.c_str());
  TTree* t = (TTree*) f->Get("SIPValidationTree");

  const int DummyNumEvents = int(t->GetEntries());
  const int NumEvents=DummyNumEvents;
  cout<<"NumEvents = "<<NumEvents <<endl;
  
  double lanWidth,lanMPV,area,GWidth;
  double err_lanWidth,err_lanMPV,err_area,err_GWidth;
  double fitChisqNdof;

  int run;

  double x[NumEvents];
  double x_err[NumEvents];
  double y_fake[NumEvents];
  double y_lanWidth[NumEvents];
  double y_lanMPV[NumEvents];
  double y_area[NumEvents];
  double y_GWidth[NumEvents];
  double y_err_lanWidth[NumEvents];
  double y_err_lanMPV[NumEvents];
  double y_err_area[NumEvents];
  double y_err_GWidth[NumEvents];

  double noDatax[10000];
  double noDatay[10000];

  t->SetBranchAddress("run"             , &run         );
  t->SetBranchAddress("lanWidth"        , &lanWidth    );                                                                         
  t->SetBranchAddress("lanMPV"          , &lanMPV      );                                                                         
  t->SetBranchAddress("area"            , &area        );                                                                         
  t->SetBranchAddress("GWitdh"          , &GWidth      );                                                                         
                                                                                                                        
  t->SetBranchAddress("err_lanWidth"    , &err_lanWidth);                                                                         
  t->SetBranchAddress("err_lanMPV	", &err_lanMPV  );                                                                         
  t->SetBranchAddress("err_area"        , &err_area    );                                                                         
  t->SetBranchAddress("err_GWitdh"      , &err_GWidth  );                                                                         
  t->SetBranchAddress("fitChisqNdof"    , &fitChisqNdof);  

  //-----------------------------------------------

  int numEventsGood=0;
  int numEventsBad=0;
  
  t->BuildIndex("run");
  TTreeIndex *index = (TTreeIndex*)t->GetTreeIndex();
  for( int i = 0; i <= index->GetN(); i++ ) {
    Long64_t local = t->LoadTree( index->GetIndex()[i] );
   
    t->GetEntry(local);

    if(fitChisqNdof<0.)
      continue;
    
    if(fitChisqNdof>100.){
      noDatay[numEventsBad] = -999.;
      noDatax[numEventsBad] = i;
      numEventsBad++;

    } else {

      y_fake   [numEventsGood]      = 0.;
      y_lanWidth[numEventsGood]     = lanWidth;
      y_lanMPV[numEventsGood]       = lanMPV;
      y_area[numEventsGood]         = area;
      y_GWidth[numEventsGood]       = GWidth;
      
      y_err_lanWidth[numEventsGood] = err_lanWidth;
      y_err_lanMPV[numEventsGood]   = err_lanMPV;
      y_err_area[numEventsGood]     = err_area;
      y_err_GWidth[numEventsGood]   = err_GWidth;
      
      //x[numEventsGood] = run;
      x_err[numEventsGood] = 0.;
      x[numEventsGood] = numEventsGood;

      cout<<" x[" << numEventsGood << "]" << x[numEventsGood] << " y_lanMPV["<<numEventsGood<<"] " <<  y_lanMPV[numEventsGood] << " +/- y_err_lanMPV["<<numEventsGood<<"] " <<   y_err_lanMPV[numEventsGood] << endl; 
      numEventsGood++;
      
    }

  }

  //////***************************************************///////

  Double_t max_scale = 0.9;
  Double_t min_scale = 0.4;

  TGraphErrors *dummyGraph = new TGraphErrors(numEventsGood,x,y_fake,0,0);

  TCanvas *c = new TCanvas("c","", 1000,500);
  c->SetFillColor(0);
  c->SetBottomMargin(0.14);
  c->SetLeftMargin(0.07);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.08);  
  c->SetGrid();

  TGraphErrors *graph_lanMPV = new TGraphErrors(numEventsGood,x,y_lanMPV,x_err,y_err_lanMPV);
  TGraph *graph_lanMPV_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_lanMPV_shade->SetPoint(i,x[i],y_lanMPV[i]+1*y_err_lanMPV[i]);
    graph_lanMPV_shade->SetPoint(numEventsGood+i,x[numEventsGood-i-1],y_lanMPV[numEventsGood-i-1]-1*y_err_lanMPV[numEventsGood-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_lanMPV_shade->SetFillStyle(3013);
  graph_lanMPV_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_lanMPV);

  graph_lanMPV->GetXaxis()->SetTitle("Run Number index");
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

  //******************************

  TCanvas *c2 = new TCanvas("c2","", 1000,500);
  c2->SetFillColor(0);
  c2->SetBottomMargin(0.14);
  c2->SetLeftMargin(0.07);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.08);  
  c2->SetGrid();

  TGraphErrors *graph_GWidth = new TGraphErrors(numEventsGood,x,y_GWidth,x_err,y_err_GWidth);
  TGraph *graph_GWidth_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_GWidth_shade->SetPoint(i,x[i],y_GWidth[i]+1*y_err_GWidth[i]);
    graph_GWidth_shade->SetPoint(numEventsGood+i,x[numEventsGood-i-1],y_GWidth[numEventsGood-i-1]-1*y_err_GWidth[numEventsGood-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_GWidth_shade->SetFillStyle(3013);
  graph_GWidth_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_GWidth);

  graph_GWidth->GetXaxis()->SetTitle("Run Number index");
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

  //******************************

  TCanvas *c3 = new TCanvas("c3","", 1000,500);
  c3->SetFillColor(0);
  c3->SetBottomMargin(0.14);
  c3->SetLeftMargin(0.07);
  c3->SetRightMargin(0.05);
  c3->SetTopMargin(0.08);  
  c3->SetGrid();

  TGraphErrors *graph_area = new TGraphErrors(numEventsGood,x,y_area,x_err,y_err_area);
  TGraph *graph_area_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_area_shade->SetPoint(i,x[i],y_area[i]+1*y_err_area[i]);
    graph_area_shade->SetPoint(numEventsGood+i,x[numEventsGood-i-1],y_area[numEventsGood-i-1]-1*y_err_area[numEventsGood-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_area_shade->SetFillStyle(3013);
  graph_area_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_area);

  graph_area->GetXaxis()->SetTitle("Run Number index");
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

  //******************************

  TCanvas *c4 = new TCanvas("c4","", 1000,500);
  c4->SetFillColor(0);
  c4->SetBottomMargin(0.14);
  c4->SetLeftMargin(0.07);
  c4->SetRightMargin(0.05);
  c4->SetTopMargin(0.08);  
  c4->SetGrid();

  TGraphErrors *graph_lanWidth = new TGraphErrors(numEventsGood,x,y_lanWidth,x_err,y_err_lanWidth);
  TGraph *graph_lanWidth_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_lanWidth_shade->SetPoint(i,x[i],y_lanWidth[i]+1*y_err_lanWidth[i]);
    graph_lanWidth_shade->SetPoint(numEventsGood+i,x[numEventsGood-i-1],y_lanWidth[numEventsGood-i-1]-1*y_err_lanWidth[numEventsGood-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_lanWidth_shade->SetFillStyle(3013);
  graph_lanWidth_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_lanWidth);

  graph_lanWidth->GetXaxis()->SetTitle("Run Number index");
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

