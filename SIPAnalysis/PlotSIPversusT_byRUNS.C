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

  string inFileName="fit_"+tag+".root";
  TFile *f = new TFile(inFileName.c_str());
  TTree* t = (TTree*) f->Get("SIPValidationTree");

  string outFileName=tag+"_byRun.root";
  TFile *fout = new TFile(outFileName.c_str(),"RECREATE");

  const int DummyNumEvents = int(t->GetEntries());
  const int NumEvents=DummyNumEvents;
  cout<<"NumEvents = "<<NumEvents <<endl;
  
  // initialize variables

  double mean,meanError;
  double rms,rmsError;
  double lanWidth,lanMPV,area,GWidth;
  double err_lanWidth,err_lanMPV,err_area,err_GWidth;
  double fitChisqNdof;

  int run;

  double x[NumEvents];
  double x_err[NumEvents];
  double y_fake[NumEvents];
  double y_mean[NumEvents];
  double y_rms[NumEvents];
  double y_lanWidth[NumEvents];
  double y_lanMPV[NumEvents];
  double y_area[NumEvents];
  double y_fitChisq[NumEvents];
  double y_GWidth[NumEvents];
  double y_err_mean[NumEvents];
  double y_err_rms[NumEvents];
  double y_err_lanWidth[NumEvents];
  double y_err_lanMPV[NumEvents];
  double y_err_area[NumEvents];
  double y_err_GWidth[NumEvents];

  double noDatax[10000];
  double noDatay[10000];

  t->SetBranchAddress("run"             , &run         );
  t->SetBranchAddress("mean"            , &mean        );
  t->SetBranchAddress("meanError"       , &meanError   );
  t->SetBranchAddress("rms"             , &rms        );
  t->SetBranchAddress("rmError"         , &rmsError   );
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

  int numEventsGood=0;
  int numEventsBad=0;
  
  t->BuildIndex("run");
  TTreeIndex *index = (TTreeIndex*)t->GetTreeIndex();
  for( int i = 0; i <= index->GetN(); i++ ) {
    Long64_t local = t->LoadTree( index->GetIndex()[i] );
   
    t->GetEntry(local);

    if(fitChisqNdof<0.)
      continue;
    
    // protect against crazy fits
    if(fitChisqNdof>300. || 
       (err_lanMPV/lanMPV > 0.5) || 
       (err_lanWidth/lanWidth > 0.5) ||
       (err_area/area > 0.5) ||
       (err_GWidth/GWidth > 0.5)
       )
      {
      noDatay[numEventsBad] = -999.;
      noDatax[numEventsBad] = i;
      numEventsBad++;

    } else {

      y_fake   [numEventsGood]      = 0.;
      y_lanWidth[numEventsGood]     = lanWidth;
      y_lanMPV[numEventsGood]       = lanMPV;
      y_area[numEventsGood]         = area;
      y_GWidth[numEventsGood]       = GWidth;
      y_fitChisq[numEventsGood]     = fitChisqNdof;
      y_mean[numEventsGood]         = mean;
      y_rms[numEventsGood]          = rms;

      y_err_lanWidth[numEventsGood] = err_lanWidth;
      y_err_lanMPV[numEventsGood]   = err_lanMPV;
      y_err_area[numEventsGood]     = err_area;
      y_err_GWidth[numEventsGood]   = err_GWidth;
      y_err_mean[numEventsGood]     = meanError;
      y_err_rms[numEventsGood]      = rmsError;

      x[numEventsGood] = run;
      x_err[numEventsGood] = 0.;
      //x[numEventsGood] = numEventsGood;

      cout<<" x[" << numEventsGood << "]" << x[numEventsGood] << " y_lanMPV["<<numEventsGood<<"] " <<  y_lanMPV[numEventsGood] << " +/- y_err_lanMPV["<<numEventsGood<<"] " <<   y_err_lanMPV[numEventsGood] << endl; 
      numEventsGood++;
      
    }

  }

  //////***************************************************///////

  Double_t max_scale = 0.9;
  Double_t min_scale = 0.4;

  TGraphErrors *dummyGraph = new TGraphErrors(numEventsGood,x,y_fake,0,0);

  // -- mean

  TCanvas *cMean = new TCanvas("cMean","", 1000,500);
  cMean->SetFillColor(0);
  cMean->SetBottomMargin(0.14);
  cMean->SetLeftMargin(0.07);
  cMean->SetRightMargin(0.05);
  cMean->SetTopMargin(0.08);  
  cMean->SetGrid();

  TGraphErrors *graph_mean = new TGraphErrors(numEventsGood,x,y_mean,x_err,y_err_mean);
  TGraph *graph_mean_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_mean_shade->SetPoint(i,x[i],y_mean[i]+1*y_err_mean[i]);
    graph_mean_shade->SetPoint(numEventsGood+i,x[numEventsGood-i-1],y_mean[numEventsGood-i-1]-1*y_err_mean[numEventsGood-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_mean_shade->SetFillStyle(3013);
  graph_mean_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_mean);

  graph_mean->SetName("graph_mean");
  graph_mean->GetXaxis()->SetTitle("Run Number index");
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

  TGraphErrors *graph_rms = new TGraphErrors(numEventsGood,x,y_rms,x_err,y_err_rms);
  TGraph *graph_rms_shade = new TGraph(2*numEventsGood);
  for (Int_t i=0;i<numEventsGood;i++) {
    graph_rms_shade->SetPoint(i,x[i],y_rms[i]+1*y_err_rms[i]);
    graph_rms_shade->SetPoint(numEventsGood+i,x[numEventsGood-i-1],y_rms[numEventsGood-i-1]-1*y_err_rms[numEventsGood-i-1]);
    //graph_dxy_phi->GetXaxis()->SetBinLabel(i+1,Form("%i",x[i]));
  }
  graph_rms_shade->SetFillStyle(3013);
  graph_rms_shade->SetFillColor(kBlue);

  makeNiceTGraph(graph_rms);

  graph_rms->SetName("graph_rms");
  graph_rms->GetXaxis()->SetTitle("Run Number index");
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

  graph_lanMPV->SetName("graph_lanMPV");
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

  graph_GWidth->SetName("graph_GWidth");
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

  graph_area->SetName("graph_area");
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

  graph_lanWidth->SetName("graph_lanWidth");
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

  TGraphErrors *graph_fitChisq = new TGraphErrors(numEventsGood,x,y_fitChisq,x_err,y_fake);
  makeNiceTGraph(graph_fitChisq);

  graph_fitChisq->SetName("graph_fitChisq");
  graph_fitChisq->GetXaxis()->SetTitle("Run Number index");
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

