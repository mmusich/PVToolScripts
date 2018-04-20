#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TObject.h"
#include <iostream>
#include "TLegend.h"
#include "TClass.h"
#include "TRatioPlot.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"

TFile *sourceFile1, *sourceFile2;

std::pair<Double_t,Double_t> getExtrema(TObjArray *array);
template<typename T> void MakeNicePlotStyle(T *hist);

//void MakeNicePlotStyle(TH1 *hist);
void plot2Histograms(TH1* h1, TH1* h2, const TString& label1,const TString& label2);
void recurseOverKeys( TDirectory *target1, const TString& label1,const TString& label2 );

/************************************************/
void loopAndPlot(TString filename1, TString filename2, TString label1,TString label2)
/************************************************/
{

  sourceFile1 = TFile::Open(filename1,"r");
  sourceFile2 = TFile::Open(filename2,"r");

  //for(auto&& keyAsObj : *file->GetListOfKeys()){
  //  auto key = (TKey*) keyAsObj;
  // std::cout << key->GetName() << " " << key->GetClassName() << std::endl;
  //}
  
  recurseOverKeys(sourceFile1,label1,label2);

  sourceFile1->Close();
  sourceFile2->Close();
}


/************************************************/
void recurseOverKeys( TDirectory *target1,const TString& label1,const TString& label2) 
/************************************************/
{
  // Figure out where we are
  TString path( (char*)strstr( target1->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  sourceFile1->cd( path );
  
  std::cout<<path<<std::endl;

  TDirectory *current_sourcedir = gDirectory;

  TKey *key;
  TIter nextkey(current_sourcedir->GetListOfKeys());

  while ( (key = (TKey*)nextkey()) ) {

    auto obj = key->ReadObj();

    // Check if this is a 1D histogram or a directory
    if (obj->IsA()->InheritsFrom("TH1F")) {

      // **************************
      // Plot & Save this Histogram
      TH1F *htemp1, *htemp2;

      htemp1 = (TH1F*)obj;
      TString histName = htemp1->GetName();

      if (path != "") {
        sourceFile2->GetObject(path+"/"+histName, htemp2);
      } else {
	sourceFile2->GetObject(histName, htemp2);
      }

      //outputFilename=histName;
      //plot2Histograms(htemp1, htemp2, outputFolder+path+"/"+outputFilename+"."+imageType);
      plot2Histograms(htemp1, htemp2,label1,label2);
      

    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory

      cout << "Found subdirectory " << obj->GetName() << endl;
      //gSystem->MakeDirectory(outputFolder+path+"/"+obj->GetName());

      // obj is now the starting point of another round of merging
      // obj still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion

      if( (TString(obj->GetName())).Contains("Residuals") ) continue;

      recurseOverKeys( (TDirectory*)obj , label1,label2);

    } // end of IF a TDriectory 
  }
}

/************************************************/
void plot2Histograms(TH1* h1, TH1* h2, const TString& label1,const TString& label2) {
/************************************************/

  TGaxis::SetMaxDigits(3);

  auto c1 = new TCanvas(Form("c1_%s",h1->GetName()), "A ratio example",800,800);
  gStyle->SetOptStat(0);

  MakeNicePlotStyle<TH1>(h1);
  MakeNicePlotStyle<TH1>(h2);

  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  TObjArray *array = new TObjArray(2); 
  array->Add(h1);
  array->Add(h2);

  std::pair<Double_t,Double_t> extrema =  getExtrema(array);
 
  delete array;
  
  float min = (extrema.first>0) ? (extrema.first)*0.7 : (extrema.first)*1.3;

  h1->GetYaxis()->SetRangeUser(min,extrema.second*1.3);
  h2->GetYaxis()->SetRangeUser(min,extrema.second*1.3);
  
  auto rp = new TRatioPlot(h1, h2);
  c1->SetTicks(0, 1);
  rp->Draw();

  //rp->GetUpperPad()->SetTopMargin(0.09);
  //rp->GetUpperPad()->SetLeftMargin(0.15);
  //rp->GetUpperPad()->SetRightMargin(0.03);
  //rp->GetLowerPad()->SetBottomMargin(0.5);
  
  rp->SetLeftMargin(0.15); 
  rp->SetRightMargin(0.03);
  rp->SetSeparationMargin(0.01);
  rp->SetLowBottomMargin(0.35); 

  rp->GetUpperPad()->cd();
  // Draw the legend
  TLegend *infoBox = new TLegend(0.75, 0.75, 0.97, 0.90, "");
  infoBox->AddEntry(h1,label1,"L");
  infoBox->AddEntry(h2,label2,"L");
  infoBox->SetShadowColor(0);  // 0 = transparent
  infoBox->SetFillColor(kWhite); 
  infoBox->Draw("same");
  
  MakeNicePlotStyle<TGraph>(rp->GetLowerRefGraph());
  rp->GetLowerRefGraph()->GetYaxis()->SetTitle("ratio");
  rp->GetLowerRefGraph()->SetMinimum(0.);
  rp->GetLowerRefGraph()->SetMaximum(2.);
  c1->Update();

  //rp->GetLowerPad()->cd();
  //c1->Update();

  c1->SaveAs(TString(h1->GetName())+".png");
  delete c1;
}

/*--------------------------------------------------------------------*/
template<typename T>
void MakeNicePlotStyle(T *hist)
/*--------------------------------------------------------------------*/
{ 
  //hist->SetStats(kFALSE);  
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetNdivisions(505);
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  hist->GetXaxis()->SetTitleFont(42); 
  hist->GetYaxis()->SetTitleFont(42);  
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  if( ((TObject*)hist)->IsA()->InheritsFrom("TGraph") ){
    hist->GetYaxis()->SetLabelSize(.025);
    //hist->GetYaxis()->SetNdivisions(505);
  } else {
    hist->GetYaxis()->SetLabelSize(.05);
  }
  hist->GetXaxis()->SetLabelSize(.05);
}

//*****************************************************//
std::pair<Double_t,Double_t> getExtrema(TObjArray *array)
//*****************************************************//
{
  Double_t theMaximum = (static_cast<TH1*>(array->At(0)))->GetMaximum();
  Double_t theMinimum = (static_cast<TH1*>(array->At(0)))->GetMinimum();
  for(Int_t i = 0; i< array->GetSize(); i++){
    if( (static_cast<TH1*>(array->At(i)))->GetMaximum() > theMaximum){
      theMaximum = (static_cast<TH1*>(array->At(i)))->GetMaximum();
    }
    if ( (static_cast<TH1*>(array->At(i)))->GetMinimum() < theMinimum){
      theMinimum = (static_cast<TH1*>(array->At(i)))->GetMinimum();
    }
  }
  return std::make_pair(theMinimum,theMaximum);
}
