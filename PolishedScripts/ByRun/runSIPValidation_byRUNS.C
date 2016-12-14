#include <Riostream.h>
#include <string>
#include "TROOT.h"
#include <vector>
#include <sstream>
#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"

struct Bands {
  int run;
  double fitChi2ndof;
  vector<double> centrals;
  vector<double> errors;
  void init(){
    fitChi2ndof = -999.;
    run = 0; 
    centrals.clear();
    errors.clear();
  }
  void set(){
    fitChi2ndof = - 999.;
    run = 0;
    for(Int_t i=0;i<4;i++){
      centrals.push_back(-999.);
      errors.push_back(-999.);
    }
  }
};

//*************************************************************
Double_t langaufun(Double_t *x, Double_t *par)  
//*************************************************************
{

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;     // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  return (par[2] * step * sum * invsq2pi / par[3]);
}


//*************************************************************
TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
//*************************************************************
{
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[4]  reasonable start values for the fit
  //   parlimitslo[4]  lower parameter limits
  //   parlimitshi[4]  upper parameter limits
  //   fitparams[4]    returns the final fit parameters
  //   fiterrors[4]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf
  
  Int_t i;
  Char_t FunName[100];
  
  sprintf(FunName,"Fitfcn_%s",his->GetName());
  
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;
  
  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma");
  
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  
  his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
  
  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  
  return (ffit);              // return fit function
  
}


//*************************************************************
Bands* langaus(TString file,TCanvas& dummyC) 
//*************************************************************
{
  
  TFile *f = TFile::Open(file);
  TH1F *hSNR =  (TH1F*)f->Get("PVValidation/ProbeTrackFeatures/h_probeRefitVSig3D");
  // Fitting SNR histo
  printf("Fitting...\n");
  
  TObjArray *tokens = file.Tokenize("_");
  TString Alignment = TString(tokens->At(1)->GetName());
  TString Run       = (TString(tokens->At(2)->GetName())).ReplaceAll(".root","");

  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  //fr[0]=0.01*hSNR->GetMean();
  //fr[1]=10.0*hSNR->GetMean();
  
  fr[0] = 0.15;
  fr[1] = 3.5;

  pllo[0]=0.1; 
  pllo[1]=0.1; 
  pllo[2]=1000.0; 
  pllo[3]=0.1;
  
  plhi[0]=10.0; 
  plhi[1]=10.0; 
  plhi[2]=100000.0; 
  plhi[3]=10.0;
  
  sv[0]=0.5; 
  sv[1]=1.0; 
  sv[2]=25000.0; 
  sv[3]=1.0;
  
  Double_t chisqr;
  Int_t    ndf;

  Bands* result = new Bands();

  // make sure there is at least 1000k tracks in the plot
  if(hSNR->GetEntries()<5000){
    std::cout<<"not enough tracks to make sensible plot"<<std::endl;
    result->set();
    return result;
  }

  TF1 *fitsnr = langaufit(hSNR,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  
  Double_t SNRPeak, SNRFWHM;
  //langaupro(fp,SNRPeak,SNRFWHM);
  
  printf("Fitting done\nPlotting results...\n");
  
  // Global style settings
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
  
  hSNR->GetXaxis()->SetRange(0,70);
  hSNR->Draw();
  fitsnr->Draw("lsame");

  TPaveText *pt = new TPaveText(0.2,0.75,0.5,0.87,"NDC");
  pt->SetFillColor(10);
  pt->SetTextColor(1);
  pt->SetTextFont(42);
  pt->SetTextAlign(11);
  TText *text1 = pt->AddText(Alignment);
  TText *text2 = pt->AddText("\nRun: "+Run);
  text1->SetTextSize(0.04);
  text1->SetTextColor(4);
  text2->SetTextSize(0.04);
  text2->SetTextColor(4);
  pt->Draw("same");

  dummyC.Print("fits.pdf");

  if(fitsnr){
    result->fitChi2ndof = chisqr/ndf;
    result->centrals.push_back(fitsnr->GetParameter(0));
    result->centrals.push_back(fitsnr->GetParameter(1));
    result->centrals.push_back(fitsnr->GetParameter(2));
    result->centrals.push_back(fitsnr->GetParameter(3));
    
    result->errors.push_back(fitsnr->GetParError(0));
    result->errors.push_back(fitsnr->GetParError(1));
    result->errors.push_back(fitsnr->GetParError(2));
    result->errors.push_back(fitsnr->GetParError(3));
  } else {
    result->set();
  }

  return result;

}

// 
// THIS IS THE MAIN FUNCTION
//

//*************************************************************
void runSIPValidation_byRUNS(string tag="PromptGT")
//*************************************************************
{
  //-----User set variables---------
  string temp;
  string directory=".";

  string stem="/PVValidation_"+tag+"_";
  string filename;

  int startRun = 272007;
  int endRun   = 274955;

  temp=tag+".root";
  TFile *file = new TFile(temp.c_str(),"RECREATE");
  TTree *t = new TTree("SIPValidationTree","tree w/ fit results");

  double lanWidth,lanMPV,area,GWidth;
  double err_lanWidth,err_lanMPV,err_area,err_GWidth;
  double fitChisqNdof;

  int run; 

  t->Branch("run"             , &run         );

  t->Branch("lanWidth"        , &lanWidth    );
  t->Branch("lanMPV"          , &lanMPV      ); 
  t->Branch("area"            , &area        );
  t->Branch("GWidth"          , &GWidth      );

  t->Branch("err_lanWidth"    , &err_lanWidth);
  t->Branch("err_lanMPV"      , &err_lanMPV  );
  t->Branch("err_area"        , &err_area    );
  t->Branch("err_GWidth"      , &err_GWidth  );
  t->Branch("fitChisqNdof"    , &fitChisqNdof);

  //-----END user set variables ----

  bool DEBUG=true;
 
  bool fileExists;
  ifstream inputFile;
  stringstream out;
  string runStr;

  TCanvas dummyC;
  dummyC.Print("fits.pdf[");

  while(startRun<=endRun){

    //------ parse string and build filename-----
    
    Bands *fitParams = new Bands();
    fitParams->init();
    
    filename=directory+stem;
    
    cout << "The run is: " << startRun << endl;
    out.str("");
    out << startRun;
    runStr = out.str();
    out.str("");
    
    filename=filename+runStr+".root";
    
    cout << "Filename: " << filename << endl;      
    
    // -----------check if file exists -----------
    
    inputFile.open(filename.c_str());
    if( inputFile.is_open() ){
      fileExists=true;
      inputFile.close();
      if(DEBUG)
	cout << "checkpoint: FILE EXISTS" << endl;
    }else{
      fileExists=false;
      if(DEBUG)
	cout << "checkpoint: FILE DOES NOT EXIST" << endl;
      inputFile.close();
    } // end of else
      // ---------------- checked ----------------
    
      // ---------- Open file and run Plotter -------
      //fileExists=true;
    if(fileExists){
      // filename=filename+"="+runStr;
      
      if(DEBUG)
	cout << "checkpoint: RUNNING PLOTTER WITH FILENAME: " << filename << endl;
      
      fitParams=langaus(filename,dummyC);
      
      cout << "check fitParams size: " << fitParams->centrals.size() << endl;
      
      if(fitParams->centrals.size()==4){
	
	if(DEBUG)
	  cout << "checkpoint: RETURN VALUE EXISTS" << endl;
	
	run             = startRun;
	fitChisqNdof    = fitParams->fitChi2ndof;
	lanWidth        = fitParams->centrals[0];
	lanMPV          = fitParams->centrals[1];
	area            = fitParams->centrals[2];
	GWidth          = fitParams->centrals[3];
	                  
	err_lanWidth    = fitParams->errors[0];
	err_lanMPV      = fitParams->errors[1];
	err_area        = fitParams->errors[2];
	err_GWidth      = fitParams->errors[3];

      } // end of check fitParams
      else
	cout << "ERROR: the fit results return null pointer..." << endl;
    }else{
      cout << "ERROR: the file for this run does not appear to exist..." << endl;
      run              = startRun;

      fitChisqNdof     = -999.9;
      lanWidth         = -999.9; 
      lanMPV           = -999.9;
      area             = -999.9;
      GWidth           = -999.9;

      err_lanWidth     = -999.9;
      err_lanMPV       = -999.9;
      err_area         = -999.9;
      err_GWidth       = -999.9;

    } // end of else
      //-------------------done plotting ------------
    
      //------------Fill tree with above data---------
    
      // -----------------Filled --------------------
    if(DEBUG)
      cout << "checkpoint: FILLING TREE" << endl;
    t->Fill();
    
    if(DEBUG)
      cout << "checkpoint: ITERATING i" << endl;
    
    startRun++;
    //delete fitParams;
    gROOT->Reset();
    
  } // end of loop over runs

  dummyC.Print("fits.pdf]");

  file->cd(); 
  t->Write();
  file->Close();
}
