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
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"

//____________________________________________________________
// Auxilliary struct to store the fit results
struct Bands 
{
  int run;
  double mean_;
  double rms_;
  double meanError_;
  double rmsError_;
  double fitChi2ndof;
  vector<double> centrals;
  vector<double> errors;
  void init(){
    mean_=0;	     
    rms_=0;	     
    meanError_=0;  
    rmsError_=0;

    fitChi2ndof = -999.;
    run = 0; 
    centrals.clear();
    errors.clear();
  }
  void set(){
    fitChi2ndof = - 999.;
    run = 0;
    
    mean_=-999.;	     
    rms_=-999.;	     
    meanError_=-999.;  
    rmsError_=-999.;

    for(Int_t i=0;i<4;i++){
      centrals.push_back(-999.);
      errors.push_back(-999.);
    }
  }
};

//*************************************************************
Double_t CBshape(Double_t *x, Double_t *par)  
//*************************************************************
{
  Double_t var;
  Double_t t = (x[0]-par[1])/par[2];
  if (t > (-par[3])){
    var = par[0]*TMath::Exp(-t*t/2.);
  } else if (t <= (-par[3])) {
    Double_t AA = TMath::Power(par[4]/TMath::Abs(par[3]),par[4])*TMath::Exp(-TMath::Abs(par[3])*TMath::Abs(par[3])/2.);
    Double_t BB = par[4]/TMath::Abs(par[3])-TMath::Abs(par[3]);
    if(TMath::Power((BB-t),par[4])!=0){
      var = par[0]*AA/TMath::Power((BB-t),par[4]);
    } else   var = 0;
  }   
  return var;
}

//*************************************************************
double landaufun(double *x,double *par){
//*************************************************************
   double t;
   t=par[0]*TMath::Landau(x[0],par[1],par[2],0); //p1 is MP,p2 is sigma
   return t;
}

//*************************************************************
Double_t moyfun(Double_t *x, Double_t *par) {
//*************************************************************
  /* Moyal distribution, i.e. "approximate Landau function" (1955)
   *
   * Function: f(x) = \frac{1}{\sqrt{2\pi}} e^{-\frac{1}{2}(x+e^{-x})}, normalized to 1 
   * MPV:      x_{mpv} = 0, f(x_{mpv}) = \frac{1}{\sqrt{2 e \pi}} = 0.241971 
   *
   * As a fitting function I substitute x with (x - mpv)/width and introduce a normalization.  
   * par[0] = Total area (integral -inf to inf, normalization constant) 
   * par[1] = Most Probable Value
   * par[2] = Width (scale) parameter
   */  
  Double_t area   = par[0];
  Double_t lambda = (x[0] - par[1])/par[2];
  return (area/sqrt(2*TMath::Pi()))*TMath::Exp(-0.5*(lambda + TMath::Exp(-lambda)));
}

//*************************************************************
Double_t moygaufun(Double_t *x, Double_t *par) {
//*************************************************************
  /* Convolution of the Moyal distribution with gaussian
   *
   * par[0] = Width (scale) parameter of Moyal distribution
   * par[1] = Most Probable Value
   * par[2] = Total area (integral -inf to inf, normalization constant)
   * par[3] = Width (sigma) of convoluted Gaussian function
   */
  // control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  // variables
  Double_t xx;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  Double_t lambda;
  // range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
  // integral
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    lambda = (xx - par[1])/par[0];
    fland = TMath::Exp(-0.5*(lambda + TMath::Exp(-lambda)));
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    lambda = (xx - par[1])/par[0];
    fland = TMath::Exp(-0.5*(lambda + TMath::Exp(-lambda)));
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  return (par[2] * step * sum / sqrt(2*TMath::Pi()));
}

//*************************************************************
Bands* myFit(TString file,TCanvas& dummyC) 
//*************************************************************
{

  Bands* result = new Bands();

  TFile *f = TFile::Open(file);
  if (f->IsZombie()) {
    //something very wrong, cannot use this file, exit
    std::cout<< "puppa" << std::endl;
    return result;
  }

  TH1F *hSNR =  (TH1F*)f->Get("PVValidation/ProbeTrackFeatures/h_probeRefitVSig3D");
  // Fitting SNR histo
  printf("Fitting...\n");
  
  TObjArray *tokens = file.Tokenize("_");
  TString Alignment = TString(tokens->At(1)->GetName());
  TString Run       = (TString(tokens->At(2)->GetName())).ReplaceAll(".root","");
              
  Double_t chisqr;
  Int_t    ndf;

  //Bands* result = new Bands();

  // make sure there is at least 1000k tracks in the plot
  if(hSNR->GetEntries()<5000){
    std::cout<<"not enough tracks to make sensible plot"<<std::endl;
    result->set();
    return result;
  }

  // Setting fit range and start values
  Double_t fr[2];
  fr[0] = 0.15;
  fr[1] = 2.5;

  /*

    // for CB fit 

    Double_t sv[5], pllo[5], plhi[5];
    pllo[0]=0.1;         // Normalization
    pllo[1]=0.01;        // peak position
    pllo[2]=0.01;        // width
    pllo[3]=0.01;        // CB tails
    pllo[4]=0.01;        // CB tails       
    
    plhi[0]=10000000.0;  // Normalization	 
    plhi[1]=10.0;        // peak position	 
    plhi[2]=10.0;        // width		 
    plhi[3]=10.0;	       // CB tails	 
    plhi[4]=10.0;        // CB tails       
    
    sv[0]=1000.0;        // Normalization	 
    sv[1]=0.7; 	       // peak position	 
    sv[2]=0.08;          // width		 
    sv[3]=1.5;	       // CB tails	 
    sv[4]=0.1;	       // CB tails       
  
    // fitsnr->SetParameter(0,10000); // normalization  
    // fitsnr->SetParameter(1,0.7);   // peak position 
    // fitsnr->SetParameter(2,0.08);  // gaussian width
    // fitsnr->FixParameter(3,0.98);  // Crystal Ball tails
    // fitsnr->FixParameter(4,5.2);   // Crystal Ball tails
    
    TF1 *fitsnr = new TF1("fitsnr",CBshape,fr[0],fr[1],5); 
    fitsnr->SetParameters(sv);
    for(int i=0;i<5;i++){
    fitsnr->SetParLimits(i,pllo[i],plhi[i]);
    }
    
  */

  /*
    // for moyal function fit

    Double_t sv[3], pllo[3], plhi[3];
    pllo[0]=1000;         // Normalization
    pllo[1]=0.3;          // peak position
    pllo[2]=0.1;          // width
    
    plhi[0]=999999999;   // Normalization	 
    plhi[1]=5.0;         // peak position	 
    plhi[2]=10.0;        // width		 
    
    sv[0]=50000.0;       // Normalization	 
    sv[1]=1.0; 	       // peak position	 
    sv[2]=1.0;           // width		 
    
    TF1 *fitsnr = new TF1("fitsnr",moyfun,fr[0],fr[1],3);  
    fitsnr->SetParameters(sv);
    for(int i=0;i<3;i++){
    fitsnr->SetParLimits(i,pllo[i],plhi[i]);
    }
  */

  Double_t sv[4], pllo[4], plhi[4];
  pllo[0]=0.1;        
  pllo[1]=0.1;         
  pllo[2]=10000;         
  pllo[3]=0.1;         
  
  plhi[0]=10.0;   
  plhi[1]=5.0;         
  plhi[2]=999999999;        
  plhi[3]=10.0;        
  
  sv[0]=1.0;       
  sv[1]=1.0; 	       
  sv[2]=50000.0;           
  sv[3]=1.0;           
  
  TF1 *fitsnr = new TF1("fitsnr",moygaufun,fr[0],fr[1],4);  
  fitsnr->SetParameters(sv);
  for(int i=0;i<4;i++){
    fitsnr->SetParLimits(i,pllo[i],plhi[i]);
  }

  fitsnr->SetLineColor(kBlue);
  fitsnr->SetLineWidth(2);

  TFitResultPtr r = hSNR->Fit(fitsnr,"RIlS");
  chisqr = fitsnr->GetChisquare();
  ndf    = fitsnr->GetNDF();

  printf("Fitting done\nPlotting results...\n");
  
  // Global style settings
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
  
  hSNR->GetXaxis()->SetRange(0,70);
  hSNR->SetMarkerStyle(20);
  hSNR->SetMarkerSize(1);
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

  result->mean_      = hSNR->GetMean();
  result->meanError_ = hSNR->GetMeanError();
  result->rms_       = hSNR->GetRMS();
  result->rmsError_  = hSNR->GetRMSError();

  if(fitsnr){
    // fit parameters
    result->fitChi2ndof = chisqr/ndf;
    result->centrals.push_back(fitsnr->GetParameter(0));
    result->centrals.push_back(fitsnr->GetParameter(1));
    result->centrals.push_back(fitsnr->GetParameter(2));
    result->centrals.push_back(fitsnr->GetParameter(3));
    
    // fit parameters
    result->errors.push_back(fitsnr->GetParError(0));
    result->errors.push_back(fitsnr->GetParError(1));
    result->errors.push_back(fitsnr->GetParError(2));
    result->errors.push_back(fitsnr->GetParError(3));

  } else {
    result->set();
  }

  return result;

}

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

  Bands* result = new Bands();

  TFile *f = TFile::Open(file);
  if (f->IsZombie()) {
    //something very wrong, cannot use this file, exit
    std::cout<< "puppa" << std::endl;
    return result;
  }

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
  plhi[2]=9999999.0; 
  plhi[3]=10.0;
  
  sv[0]=0.5; 
  sv[1]=1.0; 
  sv[2]=25000.0; 
  sv[3]=1.0;
              
  Double_t chisqr;
  Int_t    ndf;

  //Bands* result = new Bands();

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

  result->mean_      = hSNR->GetMean();
  result->meanError_ = hSNR->GetMeanError();
  result->rms_       = hSNR->GetRMS();
  result->rmsError_  = hSNR->GetRMSError();

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
  //int startRun  = 274422;
  int endRun   = 284044;

  temp="fit_"+tag+".root";
  TFile *fileOUT = new TFile(temp.c_str(),"RECREATE");
  TTree *t = new TTree("SIPValidationTree","tree w/ fit results");

  std::vector<TH1F*> inputTH1s;
  
  double mean,meanError;
  double RMS,RMSError;
  double lanWidth,lanMPV,area,GWidth;
  double err_lanWidth,err_lanMPV,err_area,err_GWidth;
  double fitChisqNdof;

  TH1F* n_chiSq = new TH1F("fitChi2","fitChi2",1000,0.,1000.);

  int run; 

  t->Branch("run"             , &run         );

  t->Branch("mean"            , &mean        );
  t->Branch("rms"             , &RMS         );

  t->Branch("lanWidth"        , &lanWidth    );
  t->Branch("lanMPV"          , &lanMPV      ); 
  t->Branch("area"            , &area        );
  t->Branch("GWidth"          , &GWidth      );

  t->Branch("err_lanWidth"    , &err_lanWidth);
  t->Branch("err_lanMPV"      , &err_lanMPV  );
  t->Branch("err_area"        , &err_area    );
  t->Branch("err_GWidth"      , &err_GWidth  );
  t->Branch("fitChisqNdof"    , &fitChisqNdof);

  t->Branch("meanError"       , &meanError   );
  t->Branch("rmsError"        , &RMSError    );


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
      
      // langaussian fit
      // fitParams=langaus(filename,dummyC);
      
      TFile *fIN = TFile::Open(filename.c_str(),"READ");
      TH1F *hSNR =  (TH1F*)fIN->Get("PVValidation/ProbeTrackFeatures/h_probeRefitVSig3D");
      TH1F* myClone = (TH1F*)hSNR->Clone();
      myClone->SetName(Form("h_%s_%s",hSNR->GetName(),runStr.c_str()));
      myClone->SetDirectory(0);
      inputTH1s.push_back(myClone);
      fIN->Close();

      // moyal function fit
      fitParams=myFit(filename,dummyC);
      cout << "check fitParams size: " << fitParams->centrals.size() << endl;
       
      // fill the moments of the distribution
      mean      = fitParams->mean_;
      meanError = fitParams->meanError_;
      RMS       = fitParams->rms_;
      RMSError  = fitParams->rmsError_;
      
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
      
    } else{
      cout << "ERROR: the file for this run does not appear to exist..." << endl;

      mean             = -999.9;
      meanError        = -999.9;
      RMS              = -999.9;
      RMSError         = -999.9;

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
    if(fitChisqNdof>0.)
      n_chiSq->Fill(fitChisqNdof);


    if(DEBUG)
      cout << "checkpoint: ITERATING i" << endl;
    
    startRun++;
    //delete fitParams;
    gROOT->Reset();
    
  } // end of loop over runs

  dummyC.Print("fits.pdf]");

  TCanvas* chi2Canv = new TCanvas("chi2Canv","chi2Canv",800,600);
  chi2Canv->cd()->SetLogy();
  n_chiSq->Draw();
  chi2Canv->SaveAs("nChi2Sq.pdf");

  fileOUT->cd(); 
  t->Write();

  fileOUT->mkdir("Histograms");
  fileOUT->cd("Histograms");
  
  TH2F* SIPvsRun = new TH2F("SIPvsRun","SIP vs run number;run number;SIP;tracks",inputTH1s.size(),0.,inputTH1s.size(),100,0.,100.);


  for ( unsigned int ih=0; ih< inputTH1s.size(); ih++ ){
    std::cout<<ih<<std::endl;
    inputTH1s[ih]->Write();
    for(int jS=1;jS<=100; jS++){
      //Float_t normBinCont = inputTH1s[ih]->GetBinContent(jS)/inputTH1s[ih]->GetEntries();
      //SIPvsRun->SetBinContent(ih,jS,normBinCont);
      SIPvsRun->SetBinContent(ih,jS,inputTH1s[ih]->GetBinContent(jS));
    }
  }

  TProfile *SIPvsRun_p = (TProfile*) SIPvsRun->ProfileX("_pfx",1,-1,"o");
  SIPvsRun->Write();
  SIPvsRun_p->Write();

  fileOUT->Close();

}
