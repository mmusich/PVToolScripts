#include <Riostream.h>
#include <string>
#include "TROOT.h"
#include <vector>
#include <sstream>
#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "FitPVResiduals_forLoop2.C"

/*
struct Bands {
  int run;
  vector<double> mean;
  vector<double> width;
  void init(){
    mean.push_back(-999.);
    width.push_back(-999.);
  }
};
*/

void runPlotPVValidation_byRUNS(string tag="PromptGT"){

  // gROOT->ProcessLine(".L FitPVResiduals_forLoop2.C++");
  //-----User set variables---------
  string temp;
  string directory=".";

  string stem="/PVValidation_"+tag+"_";
  string filename;

  int startRun = 251027;
  //int endRun   = 251040;
  int endRun   = 251883;

  temp=tag+".root";
  TFile *file = new TFile(temp.c_str(),"RECREATE");
  TTree *t = new TTree("PVValidationTree","tree w/ fit results");

  double mean_dxy_phi,width_dxy_phi,mean_dz_phi,width_dz_phi,mean_dxy_eta,width_dxy_eta,mean_dz_eta,width_dz_eta;
  double min_dxy_phi,max_dxy_phi,min_dz_phi,max_dz_phi,min_dxy_eta,max_dxy_eta,min_dz_eta,max_dz_eta;

  double mean_n_dxy_phi,width_n_dxy_phi,mean_n_dz_phi,width_n_dz_phi,mean_n_dxy_eta,width_n_dxy_eta,mean_n_dz_eta,width_n_dz_eta;
  double min_n_dxy_phi,max_n_dxy_phi,min_n_dz_phi,max_n_dz_phi,min_n_dxy_eta,max_n_dxy_eta,min_n_dz_eta,max_n_dz_eta;  
  double dz_fit,dz_fit_error;

  int run; 
  double nevents,ntracks;

  t->Branch("run"             , &run         );
  t->Branch("nevents"         , &nevents     );
  t->Branch("ntracks"         , &ntracks     );

  t->Branch("dz_fit"          , &dz_fit);
  t->Branch("dz_fit_error"    , &dz_fit_error);

  t->Branch("mean_dxy_phi"    , &mean_dxy_phi);
  t->Branch("width_dxy_phi"   , &width_dxy_phi);
  t->Branch("min_dxy_phi"     , &min_dxy_phi);
  t->Branch("max_dxy_phi"     , &max_dxy_phi);
  
  t->Branch("mean_dz_phi"     , &mean_dz_phi);
  t->Branch("width_dz_phi"    , &width_dz_phi);
  t->Branch("min_dz_phi"      , &min_dz_phi);
  t->Branch("max_dz_phi"      , &max_dz_phi);
  
  t->Branch("mean_dxy_eta"    , &mean_dxy_eta);
  t->Branch("width_dxy_eta"   , &width_dxy_eta);  
  t->Branch("min_dxy_eta"     , &min_dxy_eta);
  t->Branch("max_dxy_eta"     , &max_dxy_eta);

  t->Branch("mean_dz_eta"     , &mean_dz_eta);
  t->Branch("width_dz_eta"    , &width_dz_eta);
  t->Branch("min_dz_eta"      , &min_dz_eta);
  t->Branch("max_dz_eta"      , &max_dz_eta);

  t->Branch("mean_n_dxy_phi"  , &mean_n_dxy_phi);
  t->Branch("width_n_dxy_phi" , &width_n_dxy_phi);
  t->Branch("min_n_dxy_phi"   , &min_n_dxy_phi);
  t->Branch("max_n_dxy_phi"   , &max_n_dxy_phi);
   
  t->Branch("mean_n_dz_phi"   , &mean_n_dz_phi);
  t->Branch("width_n_dz_phi"  , &width_n_dz_phi);
  t->Branch("min_n_dz_phi"    , &min_n_dz_phi);
  t->Branch("max_n_dz_phi"    , &max_n_dz_phi);

  t->Branch("mean_n_dxy_eta"  , &mean_n_dxy_eta);
  t->Branch("width_n_dxy_eta" , &width_n_dxy_eta);  
  t->Branch("min_n_dxy_eta"   , &min_n_dxy_eta);
  t->Branch("max_n_dxy_eta"   , &max_n_dxy_eta); 

  t->Branch("mean_n_dz_eta"   , &mean_n_dz_eta);
  t->Branch("width_n_dz_eta"  , &width_n_dz_eta);    
  t->Branch("min_n_dz_eta"    , &min_n_dz_eta);
  t->Branch("max_n_dz_eta"    , &max_n_dz_eta);   

  //-----END user set variables ----

  bool DEBUG=true;
 
  bool fileExists;
  ifstream inputFile;
  stringstream out;
  string runStr;

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
      filename=filename+"="+runStr;
      
      if(DEBUG)
	cout << "checkpoint: RUNNING PLOTTER WITH FILENAME: " << filename << endl;
      
      fitParams=FitPVResiduals_forLoop(filename,true,false,runStr);
      
      cout << "check fitParams size: " << fitParams->mean.size() << endl;
      
      if(fitParams->mean.size()==8){
	
	if(DEBUG)
	  cout << "checkpoint: RETURN VALUE EXISTS" << endl;
	
	run             = startRun;
	ntracks         = fitParams->ntracks;
	nevents         = fitParams->nevents;
	dz_fit          = fitParams->dzfit;
	dz_fit_error    = fitParams->dzfitError;

	mean_dxy_phi    = fitParams->mean[0];  	     
	width_dxy_phi   = fitParams->width[0]; 
	mean_dz_phi     = fitParams->mean[1];  	     
	width_dz_phi    = fitParams->width[1]; 	     
	mean_dxy_eta    = fitParams->mean[2];  	     
	width_dxy_eta   = fitParams->width[2]; 
	mean_dz_eta     = fitParams->mean[3];  	     
	width_dz_eta    = fitParams->width[3]; 	     
	mean_n_dxy_phi  = fitParams->mean[4];  	     
	width_n_dxy_phi = fitParams->width[4]; 
	mean_n_dz_phi   = fitParams->mean[5];  
	width_n_dz_phi  = fitParams->width[5]; 
	mean_n_dxy_eta  = fitParams->mean[6];  
	width_n_dxy_eta = fitParams->width[6]; 
	mean_n_dz_eta   = fitParams->mean[7];  
	width_n_dz_eta  = fitParams->width[7]; 
	
	min_dxy_phi    = fitParams->min[0];  	     
	max_dxy_phi    = fitParams->max[0]; 
	min_dz_phi     = fitParams->min[1];  	     
	max_dz_phi     = fitParams->max[1]; 	     
	min_dxy_eta    = fitParams->min[2];  	     
	max_dxy_eta    = fitParams->max[2]; 
	min_dz_eta     = fitParams->min[3];  	     
	max_dz_eta     = fitParams->max[3]; 	     
	min_n_dxy_phi  = fitParams->min[4];  	     
	max_n_dxy_phi  = fitParams->max[4]; 
	min_n_dz_phi   = fitParams->min[5];  
	max_n_dz_phi   = fitParams->max[5]; 
	min_n_dxy_eta  = fitParams->min[6];  
	max_n_dxy_eta  = fitParams->max[6]; 
	min_n_dz_eta   = fitParams->min[7];  
	max_n_dz_eta   = fitParams->max[7]; 
	
      } // end of check fitParams
      else
	cout << "ERROR: the fit results return null pointer..." << endl;
    }else{
      cout << "ERROR: the file for this run does not appear to exist..." << endl;
      run             = startRun;
      ntracks         = -999.9;
      nevents         = -999.9;
      dz_fit          = -999.9;
      dz_fit_error    = -999.9;
      mean_dxy_phi    = -999.9;	     
      width_dxy_phi   = -999.9;    
      mean_dz_phi     = -999.9;	     
      width_dz_phi    = -999.9;	     
      mean_dxy_eta    = -999.9;	     
      width_dxy_eta   = -999.9;    
      mean_dz_eta     = -999.9;	     
      width_dz_eta    = -999.9;	     
      mean_n_dxy_phi  = -999.9;	     
      width_n_dxy_phi = -999.9;    
      mean_n_dz_phi   = -999.9;    
      width_n_dz_phi  = -999.9;   
      mean_n_dxy_eta  = -999.9;   
      width_n_dxy_eta = -999.9;  
      mean_n_dz_eta   = -999.9;    
      width_n_dz_eta  = -999.9;   

      min_dxy_phi   = -999.9;	     
      max_dxy_phi   = -999.9;    
      min_dz_phi    = -999.9;	     
      max_dz_phi    = -999.9;	     
      min_dxy_eta   = -999.9;	     
      max_dxy_eta   = -999.9;    
      min_dz_eta    = -999.9;	     
      max_dz_eta    = -999.9;	     
      min_n_dxy_phi = -999.9;	     
      max_n_dxy_phi = -999.9;    
      min_n_dz_phi  = -999.9;    
      max_n_dz_phi  = -999.9;   
      min_n_dxy_eta = -999.9;   
      max_n_dxy_eta = -999.9;  
      min_n_dz_eta  = -999.9;    
      max_n_dz_eta  = -999.9;

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
    delete fitParams;
    gROOT->Reset();
    
  } // end of loop over runs
  
  file->cd(); 
  t->Write();
  file->Close();
}
