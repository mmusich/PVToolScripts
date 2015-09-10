#include "TRandom.h"
#include "TMath.h"
#include "TGaxis.h"
#include "time.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "TStyle.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TObject.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH1.h"

void checkMissing(TString filename1_,TString filename2_){
  using namespace std;
  
  ifstream infile1,infile2;
  infile1.open(filename1_.Data());
  infile2.open(filename2_.Data());

  //variables to be read

  Int_t date1, nevents1, date2, nevents2;
  Double_t dz1, sigmadz1, dz2, sigmadz2;
  TString flag1, flag2;
  
  const Int_t nlines=1000;

  Double_t issues1[nlines], issues2[nlines];
  
  //counter
  UInt_t theline1(0),theline2(0);

  while (!infile1.eof()) {
    infile1 >> date1  >> dz1 >> sigmadz1 >> nevents1 >> flag1;
    issues1[theline1]=date1;
    theline1++;
  }

  while (!infile2.eof()) {
    infile2 >> date2  >> dz2 >> sigmadz2 >> nevents2 >> flag2;
    issues2[theline2]=date2;
    theline2++;
  }
  
  std::cout<<filename1_<<" has "<<theline1<<" lines.  "<<filename2_<<" has "<<theline2<<" lines."<<std::endl;


  if(theline1>theline2){
    // first this way
    for(UInt_t i=0;i<theline1;i++){
      Bool_t isFound_(false); 
      for(UInt_t j=0;j<theline2;j++){
	if(issues1[i]==issues2[j]){
	  //cout<<"i:"<<i<<" j:"<<j<<" "<<issues1[i]<<"=="<<issues2[j]<<std::endl;
	  isFound_=true;
	}
      }
      if(isFound_==false){
	std::cout<<"Houston we got a problem! "
		 <<"Entry "<<i<<" \""<<issues1[i]<<"\" of file: "<<filename1_<<" is not found anywhere in file: "<<filename2_<<std::endl;
      }
    }
  } else {
    // then the other way around
    for(UInt_t i=0;i<theline2;i++){
      Bool_t isFound_(false); 
      for(UInt_t j=0;j<theline1;j++){
	if(issues2[i]==issues1[j]){
	  //cout<<"i:"<<i<<" j:"<<j<<" "<<issues2[i]<<"=="<<issues1[j]<<std::endl;
	  isFound_=true;
	}
      }
      if(isFound_==false){
	std::cout<<"Houston we got a problem! "
		 <<"Entry "<<i<<" \""<<issues2[i]<<"\" of file: "<<filename2_<<" is not found anywhere in file: "<<filename1_<<std::endl;
      }
    }
  }

}
