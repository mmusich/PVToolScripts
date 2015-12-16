#include "TString.h"
#include "TROOT.h"

void runMe(const char* theString,const char* theDate){
  gROOT->ProcessLine(".L FitPVResiduals_forLoop_C.so");
  gROOT->ProcessLine(TString::Format("FitPVResiduals_forLoop(\"%s\",true,false,\"%s\")",theString,theDate));
}
