#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TClass.h"
#include <iostream>
#include <string>
#include "TPaveText.h"

void compareAll(TString file1,TString file2){

  TFile *f1 = TFile::Open(file1);
  TFile *f2 = TFile::Open(file2);

  f1->cd("PVValidation");

  int i=0;
  int j=0;

  TCanvas dummyC;
  dummyC.Print("diff.pdf[");

  TIter nextkey(gDirectory->GetListOfKeys());
  while (TKey* key = (TKey*)nextkey()) {
    //std::cout << "i: " << i << std::endl;
    ++i;
    TObject* obj = key->ReadObj();
    std::string name = obj->GetName();
    std::cout << "name: " << name << std::endl;
    if ( obj->IsA()->InheritsFrom("TDirectory") ){
      f1->cd(("PVValidation/"+name).c_str());
      TIter nextkey(gDirectory->GetListOfKeys());                                                                             
      while (key = (TKey*)nextkey()) { 
	obj = key->ReadObj();                                                                                         
	if (obj->IsA()->InheritsFrom("TH1")) {                                                                                
	  TH1* h = (TH1*)obj;                                                                              
	  TString fullpath = "PVValidation/"+name+"/"+h->GetName();
	  //std::cout << "j: " << j << " "<< h->GetName() <<" "<< fullpath << std::endl;
	  ++j;
	  if (obj->IsA()->InheritsFrom("TH2")) continue;
	  TCanvas *c1 = new TCanvas(h->GetName(),h->GetName(),800,600);
	  c1->cd();
	  h->SetMarkerColor(kRed);
	  h->SetLineColor(kRed);
	  h->SetMarkerStyle(kOpenSquare);
	  h->Draw();
	  TH1 *h2 = (TH1*)f2->Get(fullpath.Data());

	  if(h2==nullptr){
	    std::cout<<"WARNING!!!!!! "<<fullpath<<" does NOT exist in second file!!!!!"<<std::endl;
	    continue;
	  }

	  h2->SetMarkerColor(kBlue);
	  h2->SetLineColor(kBlue);
	  h2->SetMarkerStyle(kOpenCircle);
	  h2->Draw("same");
	  TString savename = fullpath.ReplaceAll("/","_");
	  double ksProb = 0;
	  ksProb = h->KolmogorovTest(h2);
	  if(ksProb!=1.){
	    //c1->SaveAs(savename+".pdf");
	    TPaveText ksPt(0,0, 0.35, 0.04, "NDC"); ksPt.SetBorderSize(0); ksPt.SetFillColor(0);
	    ksPt.AddText(Form("P(KS)=%g, ered %g eblue %g",ksProb, h->GetEntries(), h2->GetEntries()));
	    ksPt.Draw();
	    c1->Print("diff.pdf");
	    std::cout<<"histogram # "<<j<<": "<<fullpath<<" |kolmogorov: "<<ksProb<<std::endl;
	  }
	  
	  delete c1;
	  delete h;
	  delete h2;
	}                                                                                                                
      }
    }
  }//while
  f1->Close();
  f2->Close();
  dummyC.Print("diff.pdf]");
}
