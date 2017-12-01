#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include "TDatime.h"
#include <vector>
#include <string>
#include <map>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

std::vector<std::string> split(const std::string& s,char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)){
    tokens.push_back(token);
  }
  return tokens;
}

void testFileParsing(){

  std::map<int,TDatime> times;
   
  std::ifstream infile("times.txt");
  std::string line;
  while (std::getline(infile, line))
    {
      std::istringstream iss(line);
      std::string a,b,c;
      if (!(iss >> a >> b >> c)) { break; } // error
      
      //std::cout<<a<<"  "<<b<<"   "<<c<<"   "<<std::endl;
      
      int run  = std::stoi(a);
      auto tokens_b = split(b,'-');
      int year  = std::stoi(tokens_b[0]);
      int month = std::stoi(tokens_b[1]);
      int day   = std::stoi(tokens_b[2]);
      
      auto tokens_c  = split(c,'.');
      auto tokens_c1 = split(tokens_c[0],':');

      int hour   = std::stoi(tokens_c1[0]);
      int minute = std::stoi(tokens_c1[2]);
      int second = std::stoi(tokens_c1[2]);

      //std::cout<<run<<" "<<year<<" "<<month<<" "<<day<<" "<<hour<<" "<<minute<<" "<<second<<" "<<std::endl;
	
      TDatime da(year,month,day,hour,minute,second);
      times[run]=da;
    }

  float x[times.size()],y[times.size()];
  
  std::cout<<times.size()<<std::endl;

  int i=0;
  for(const auto &element : times){
    x[i]=(element.second).Convert();
    y[i]=(element.first);
    (element.second).Print();
    std::cout<<x[i]<<" "<<y[i]<<std::endl;
    i++;

  }

  TCanvas *c_out = new TCanvas("cout","cout",1600,800);
  c_out->cd();

  TGraph *mgr = new TGraph(times.size(),x,y);
  mgr->SetMarkerStyle(20);
  mgr->Draw("ap");

  mgr->GetXaxis()->SetTimeDisplay(1);
  mgr->GetXaxis()->SetNdivisions(510);
  mgr->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  mgr->GetXaxis()->SetTimeOffset(0,"gmt");



}




