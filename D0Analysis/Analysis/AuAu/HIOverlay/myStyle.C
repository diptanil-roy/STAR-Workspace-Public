#ifndef __myStyle_C__
#define __myStyle_C__

#include <string>
#include <iostream>
#include <sstream>

#include "TPaveText.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TText.h"

using namespace std;

void myStyle() {

  // define names for colours
  //  Int_t black  = 1;
  //Int_t red    = 2;
  //Int_t green  = 3;
  //Int_t blue   = 4;
  //Int_t yellow = 5; 
  //Int_t magenta= 6;
  //Int_t cyan   = 7;
  //Int_t purple = 9;
  

////////////////////////////////////////////////////////////////////


  // Use times new roman, precision 2 
  Int_t myFont        = 132;
  // Line thickness
  Double_t myWidth    = 1.0; 
  // Text size
  Double_t myTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *myStyle= new TStyle("myStyle","My plots style");
    
  myStyle->SetFillColor(1);
  myStyle->SetFillStyle(1001);   // solid
  myStyle->SetFrameFillColor(0);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetStatColor(0);
  myStyle->SetLegendBorderSize(0);

  // If you want the usual gradient palette (blue -> red)
  //myStyle->SetPalette(56);
  // gROOT->ProcessLine(".x ../stpal.C");
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
//  myStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  myStyle->SetPaperSize(20,26);
  myStyle->SetPadTopMargin(0.05);
  myStyle->SetPadRightMargin(0.05); // increase for colz plots
  myStyle->SetPadBottomMargin(0.16);
  myStyle->SetPadLeftMargin(0.14);
  
  // use large fonts
  myStyle->SetTextFont(myFont);
  myStyle->SetTextSize(myTSize);
  myStyle->SetLabelFont(myFont,"x");
  myStyle->SetLabelFont(myFont,"y");
  myStyle->SetLabelFont(myFont,"z");
  myStyle->SetLabelSize(myTSize,"x");
  myStyle->SetLabelSize(myTSize,"y");
  myStyle->SetLabelSize(myTSize,"z");
  myStyle->SetTitleFont(myFont);
  myStyle->SetTitleFont(myFont,"x");
  myStyle->SetTitleFont(myFont,"y");
  myStyle->SetTitleFont(myFont,"z");
  myStyle->SetTitleSize(1.2*myTSize,"x");
  myStyle->SetTitleSize(1.2*myTSize,"y");
  myStyle->SetTitleSize(1.2*myTSize,"z");

  // use medium bold lines and thick markers
  myStyle->SetLineWidth(2);
  myStyle->SetFrameLineWidth(myWidth);
  myStyle->SetHistLineWidth(myWidth);
  myStyle->SetFuncWidth(myWidth);
  myStyle->SetGridWidth(myWidth);
  myStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  myStyle->SetMarkerStyle(20);
  myStyle->SetMarkerSize(1.);

  // label offsets
  myStyle->SetLabelOffset(0.010,"X");
  myStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  myStyle->SetOptStat(0);  
  //myStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  myStyle->SetStatFormat("6.3g"); // specified as c printf options
  myStyle->SetOptTitle(0);
  myStyle->SetOptFit(0);
  //myStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  myStyle->SetTitleOffset(0.95,"X");
  myStyle->SetTitleOffset(0.95,"Y");
  myStyle->SetTitleOffset(1.2,"Z");
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleStyle(0);
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleFont(myFont,"title");
  myStyle->SetTitleX(0.0);
  myStyle->SetTitleY(1.0); 
  myStyle->SetTitleW(1.0);
  myStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  myStyle->SetStatBorderSize(0);
  myStyle->SetStatFont(myFont);
  myStyle->SetStatFontSize(0.05);
  myStyle->SetStatX(0.9);
  myStyle->SetStatY(0.9);
  myStyle->SetStatW(0.25);
  myStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  myStyle->SetNdivisions(505,"x");
  myStyle->SetNdivisions(510,"y");
  
  gROOT->SetStyle("myStyle");
  gROOT->ForceStyle();

  // add My label
  TPaveText *myName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                           0.87 - gStyle->GetPadTopMargin(),
                           gStyle->GetPadLeftMargin() + 0.20,
                           0.95 - gStyle->GetPadTopMargin(),
                           "BRNDC");
  myName->AddText("STAR");
  myName->SetFillColor(0);
  myName->SetTextAlign(12);
  myName->SetBorderSize(0);

  TText *myLabel = new TText();
  myLabel->SetTextFont(myFont);
  myLabel->SetTextColor(1);
  myLabel->SetTextSize(myTSize);
  myLabel->SetTextAlign(12);

  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(myFont);
  myLatex->SetTextColor(1);
  myLatex->SetTextSize(myTSize);
  myLatex->SetTextAlign(12);

  cout << "-------------------------" << endl;  
  cout << "Set My Style - July 2019" << endl;
  cout << "-------------------------" << endl;  
  
}


TPaveText* printMyMy(string optLR="L", string optPrelim="Final", string optText="") {
//////////////////////////////////////////////////////////////////////////
// routine to print 'My', 'My Preliminary' on plots 
// options: optLR=L (top left) / R (top right) of plots
//          optPrelim= Final (My), Prelim (My Preliminary), Other
//          optText= text printed if 'Other' specified
////////////////////////////////////////////////////////////////////
  TPaveText* myName = 0;
  const double rightmargin = gStyle->GetPadRightMargin();
  gStyle->SetPadTopMargin(0.05);
  const double topmargin = gStyle->GetPadTopMargin();
  
  if (optLR==string("R") && optPrelim=="Prelim"){    
    myName = new TPaveText(0.70 - rightmargin,
			     0.80 - topmargin,
			     0.95 - rightmargin,
			     0.90 - topmargin,
			     "BRNDC");
  }
  else if (optLR==string("R")){    
    myName = new TPaveText(0.70 - rightmargin,
			     0.85 - topmargin,
			     0.95 - rightmargin,
			     0.95 - topmargin,
			     "BRNDC");
  }
  else if (optLR==string("L")){
    myName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
			     0.87 - topmargin,
			     gStyle->GetPadLeftMargin() + 0.30,
			     0.95 - gStyle->GetPadTopMargin(),
			     "BRNDC");
  }
  else{
   cout << "printMy: option unknown" << optLR << endl;   
  }
  if (optPrelim==string("Final")){
    myName->AddText("#splitline{My}{#scale[1.0]{#sqrt{s} = 7 TeV data}}");
  }
  else if (optPrelim==string("Prelim")){
    myName->AddText("#splitline{My}{#scale[1.0]{Preliminary}}");  
  }
  else if (optPrelim==string("Other")){
    myName->AddText(optText.c_str());
  }
  else{
    cout << "printMy: option unknown " << optPrelim << endl;   
  }

  myName->SetFillColor(0);
  myName->SetTextAlign(12);
  myName->SetBorderSize(0);
  return myName;
}
#endif
