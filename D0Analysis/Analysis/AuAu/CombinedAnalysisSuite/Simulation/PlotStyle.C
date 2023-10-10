#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

void PlotStyle(){
    // Use times new roman, precision 2 
    Int_t myFont        = 132;
    // Line thickness
    Double_t myWidth    = 1.0; 
    // Text size
    Double_t myTSize    = 0.05; 

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
    // gROOT->ProcessLine(".x ./stpal.C");
    // If you want colors that correspond to gray scale in black and white:
    int colors[8] = {0,5,7,3,6,2,4,1};
    //  myStyle->SetPalette(8,colors);

    // set the paper & margin sizes
    myStyle->SetPaperSize(20,26);
    myStyle->SetPadTopMargin(0.07);
    myStyle->SetPadRightMargin(0.05); // increase for colz plots
    myStyle->SetPadBottomMargin(0.14);
    myStyle->SetPadLeftMargin(0.22);

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
    myStyle->SetTitleSize(myTSize,"x");
    myStyle->SetTitleSize(myTSize,"y");
    myStyle->SetTitleSize(myTSize,"z");

    // use medium bold lines and thick markers
    // myStyle->SetLineWidth(2);
    // myStyle->SetFrameLineWidth(myWidth);
    // myStyle->SetHistLineWidth(myWidth);
    // myStyle->SetFuncWidth(myWidth);
    // myStyle->SetGridWidth(myWidth);
    // myStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    // myStyle->SetMarkerStyle(20);
    // myStyle->SetMarkerSize(1.);

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
    // myStyle->SetTitleOffset(1.2,"X");
    // myStyle->SetTitleOffset(1.9,"Y");
    // myStyle->SetTitleOffset(1.2,"Z");
    // myStyle->SetTitleFillColor(0);
    // myStyle->SetTitleStyle(0);
    // myStyle->SetTitleBorderSize(0);
    // myStyle->SetTitleFont(myFont,"title");
    // myStyle->SetTitleX(0.0);
    // myStyle->SetTitleY(1.0); 
    // myStyle->SetTitleW(1.0);
    // myStyle->SetTitleH(0.05);

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
    myStyle->SetNdivisions(504,"y");

    gROOT->SetStyle("myStyle");
    // gROOT->ForceStyle();
}