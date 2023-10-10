#ifndef MYFUNCTION_H_INCLUDED
#define MYFUNCTION_H_INCLUDED

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TString.h"
#include "TClonesArray.h"

#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"

#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraphAsymmErrors.h"

#include "TStyle.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TLatex.h"
#include "TPaveStats.h"


#include "TStopwatch.h"
// #include "TSystem.h"
// #include "TMemStat.h"

void globalSetting()
{
    gStyle->Reset("plain");
    gStyle->SetPalette(1,0);  //1--for colored TH2F
    gStyle->SetFillColor(0); //
    gStyle->SetOptStat(0);   //0--not show histograph para, 1--show histo para
    gStyle->SetOptFit(000);  //000--not show fit para, 111--show fit para.
    gStyle->SetEndErrorSize(0.01);
    TGaxis::SetMaxDigits(4);
    
    gStyle->SetTitleX(0.55); //center title
    gStyle->SetTitleW(0.6);
    gStyle->SetTitleAlign(23);

	gRandom->SetSeed(0);
}

//__________________________________________________________________________________________________________
// two TH1 geometry mean
template <typename HIST_CaclGeoMean>
HIST_CaclGeoMean* caclGeoMean(HIST_CaclGeoMean* h1, HIST_CaclGeoMean* h2) {
    if(!h1 || !h2) return NULL;
    if(h1->GetNbinsX() != h2->GetNbinsX()) return NULL;
    
    HIST_CaclGeoMean* hh = (HIST_CaclGeoMean*)h1->Clone("HIST_CaclGeoMean");
    hh->Sumw2();
    for(int ibin=0; ibin<h1->GetNbinsX(); ibin++) {
        double num1 = h1->GetBinContent(ibin+1);
        double err1 = h1->GetBinError(ibin+1);
        double num2 = h2->GetBinContent(ibin+1);
        double err2 = h2->GetBinError(ibin+1);
        
        double num, err;
        if(num1<0 || num2<0) return NULL;
        else if(num1==0 || num2==0) {
            num = (num1+num2)/2.;
            err = (err1+err2)/2.;
        }
        else {
            num = sqrt(num1*num2);
            err = num*sqrt(pow(err1/num1,2)+pow(err2/num2,2))/2.;
        }
        hh->SetBinContent(ibin+1,num);
        hh->SetBinError(ibin+1,err);
    }
    cout << h1->Integral() << endl;
    cout << h2->Integral() << endl;
    cout << (h1->Integral() + h2->Integral())/2. << "\t" << hh->Integral() << endl;
    return hh;
}


//__________________________________________________________________________________________________________
// TH1 - TF1 integral and error
template <typename HISTtoTF1_Integral>
void hIntegralAndError(HISTtoTF1_Integral* hh, const float x_lw, const float x_up, float& sum, float& err) {
    int bin_lw = hh->FindBin(x_lw);
    int bin_up = hh->FindBin(x_up);
    const float binW_lw = hh->GetBinWidth(bin_lw);
    const float binW_up = hh->GetBinWidth(bin_up);
    const float edge_lw = hh->GetBinLowEdge(bin_lw);
    const float edge_up = hh->GetBinLowEdge(bin_up);
    const float ratio_first = (edge_lw+binW_lw-x_lw)/binW_lw;
    const float ratio_last = (x_up-edge_up)/binW_up;
    sum = 0;
    err = 0;
    for(int i=bin_lw; i<=bin_up; i++) {
        float ratio = 1.0;
        if(i==bin_lw) ratio = ratio_first;
        if(i==bin_up) ratio = ratio_last;
        sum += ratio*hh->GetBinContent(i);
        err += pow(ratio*hh->GetBinError(i),2);
    }
    err = sqrt(err);
    
    /*cout << endl << endl;
     cout << x_lw << "\t" << x_up << endl;
     cout << bin_lw << "\t" << bin_up << endl;
     cout << hh->FindBin(x_lw) << "\t" << hh->FindBin(x_up) << endl;
     cout << sum << endl << endl;*/
}
//__________________________________________________________________________________________________________________
// reset each bin by poisson sample
template <typename HISTPOISSON>    //TH1, TH1F, TH1D, TH1I
void poissonHist(HISTPOISSON* hh) {
    int Nbin = hh->GetNbinsX();
    for(int ii=1; ii<=Nbin; ii++) {
        float mean = hh->GetBinContent(ii);
        float new_mean = gRandom->Poisson(mean);
        hh->SetBinContent(ii, new_mean);
    }
}

//__________________________________________________________________________________________________________________
// reset each bin error == sqrt(bin content)
template <typename HISTBinError>    //TH1, TH1F, TH1D, TH1I
void resetHistError(HISTBinError* hh) {
    int Nbin = hh->GetNbinsX();
    for(int i=1; i<=Nbin; i++) {
        float mean = hh->GetBinContent(i);
        hh->SetBinError(i, sqrt(mean));
    }
}

//__________________________________________________________________________________________________________________
// divide each bin by bin width
template <typename HIST>    //TH1, TH1F, TH1D, TH1I
void densityHist(HIST* h) {
    if(!h) return;
    int nbinsX = h->GetNbinsX();
    if(nbinsX<=0) return;
    for(int ibinX=1; ibinX<=nbinsX; ibinX++) {
        float num = h->GetBinContent(ibinX);
        float err = h->GetBinError(ibinX);
        float binWidth = h->GetBinWidth(ibinX);
        //cout << "before: " << "num = " << num << "\t" << "err = " << err << endl;
        num /= binWidth;
        err /= binWidth;
        //cout << "after: " << "num = " << num << "\t" << "err = " << err << endl;
        h->SetBinContent(ibinX,num);
        h->SetBinError(ibinX,err);
    }
}
template <typename HIST2D>    //TH2, TH2F, TH2D, TH2I
void densityHist2D(HIST2D* h) {
    if(!h) return;
    int nbinsX = h->GetNbinsX();
    int nbinsY = h->GetNbinsY();
    if(nbinsX<=0) return;
    for(int ibinX=1; ibinX<=nbinsX; ibinX++) {
        for(int ibinY=1; ibinY<=nbinsY; ibinY++) {
            float num = h->GetBinContent(ibinX,ibinY);
            float err = h->GetBinError(ibinX,ibinY);
            float binWidthX = h->GetXaxis()->GetBinWidth(ibinX);
            float binWidthY = h->GetYaxis()->GetBinWidth(ibinY);
            //cout << "before: " << "num = " << num << "\t" << "err = " << err << endl;
            num = num / binWidthX / binWidthY;
            err = err / binWidthX / binWidthY;
            //cout << "after: " << "num = " << num << "\t" << "err = " << err << endl;
            h->SetBinContent(ibinX,ibinY,num);
            h->SetBinError(ibinX,ibinY,err);
        }
    }
}

//__________________________________________________________________________________________________________________
// Scale TGraph, TGraphErrors, TGraphAsymmErrors
TGraph* ScaleGraph(TGraph* gr, float scale) {
    for(int i=0; i<gr->GetN(); i++) {
        gr->GetY()[i] *= scale;
    }
    return gr;
}
TGraphErrors* ScaleGraph(TGraphErrors* gr, float scale) {
    for(int i=0; i<gr->GetN(); i++) {
        gr->GetY()[i] *= scale;
        gr->GetEY()[i] *= scale;
    }
    return gr;
}
TGraphAsymmErrors* ScaleGraph(TGraphAsymmErrors* gr, float scale) {
    for(int i=0; i<gr->GetN(); i++) {
        gr->GetY()[i] *= scale;
        gr->GetEYlow()[i] *= scale;
        gr->GetEYhigh()[i] *= scale;
    }
    return gr;
}

//__________________________________________________________________________________________________________________
// two correlated histogram divide
template <typename HISTDivide>
HISTDivide* histDivide(HISTDivide* h1, HISTDivide* h2, bool isWg = false) {
    HISTDivide* heff = new HISTDivide(*h1);
    if(!h1 || !h2 || h1->GetNbinsX() != h2->GetNbinsX()) {
        cout << "error in hist input for histDivide function!!!" << endl;
        return NULL;
    }
    for(int i=0; i<h1->GetNbinsX(); i++) {
        float m = h1->GetBinContent(i+1);
        float merr = h1->GetBinError(i+1);
        float n = h2->GetBinContent(i+1);
        float nerr = h2->GetBinError(i+1);
        
        float eff = (m+1.)/(n+2.);  // not weighted h1 & h2
        float err = sqrt(eff*((m+2.)/(n+3.)-eff));
        
        if(isWg) {
            float meff = m*m/merr/merr;  // calculate the effective entries
            float neff = n*n/nerr/nerr;
            
            eff = n>0 ? m/n : 0;
            err = eff>0 ? sqrt(eff*(1-eff)/neff) : 0;
        }
        
        heff->SetBinContent(i+1,eff);
        heff->SetBinError(i+1,err);
    }
    return heff;
}

//__________________________________________________________________________________________________________________
template <typename HIST3>
HIST3* setHisto(HIST3* h, const char* titlename, const char* xtitle, const char* ytitle,
               int optSize=1,
               double titleSize=0.05, double labelSize=0.05,
               int optOffset=1,
               double x_Offset=1.0, double y_Offset=1.73,
               int optRange=0,
               double x_low=0, double x_up=1, double y_low=0, double y_up=1)
{
    if(h==0) h = new HIST3();
    h->SetTitle(titlename);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitleFont(22);
    h->GetYaxis()->SetTitleFont(22);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->GetZaxis()->SetTitleFont(22);
    h->GetZaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelOffset(0.021);

    h->SetLineWidth(2);

    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
    h->GetZaxis()->SetNdivisions(505);

    if(optSize) {
        h->GetXaxis()->SetTitleSize(titleSize);
        h->GetXaxis()->SetLabelSize(labelSize);  //axis number size
        h->GetYaxis()->SetTitleSize(titleSize);
        h->GetYaxis()->SetLabelSize(labelSize);
        //h->GetZaxis()->SetTitleSize(titleSize);
        //h->GetZaxis()->SetLabelSize(labelSize);
        }

    if(optOffset) {
        h->GetXaxis()->SetTitleOffset(x_Offset);
        h->GetYaxis()->SetTitleOffset(y_Offset);}

    if(optRange) {
        h->GetXaxis()->SetRangeUser(x_low, x_up);
        h->GetYaxis()->SetRangeUser(y_low, y_up);}
    return h;
}
//__________________________________________________________________________________________________________________
template <typename HIST1>
TH1F* getHisto(char name[], Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, char xTitle[], char yTitle[])
{
    TH1F* dd = new TH1F(name,"",100,xlow,xup);
    dd->SetMinimum(ylow);
    dd->SetMaximum(yup);
    dd->GetXaxis()->SetTitle(xTitle);
    dd->GetYaxis()->SetTitle(yTitle);
    dd->GetXaxis()->SetTitleFont(22);
    dd->GetYaxis()->SetTitleFont(22);
    dd->GetXaxis()->SetLabelFont(62);
    dd->GetYaxis()->SetLabelFont(62);

    dd->GetXaxis()->SetTitleSize(0.05);
    dd->GetXaxis()->SetTitleOffset(0.9);
    dd->GetXaxis()->SetLabelSize(0.05);
    dd->GetYaxis()->SetTitleSize(0.05);
    dd->GetYaxis()->SetTitleOffset(1);
    dd->GetYaxis()->SetLabelSize(0.05);
    //dd->GetXaxis()->CenterTitle(kTRUE);
    //dd->GetYaxis()->CenterTitle(kTRUE);
    dd->GetXaxis()->SetNdivisions(512);
    return dd;
}
//__________________________________________________________________________________________________________________
void setPad(TPad* pad, float left=0.17, float right=0.05, float top=0.05, float bottom=0.12, int nlogy=0, int nlogz=1)
{
    pad->SetLogy(nlogy);
    pad->SetLogz(nlogz);
    pad->SetFillColor(10);
    pad->SetBorderMode(0);
    pad->SetBorderSize(0);
    pad->SetFrameFillColor(10);
    pad->SetFrameBorderMode(0);
    pad->SetFrameBorderSize(0);

    #if 1
    pad->SetGridx(0);
    pad->SetGridy(0);
    pad->SetTickx(1);
    pad->SetTicky(1);
    #endif

    #if 1
    pad->SetLeftMargin(left);
    pad->SetRightMargin(right);
    pad->SetTopMargin(top);
    pad->SetBottomMargin(bottom);
    #endif
}

void setPad(float left=0.17, float right=0.05, float top=0.05, float bottom=0.12, int nlogy=0, int nlogz=1)
{
    gPad->SetLogy(nlogy);
    gPad->SetLogz(nlogz);
    gPad->SetFillColor(10);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFrameFillColor(10);
    gPad->SetFrameBorderMode(0);
    gPad->SetFrameBorderSize(0);
    
#if 1
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
#endif
    
#if 1
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetTopMargin(top);
    gPad->SetBottomMargin(bottom);
#endif
}
//__________________________________________________________________________________________________________________

TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex, bool useHistAxis=false)
{
    TLatex *latex = new TLatex(x,y,text);
    if(!useHistAxis) latex->SetNDC();
    latex->SetTextFont(textFont);
    latex->SetTextSize(textSize);
    latex->SetTextColor(colorIndex);
    latex->Draw("same");
    return latex;
}
//__________________________________________________________________________________________________________________
template <typename HIST2>
TPaveStats* setPaveStats(HIST2* h, float x_low = 0.15, float y_low = 0.74, float x_up = 0.55, float y_up = 0.94, int n = 0, char name[][250] = NULL)
{
    TPaveStats* ptstats = new TPaveStats(x_low,y_low,x_up,y_up,"brNDC");
    ptstats->SetName("stats");
    ptstats->SetBorderSize(2);
    //ptstats->SetFillColor(10);
    ptstats->SetTextAlign(12);
    ptstats->SetTextSize(0.04);
    ptstats->SetTextFont(22);
    ptstats->SetShadowColor(0);
    //ptstats->SetLineColor(0);
    ptstats->SetFitFormat("5.3g");
    ptstats->SetStatFormat("6.3g");

    if(name)
    {
        gStyle->SetOptStat(0);   //0--not show histograph para, 1--show histo para
        gStyle->SetOptFit(000);  //000--not show fit para, 111--show fit para.
    }
    else
    {
        gStyle->SetOptStat(0);   //0--not show histograph para, 1--show histo para
        gStyle->SetOptFit(111);  //000--not show fit para, 111--show fit para.
    }

    TText* text;
    for(int i=0; i<n; i++) text = ptstats->AddText(name[i]);

    ptstats->SetOptStat(0);
    ptstats->SetOptFit(000);
    h->GetListOfFunctions()->Add(ptstats);
    ptstats->SetParent(h);

    /*#if 0
    TList *list = ps->GetListOfLines();
    // Remove the RMS line
    TText *tconst = ps->GetLineWith("RMS");
    list->Remove(tconst);
    // Add a new line in the stat box.
    // Note that "=" is a control character
    TLatex *myt = new TLatex(0,0,"Test = 10");
    myt ->SetTextFont(42);
    myt ->SetTextSize(0.04);
    myt ->SetTextColor(kRed);
    list->Add(myt);
    #endif 0*/

    ptstats->Draw();

    //the following line is needed to avoid that the automatic redrawing of stats
    //h->SetStats(0);
    return ptstats;
}

//__________________________________________________________________________________________________________________
// e purity------------one two five Gauss-------------------//

Double_t oneGaus(Double_t* x, Double_t* par)
{
    return par[0]*TMath::Gaus(x[0], par[1], par[2], kTRUE)*par[3];
}

Double_t twoGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<2; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[6];
    return sum;
}

Double_t threeGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<3; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[9];
    return sum;
}

Double_t fourGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<4; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[12];
    return sum;
}

Double_t fiveGaus(Double_t* x, Double_t* par)
{
    Double_t sum = 0;
    for(Int_t i=0; i<5; i++) sum += par[3*i+0]*TMath::Gaus(x[0], par[3*i+1], par[3*i+2], kTRUE)*par[15];
    return sum;
}

Double_t fInvBeta(Double_t* x, Double_t* par)  //par[0]---m, par[1]---per
{
    Double_t p = x[0];
    Double_t mass = par[0];
    Double_t percent = par[1];
    return sqrt(pow(mass,2)/pow(p,2)+1) * (1+percent);
}
// ------------possion function-------------------//
//Double_t myPossion(Double_t* x, Double_t* par) {
    
//}

/*Double_t fMass2(Double_t* x, Double* par)
 {
 Double_t p = x[0];
 Double_t mass = par[0];
 Double_t percent = par[1];
 Double_t mass_c = mass * (1+percent);
 return pow(p,2) * sqrt(pow)
 }*/

Double_t oneExp(Double_t* x, Double_t* par)
{
    return par[0]*pow(10,par[1]*x[0])+par[2];
}

//__________________________________________________________________________________________________________________
double fk_lw(double* x, double* par)
{
    if(x[0]<1.372653) return par[0] + par[1]*x[0] +par[2]*x[0]*x[0];
    else return par[0] + par[1]*1.372653 +par[2]*1.372653*1.372653;
}
double fk_up(double* x, double* par)
{
    if(x[0]<1.899554) return par[0] + par[1]*x[0] +par[2]*x[0]*x[0];
    else return par[0] + par[1]*1.899554 +par[2]*1.899554*1.899554;
}
double fpi_lw(double* x, double* par)
{
    return -2.0;
}
double fpi_up(double* x, double* par)
{
    if(x[0]<1.600093) return par[0] + par[1]*x[0];
    else return par[0] + par[1]*1.600093;
}

//__________________________________________________________________________________________________________________
//double HtoF(double x, TH1F* h)
//{ 
//  int bin = (int)h->FindBin(x);
//  return h->GetBinContent(bin);
//} 
  

//__________________________________________________________________________________________________________________
double round(double val)
{
    return (val> 0.0) ? floor(val+ 0.5) : ceil(val- 0.5);
}
//__________________________________________________________________________________________________________________
//function used in centrality define
double gRefTail(double* x, double* par)
{
    double A = par[0];
    double h = par[1];
    double sigma = par[2];
    return A*TMath::Erf(-sigma*(x[0]-h))+A;
}

double fVzDepand(double* x, double* par)  //polN
{
    double sum = par[0];
    double pol = 1;
    int n=1;
    for(int i=0; i<n; i++) {
        pol *= x[0];
        sum += par[i+1]*pol;
    }
    return sum;
}
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
#endif // MYFUNCTION_H_INCLUDED
