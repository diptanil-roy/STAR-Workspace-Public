/////////////////////////////////////////////////////////////////////////
// SPlot for HF Jet
// author: Matthew Kelsey
// date Mar. 1 2021 
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"

// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;

// see below for implementation

void AddModel(RooWorkspace*);
void AddData(RooWorkspace*, TTree*, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh);
void DoSPlot(RooWorkspace*, TTree*, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh);
void MakePlots(RooWorkspace*, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh, TString PlotDir = "");
void Write(RooWorkspace*, TTree*, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh, TString OutDir = "");
void makeweights(TString PlotDir = "", TString OutDir = "", double _pt_low = 1., double _pt_high=10.,  int _cen_low=0, int _cen_high=10, double _jetpt_low=-1000, double _jetpt_high=1000., TString mode = "CS", bool widebin = false){

    if (mode == "CS"){
        _jetpt_low=_pt_low;
    }
    if (mode == "AreaBased" && widebin){
        _jetpt_low=-1000;
    }

    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    //Get Signal reduced tree
    TFile *f_D;
    if (mode == "AreaBased")f_D = new TFile("./out_final_AreaBased_Jun20.root");
    // else if (mode == "CS") f_D = new TFile("./out_final_CS1_May27.root");
    else if (mode == "CS") f_D = new TFile("./out_final_CS1_Jun20.root");
    // TFile *f_D = new TFile("./out_final_CS1_May2.root");
    // TFile *f_D = new TFile("../SPlot_Mar2_2023/out_final.root");
    TTree *t = (TTree*)f_D->Get("t");
    // Create a new workspace to manage the project.
    RooWorkspace* wspace = new RooWorkspace("myWS");
    // add the signal and background models to the workspace.
    // Inside this function you will find a discription of the model.
    AddModel(wspace);//worspace,ptbin,ht trigger (0=MB,1=HT),systematic iteration
    // add data to the workspace
    AddData(wspace,t, _pt_low, _pt_high, _cen_low, _cen_high, _jetpt_low, _jetpt_high);
    // do sPlot.  
    //This wil make a new dataset with sWeights added for every event.
    DoSPlot(wspace,t, _pt_low, _pt_high, _cen_low, _cen_high, _jetpt_low, _jetpt_high);
    // Make some plots showing the discriminating variable
    MakePlots(wspace, _pt_low, _pt_high, _cen_low, _cen_high, _jetpt_low, _jetpt_high, PlotDir);
    //Now write everything to new tree
    Write(wspace,t, _pt_low, _pt_high, _cen_low, _cen_high, _jetpt_low, _jetpt_high, OutDir);
    // cleanup
    delete wspace;
}
//____________________________________
void AddModel(RooWorkspace* ws){
    Double_t lowRange = 1.75, highRange = 2.02;
    // make a RooRealVar for the observables
    RooRealVar mM("mM", "m(K#pi) [GeV/#it{c^{2}}]", lowRange, highRange,"");
    /////////////////////////////////////////////
    std::cout << ">>> Make model" << std::endl;
    // model for signal

    RooRealVar m1("m1","D0 Mass", 1.868, 1.855, 1.90);
    RooRealVar sigma1("sigma1", "sigma1",0.02 ,0.005,0.024);
    RooGaussian myModel("myModel", "myModel", mM, m1, sigma1);
      
    RooRealVar bg1 ("bg1","bg1",-5.08148e-01,-10,10);
    RooRealVar bg2("bg2","bg2",-1.46359e-01,-10,10);
    RooRealVar bg3("bg3","bg3",7.85434e-02,-10,10);
    //RooRealVar bg3("bg3","bg3",7.85434e-03,-1,1);
    RooPolynomial bkgModel("bkgModel","bkgModel",mM,RooArgList(bg1,bg2,bg3));
  
    RooRealVar Yield("Yield","fitted yield",5000 ,100.,10000000) ;
    RooRealVar bkgYield("bkgYield","bkg fitted yield",50000 ,0.,10000000) ;
    // now make the combined model                                                                                                                                                                       
    std::cout << ">>> Make full model" << std::endl;
    RooAddPdf model("model","model",
        RooArgList(myModel,bkgModel),
        RooArgList(Yield,bkgYield));


    std::cout << ">>> Import model" << std::endl;
    ws->import(model);
}

//____________________________________
void AddData(RooWorkspace* ws, TTree* t, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh){
  std::cout << ">>> Make data set and import to workspace" << std::endl;
    RooRealVar* mM = ws->var("mM");
    RooRealVar mPt("mPt", "mPt", -1e10, 1e10);
    RooRealVar mEta("mEta", "mEta", -1e10, 1e10);
    RooRealVar mJEta("mJEta", "mJEta", -1e10, 1e10);
    RooRealVar mJPt("mJPt", "mJPt", -1e10, 1e10);
    RooRealVar mJPtArea("mJPtArea", "mJPtArea", -1e10, 1e10);
    // RooRealVar mHPt("mHPt", "mHPt", -1e10, 1e10);
    RooRealVar mZ("mZ", "mZ", -1e10, 1e10);
    RooRealVar mZArea("mZArea", "mZArea", -1e10, 1e10);
    RooRealVar mR("mR", "mR", -1e10, 1e10);
    RooRealVar mCen("mCen", "mCen", -1e10, 1e10);
    RooRealVar mgRefMultCorr("mgRefMultCorr", "mgRefMultCorr", -1e10, 1e10);
    RooRealVar mWeight("mWeight", "mWeight", -1e10, 1e10);
    RooRealVar mJetArea("mJetArea", "mJetArea", -1e10, 1e10);
    RooRealVar mJetNConstArea("mConstArea", "mConstArea", -1e10, 1e10);
    RooRealVar mJetNConstCS("mConstCS", "mConstCS", -1e10, 1e10);
    RooArgSet* variables = new RooArgSet();
    variables->add(*mM);
    variables->add(mPt);
    variables->add(mEta);
    variables->add(mJEta);
    variables->add(mJPt);
    variables->add(mJPtArea);
    // variables->add(mHPt);
    variables->add(mZ);
    variables->add(mZArea);
    variables->add(mR);
    variables->add(mCen);
    variables->add(mgRefMultCorr);
    variables->add(mWeight);
    variables->add(mJetArea);
    variables->add(mJetNConstArea);
    variables->add(mJetNConstCS);

    char cuts[200];
    //sprintf(cuts,"mJPt>%1.2f",5.0);
    // sprintf(cuts,"mEta>%1.2f && mEta<%1.2f &&  mJPt>%1.2f && mJPt<%1.2f && mPt>%1.2f && mPt<%1.2f && mCen>=%i && mCen<%i",-10.,10.,_jetpt_low,_jetpt_high,_pt_low,_pt_high,_cen_low,_cen_high);
    sprintf(cuts,"mEta>%1.2f && mEta<%1.2f &&  mJPt>%1.2f && mJPt<%1.2f && mPt>%1.2f && mPt<%1.2f && mCen>=%i && mCen<%i",-10.,10.,jetptlow,jetpthigh,ptlow,pthigh,cenlow,cenhigh);
    cout << "Cuts used " << cuts << endl;
    RooDataSet* data = new RooDataSet("data","data",*variables,Import(*t),Cut(cuts));
    // import data into workspace
    ws->import(*data, Rename("data"));
}

//____________________________________
void DoSPlot(RooWorkspace* ws,TTree *tree, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh){
  std::cout << ">>> Calculate sWeights" << std::endl;
  
    // get what we need out of the workspace to do the fit                                                                                                                                            
    RooAbsPdf* model = ws->pdf("model");
    RooRealVar* Yield = ws->var("Yield");
    RooRealVar* bkgYield = ws->var("bkgYield");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    // fit the model to the data.                                                                                                                                                                    
    model->fitTo(*data, Extended() );
    // The sPlot technique requires that we fix the parameters                                                                                                                                    
    // of the model that are not yields after doing the fit.                                                                                                                                         
    RooRealVar* sigma1 = ws->var("sigma1");
    RooRealVar* m1 = ws->var("m1");
    RooRealVar* bg1 = ws->var("bg1");
    RooRealVar* bg2 = ws->var("bg2");
    RooRealVar* bg3 = ws->var("bg3");
    sigma1->setConstant();
    m1->setConstant();
    bg1->setConstant();
    bg2->setConstant();
    bg3->setConstant();
    RooMsgService::instance().setSilentMode(true);


    RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
             *data, model, RooArgList(*Yield,*bkgYield) );
    // Check that our weights have the desired properties                                                                                                                                             
    std::cout << ">>> Check SWeights:" << std::endl;


    std::cout << std::endl <<  "Yield of D0 is "
        << Yield->getVal() << ".  From sWeights it is "
        << sData->GetYieldFromSWeight("Yield") << std::endl;


    std::cout << "Yield of BKG is "
        << bkgYield->getVal() << ".  From sWeights it is "
        << sData->GetYieldFromSWeight("bkgYield") << std::endl
        << std::endl;
    std::cout << "import new dataset with sWeights" << std::endl;
    ws->import(*data, Rename("dataWithSWeights"));
 
}

void MakePlots(RooWorkspace* ws, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh, TString PlotDir = "plots"){

    float lmargin = 0.15;
    float rmargin = 0.05;
    float bmargin = 0.14;
    float tmargin = 0.07;

    gStyle->SetPaperSize(20,26);
    gStyle->SetPadTopMargin(tmargin);
    gStyle->SetPadRightMargin(rmargin); // increase for colz plots
    gStyle->SetPadBottomMargin(bmargin);
    gStyle->SetPadLeftMargin(lmargin);

    gStyle->SetLabelOffset(0.010,"X");
    gStyle->SetLabelOffset(0.010,"Y");

    gStyle->SetTitleOffset(1.1,"X");
    gStyle->SetLabelOffset(1.1,"Y");

    // put tick marks on top and RHS of plots
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // histogram divisions: only 5 in x to avoid label overlaps
    gStyle->SetNdivisions(505,"x");
    gStyle->SetNdivisions(505,"y");
    gStyle->SetOptStat(0);  
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);

    // gROOT->ProcessLine(".x ~/myStyle.C"); 
    std::cout << ">>> Make plots" << std::endl;
    RooAbsPdf* model = ws->pdf("model");
    RooAbsPdf* myModel = ws->pdf("myModel");
    RooAbsPdf* bkgModel = ws->pdf("bkgModel");
    RooRealVar* mM = ws->var("mM");
    RooRealVar* Yield = ws->var("Yield");
    RooRealVar* bkgYield = ws->var("bkgYield");
    // note, we get the dataset with sWeights                                                                                                                                                          
    RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
    //RooDataSet* data_WS = (RooDataSet*) ws->data("data_WS");                                                                                                                                         
    model->fitTo(*data, Extended() );

    TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 800);
    //cdata->Divide(1,3);                                                                                                                                                                              
    cdata->cd(1);
    RooPlot* frame = mM->frame(Bins(100)) ;
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("m(K#pi) [GeV/#it{c^{2}}]");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(0.9);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitle("Counts(/27 [MeV/#it{c^{2}}])");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleOffset(0.9);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelOffset(0.012);
    frame->GetYaxis()->SetRangeUser(0,9900);
    frame->GetYaxis()->SetMaxDigits(2);

    data->plotOn(frame, Name("Data") ) ;
    //data_WS->plotOn(frame,Bins(50),Name("Histws"),LineStyle(1),LineColor(15),FillColor(17),FillStyle(3007));                                                                                         
    model->plotOn(frame, Name("Model")) ;
    model->plotOn(frame,Components(*myModel),Name("Signal"),LineStyle(kDashed),LineColor(kGreen-2),Normalization(1.0,RooAbsReal::RelativeExpected)) ;
    model->plotOn(frame,Components(*bkgModel),Name("Background"),LineStyle(kDashed),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected)) ;
    // frame->SetTitle("Fit of model to discriminating variable");
    frame->Draw() ;
    char save[100];

    TH1D *DataDummy = new TH1D("DataDummy","",1,1,2);
    DataDummy->SetMarkerStyle(20);
    DataDummy->SetMarkerSize(1.5);
    DataDummy->SetMarkerColor(kBlack);
    DataDummy->SetLineColor(kBlack);
    // DataDummy->SetMinimum(0);
    // DataDummy->SetMaximum(0.05*data->sumEntries());
    // DataDummy->GetXaxis()->SetTitle("m(K#pi) [GeV/#it{c^{2}}]");
    // DataDummy->GetXaxis()->SetTitleSize(0.05);
    // DataDummy->GetXaxis()->SetTitleOffset(0.9);
    // DataDummy->GetXaxis()->SetLabelSize(0.05);
    // DataDummy->GetYaxis()->SetTitle("Counts(/27 [MeV/#it{c^{2}}])");
    // DataDummy->GetYaxis()->SetTitleSize(0.05);
    // DataDummy->GetYaxis()->SetTitleOffset(0.9);
    // DataDummy->GetYaxis()->SetLabelSize(0.05);
    // DataDummy->GetYaxis()->SetLabelOffset(0.12);

    
    TH1D *ModelDummy = new TH1D("ModelDummy","",1,1,2);
    ModelDummy->SetLineColor(kBlue);
    ModelDummy->SetLineWidth(3);

    TH1D *SignalDummy = new TH1D("SignalDummy","",1,1,2);
    SignalDummy->SetLineColor(kGreen-2);
    SignalDummy->SetLineStyle(kDashed);
    SignalDummy->SetLineWidth(3);

    TH1D *BackgroundDummy = new TH1D("BackgroundDummy","",1,1,2);
    BackgroundDummy->SetLineColor(kRed);
    BackgroundDummy->SetLineStyle(kDashed);
    BackgroundDummy->SetLineWidth(3);


    TLegend *Legend;
    Legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    Legend->SetBorderSize(0);
    Legend->SetFillStyle(0);
    Legend->SetTextFont(43);
    Legend->SetTextSize(24);
    Legend->AddEntry(DataDummy, "Data", "lp");
    Legend->AddEntry(ModelDummy, "Gaus + pol3", "l");
    Legend->AddEntry(SignalDummy, "Signal", "l");
    Legend->AddEntry(BackgroundDummy, "Background", "l");
    Legend->Draw("SAME");
    // gPad->BuildLegend();


    // TPaveText with Labels
    TPaveText *Label;
    Label = new TPaveText(0.2, 0.7, 0.42, 0.9, "NDC");
    Label->SetBorderSize(0);
    Label->SetFillStyle(0);
    Label->SetTextFont(43);
    Label->SetTextSize(24);
    TText *t1 = Label->AddText(Form("#bf{#it{STAR}}"));
    Label->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
    Label->AddText(Form("0 - 80 %%"));
    Label->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", 1, 10));
    Label->Draw("SAME");

    // TPaveText with fit parameters

    TPaveText* pt = new TPaveText(0.58,0.6,0.88,0.7,"brNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.03);
    pt->SetTextColor(kBlack);
    pt->AddText(Form("Yield = %1.2f #pm %1.2f",Yield->getVal(),Yield->getError()));
    // pt->AddText(Form("Bkg = %1.2f #pm %1.2f",bkgYield->getVal(),bkgYield->getError()));
    pt->Draw("SAME");

    sprintf(save,"%s/D0_Mass_sPlot_PT_%1.f_%1.f_Cen_%i_%i.pdf",PlotDir.Data(),ptlow,pthigh,cenlow,cenhigh);
    cout << save << endl;
    cdata->SaveAs(save);


    double bkg_y = bkgYield->getVal();
    

    double bkg_sb = 0;
    mM->setRange("range",1.75,2.02);
    RooAbsReal* bkg = bkgModel->createIntegral(*mM,*mM,"range");

    mM->setRange("range1",1.75,1.8);
    RooAbsReal* bkg1 = bkgModel->createIntegral(*mM,*mM,"range1");
    mM->setRange("range2",1.8,1.92);
    RooAbsReal* bkg2 = bkgModel->createIntegral(*mM,*mM,"range2");
    mM->setRange("range3",1.92,2.02);
    RooAbsReal* bkg3 = bkgModel->createIntegral(*mM,*mM,"range3");

    double frac1 = bkg1->getVal();///bkg->getVal();
    double frac2 = bkg2->getVal();///bkg->getVal();
    double frac3 = bkg3->getVal();///bkg->getVal();

    cout << "Now looking at BKG fractions: " << endl;
    cout << "Low SB " << frac1 <<endl;
    cout << "High SB " << frac3 <<endl;
    cout << "Signal Region " << frac2 <<endl;

    cout <<"Ratio for SB subtractions " << endl;
    cout << frac3/(frac1+frac2) << endl;

}

void Write(RooWorkspace* ws, TTree* tree, double ptlow, double pthigh, int cenlow, int cenhigh, double jetptlow, double jetpthigh, TString OutputDirName = ""){
    char outfile[100];
    sprintf(outfile,"%s/NewsWeight_PT_%1.f_%1.f_Cen_%i_%i.root",OutputDirName.Data(),ptlow,pthigh,cenlow,cenhigh);
    TFile out(outfile,"RECREATE");
    RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
    TTree* Signal_sw = new TTree("Signal_sw","Signal_sw");//tree->CloneTree(0);
    float weight;
    float MM;
    float MPt;
    float MJPt;
    float MJPtArea;
    float MHPt;
    float MZ;
    float MZArea;
    float MR;
    int MCEN;
    float MGREFMULT;
    float MWEIGHT;
    float MJETAREA;
    float META;
    float MJETA;
    int MCONSTAREA;
    int MCONSTCS;
    Signal_sw->Branch("sWeight",&weight,"sWeight/F") ;
    Signal_sw->Branch("mM",&MM,"mM/F") ;
    Signal_sw->Branch("mPt",&MPt,"mPt/F") ;
    Signal_sw->Branch("mJPt",&MJPt,"mJPt/F") ;
    Signal_sw->Branch("mJPtArea",&MJPtArea,"mJPtArea/F") ;
    // Signal_sw->Branch("mHPt",&MHPt,"mHPt/F") ;
    Signal_sw->Branch("mZ",&MZ,"mZ/F") ;
    Signal_sw->Branch("mZArea",&MZArea,"mZArea/F") ;
    Signal_sw->Branch("mR",&MR,"mR/F") ;
    Signal_sw->Branch("mCen",&MCEN,"mCen/I") ;
    Signal_sw->Branch("mGRefMult",&MGREFMULT,"mGRefMult/F") ;
    Signal_sw->Branch("mWeight",&MWEIGHT,"mWeight/F") ;
    Signal_sw->Branch("mJetArea",&MJETAREA,"mJetArea/F") ;
    Signal_sw->Branch("mEta",&META,"mEta/F") ;
    Signal_sw->Branch("mJEta",&MJETA,"mJEta/F") ;
    Signal_sw->Branch("mConstArea",&MCONSTAREA,"mConstArea/I") ;
    Signal_sw->Branch("mConstCS",&MCONSTCS,"mConstCS/I") ;
    data->get()->Print();
    RooRealVar* Yield_sw = ws->var("Yield_sw");
    RooRealVar* mM = ws->var("mM");
    RooRealVar* mPt = ws->var("mPt");
    RooRealVar* mJPt = ws->var("mJPt");
    RooRealVar* mJPtArea = ws->var("mJPtArea");
    // RooRealVar* mHPt = ws->var("mHPt");
    RooRealVar* mZ = ws->var("mZ");
    RooRealVar *mZArea = ws->var("mZArea");
    RooRealVar* mR = ws->var("mR");
    RooRealVar* mCen = ws->var("mCen");
    RooRealVar* mgRefMult = ws->var("mgRefMultCorr");
    RooRealVar* mWeight = ws->var("mWeight");
    RooRealVar* mJetArea = ws->var("mJetArea");
    RooRealVar* mJEta = ws->var("mJEta");
    RooRealVar* mEta = ws->var("mEta");
    RooRealVar* mJetNConstArea = ws->var("mConstArea");
    RooRealVar* mJetNConstCS = ws->var("mConstCS");
    RooArgSet* set;
    for (int i=0 ; i< data->sumEntries() ; i++) {
        RooArgSet* row = (RooArgSet*)data->get(i);  
  
        Yield_sw = (RooRealVar*)row->find(Yield_sw->GetName()); 
        weight = Yield_sw->getVal();
  
        mM = (RooRealVar*)row->find(mM->GetName());
        MM = mM->getVal();

        mPt = (RooRealVar*)row->find(mPt->GetName());
        MPt = mPt->getVal();
  
        mJPt = (RooRealVar*)row->find(mJPt->GetName());
        MJPt = mJPt->getVal();

        mJPtArea = (RooRealVar*)row->find(mJPtArea->GetName());
        MJPtArea = mJPtArea->getVal();

        // mHPt = (RooRealVar*)row->find(mHPt->GetName());
        // MHPt = mHPt->getVal();

        mZ = (RooRealVar*)row->find(mZ->GetName());
        MZ = mZ->getVal();

        mZArea = (RooRealVar*)row->find(mZArea->GetName());
        MZArea = mZArea->getVal();

        mR = (RooRealVar*)row->find(mR->GetName());
        MR = mR->getVal();

        mCen = (RooRealVar*)row->find(mCen->GetName());
        MCEN = mCen->getVal();

        mgRefMult = (RooRealVar*)row->find(mgRefMult->GetName());
        MGREFMULT = mgRefMult->getVal();

        mWeight = (RooRealVar*)row->find(mWeight->GetName());
        MWEIGHT = mWeight->getVal();

        mJetArea = (RooRealVar*)row->find(mJetArea->GetName());
        MJETAREA = mJetArea->getVal();

        mJEta = (RooRealVar*)row->find(mJEta->GetName());
        MJETA = mJEta->getVal();
        
        mEta = (RooRealVar*)row->find(mEta->GetName());
        META = mEta->getVal();

        mJetNConstArea = (RooRealVar*)row->find(mJetNConstArea->GetName());
        MCONSTAREA = mJetNConstArea->getVal();

        mJetNConstCS = (RooRealVar*)row->find(mJetNConstCS->GetName());
        MCONSTCS = mJetNConstCS->getVal();

        Signal_sw->Fill();
    }
    
    Signal_sw->Write();
    // cout << Form("FileName = sWeights/NewsWeight_PT_%1.f_%1.f_Cen_%i_%i.root",_pt_low,_pt_high,_cen_low,_cen_high) << endl;
    out.Close();
    out.Delete(); 
}
