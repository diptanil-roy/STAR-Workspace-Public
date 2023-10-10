#include "BinDef.h"
#include "NewBinDef.h"

using namespace std;

void MakeFONLLHistogram(double lowptcutoff = 1.0, double highptcutoff = 20.0){

	int nptbins = (highptcutoff - lowptcutoff)*4;

	TString FileName = "FONLL.txt";
	TH1D *FONLL = new TH1D(Form("FONLL"), Form("FONLL"), nptbins, lowptcutoff, highptcutoff);

	double pt, ptweight;
	ifstream myfile (FileName.Data());
	if (myfile.is_open())
	{
		while ( !myfile.eof() )
		{
			myfile >> pt >> ptweight;
			// cout << x << "\t" << y << endl;
			if (pt < lowptcutoff || pt > highptcutoff) continue;
			int bin = FONLL->FindBin(pt+0.0001);
			FONLL->SetBinContent(bin, ptweight);
		}
	}

	FONLL->Scale(1./FONLL->Integral());
	FONLL->Draw();
	gPad->SetLogy();

	TFile *f = new TFile(Form("FONLL_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "RECREATE");
	f->cd();
	FONLL->Write();
	f->Close();
}
