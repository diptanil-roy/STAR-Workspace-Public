R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int SUPERITERATION = 1, int iteration = 4, bool mimicDataUncertainty = kFALSE){
	// for (int i = 0; i < njpt_gen_bins_var+1; i++){
	// 	jetpt_var_bin[i] = 5 + (20. - 5.)/njpt_gen_bins_var * i;
	// }

	TString RespFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins_var, nz_gen_bins, 1);
	TFile *RespFile = new TFile(RespFileName.Data());
	RespFile->cd();

	RooUnfoldResponse *resp[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
	}

	TString FoldedFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins_var, nz_gen_bins, 2); // Step 0 is always where I store the folded distro to be unfolded
	TFile *FoldedFile = new TFile(FoldedFileName.Data());

	THnSparseF *hJet[3];

	for (int i = 0; i < 3; i++){
		hJet[i] = (THnSparseF *)gDirectory->Get(Form("hJetZNormPtNorm_%i", i));
	}

	TH2D *Measured[3];
	TH2D *Truth[3];

	for (int i = 0; i < 3; i++){

		TH2D *tmpMeas = (TH2D *)hJet[i]->Projection(3, 2, "E");
		TH2D *tmpTrue = (TH2D *)hJet[i]->Projection(1, 0, "E");

		cout << i << "\t" << tmpMeas->Integral() << "\t" << tmpTrue->Integral() << endl;

		Truth[i] = new TH2D(Form("Truth_%i", i), Form("Truth_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
		Measured[i] = new TH2D(Form("Measured_%i", i), Form("Measured_%i", i), njpt_bins, nbinsjetpt, nz_bins, nbinsz);

		// Truth Distribution
		for (int binx = 1; binx <= Truth[i]->GetNbinsX(); binx++){
			for (int biny = 1; biny <= Truth[i]->GetNbinsY(); biny++){
				double xlow = Truth[i]->GetXaxis()->GetBinLowEdge(binx);
				double xup = Truth[i]->GetXaxis()->GetBinUpEdge(binx);
				double ylow = Truth[i]->GetYaxis()->GetBinLowEdge(biny);
				double yup = Truth[i]->GetYaxis()->GetBinUpEdge(biny);

				int binxlow = tmpTrue->GetXaxis()->FindBin(xlow + 0.00001);
				int binxup = tmpTrue->GetXaxis()->FindBin(xup - 0.00001);
				int binylow = tmpTrue->GetYaxis()->FindBin(ylow + 0.00001);
				int binyup = tmpTrue->GetYaxis()->FindBin(yup - 0.00001);

				// cout << xlow << "\t" << xup << "\t" << ylow << "\t" << yup << endl;
				// cout << binxlow << "\t" << binxup << "\t" << binylow << "\t" << binyup << endl;

				double error = 0.;
				double content = tmpTrue->IntegralAndError(binxlow, binxup, binylow, binyup, error, "");

				// cout << content << "\t" << error << endl;
				
				Truth[i]->SetBinContent(binx, biny, content);
				Truth[i]->SetBinError(binx, biny, error);
			}
		}

		// Measured Distribution
		for (int binx = 1; binx <= Measured[i]->GetNbinsX(); binx++){
			for (int biny = 0; biny <= Measured[i]->GetNbinsY()+1; biny++){
				double xlow = Measured[i]->GetXaxis()->GetBinLowEdge(binx);
				double xup = Measured[i]->GetXaxis()->GetBinUpEdge(binx);
				double ylow = Measured[i]->GetYaxis()->GetBinLowEdge(biny);
				double yup = Measured[i]->GetYaxis()->GetBinUpEdge(biny);

				int binxlow = tmpMeas->GetXaxis()->FindBin(xlow + 0.00001);
				int binxup = tmpMeas->GetXaxis()->FindBin(xup - 0.00001);
				int binylow = tmpMeas->GetYaxis()->FindBin(ylow + 0.00001);
				int binyup = tmpMeas->GetYaxis()->FindBin(yup - 0.00001);

				double error = 0.;
				double content = tmpMeas->IntegralAndError(binxlow, binxup, binylow, binyup, error);
				
				Measured[i]->SetBinContent(binx, biny, content);
				Measured[i]->SetBinError(binx, biny, error);
			}
		}
	}

	if (mimicDataUncertainty){
		TFile *DataFile = new TFile("../SPlot_Mar2_2023/Histograms3_D01_10GeV_NewBinningOldCuts_Apr6_WiderJpTRange.root", "READ");
		DataFile->cd();

		for (int cent = 0; cent < 3; cent++){
			TH2D *tmp;
			if (cent == 0) tmp = (TH2D *)gDirectory->Get("ZPt_0_10");
			else if (cent == 1) tmp = (TH2D *)gDirectory->Get("ZPt_10_40");
			else if (cent == 2) tmp = (TH2D *)gDirectory->Get("ZPt_40_80");

			TF1* fRand = new TF1("fRand","TMath::Gaus(x,[0],[1],1)",-100,100);
			fRand->SetNpx(10000);
			for (int binx = 0; binx <= Measured[cent]->GetNbinsX()+1; binx++ ){
				for (int biny = 0; biny <= Measured[cent]->GetNbinsY()+1; biny++ ){
					double binerrorpercent = 0.;
					if (tmp->GetBinContent(binx, biny) != 0) binerrorpercent = tmp->GetBinError(binx, biny)/(1.0*tmp->GetBinContent(binx, biny));
					double bincontent = Measured[cent]->GetBinContent(binx, biny);
					double binerror = bincontent*binerrorpercent;

					double newbincontent = bincontent;
					if (binerrorpercent == 0.) newbincontent = 0.;

					// cout << binx << "\t" << biny << "\t" << bincontent << "\t" << binerror << endl;
					if (binerrorpercent != 0.){
						fRand->SetRange(bincontent-6.*binerror, bincontent+6.*binerror);
					    fRand->SetParameter(0,bincontent);
					    fRand->SetParameter(1,binerror);
					    newbincontent = fRand->GetRandom();
						while (newbincontent<=0.) {newbincontent = fRand->GetRandom();}
					}

					Measured[cent]->SetBinContent(binx, biny, newbincontent);
					Measured[cent]->SetBinError(binx, biny, newbincontent*binerrorpercent);
				}
			}	
		}
	}


	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

	RooUnfoldBayes unfoldcent (resp[0], Measured[0], iteration);
	RooUnfoldBayes unfoldmid (resp[1], Measured[1], iteration);
	RooUnfoldBayes unfoldperi (resp[2], Measured[2], iteration);

	TH2D *Unfolded[3];

	Unfolded[0] = (TH2D *)unfoldcent.Hreco(errorTreatment);
	Unfolded[1] = (TH2D *)unfoldmid.Hreco(errorTreatment);
	Unfolded[2] = (TH2D *)unfoldperi.Hreco(errorTreatment);

	TH1D *MeasuredPt[3];
	TH1D *MeasuredZ[3];

	TH1D *TruthPt[3];
	TH1D *TruthZ[3];

	TH1D *UnfoldedPt[3];
	TH1D *UnfoldedZ[3];

	TH1D *UnfoldedZ_JPtBins[3][njpt_gen_bins_var];

	TCanvas *c[3];

	for (int i = 0; i < 3; i++){
		MeasuredPt[i] = (TH1D *)Measured[i]->ProjectionX();
		MeasuredZ[i] = (TH1D *)Measured[i]->ProjectionY();

		TruthPt[i] = (TH1D *)Truth[i]->ProjectionX();
		TruthZ[i] = (TH1D *)Truth[i]->ProjectionY();

		UnfoldedPt[i] = (TH1D *)Unfolded[i]->ProjectionX();
		UnfoldedZ[i] = (TH1D *)Unfolded[i]->ProjectionY();

		for (int jptbin = 1; jptbin <= njpt_gen_bins_var; jptbin++){
			TH2D *h = (TH2D *)Unfolded[i]->Clone("tmp");
			h->GetXaxis()->SetRange(jptbin, jptbin);
			UnfoldedZ_JPtBins[i][jptbin-1] = (TH1D *)h->ProjectionY(Form("ZUnf_%i_%i", i, jptbin));	
			// UnfoldedZ_JPtBins[i][jptbin-1] = NULL;
		}

		c[i] = new TCanvas(Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 800, 800);
		c[i]->Divide(2);

		SetName(MeasuredPt[i], Form("MeasuredPt_%i", i));
		SetName(TruthPt[i], Form("TruthPt_%i", i));
		SetName(UnfoldedPt[i], Form("UnfoldedPt_%i", i));
		SetColor(MeasuredPt[i], kBlue, 24);
		SetColor(TruthPt[i], kRed, 24);
		SetColor(UnfoldedPt[i], kBlack, 20);

		c[i]->cd(1);
		MeasuredPt[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		MeasuredPt[i]->Draw("EP");
		TruthPt[i]->Draw("EP SAME");
		UnfoldedPt[i]->Draw("EP SAME");

		gPad->SetLogy();

		SetName(MeasuredZ[i], Form("MeasuredZ_%i", i));
		SetName(TruthZ[i], Form("TruthZ_%i", i));
		SetName(UnfoldedZ[i], Form("UnfoldedZ_%i", i));
		SetColor(MeasuredZ[i], kBlue, 24);
		SetColor(TruthZ[i], kRed, 24);
		SetColor(UnfoldedZ[i], kBlack, 20);

		c[i]->cd(2);
		MeasuredZ[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		MeasuredZ[i]->Draw("EP");
		TruthZ[i]->Draw("EP SAME");
		UnfoldedZ[i]->Draw("EP SAME");

		gPad->SetLogy();
	}

	TString FinalOutName = RespFileName;
	FinalOutName.ReplaceAll("Response_", "Output_");
	FinalOutName.ReplaceAll("_Split_1", "");

	TH1D *RatioPt[3];
	TH1D *RatioZ[3];

	TFile *out = new TFile(FinalOutName.Data(), "RECREATE");
	out->cd();

	for (int i = 0; i < 3; i++){
		SetName(Measured[i], Form("Measured_%i", i));
		Measured[i]->Write();
		SetName(Truth[i], Form("Truth_%i", i));
		Truth[i]->Write();
    	SetName(Unfolded[i], Form("Unfolded_%i", i));
		Unfolded[i]->Write();

		MeasuredPt[i]->Write();
		MeasuredZ[i]->Write();

		TruthPt[i]->Write();
		TruthZ[i]->Write();

		UnfoldedPt[i]->Write();
		UnfoldedZ[i]->Write();

		RatioPt[i] = (TH1D *)UnfoldedPt[i]->Clone(Form("RatioPt_%i", i));
		RatioPt[i]->Divide(TruthPt[i]);
		RatioPt[i]->Write();
		RatioZ[i] = (TH1D *)UnfoldedZ[i]->Clone(Form("RatioZ_%i", i));
		RatioZ[i]->Divide(TruthZ[i]);
		RatioZ[i]->Write();

		c[i]->Write();
	}
	for (int i = 0; i < 3; i++){
		for (int jptbin = 1; jptbin <= njpt_gen_bins_var; jptbin++){
			SetName(UnfoldedZ_JPtBins[i][jptbin-1], Form("ZUnf_%i_%i", i, jptbin));
			UnfoldedZ_JPtBins[i][jptbin-1]->Write();
		}
	}
	out->Close();

}

void UnfoldMC(TString DirName = "MCMCUnf", int SUPERITERATION = 1, int iteration = 4, bool mimicDataUncertainty = kFALSE){
	Method(DirName, SUPERITERATION, iteration, mimicDataUncertainty);
}