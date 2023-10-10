R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, int mode = 0){

	cout << "File Reader Method" << endl;

	// for (int i = 0; i < njpt_gen_bins_var+1; i++){
	// 	jetpt_var_bin[i] = 5 + (20. - 5.)/njpt_gen_bins_var * i;
	// }

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	THnSparseF *hJ[3];
	THnSparseF *hJet[3];
	THnSparseF *hJetZNorm[3];
	THnSparseF *hJetZNormPtNorm[3];

	TH2D *Weight[3];
	TH1D *WeightZ[3];

	RooUnfoldResponse *resp[3]; //SuperIteration Response Matrices are defined here

	TH2D *tmpMeas[3];
	TH2D *tmpTrue[3];

	for (int i = 0; i < 3; i++){
		tmpMeas[i] = new TH2D(Form("Meas_%i", i), Form("Meas_%i", i), njpt_bins, nbinsjetpt, nz_bins, nbinsz);
		tmpTrue[i] = new TH2D(Form("True_%i", i), Form("True_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);

		hJetZNorm[i] = new THnSparseF(Form("hJetZNorm_%i", i), Form("hJetZNorm_%i", i), ndim, nbins, NULL, NULL);
		hJetZNormPtNorm[i] = new THnSparseF(Form("hJetZNormPtNorm_%i", i), Form("hJetZNormPtNorm_%i", i), ndim, nbins, NULL, NULL);

		hJetZNorm[i]->SetBinEdges(0, jetpt_var_bin);
		hJetZNorm[i]->SetBinEdges(1, z_gen_bin);
		hJetZNorm[i]->SetBinEdges(2, nbinsjetpt);
		hJetZNorm[i]->SetBinEdges(3, nbinsz);

		hJetZNormPtNorm[i]->SetBinEdges(0, jetpt_var_bin);
		hJetZNormPtNorm[i]->SetBinEdges(1, z_gen_bin);
		hJetZNormPtNorm[i]->SetBinEdges(2, nbinsjetpt);
		hJetZNormPtNorm[i]->SetBinEdges(3, nbinsz);

		hJetZNorm[i]->Sumw2();
		hJetZNormPtNorm[i]->Sumw2();
		hJetZNorm[i]->CalculateErrors();
		hJetZNormPtNorm[i]->CalculateErrors();

		resp[i] = new RooUnfoldResponse(Form("Resp_%i", i), Form("Resp_%i", i));
		resp[i]->Setup(tmpMeas[i], tmpTrue[i]); //Setup Response Matrix Definition
	}

	TFile *RespHistogramFile;

	if (mode == 0) RespHistogramFile = new TFile(Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
	else RespHistogramFile = new TFile(Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins, mode), "READ"); // Split 1 is used for Response Matrix. // Split 2 is the file unfolded.

	cout << "Response Matrix File = " << RespHistogramFile->GetName() << endl;

	RespHistogramFile->cd();

	for (int i = 0; i < 3; i++){
		hJ[i] = (THnSparseF *)gDirectory->Get(Form("hJ_%i", i));
		hJet[i] = (THnSparseF *)gDirectory->Get(Form("hJet_%i", i));
		// hJet[i]->SetDirectory(0);
	}

	// RespHistogramFile->Close();

	TFile *DataFile = new TFile(Form("../SPlot_Mar2_2023/Histograms3_D0%i_10GeV_NewBinningOldCuts_Apr6_WiderJpTRange.root", D0pTLow), "READ");
	DataFile->cd();

	TH2D *DataMeasured[3];

	DataMeasured[0] = (TH2D *)gDirectory->Get("ZPt_0_10");
	DataMeasured[1] = (TH2D *)gDirectory->Get("ZPt_10_40");
	DataMeasured[2] = (TH2D *)gDirectory->Get("ZPt_40_80");

	TH1D *DataMeasuredZ[3];

	for (int i = 0; i < 3; i++){
		DataMeasuredZ[i] = (TH1D *)DataMeasured[i]->ProjectionY();
	}
	
	TH1D *TruthZ[3];
	TH1D *TruthPt[3];
	TH1D *MeasuredPt[3];
	TH1D *MeasuredZ[3];

	for (int cent = 0; cent < 3; cent++){
		TruthZ[cent] = (TH1D *)hJet[cent]->Projection(1, "E");
		TruthPt[cent] = (TH1D *)hJet[cent]->Projection(0, "E");
		MeasuredPt[cent] = (TH1D *)hJet[cent]->Projection(2, "E");
		MeasuredZ[cent] = (TH1D *)hJet[cent]->Projection(3, "E");

		WeightZ[cent] = (TH1D *)DataMeasuredZ[cent]->Clone(Form("WeightZ_%i", cent));
		WeightZ[cent]->Divide(MeasuredZ[cent]);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Changing the Reco pT and Z distribution to match the data
	for (int cent = 0; cent < 3; cent++){
		for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){ // Reweighing on the basis of the z distribution from data
			for (int gbinx = 1; gbinx <= njpt_gen_bins_var; gbinx++ ){
				for (int gbiny = 1; gbiny <= nz_gen_bins; gbiny++ ){
					for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
						double gptmid = hJet[cent]->GetAxis(0)->GetBinCenter(gbinx);
						double gzmid = hJet[cent]->GetAxis(1)->GetBinCenter(gbiny);
						double rptmid = hJet[cent]->GetAxis(2)->GetBinCenter(rbinx);
						double rzmid = hJet[cent]->GetAxis(3)->GetBinCenter(rbiny);

						double xbin[ndim] = {gptmid, gzmid, rptmid, rzmid};

						int bin = hJet[cent]->GetBin(xbin);
						int binZNorm = hJetZNorm[cent]->GetBin(xbin);

						double bincontent = hJet[cent]->GetBinContent(bin);
						double binerror = hJet[cent]->GetBinError(bin);

						double scalefactor = 1.0;
						if (SUPERITERATION == 1) scalefactor = WeightZ[cent]->GetBinContent(rbiny);

						double newbincontent = bincontent * scalefactor;
						double newbinerror = binerror * scalefactor;

						hJetZNorm[cent]->SetBinContent(binZNorm, newbincontent);
						hJetZNorm[cent]->SetBinError(binZNorm, newbinerror);
					}
				}
			}
		}
	} //End of centrality loop

	TH1D *TruthZZNorm[3];
	TH1D *TruthPtZNorm[3];

	for (int cent = 0; cent < 3; cent++){
		TruthZZNorm[cent] = (TH1D *)hJetZNorm[cent]->Projection(1, "E");
		TruthPtZNorm[cent] = (TH1D *)hJetZNorm[cent]->Projection(0, "E");

		cout << "Integral of Truth Pt = " << TruthPtZNorm[cent]->Integral() << endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Changing the pT distribution to go back to the original distribution
	for (int cent = 0; cent < 3; cent++){
		for (int gbinx = 1; gbinx <= njpt_gen_bins_var; gbinx++ ){
			double scalefactor = 0;
			double origIntegral = 0;
			double scalIntegral = 0;
			for (int gbiny = 1; gbiny <= nz_gen_bins; gbiny++ ){
				for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
					for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){
						double gptmid = hJetZNorm[cent]->GetAxis(0)->GetBinCenter(gbinx);
						double gzmid = hJetZNorm[cent]->GetAxis(1)->GetBinCenter(gbiny);
						double rptmid = hJetZNorm[cent]->GetAxis(2)->GetBinCenter(rbinx);
						double rzmid = hJetZNorm[cent]->GetAxis(3)->GetBinCenter(rbiny);
						// double gptmid = jetpt_var_bin[0] + (gbinx - 1)*jetpt_gen_binwidth + jetpt_gen_binwidth*0.5;
						// double gzmid = z_gen_low + (gbiny - 1)*z_gen_binwidth + z_gen_binwidth*0.5; 
						// double rptmid = jetpt_low + (rbinx - 1)*jetpt_binwidth + jetpt_binwidth*0.5;
						// double rzmid = z_low + (rbiny - 1)*z_binwidth + z_binwidth*0.5;
						double xbin[ndim] = {gptmid, gzmid, rptmid, rzmid};
						int binZnorm = hJetZNorm[cent]->GetBin(xbin);
						int bin = hJet[cent]->GetBin(xbin);
						origIntegral+=hJet[cent]->GetBinContent(bin);
						scalIntegral+=hJetZNorm[cent]->GetBinContent(binZnorm);
					}
				}
				// origIntegral+=fMiss[cent]->GetBinContent(gbinx, gbiny);
				// scalIntegral+=fMissZNorm[cent]->GetBinContent(gbinx, gbiny);
			}
			if (scalIntegral > 0) {
				if (SUPERITERATION == 0) scalefactor = 1.;
				else scalefactor = 1.0*origIntegral/scalIntegral;
			}
			for (int gbiny = 1; gbiny <= nz_gen_bins; gbiny++ ){
				for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
					for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){
						double gptmid = hJetZNorm[cent]->GetAxis(0)->GetBinCenter(gbinx);
						double gzmid = hJetZNorm[cent]->GetAxis(1)->GetBinCenter(gbiny);
						double rptmid = hJetZNorm[cent]->GetAxis(2)->GetBinCenter(rbinx);
						double rzmid = hJetZNorm[cent]->GetAxis(3)->GetBinCenter(rbiny);
						// double gptmid = jetpt_var_bin[0] + (gbinx - 1)*jetpt_gen_binwidth + jetpt_gen_binwidth*0.5;
						// double gzmid = z_gen_low + (gbiny - 1)*z_gen_binwidth + z_gen_binwidth*0.5; 
						// double rptmid = jetpt_low + (rbinx - 1)*jetpt_binwidth + jetpt_binwidth*0.5;
						// double rzmid = z_low + (rbiny - 1)*z_binwidth + z_binwidth*0.5; 
						double xbin[ndim] = {gptmid, gzmid, rptmid, rzmid};
						// if (cent == 0) cout << gptmid << "\t" << gzmid << "\t" << rptmid << "\t" << rzmid << endl;
						int binZnorm = hJetZNorm[cent]->GetBin(xbin);
						int binZnormPtnorm = hJetZNormPtNorm[cent]->GetBin(xbin);
						double content = hJetZNormPtNorm[cent]->GetBinContent(binZnormPtnorm);
						double err = hJetZNormPtNorm[cent]->GetBinError(binZnormPtnorm);
						double newbincontent = content + hJetZNorm[cent]->GetBinContent(binZnorm)*scalefactor;
						double newbinerr = err + hJetZNorm[cent]->GetBinError(binZnorm)*scalefactor;
						hJetZNormPtNorm[cent]->SetBinContent(binZnormPtnorm, newbincontent);
						hJetZNormPtNorm[cent]->SetBinError(binZnormPtnorm, newbinerr);
					}
				}
			}
		}
	} //End of centrality loop


	TH1D *TruthZZNormPtNorm[3];
	TH1D *TruthPtZNormPtNorm[3];
	TH2D *MeasuredZNormPtNorm[3];

	for (int cent = 0; cent < 3; cent++){
		TruthZZNormPtNorm[cent] = (TH1D *)hJetZNormPtNorm[cent]->Projection(1, "E");
		TruthPtZNormPtNorm[cent] = (TH1D *)hJetZNormPtNorm[cent]->Projection(0, "E");
		MeasuredZNormPtNorm[cent] = (TH2D *)hJetZNormPtNorm[cent]->Projection(1, 0, "E");

		cout << "Integral of Truth Pt = " << TruthPtZNormPtNorm[cent]->Integral() << endl;
	}

	// Fill RooUnfoldResponse Object

	int* coord = new int[ndim];

	double integralinbinresp = 0.;
	double integralinresp = 0.;

	cout << "================ Limits =================" << endl;
	// cout << Form("%i < pT,D0 [GeV/c] < %i", D0pTLow, 10) << endl;
	cout << Form("%.1f < MC pT,Jet [GeV/c] < %.1f", jetpt_var_bin[0], jetpt_var_bin[njpt_gen_bins_var]) << endl;
	cout << Form("%.1f < Reco pT,Jet [GeV/c] < %.1f", nbinsjetpt[0], nbinsjetpt[njpt_bins]) << endl;
	cout << Form("%.1f < MC Z,Jet [GeV/c] < %.1f", z_gen_bin[0], z_gen_bin[nz_gen_bins]) << endl;
	cout << Form("%.1f < Reco Z,Jet [GeV/c] < %.1f", nbinsz[0], nbinsz[nz_bins]) << endl;
	cout << "=========================================" << endl;


	for (int cent = 0; cent < 3; cent++){
		int nbin = hJetZNormPtNorm[cent]->GetNbins();

		cout << hJetZNormPtNorm[cent]->GetName() << endl;
		cout << resp[cent]->GetName() << endl;

		cout << "Bins = " << nbin << endl;

		for (int bin = 0; bin < nbin; bin++){
			Double_t w = hJetZNormPtNorm[cent]->GetBinContent(bin, coord);
			Double_t pttrue = hJetZNormPtNorm[cent]->GetAxis(0)->GetBinCenter(coord[0]);
			Double_t ztrue = hJetZNormPtNorm[cent]->GetAxis(1)->GetBinCenter(coord[1]);
			Double_t ptdet =  hJetZNormPtNorm[cent]->GetAxis(2)->GetBinCenter(coord[2]); 
			Double_t zdet =  hJetZNormPtNorm[cent]->GetAxis(3)->GetBinCenter(coord[3]);

			// if (cent == 0) cout << bin << "\t" << pttrue << "\t" << ztrue << "\t" << ptdet << "\t" << zdet << "\t" << w << endl;

			if (cent == 0) integralinresp+=w;

			// resp[cent]->Fill(ptdet, zdet, pttrue, ztrue, w); 	

			if (pttrue >= jetpt_var_bin[0] && pttrue <= jetpt_var_bin[njpt_gen_bins_var] &&
				// ztrue >= z_gen_bin[0] && ztrue <= z_gen_bin[nz_gen_bins] &&
				ptdet >= nbinsjetpt[0] && ptdet <= nbinsjetpt[njpt_bins] 
				// && zdet >= nbinsz[0] && zdet <= nbinsz[nz_bins]
			)
			{
				resp[cent]->Fill(ptdet, zdet, pttrue, ztrue, w); 
				if (cent == 0) integralinbinresp+=w;
			}

			else if (ptdet >= nbinsjetpt[0] && ptdet <= nbinsjetpt[njpt_bins]
				// && zdet >= nbinsz[0] && zdet <= nbinsz[nz_bins]
			)
			{
				resp[cent]->Fake(ptdet, zdet, w);
			}

			else if (pttrue >= jetpt_var_bin[0] && pttrue <= jetpt_var_bin[njpt_gen_bins_var] 
			// && ztrue >= z_gen_bin[0] && ztrue <= z_gen_bin[nz_gen_bins]
			)
			{
				resp[cent]->Miss(pttrue, ztrue, w); 
			}	
			// resp[cent]->Fill(ptdet, zdet, pttrue, ztrue, w); 
		}
	}

	delete [] coord;

	cout << "RooUnfoldResponse filled" << endl;

	TH2D *TruthFromResp[3];
	TH1D *TruthZFromResp[3];
	TH1D *TruthPtFromResp[3];

	for (int cent = 0; cent < 3; cent++){
	
		TruthFromResp[cent] = (TH2D *)resp[cent]->Htruth();
		SetName(TruthFromResp[cent], Form("Truth From Resp %i", cent));

		TruthZFromResp[cent] = (TH1D *)TruthFromResp[cent]->ProjectionY();
		// TruthZFromResp[cent]->Rebin(4);
		TruthPtFromResp[cent] = (TH1D *)TruthFromResp[cent]->ProjectionX();
		// TruthPtFromResp[cent]->Rebin(4);
	}

	// cout << TruthZZNormPtNorm[2]->Integral() << "\t" << integralinresp << "\t" << integralinbinresp << "\t" << integralinresp - integralinbinresp << "\t" << TruthFromResp[2]->Integral() << endl;
	
	// This might not be necessary all the time. For now, I want to check all the steps.

	TH1D *RatioPt[3];
	TH1D *RatioZ[3];

	TH1D *TargetPt[3];
	TH1D *TargetZ[3];

	for (int cent = 0; cent < 3; cent++){
		TargetPt[cent] = (TH1D *)WeightZ[cent]->Clone();
		TargetZ[cent] = (TH1D *)DataMeasuredZ[cent]->Clone();
	}

	TCanvas *d[6];

	TString CentBinAppend[3] = {"Central", "MidCentral", "Peripheral"};

	for (int cent = 0; cent < 3; cent++){
		d[cent] = new TCanvas(Form("Testing the weighting mechanism %i", cent), Form("Testing the weighting mechanism %i", cent), 1200, 600);
		d[cent]->Divide(2);
		d[cent]->cd(1);
		gPad->SetLogy();
		SetAxisTitles(TruthPt[cent], "p_{T,Jet} [GeV/#it{c}]", "Normalised Yield");
		SetColor(TruthPt[cent], kBlue, 24);
		SetName(TruthPt[cent], Form("FONLL_%i", cent));
		SetColor(TruthPtZNorm[cent], kRed, 24);
		SetName(TruthPtZNorm[cent], Form("Z norm (2D) %i", cent));
		SetColor(TruthPtZNormPtNorm[cent], kBlack, 20);
		SetName(TruthPtZNormPtNorm[cent], Form("Z norm + FONLL pT %i", cent));
		SetColor(TruthPtFromResp[cent], kGreen-2, 29, 1.4);
		SetName(TruthPtFromResp[cent], Form("From Response %i", cent));
		SetColor(TargetPt[cent], kViolet, 47, 1.4);
		SetName(TargetPt[cent], Form("Target p_{T} %i", cent));
		double scale = TruthPt[cent]->Integral();
		TruthPt[cent]->Scale(1./TruthPt[cent]->Integral());
		TruthPtZNorm[cent]->Scale(1./TruthPtZNorm[cent]->Integral());
		TruthPtZNormPtNorm[cent]->Scale(1./TruthPtZNormPtNorm[cent]->Integral());
		TruthPtFromResp[cent]->Scale(1./TruthPtFromResp[cent]->Integral());
		// TargetPt[cent]->Scale(1./TargetPt[cent]->Integral());

		RatioPt[cent] = (TH1D *)TruthPtFromResp[cent]->Clone(Form("RatioPt_%i", cent));
		RatioPt[cent]->Divide(TruthPtZNormPtNorm[cent]);

		TruthPt[cent]->Draw();
		TruthPtZNorm[cent]->Draw("SAME");
		TruthPtZNormPtNorm[cent]->Draw("SAME");
		TruthPtFromResp[cent]->Draw("SAME");
		// TargetPt[cent]->Draw("HIST P SAME");

		auto legend = new TLegend(0.4,0.85,0.78,0.9);
   		legend->SetHeader(Form("Superiteration %i", SUPERITERATION),"C"); // option "C" allows to center the header
   		legend->Draw("SAME");

		gPad->BuildLegend();
	// d[cent]->Update();

		d[cent]->cd(2);
		gPad->SetLogy();
		SetAxisTitles(TruthZ[cent], "Z", "Normalised Yield");
		SetColor(TruthZ[cent], kBlue, 24);
		SetName(TruthZ[cent], Form("FONLL_%i", cent));
		SetColor(TruthZZNorm[cent], kRed, 24);
		SetName(TruthZZNorm[cent], Form("Z norm (2D) %i", cent));
		SetColor(TruthZZNormPtNorm[cent], kBlack, 20);
		SetName(TruthZZNormPtNorm[cent], Form("Z norm + FONLL pT %i", cent));
		SetColor(TruthZFromResp[cent], kGreen-2, 29, 1.4);
		SetName(TruthZFromResp[cent], Form("From Response %i", cent));
		SetColor(TargetZ[cent], kViolet, 47, 1.4);
		SetName(TargetZ[cent], Form("Target Z %i", cent));
		scale = TruthZ[cent]->Integral();
		TruthZ[cent]->Scale(1./TruthZ[cent]->Integral());
		TruthZZNorm[cent]->Scale(1./TruthZZNorm[cent]->Integral());
		TruthZZNormPtNorm[cent]->Scale(1./TruthZZNormPtNorm[cent]->Integral());
		TruthZFromResp[cent]->Scale(1./TruthZFromResp[cent]->Integral());
		// TargetZ[cent]->Scale(1./TargetZ[cent]->Integral());

		RatioZ[cent] = (TH1D *)TruthZFromResp[cent]->Clone(Form("RatioZ_%i", cent));
		RatioZ[cent]->Divide(TruthZZNormPtNorm[cent]);

		TruthZ[cent]->Draw();
		TruthZZNorm[cent]->Draw("SAME");
		TruthZZNormPtNorm[cent]->Draw("SAME");
		TruthZFromResp[cent]->Draw("SAME");
		// TargetZ[cent]->Draw("HIST P SAME");
		gPad->BuildLegend();

		if (mode == 2) d[cent]->SaveAs(Form("%s/%s/Mode_%i_Step_%i.pdf", DirName.Data(), CentBinAppend[cent].Data(), mode, SUPERITERATION));
		else d[cent]->SaveAs(Form("%s/%s/Step_%i.pdf", DirName.Data(), CentBinAppend[cent].Data(), SUPERITERATION));
	}

	d[3] = new TCanvas(Form("Testing the weighting mechanism %i", 3), Form("Testing the weighting mechanism %i", 3), 1200, 600);
	d[3]->Divide(2);
	d[3]->cd(1);
	for (int cent = 0; cent < 2; cent++){
		SetColor(RatioPt[cent], col[cent]);
		SetName(RatioPt[cent], Form("RatioPt_%i", cent));
		RatioPt[cent]->Draw("E0 SAME");
	}
	gPad->BuildLegend();

	d[3]->cd(2);
	for (int cent = 0; cent < 2; cent++){
		SetColor(RatioZ[cent], col[cent]);
		SetName(RatioZ[cent], Form("RatioZ_%i", cent));
		RatioZ[cent]->Draw("E0 SAME");
	}
	gPad->BuildLegend();

	d[4] = new TCanvas(Form("Testing the weighting mechanism %i", 4), Form("Testing the weighting mechanism %i", 4), 1200, 600);
	d[4]->Divide(2);
	d[4]->cd(1);
	TH2D *TrueFromResp = (TH2D *)resp[1]->Htruth();
	// MeasuredZNormPtNorm[1]->Divide(TrueFromResp);
	MeasuredZNormPtNorm[1]->Draw("COLZ");
	// MeasuredZNormPtNorm[1]->GetZaxis()->SetRangeUser(0.99, 1.01);
	d[4]->cd(2);
	resp[0]->Htruth()->Draw("COLZ");

	for (int binx = 1; binx <= MeasuredZNormPtNorm[1]->GetNbinsX(); binx++){
		for (int biny = 1; biny <= MeasuredZNormPtNorm[1]->GetNbinsY(); biny++){
			double x = MeasuredZNormPtNorm[1]->GetXaxis()->GetBinCenter(binx);
			double y = MeasuredZNormPtNorm[1]->GetYaxis()->GetBinCenter(biny);
			double w = MeasuredZNormPtNorm[1]->GetBinContent(binx, biny);
			double w2 = TrueFromResp->GetBinContent(binx, biny);

			// if (w != w2)  cout << x << "\t" << y << "\t" << w << "\t" << w2 << endl;
		}
	}

	d[5] = new TCanvas(Form("Testing the weighting mechanism %i", 5), Form("Testing the weighting mechanism %i", 5), 1200, 600);
	d[5]->cd();

	TH1D *tmp = (TH1D *)TruthZFromResp[0]->Clone();
	tmp->Divide(TruthZFromResp[2]);

	TH1D *tmp2 = (TH1D *)TruthZZNormPtNorm[0]->Clone();
	tmp2->Divide(TruthZZNormPtNorm[2]);

	tmp->Draw("P SAME");
	tmp2->Draw("P SAME");


	TFile *OutFile;

	if (mode == 0) OutFile = new TFile(Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins_var, nz_gen_bins), "RECREATE");
	else OutFile = new TFile(Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins_var, nz_gen_bins, mode), "RECREATE");
	cout << OutFile->GetName() << endl;
	OutFile->cd();

	for (int i = 0; i < 3; i++){
		hJet[i]->Write();
		// hJetZNorm[i]->Write();
		hJetZNormPtNorm[i]->Write();

		d[i]->Write();

		resp[i]->Write();
	}

	OutFile->Close();
}

void CreateResponseMatrixWithDataZ(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, int mode = 0){ // Mode 0 = Full MC; Mode 1 = Split 1; Mode 2 = Split 2
	Method(DirName, D0pTLow, SUPERITERATION, iteration, mode);
}