void Flattener(){

    THnSparse *Step1;
    THnSparse *Step2;
    THnSparse *Step3;

    const int nDim = 4;
    int nBins[nDim] = {10, 10, 10, 10};
    double xlow[nDim] = {5, 0, 5, 0};
    double xup[nDim] = {10, 1, 10, 1};

    Step1 = new THnSparseD("Step1", "Step1", nDim, nBins, xlow, xup);
    Step2 = new THnSparseD("Step2", "Step2", nDim, nBins, xlow, xup);
    Step3 = new THnSparseD("Step3", "Step3", nDim, nBins, xlow, xup);

    Step1->Sumw2();
    Step1->CalculateErrors();

    Step2->Sumw2();
    Step2->CalculateErrors();

    Step3->Sumw2();
    Step3->CalculateErrors();

    TRandom3 *r = new TRandom3(0);

    int* coord = new int[nDim];

    for (int i = 0; i < 1000000; i++){
        double value[nDim] = {0};
        value[0] = r->Uniform(5, 10);
        value[1] = r->Gaus(0.5, 0.5);
        value[2] = r->Uniform(5, 10);
        value[3] = r->Uniform(0, 1);

        if (value[1] < 0 || value[1] > 1) continue;
        
        int bin = Step1->GetBin(value);
        Double_t w = Step1->GetBinContent(bin, coord);
        // cout << coord[0] << " " << coord[1] << " " << coord[2] << " " << coord[3] << endl;
        if ((coord[0] == 1 && coord[1] == 3) ||
            (coord[0] == 1 && coord[1] == 4) ||
            (coord[0] == 3 && coord[1] == 2) ||
            (coord[0] == 6 && coord[1] == 5) ||
            (coord[0] == 6 && coord[1] == 7) ||
            (coord[0] == 6 && coord[1] == 9) ||
            (coord[0] == 7 && coord[1] == 2) ||
            (coord[0] == 7 && coord[1] == 3) ||
            (coord[0] == 8 && coord[1] == 2) ||
            (coord[0] == 9 && coord[1] == 2) ||
            (coord[1] == 1))
        {
            continue;
        }

        else{
            Step1->Fill(value);
        }

        // Step1->Fill(value);
    }

    TH2D *Meas = (TH2D *)Step1->Projection(3, 2, "E");
    TH2D *True = (TH2D *)Step1->Projection(1, 0, "E");

    TH2D *Weight = (TH2D *)True->Clone("Weight");
    for (int gbinx = 1; gbinx <= Weight->GetNbinsX(); gbinx++){
        for (int gbiny = 1; gbiny <= Weight->GetNbinsY(); gbiny++){
            Weight->SetBinContent(gbinx, gbiny, 1);
        }
    }

    TH1D *EmptyBinsVsPt = new TH1D("EmptyBinsVsPt", "EmptyBinsVsPt", nBins[1], xlow[1], xup[1]);

    for (int gbiny = 1; gbiny <= EmptyBinsVsPt->GetNbinsX(); gbiny++){
        int emptybin = 0;
        for (int gbinx = 1; gbinx <= True->GetNbinsX(); gbinx++){
            if (True->GetBinContent(gbinx, gbiny) == 0){
                emptybin++;
            }
        }
        EmptyBinsVsPt->SetBinContent(gbiny, (1.0*nBins[0])/(1.0*nBins[0]-1.0*emptybin));
    }

    TH2D *WeightEmptyAcc = (TH2D *)Weight->Clone("WeightEmptyAcc");
    for (int gbinx = 1; gbinx <= WeightEmptyAcc->GetNbinsX(); gbinx++){
        for (int gbiny = 1; gbiny <= WeightEmptyAcc->GetNbinsY(); gbiny++){
            WeightEmptyAcc->SetBinContent(gbinx, gbiny, Weight->GetBinContent(gbinx, gbiny) * EmptyBinsVsPt->GetBinContent(gbiny));
        }
    }

    for (int gbinx = 1; gbinx <= nBins[0]; gbinx++){
        for (int gbiny = 1; gbiny <= nBins[1]; gbiny++){
            double scalefactor = 0;
			double valintegral = 0;
            for (int rbinx = 1; rbinx <= nBins[2]; rbinx++ ){
			    for (int rbiny = 1; rbiny <= nBins[3]; rbiny++ ){
                    double gptmid = Step1->GetAxis(0)->GetBinCenter(gbinx);
					double gzmid = Step1->GetAxis(1)->GetBinCenter(gbiny);
					double rptmid = Step1->GetAxis(2)->GetBinCenter(rbinx);
					double rzmid = Step1->GetAxis(3)->GetBinCenter(rbiny);

                    double xbin[nDim] = {gptmid, gzmid, rptmid, rzmid};

					int bin = Step1->GetBin(xbin);

					valintegral+=Step1->GetBinContent(bin);
                }
            }

            if (valintegral > 0){
                scalefactor = WeightEmptyAcc->GetBinContent(gbinx, gbiny)/(1.0*valintegral);
            }

            for (int rbinx = 1; rbinx <= nBins[2]; rbinx++ ){
			    for (int rbiny = 1; rbiny <= nBins[2]; rbiny++ ){
                    double gptmid = Step1->GetAxis(0)->GetBinCenter(gbinx);
					double gzmid = Step1->GetAxis(1)->GetBinCenter(gbiny);
					double rptmid = Step1->GetAxis(2)->GetBinCenter(rbinx);
					double rzmid = Step1->GetAxis(3)->GetBinCenter(rbiny);

                    double xbin[nDim] = {gptmid, gzmid, rptmid, rzmid};

                    int bin = Step1->GetBin(xbin);
                    int bin2 = Step2->GetBin(xbin);

                    double content = Step2->GetBinContent(bin2);
                    double err = Step2->GetBinError(bin2);

                    double newbincontent = content + Step1->GetBinContent(bin)*scalefactor;
                    double newbinerr = err + Step1->GetBinError(bin)*scalefactor;

                    Step2->SetBinContent(bin2, newbincontent);
                    Step2->SetBinError(bin2, newbinerr);
                }
            }
        }
    }

    TH2D *Meas2 = (TH2D *)Step2->Projection(3, 2, "E");
    TH2D *True2 = (TH2D *)Step2->Projection(1, 0, "E");

    for (int gbinx = 1; gbinx <= nBins[0]; gbinx++){
        double scalefactor = 0;
		double origIntegral = 0;
		double scalIntegral = 0;
        for (int gbiny = 1; gbiny <= nBins[1]; gbiny++){
            for (int rbinx = 1; rbinx <= nBins[2]; rbinx++ ){
			    for (int rbiny = 1; rbiny <= nBins[3]; rbiny++ ){
                    double gptmid = Step1->GetAxis(0)->GetBinCenter(gbinx);
					double gzmid = Step1->GetAxis(1)->GetBinCenter(gbiny);
					double rptmid = Step1->GetAxis(2)->GetBinCenter(rbinx);
					double rzmid = Step1->GetAxis(3)->GetBinCenter(rbiny);

                    double xbin[nDim] = {gptmid, gzmid, rptmid, rzmid};

                    int bin = Step1->GetBin(xbin);
                    int bin2 = Step2->GetBin(xbin);

                    origIntegral+=Step1->GetBinContent(bin);
                    scalIntegral+=Step2->GetBinContent(bin2);
                }
            }
        }
        if (scalIntegral > 0){
            scalefactor = 1.0*origIntegral/scalIntegral;
        }

        for (int gbiny = 1; gbiny <= nBins[1]; gbiny++){
            for (int rbinx = 1; rbinx <= nBins[2]; rbinx++ ){
			    for (int rbiny = 1; rbiny <= nBins[3]; rbiny++ ){
                    double gptmid = Step1->GetAxis(0)->GetBinCenter(gbinx);
					double gzmid = Step1->GetAxis(1)->GetBinCenter(gbiny);
					double rptmid = Step1->GetAxis(2)->GetBinCenter(rbinx);
					double rzmid = Step1->GetAxis(3)->GetBinCenter(rbiny);

                    double xbin[nDim] = {gptmid, gzmid, rptmid, rzmid};

                    int bin2 = Step2->GetBin(xbin);
                    int bin3 = Step3->GetBin(xbin);

                    double content = Step3->GetBinContent(bin3);
                    double err = Step3->GetBinError(bin3);

                    double newbincontent = content + Step2->GetBinContent(bin2)*scalefactor;
                    double newbinerr = err + Step2->GetBinError(bin2)*scalefactor;
                    Step3->SetBinContent(bin3, newbincontent);
                    Step3->SetBinError(bin3, newbinerr);
                }
            }
        }
    }

    TH2D *Meas3 = (TH2D *)Step3->Projection(3, 2, "E");
    TH2D *True3 = (TH2D *)Step3->Projection(1, 0, "E");

    gSystem->Exec("mkdir -p FlattenerTest");
    TFile *f = new TFile("FlattenerTest/Test.root", "RECREATE");
    f->cd();
    Step1->Write();
    Meas->Write();
    True->Write();
    EmptyBinsVsPt->Write();
    Weight->Write();
    WeightEmptyAcc->Write();
    Step2->Write();
    Meas2->Write();
    True2->Write();
    Step3->Write();
    Meas3->Write();
    True3->Write();
    f->Close();
}