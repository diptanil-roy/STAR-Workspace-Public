void FlattenerWithTH2(){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TH2D *Step1;
    TH2D *Step2;
    TH2D *Step3;

    const int nDim = 2;
    int nBins[nDim] = {10, 10};
    double xlow[nDim] = {5, 0};
    double xup[nDim] = {10, 1};

    Step1 = new TH2D("Step1", "Step1", nBins[0], xlow[0], xup[0], nBins[1], xlow[1], xup[1]);
    Step2 = new TH2D("Step2", "Step2", nBins[0], xlow[0], xup[0], nBins[1], xlow[1], xup[1]);
    Step3 = new TH2D("Step3", "Step3", nBins[0], xlow[0], xup[0], nBins[1], xlow[1], xup[1]);
    
    TRandom3 *r = new TRandom3(0);

    int* coord = new int[nDim];

    for (int i = 0; i < 1000000; i++){
        double value[nDim] = {0};
        value[0] = r->Uniform(5, 10);
        value[1] = r->Gaus(0.5, 0.5);
        
        if (value[1] < 0 || value[1] > 1) continue;

        coord[0] = Step1->GetXaxis()->FindBin(value[0]);
        coord[1] = Step1->GetYaxis()->FindBin(value[1]);
        
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
            Step1->Fill(value[0], value[1]);
        }

        // Step1->Fill(value[0], value[1]);
    }

    TH1D *Step1X = (TH1D *)Step1->ProjectionX("Step1X");
    TH1D *Step1Y = (TH1D *)Step1->ProjectionY("Step1Y");


    TH2D *Weight = (TH2D *)Step1->Clone("Weight");
    for (int gbinx = 1; gbinx <= Weight->GetNbinsX(); gbinx++){
        for (int gbiny = 1; gbiny <= Weight->GetNbinsY(); gbiny++){
            Weight->SetBinContent(gbinx, gbiny, 1);
        }
    }

    TH1D *EmptyBinsVsPt = new TH1D("EmptyBins", "EmptyBins", nBins[1], xlow[1], xup[1]);

    for (int gbiny = 1; gbiny <= EmptyBinsVsPt->GetNbinsX(); gbiny++){
        int emptybin = 0;
        for (int gbinx = 1; gbinx <= Step1->GetNbinsX(); gbinx++){
            if (Step1->GetBinContent(gbinx, gbiny) == 0){
                emptybin++;
            }
        }
        if (emptybin != nBins[0]) EmptyBinsVsPt->SetBinContent(gbiny, (1.0*nBins[0])/(1.0*nBins[0]-1.0*emptybin));
        else EmptyBinsVsPt->SetBinContent(gbiny, 0);
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
			double valintegral = Step1->GetBinContent(gbinx, gbiny);
            
            if (valintegral > 0){
                scalefactor = WeightEmptyAcc->GetBinContent(gbinx, gbiny)/(1.0*valintegral);
            }

            Step2->SetBinContent(gbinx, gbiny, scalefactor*Step1->GetBinContent(gbinx, gbiny));
            Step2->SetBinError(gbinx, gbiny, scalefactor*Step1->GetBinError(gbinx, gbiny));
        }
    }

    TH1D *Step2X = (TH1D *)Step2->ProjectionX("Step2X");
    TH1D *Step2Y = (TH1D *)Step2->ProjectionY("Step2Y");


    for (int gbinx = 1; gbinx <= nBins[0]; gbinx++){
        double scalefactor = 0;
		double origIntegral = 0;
		double scalIntegral = 0;
        for (int gbiny = 1; gbiny <= nBins[1]; gbiny++){
            origIntegral+=Step1->GetBinContent(gbinx, gbiny);
            scalIntegral+=Step2->GetBinContent(gbinx, gbiny);
        }
        if (scalIntegral > 0){
            scalefactor = 1.0*origIntegral/scalIntegral;
        }

        for (int gbiny = 1; gbiny <= nBins[1]; gbiny++){
            Step3->SetBinContent(gbinx, gbiny, scalefactor*Step2->GetBinContent(gbinx, gbiny));
            Step3->SetBinError(gbinx, gbiny, scalefactor*Step2->GetBinError(gbinx, gbiny));
        }
    }

    TH1D *Step3X = (TH1D *)Step3->ProjectionX("Step3X");
    TH1D *Step3Y = (TH1D *)Step3->ProjectionY("Step3Y");

    gSystem->Exec("mkdir -p FlattenerTest");
    TFile *f = new TFile("FlattenerTest/Test.root", "RECREATE");
    f->cd();
    Step1->Write();
    Step1X->Write();
    Step1Y->Write();
    EmptyBinsVsPt->Write();
    Weight->Write();
    WeightEmptyAcc->Write();
    Step2->Write();
    Step2X->Write();
    Step2Y->Write();
    Step3->Write();
    Step3X->Write();
    Step3Y->Write();
    f->Close();
}