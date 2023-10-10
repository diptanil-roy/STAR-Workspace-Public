void test(){

    // gStyle->SetPalette(kVisibleSpectrum);
    const int numcolors = 200;
    TF1 *f[numcolors+1];

    static int colors[numcolors+1];
    Double_t Red[3]    = { 1.00, 0.00, 0.00};
    Double_t Green[3]  = { 0.00, 1.00, 0.00};
    Double_t Blue[3]   = { 1.00, 0.00, 1.00};
    Double_t Length[3] = { 0.00, 0.50, 1.00 };

    // Double_t Red[2]    = { 1.00, 0.00};
    // Double_t Green[2]  = { 0.00, 0.00};
    // Double_t Blue[2]   = { 0.00, 0.00};
    // Double_t Length[2] = { 0.00, 1.00 };

    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,numcolors+1);
    // Int_t FI = TColor::CreateGradientColorTable(2,Length,Red,Green,Blue,numcolors+1);
    // Int_t FI = TColor::SetPalette(kVisibleSpectrum, 0);
    for (int i=0; i<=numcolors; i++) colors[i] = FI+i;

    for (int i = 0; i <= numcolors; i++){
        f[i] = new TF1(Form("f%i",i), Form("%i", i), 0,1);
        f[i]->SetLineColor(colors[i]);
        if (i == 0) {f[i]->Draw(); f[i]->SetMaximum(250); f[i]->SetMinimum(-1);}
        else f[i]->Draw("lsame");
    }

}