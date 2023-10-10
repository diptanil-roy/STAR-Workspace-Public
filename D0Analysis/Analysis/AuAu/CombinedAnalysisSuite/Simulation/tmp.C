void tmp() {
    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Plot Drawing", 800, 600);

    // Set the user ranges for the axes
    double xmin = 4.7;
    double xmax = 20.1;
    double ymin = 0.0001;
    double ymax = 1.9999;

    // Create a TH2D histogram
    TH2D *hist = new TH2D("hist", "Plot Title;X Axis;Y Axis", 10, xmin, xmax, 10, ymin, ymax);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetYaxis()->SetTitleSize(0.08);
    hist->Draw("AXIS"); // Draw the histogram to create the axes with labels


    // Draw the canvas
    canvas->Update();
    canvas->Draw();
}