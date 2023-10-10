void HIOverlayMacro(int iteration = 3){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L HIOverlay.C++");
    gROOT->ProcessLine(Form(".x HIOverlay.C+(%i)", iteration));

    // gROOT->ProcessLine(".L HIOverlayData.C++");
    // gROOT->ProcessLine(Form(".x HIOverlayData.C+(%i)", iteration));

    // gROOT->ProcessLine(".L HIIter.C++");
    // gROOT->ProcessLine(Form(".x HIIter.C+()"));
    // gROOT->ProcessLine(".L HIOverlayNormalised.C++");
    // gROOT->ProcessLine(Form(".x HIOverlayNormalised.C+(%i)", iteration));
}