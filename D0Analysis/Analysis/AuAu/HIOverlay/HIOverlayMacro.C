void HIOverlayMacro(){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L HIOverlay.C++");
    gROOT->ProcessLine(Form(".x HIOverlay.C+()"));
}