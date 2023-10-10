void HIOverlayMacro(int iteration = 7){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L HIOverlay.C++");
    gROOT->ProcessLine(Form(".x HIOverlay.C+(%i)", iteration));
}