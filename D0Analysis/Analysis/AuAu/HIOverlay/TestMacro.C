void TestMacro(){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L TestForSampling.C++");
    gROOT->ProcessLine(Form(".x TestForSampling.C+()"));
}