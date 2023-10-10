void UnfoldMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "MCMCUnf"){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L UnfoldTheMCFoldedDistribution.C++");
    gROOT->ProcessLine(Form(".x UnfoldTheMCFoldedDistribution.C+(%i, %i, \"%s\")", SUPERITERATION, iteration, dir.c_str()));
}