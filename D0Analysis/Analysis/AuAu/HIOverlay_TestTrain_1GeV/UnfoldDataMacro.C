void UnfoldDataMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "DataMCUnf", std::string oldresponse = "false"){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L UnfoldTheDataFoldedDistribution.C++");
    gROOT->ProcessLine(Form(".x UnfoldTheDataFoldedDistribution.C+(%i, %i, \"%s\", %s)", SUPERITERATION, iteration, dir.c_str(), oldresponse.c_str()));
}