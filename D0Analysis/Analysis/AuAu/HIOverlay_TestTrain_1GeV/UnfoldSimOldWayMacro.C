void UnfoldSimOldWayMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "DataMCUnf", std::string oldresponse = "true"){
    gSystem->Load("/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L UnfoldTheSimFoldedDistributionOldWay.C++");
    gROOT->ProcessLine(Form(".x UnfoldTheSimFoldedDistributionOldWay.C+(%i, %i, \"%s\", %s)", SUPERITERATION, iteration, dir.c_str(), oldresponse.c_str()));
}