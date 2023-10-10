void UnfoldDataOldWayMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "DataMCUnf", std::string oldresponse = "true", std::string respname = "Resp1D"){
    gSystem->Load("/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L UnfoldTheDataFoldedDistributionOldWay.C++");
    gROOT->ProcessLine(Form(".x UnfoldTheDataFoldedDistributionOldWay.C+(%i, %i, \"%s\", %s, \"%s\")", SUPERITERATION, iteration, dir.c_str(), oldresponse.c_str(), respname.c_str()));
}