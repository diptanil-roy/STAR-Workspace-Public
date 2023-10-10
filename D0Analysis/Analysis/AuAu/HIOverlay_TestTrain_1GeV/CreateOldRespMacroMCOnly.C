void CreateOldRespMacroMCOnly(int SUPERITERATION = 1, int iteration = 3, std::string dir = "MCMCUnf", int mode = 0){
    gSystem->Load("/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L PrepareResponseMatrixUsingMCOnly.C++");
    cout << dir.c_str() << endl;
    gROOT->ProcessLine(Form(".x PrepareResponseMatrixUsingMCOnly.C+(%i, %i, \"%s\", %i)", SUPERITERATION, iteration, dir.c_str(), mode));
}