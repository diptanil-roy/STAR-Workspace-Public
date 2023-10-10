// void CreateOldRespMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "MCMCUnf", int mode = 0){
//     gSystem->Load("/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so");
//     gROOT->ProcessLine(".L PrepareResponseMatrixOldWayUsingHI.C++");
//     cout << dir.c_str() << endl;
//     gROOT->ProcessLine(Form(".x PrepareResponseMatrixOldWayUsingHI.C+(%i, %i, \"%s\", false, %i)", SUPERITERATION, iteration, dir.c_str(), mode));
// }

// void CreateOldRespMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "MCMCUnf", int mode = 0){
//     gSystem->Load("/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so");
//     gROOT->ProcessLine(".L PrepareResponseMatrixOldWayUsingHI_2D.C++");
//     cout << dir.c_str() << endl;
//     gROOT->ProcessLine(Form(".x PrepareResponseMatrixOldWayUsingHI_2D.C+(%i, %i, \"%s\", false, %i)", SUPERITERATION, iteration, dir.c_str(), mode));
// }

void CreateOldRespMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "MCMCUnf", int mode = 0){
    gSystem->Load("/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L PrepareResponseMatrixOldWayUsingHI_ChangedCriterion.C++");
    cout << dir.c_str() << endl;
    gROOT->ProcessLine(Form(".x PrepareResponseMatrixOldWayUsingHI_ChangedCriterion.C+(%i, %i, \"%s\", false, %i)", SUPERITERATION, iteration, dir.c_str(), mode));
}