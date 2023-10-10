void CreateRespMacro(int SUPERITERATION = 1, int iteration = 3, std::string dir = "MCMCUnf"){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L CreateResponseMatrix.C++");
    cout << dir.c_str() << endl;
    gROOT->ProcessLine(Form(".x CreateResponseMatrix.C+(%i, %i, \"%s\")", SUPERITERATION, iteration, dir.c_str()));
}