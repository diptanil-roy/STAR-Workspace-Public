void Chi2Macro(){
	gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L Chi2ForSuperIteration.C++");
    gROOT->ProcessLine(Form(".x Chi2ForSuperIteration.C+()"));
}