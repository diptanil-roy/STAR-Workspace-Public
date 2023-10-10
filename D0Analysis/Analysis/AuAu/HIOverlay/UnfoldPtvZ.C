void UnfoldPtvZ(int mode = 1){
	gSystem->Load("./RooUnfold/libRooUnfold.so");
    if (mode == 1){
        gROOT->ProcessLine(".L ResponseMatrix.C++");
        gROOT->ProcessLine(".x ResponseMatrix.C+()");
    }
	else if (mode == 2){
        gROOT->ProcessLine(".L JetPtvsZ.C++");
        gROOT->ProcessLine(".x JetPtvsZ.C+()");
    }
}