void ClosureMacro(int superiter = 100, int iter = 4){
	gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L Closure.C++");
    gROOT->ProcessLine(Form(".x Closure.C+(%i, %i)", superiter, iter));
}