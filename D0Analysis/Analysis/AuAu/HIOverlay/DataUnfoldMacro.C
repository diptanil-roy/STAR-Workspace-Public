void DataUnfoldMacro(int superiter = 2, int iter = 4){
	gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".x myStyle.C");
    gROOT->ProcessLine(".L DataUnfold.C++");
    gROOT->ProcessLine(Form(".x DataUnfold.C+(%i, %i)", superiter, iter));
}