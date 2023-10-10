void CrossExamMacro(){
    gSystem->Load("./RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(".L UnfoldTheSimFoldedDistributionForCrossExam.C++");
    gROOT->ProcessLine(Form(".x UnfoldTheSimFoldedDistributionForCrossExam.C+"));
}