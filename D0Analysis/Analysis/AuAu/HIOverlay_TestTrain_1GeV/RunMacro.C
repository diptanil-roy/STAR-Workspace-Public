void RunMacro(TString MacroName = "Compare2DHistogramsForClosure.C"){
    gSystem->Load("/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so");
    gROOT->ProcessLine(Form(".L %s++", MacroName.Data()));
    gROOT->ProcessLine(Form(".x %s+()", MacroName.Data()));
}