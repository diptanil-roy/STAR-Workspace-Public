void MuMacro(){

    gROOT->ProcessLine(".L MuHijingMaker.C");
    gROOT->ProcessLine("Test t;");
    gROOT->ProcessLine("t.Loop()");
}

