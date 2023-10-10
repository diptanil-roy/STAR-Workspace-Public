void RunMacro(Int_t nEvents = 10, const int random_seed = 99999, const char* VertexFileName = "\"./VERTEXFULL/Vertex_st_physics_15094070_raw_1000002.txt\""){
    // Load Pythia 8
    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");
    // // Load fastjet libraries and includes
    // gSystem->AddIncludePath( "-I$FASTJET/include");
    // gSystem->Load("$FASTJET/lib/libfastjet");

    // Compile Macro
    gROOT->ProcessLine(".L PYTHIA8ForD0Decay.C++");
    gROOT->ProcessLine(Form(".x PYTHIA8ForD0Decay.C+(%i, %i, %s)", nEvents, random_seed, VertexFileName));
}