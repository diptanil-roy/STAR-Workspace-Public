void run(){
	gSystem->Load("libEG");
  	gSystem->Load("libEGPythia8");

  	gROOT->ProcessLine(".L PythiaCrossSections.C++");
  	gROOT->ProcessLine(".x PythiaCrossSections.C+");
}