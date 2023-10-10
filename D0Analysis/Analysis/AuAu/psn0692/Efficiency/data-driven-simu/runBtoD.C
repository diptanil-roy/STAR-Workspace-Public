void runBtoD(int npart = 10, TString output = "D0.toyMc.root", TString particleName = "D0", bool isCombinB = false)
{
  gSystem->Load("toyMcBtoD_C");
  toyMcBtoD(npart,output,particleName,isCombinB);
}
