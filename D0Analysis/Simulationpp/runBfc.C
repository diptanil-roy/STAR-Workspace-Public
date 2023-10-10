void runBfc( int n = 2, const char* filename = "./rcf22000_13039122_0_10evts.fzd" )
{

  TString chainOpts;//("fzin gen_T geomT sim_T TpcRS -ittf -tpc_daq nodefault y2014x usexgeom misalign sti   DbV20160418 P2014a pxlHit istHit btof mtd mtdCalib BEmcChkStat CorrX OSpaceZ2 OGridLeak3D ODistoSmear -hitfilt  -evout  vfminuit -vfmce tpxclu pxlslowsim istslowsim nosvtit nossdit " ); // picoWrite PicoVtxDefault            DbV20150316 P2014a pxlHit istHit btof mtd mtdCalib BEmcChkStat CorrX OSpaceZ2 OGridLeak3D -hitfilt

  // chainOpts += " DbV20130212 DbV20160506_EMC_Calibrations fzin Sti pp2012b agml mtdDat btof fmsDat VFPPVnoCTB useBTOF4Vtx beamline BEmcChkStat Corr4 OSpaceZ2 OGridLeak3D -hitfilt ";
  // // chainOpts += " DbV20130212 DbV20160506_EMC_Calibrations fzin Sti pp2012b agml mtdDat btof fmsDat -VFMinuit vfmce useBTOF4Vtx BEmcChkStat Corr4 OSpaceZ2 OGridLeak3D -hitfilt ";
  // chainOpts += " tpcdb sdt20120423 ";
  // chainOpts += " bana IAna TPCRS bbcsim btofdat btofsim btofmatch MakeEvent ";
  // chainOpts += " simu dEdxY2 ";                    // Runs TPC fast simulation
  // chainOpts += " -emcDY2 ";
  // chainOpts += " sti ittf ";                      // Runs track finding and reconstruction using the "sti" tracker
  // chainOpts += " gen_T,geomT,sim_T,AgML, emcY2 eess "; // Remove this later.
  // chainOpts += " GeantOut,MiniMcMk,-in,useInTracker, EEss, cmudst";
  //  chainOpts += ",picowrite,picovtxvdefault";

  // chainOpts += " fzin P2012a -emcDY2 emcSim EEss btofSim btofMatch bbcSim sti VFPPVnoCTB ittf cmudst -evout DbV20160418 -hitfilt tpxclu tpcdb tpcrs simu dEdxY2 -in idtruth GeantOut MiniMcMk ";
  // chainOpts += " fzin pp2012b -emcDY2 emcSim EEss btofSim btofMatch bbcSim sti VFPPVnoCTB beamline useBTOF4Vtx ittf cmudst -evout DbV20130212 -hitfilt tpxclu tpcdb tpcrs simu dEdxY2 -in idtruth GeantOut MiniMcMk ";

  chainOpts += "fzin DbV20130212 sdt20120304 PP2012b AgML mtdDat btof fmsDat VFPPVnoCTB useBTOF4Vtx beamline BEmcChkStat Corr4 OSpaceZ2 OGridLeak3D -hitfilt ";
  chainOpts += "btofSim bbcSim btofMatch btofutil db ";
  chainOpts += "emcSim EEss sti cmudst -evout -hitfilt tpcrs tpxclu tpcdb simu dEdxY2 -in idtruth GeantOut MiniMcMk ";

  gSystem->Load("StarRoot.so");
  gROOT->LoadMacro("bfc.C");

  bfc(2, 3,chainOpts,filename);

}
