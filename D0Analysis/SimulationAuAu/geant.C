// $Id: geant.C,v 1.14 2001/05/30 18:24:59 perev Exp $
//=======================================================================
// owner: Pavel Nevski
// what it does: 
//=======================================================================
//
//#define gtrack
#include "TH1F.h"

TBrowser *b = 0;
class St_geant_Maker;
St_geant_Maker *geant=0;
void geant(const Int_t Nevents=10,
     // const Char_t *fzfile ="TESTSAMPLE/st_physics_15094070_raw_2000044_3_0_cc_MONASH_TXFile_NoD0Decay_15094070_fieldon_misalign_sdt20140216_10000evts.fzd")

      // const Char_t *fzfile ="/gpfs01/star/pwg_tasks/jetcorr03/HFJets_Mar17_2022_SimulatedFiles/out_production/rcf22000_15165031_0014_1_100evts.fzd")
  const Char_t *fzfile ="/gpfs01/star/pwg_tasks/jetcorr03/zweighted_sample/test/out/rcf22000_15147028_001001_1_10evts.fzd")
  // const Char_t *fzfile = "/gpfs01/star/pwg_tasks/jetcorr03/HFJets_AllHardQCD_3_inf/out/rcf22000_15166041_0010_1_200evts.fzd")
{

  // TH1F *hD0Pt =  new TH1F("hD0Pt", "hD0Pt", 20, 0, 10);
  gROOT->LoadMacro("bfc.C");
  bfc(0,"fzin sim_T globT gen_T",fzfile);
  Int_t i=0;
  for (Int_t i =1; i <= 10; i++){
    chain->Clear();
    cout << "============================================================" << endl;
    cout << "Event # " << i << endl;
    if (chain->Make(i)>=kStEOF){cout << "End of File \n"; break;}
    St_g2t_track *track = (St_g2t_track *) chain->FindObject("g2t_track");
     g2t_track_st *trk = track->GetTable();
     for (Int_t j = 0; j < track->GetNRows(); j++,trk++) {

      // printf("PID: %d ", trk->ge_pid);
       // if (TMath::Abs(trk->charge) < 0.5) continue;
       // if (! trk->eg_label) continue;
       // if (trk->n_tpc_hit < 26) continue;
       // hyp = -1;
       // if (trk->ge_pid == 2 || trk->ge_pid == 3) hyp = 3;
       // if (trk->ge_pid == 8 || trk->ge_pid == 9) hyp = 2;
       // if (trk->ge_pid ==11 || trk->ge_pid ==12) hyp = 1;
       // if (trk->ge_pid ==14 || trk->ge_pid ==15) hyp = 0;

       Double_t pT =  trk->pt;
       Double_t Eta = trk->eta;
       Double_t Y = trk->rapidity;

       // if (trk->ge_pid == 37 || trk->ge_pid == 38){
        printf("eg: %d ge: %d  eg: %d start_vertex: %d stop_vertex: %d pT %f eta %f  Nvpd %d  NBEMC %d\n",
       trk->eg_label,trk->ge_pid,trk->eg_pid,trk->start_vertex_p,trk->stop_vertex_p,pT,Eta,trk->n_vpd_hit, trk->n_emc_hit);
         // if (abs(Y) < 1)hD0Pt->Fill(pT, 1./pT);
       // }
     }
    printf ("=========================================== Done with Event no. %d\n",i);
  }
  // TCanvas *c1 = new TCanvas("c1", "c1", 5,5, 600, 600);
  // c1->cd();
  // hD0Eta->Draw();
  // c1->SaveAs("D0Eta.pdf");

  // TFile *f = new TFile("D0Pt.root", "RECREATE");
  //     f->cd();
  //     hD0Pt->Write();
  //     f->Close();
}


/*
 input >> track.eg_label;
 input >> track.eg_pid;
 input >> track.ge_pid;
 input >> track.hit_ctb_p;
 input >> track.hit_eem_p;
 input >> track.hit_emc_p;
 input >> track.hit_esm_p;
 input >> track.hit_ftp_p;
 input >> track.hit_mwc_p;
 input >> track.hit_pgc_p;
 input >> track.hit_psc_p;
 input >> track.hit_smd_p;
 input >> track.hit_svt_p;
 input >> track.hit_tof_p;
 input >> track.hit_tpc_p;
 input >> track.hit_vpd_p;
 input >> track.id;
 input >> track.is_shower;
 input >> track.itrmd_vertex_p;
 input >> track.n_ctb_hit;
 input >> track.n_eem_hit;
 input >> track.n_emc_hit;
 input >> track.n_esm_hit;
 input >> track.n_ftp_hit;
 input >> track.n_mwc_hit;
 input >> track.n_pgc_hit;
 input >> track.n_psc_hit;
 input >> track.n_smd_hit;
 input >> track.n_svt_hit;
 input >> track.n_tof_hit;
 input >> track.n_tpc_hit;
 input >> track.n_vpd_hit;
 input >> track.next_parent_p;
 input >> track.next_vtx_trk_p;
 input >> track.start_vertex_p;
 input >> track.stop_vertex_p;
 input >> track.charge;
 input >> track.e;
 input >> track.eta;
 input >> track.p[0] >> track.p[1] >> track.p[2];
 input >> track.pt;
 input >> track.ptot;
 input >> track.rapidity;
*/