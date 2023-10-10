#include "../anaCuts.h"
#include "../myConst.h"
void write_CrossS_Rcp_err_WRONG()  // this one calculate the correlated sys for Rcp vs Nbin
{
   char CMD[250];
   char dir[250];
   sprintf(dir, "data");
   sprintf(CMD, "[ -d %s ] || mkdir -p %s", dir, dir);
   // gSystem->Exec(CMD);

   float ptmean[npt], pterr[npt];
   for (int ipt = 0; ipt < npt; ipt++)
   {
      ptmean[ipt] = 0.5 * (nptbin[ipt] + nptbin[ipt + 1]);
      pterr[ipt] = 0.5 * (-nptbin[ipt] + nptbin[ipt + 1]);
   }

   //caculate err
   float yCS1[ncent], yCS1sys[ncent];
   float yCS2[ncent], yCS2sys[ncent];
   float tmp[ncent];
   float tmp1[ncent], tmp11[ncent];
   float tmp2[ncent], tmp22[ncent];
   ifstream in;
   ofstream out;

   float unit[npt];
   for (int ipt = 0; ipt < npt; ipt++)  unit[ipt] = 1.;

   //for Rcp vs Nbin
   for (int icent = 0; icent < ncent; icent++)
   {
      //init

      //statistics error
      in.open(Form("../default/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No default yield error file!!!" << endl;
         exit(1);
      }
      in >> yCS1[icent] >> yCS2[icent];
      in.close();


      //sys 3, Daughter pt Cut scan // choose the maximum difference
      //sys 3.1 -- daughter pt cut1
      in.open(Form("../ptCut1/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No pt = 0.3 sys error file!!!" << endl;
         exit(1);
      }
      in >> tmp1[icent] >> tmp2[icent];
      in.close();

      //sys 3.2 -- daughter pt cut2
      in.open(Form("../ptCut2/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No pt = 0.5 sys error file!!!" << endl;
         exit(1);
      }
      in >> tmp11[icent] >> tmp22[icent];
      in.close();
   }

   for (int icent = 0; icent < ncent; icent++)
   {
      tmp[icent] = fabs(tmp1[icent]/tmp1[4] - yCS1[icent]/yCS1[4]) > fabs(tmp11[icent]/tmp11[4] - yCS1[icent]/yCS1[4]) ? tmp1[icent] : tmp11[icent];
      cout << "tmp["<< icent << "]=" << tmp[icent]  << endl;
      yCS1sys[icent] = ((tmp[icent]/tmp[4]) *1. / (yCS1[icent]/yCS1[4])) - 1.;
      cout << "yCS1["<< icent << "]=" << yCS1[icent]  << endl;
      double tttt = tmp[icent]/1./tmp[4];
      cout << tttt <<endl;
      // cout << yCS1sys[icent] << endl;

      tmp[icent] = fabs(tmp2[icent]/tmp2[4] - yCS2[icent]/yCS2[4]) > fabs(tmp22[icent]/tmp22[4] - yCS2[icent]/yCS2[4]) ? tmp2[icent] : tmp22[icent];
      yCS2sys[icent] = ((tmp[icent]/tmp[4]) / (yCS2[icent]/yCS2[4])) - 1.;
      cout << yCS2sys[icent] << endl;
   }

   for (int icent = 0; icent < ncent; icent++)
   {
      //sys 4 topological cut
      //4.1 -- tight topo cuts
      in.open(Form("../topoCut1/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No tight topo sys error file!!!" << endl;
         exit(1);
      }
      in >> tmp1[icent] >> tmp2[icent];
      in.close();

      //4.2 -- loose topo cuts
      in.open(Form("../topoCut2/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No loose topo sys error file!!!" << endl;
         exit(1);
      }
      in >> tmp11[icent] >> tmp22[icent];
      in.close();
   }

   for (int icent = 0; icent < ncent; icent++)
   {
      tmp[icent] = ((fabs((tmp1[icent]/tmp1[4]) / (yCS1[icent]/yCS1[4])) - 1.) + (fabs((tmp11[icent]/tmp11[4]) / (yCS1[icent]/yCS1[4])) - 1.)) * 0.5;
      yCS1sys[icent] = sqrt(pow(tmp[icent], 2) + pow(yCS1sys[icent], 2));

      tmp[icent] = ((fabs((tmp2[icent]/tmp2[4]) / (yCS2[icent]/yCS2[4])) - 1.) + (fabs((tmp22[icent]/tmp22[4]) / (yCS2[icent]/yCS2[4])) - 1.)) * 0.5;
      yCS2sys[icent] = sqrt(pow(tmp[icent], 2) + pow(yCS2sys[icent], 2));

   }

   out.open(Form("data/CrossX_Rcp1_sys_Summary.txt"));
   out << "sys1(pt>0)" << "\t\t" << "sys2.(pt>4) relative" << endl;
   for (int icent = 0; icent < ncent; icent++)
   {
      out <<nameCent1[icent] << "\t" << yCS1sys[icent] << "\t" << yCS2sys[icent] << endl;
   }
   out.close();

}
