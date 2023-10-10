#include "../anaCuts.h"
#include "../myConst.h"
void write_CrossS_err()  // this one calculate the correlated sys for Rcp vs Nbin
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

   // yield
   // for yield extraction and double count and Vtxsys
   for (int icent = 0; icent < ncent; icent++)
   {
      //sys 2
      //
      //-- 2.0 count with side band
      in.open(Form("../default/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No cout sys error file!!!" << endl;
         exit(1);
      }
      // in >> yCS1[icent] >> yCS2[icent];
      in.close();

      //-- 2.1 count with side band
      in.open(Form("../count/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No cout sys error file!!!" << endl;
         exit(1);
      }
      // in >> yCS1[icent] >> yCS2[icent];
      in.close();

      //-- 2.2 count with fit change fit range
      in.open(Form("../fitRange/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No change fit range sys error file!!!" << endl;
         exit(1);
      }
      // in >> yCS1[icent] >> yCS2[icent];
      in.close();

      //-- 2.3 count with likesign bkg subtraction
      in.open(Form("../likeSign/data/crossX_%s.txt", nameCent1[icent]));
      if (in.eof())
      {
         cout << "No like-sign sys error file!!!" << endl;
         exit(1);
      }
      // in >> yCS1[icent] >> yCS2[icent];
      in.close();

   }


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

      tmp[icent] = fabs(tmp1[icent] - yCS1[icent]) > fabs(tmp11[icent] - yCS1[icent]) ? tmp1[icent] : tmp11[icent];
      yCS1sys[icent] = sqrt(pow(tmp[icent] / yCS1[icent] - 1., 2));

      tmp[icent] = fabs(tmp2[icent] - yCS2[icent]) > fabs(tmp22[icent] - yCS2[icent]) ? tmp2[icent] : tmp22[icent];
      yCS2sys[icent] = sqrt(pow(tmp[icent] / yCS2[icent] - 1., 2));

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

      tmp[icent] = (fabs(tmp1[icent] / yCS1[icent] - 1.) + fabs(tmp11[icent] / yCS1[icent] - 1.)) * 0.5;
      yCS1sys[icent] = sqrt(pow(tmp[icent], 2) + pow(yCS1sys[icent], 2));

      tmp[icent] = (fabs(tmp2[icent] / yCS2[icent] - 1.) + fabs(tmp22[icent] / yCS2[icent] - 1.)) * 0.5;
      yCS2sys[icent] = sqrt(pow(tmp[icent], 2) + pow(yCS2sys[icent], 2));


      //out put
      out.open(Form("data/CrossX_Rcp_sys_%s.txt", nameCent1[icent]));
      out << "sys1(pt>0)" << "\t\t" << "sys2.(pt>4) relative" << endl;
      out << yCS1sys[icent] << "\t" << yCS2sys[icent] << endl;
      out.close();

   }

   out.open(Form("data/CrossX_sys_Summary.txt"));
   out << "sys1(pt>0)" << "\t\t" << "sys2.(pt>4) relative" << endl;
   for (int icent = 0; icent < ncent; icent++)
   {
      out <<nameCent1[icent] << "\t" << yCS1sys[icent] << "\t" << yCS2sys[icent] << endl;
   }
   out.close();

}
