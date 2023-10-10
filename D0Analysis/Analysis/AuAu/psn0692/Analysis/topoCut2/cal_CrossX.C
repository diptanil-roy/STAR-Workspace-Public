#include "../anaCuts.h"
void cal_CrossX()  //re calculate the cross section for Rcp vs Nbin
{
   ifstream in;
   ofstream out;

   float y1[npt], y1err[npt];//
   float y[npt], yerr[npt];//
   
   float XC1[ncent];
   float XC2[ncent];
   for (int icent = 0; icent < ncent; icent++)
   {
     XC1[icent] = 0;
     XC2[icent] = 0;
   }

   float pt_mean[npt], pt_err[npt], pt_width[npt];
   for (int ipt = 0; ipt < npt; ipt++)
   {
     pt_mean[ipt] = 0.5*(nptbin[ipt] + nptbin[ipt+1]);
     pt_err[ipt] = 0.5*(-nptbin[ipt] + nptbin[ipt+1]);
     pt_width[ipt] = nptbin[ipt+1] - nptbin[ipt];
   }

   for (int icent = 0; icent < ncent; icent++)
   {
      in.open(Form("../topoCut2/data/yield_%s.txt", nameCent1[icent]));
      for (int ipt = 0; ipt < npt; ipt++)
         in >> y1[ipt] >> y1err[ipt];
      in.close();

      //re sign some specific data point
      for (int ipt = 0; ipt < npt; ipt++)
      {
      if(pt_mean[ipt]> 0.)  XC1[icent] += y1[ipt];
      if(pt_mean[ipt]> 4.)  XC2[icent] += y1[ipt];
      }

      out.open(Form("../topoCut2/data/crossX_%s.txt", nameCent1[icent]));
         out << XC1[icent] << "\t" << XC2[icent] << endl;
      out.close();

   } 

}
