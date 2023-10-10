#include <math.h>

/////////////////////////////////////////////////////////////////////////////
// Levy function for dN/ptdpt
/////////////////////////////////////////////////////////////////////////////
Double_t LevyFcn(Double_t *x, Double_t *par)
{
 Double_t A  = par[0];
 Double_t n  = par[1];
 Double_t T  = par[2];
 Double_t m0 = par[3];
 Double_t mT = sqrt(x[0]*x[0]+m0*m0);

 Double_t a1 = A*(n-1)*(n-2);
 Double_t a2 = n*T*(n*T+m0*(n-2));
 Double_t a3 = pow(1+(mT-m0)/n/T,-n);

 return  a1/a2*a3;       
}

/////////////////////////////////////////////////////////////////////////////
// Levy function for dN/dpt
/////////////////////////////////////////////////////////////////////////////
Double_t LevyFcnPt(Double_t *x, Double_t *par)
{
 Double_t A  = par[0];
 Double_t n  = par[1];
 Double_t T  = par[2];
 Double_t m0 = par[3];
 Double_t mT = sqrt(x[0]*x[0]+m0*m0);

 Double_t a1 = A*(n-1)*(n-2);
 Double_t a2 = n*T*(n*T+m0*(n-2));
 Double_t a3 = pow(1+(mT-m0)/n/T,-n);

 return  a1/a2*a3*x[0];       
}
