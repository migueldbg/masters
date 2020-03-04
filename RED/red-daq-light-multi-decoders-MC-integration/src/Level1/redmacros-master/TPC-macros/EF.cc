#include <TMath.h>
void GetFields(double s /*cm*/, double d /*cm*/, double E1 /*V/cm*/, double E2 /*V/cm*/, bool Crossed=true )
{
  double PI=TMath::Pi();

  double Cm,Cd;
  // Crossed wires
  if ( Crossed )
    {
      Cm=log(s/PI/d)/4./PI +  0.07*d/s  +  0.4*pow(d/s,2.);
      Cd=0.25*(d/s)-0.5*pow(d/s,2.);
    }
  // Paralell wires
  else
    {
      Cm=log(s/PI/d)/2./PI +  0.25*d/s;
      Cd=0.25*(d/s)-0.36*pow(d/s,2.);
    }

  double Phi1=(E1-E2)*s*Cm+(E1+E2)*s*Cd;
  double Phi2=(E1-E2)*s*Cm-(E1+E2)*s*Cd;
  printf("Phi=%5.2lf %5.2lf |  d=%5.2lf mm s=%5.2lf mm E1=%5.2lf V/cm E2=%5.2lf V/cm \n",Phi1,Phi2,d*10.,s*10.,E1,E2);
  printf("Cm=%5.4lf Cd=%5.4lf d/s=%5.2lf \n",Cm,Cd,d/s);
}
