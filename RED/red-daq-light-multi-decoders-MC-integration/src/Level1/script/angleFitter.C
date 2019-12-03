#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm> 
#include <vector> 
#include <sstream> 
#include <fstream>
#include <map>

#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLegend.h>
#include <TCut.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TMath.h>

#include "TString.h"
#include "Math/SpecFunc.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"

using namespace std;

double CalculateEnergies(Double_t beamEne, Double_t theta, TString peakcode);

struct fitresult{
  Double_t chi2;
  Double_t a;
  Double_t b;
  Double_t c;
  Double_t theta;
};

void angleFitter(TString filename="par_fit.out")
{    
  ifstream infile(filename.Data());
  
  if (!infile.good()) {
    
    cout << "file no buono" << endl;
    return;    
  }
    
  //Read data: beam energy, E, DeltaE (w/ errors), peak code    
  vector<Double_t> beamEne;
  vector<Double_t> DeltaE;
  vector<Double_t> E;
  vector<Double_t> err1;
  vector<Double_t> err2;
  vector<Double_t> rho;
  vector<TString> peakcode;
    
  double deltaE,En,Err1,Err2,Rho,enebeam;
  TString Peakcode;
  double deg=180./TMath::Pi();
  double rad=1./deg;
  double theta=5.*rad;
  double thetamin=0;
  double minf = 1e20;
  fitresult bestResult;

  while(infile >> enebeam >> En >> Err1 >> deltaE >> Err2 >> Rho >> Peakcode) 
    {
      beamEne.push_back(enebeam-0.13);
      DeltaE.push_back(deltaE);
      E.push_back(En);
      err1.push_back(Err1);
      err2.push_back(Err2);
      rho.push_back(Rho);
      peakcode.push_back(Peakcode);    
    }
  
  infile.close();
    
  for(size_t i=0;i<E.size();i++){
    
    cout << "" << beamEne[i] << " " << E[i] << " " << err1[i] << " " << DeltaE[i] << " " << err2[i]
	 << " " << rho[i] << " " << peakcode[i] << endl;
  }
  //cout << peakcode.size() << endl;
    std::ofstream outfile("data.out");
  for (theta=4;theta<5.8;theta+=0.005)
    {
        Double_t thetarad = theta*rad;
        std::pair<Double_t,TString> key;
        std::map< std::pair<Double_t,TString> ,Double_t> refEnergies;
    
  for (size_t i=0;i<peakcode.size();i++)
    {
      
      Double_t Erec =  CalculateEnergies(beamEne[i], thetarad, peakcode[i]);
      //cout << "Sono il picco: " << peakcode[i] << " e la mia energia e " << Erec << " MeV" << endl;
      key = std::make_pair(beamEne[i],peakcode[i]);
      refEnergies.insert(std::pair< std::pair<Double_t,TString> ,Double_t>(key,Erec) );
    }
  
  //Function to minimize
  auto chi2Function = [&](const Double_t *par)
    {
      //minimisation function computing the sum of squares of residuals
      // looping at the data points
      Double_t f=0;
      for (size_t i=0;i<DeltaE.size();i++) {
	key = std::make_pair(beamEne[i],peakcode[i]);
	Double_t expectedEnergy = refEnergies.find(key)->second;
	Double_t p1 = (par[0] + par[1]*DeltaE[i]+par[2]*E[i] -
		       expectedEnergy);
    Double_t errorsquare = ((err1[i]/DeltaE[i])*(err1[i]/DeltaE[i])+(err2[i]/E[i])*(err2[i]/E[i]))
          *expectedEnergy*expectedEnergy;
	f += p1*p1/errorsquare; //minimization!
      }   
      return f;
    };
  // wrap chi2 funciton in a function object for the fit
  // 3 is the number of fit parameters (size of array par)
  ROOT::Math::Functor fcn(chi2Function,3);
  ROOT::Fit::Fitter  fitter;
  double pStart[3] = {0,1,1};
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("a");
  fitter.Config().ParSettings(1).SetName("b");
  fitter.Config().ParSettings(2).SetName("c");
  
  //fitter.Config().ParSettings(2).SetValue(0.00243);
  //fitter.Config().ParSettings(2).Fix();
  // do the fit
  bool ok = fitter.FitFCN();
  if (!ok) {
    Error("line3Dfit","Line3D Fit failed");
  }   
  const ROOT::Fit::FitResult & result = fitter.Result();
  //result.Print(std::cout);
  
  Double_t a = result.Parameter(0);
  Double_t b = result.Parameter(1);
  Double_t c = result.Parameter(2);

  for (size_t i=0;i<DeltaE.size();i++) {
    key = std::make_pair(beamEne[i],peakcode[i]);
    Double_t expectedEnergy = refEnergies.find(key)->second;
    Double_t p1 = a+b*DeltaE[i]+c*E[i];
    //cout << "Peak: " << peakcode[i] << ", calibrated energy: " << p1 << ", tabular energy:"  << expectedEnergy << endl;
  }
  fitresult aResult;
  Double_t f = chi2Function(result.GetParams());
  aResult.chi2 = 0;
  aResult.a = a;
  aResult.b = b;
  aResult.c = c;
  aResult.theta = theta;
  cout << "Value of f: " << theta << " " << f <<  endl;
        outfile << theta << " " << f << endl;
        if (f < minf){
            minf=f;
            thetamin=theta;
	    bestResult = aResult;
        }
    }
  cout << "Minimum for theta= " << thetamin << endl;
  cout << "a= " << bestResult.a << endl;
  cout << "b= " << bestResult.b << endl;
  cout << "c= " << bestResult.c << endl;
  outfile.close();

  for (size_t i=0;i<DeltaE.size();i++) {
    Double_t expectedEnergy =  CalculateEnergies(beamEne[i], 
						 bestResult.theta*rad, peakcode[i]);
    Double_t p1 = bestResult.a+bestResult.b*DeltaE[i]+bestResult.c*E[i];
    cout << "Peak: " << peakcode[i] << " at " << beamEne[i] << " MeV: calibrated energy: " << p1 << ", tabular energy:"  << expectedEnergy << endl;
  }

  
  //return 1;
}

//Calculate the expected energy of each peak, given
// beam emergy, theta

double CalculateEnergies(Double_t beamEne, Double_t theta, TString peakcode)
{  
  double xmli=7.016004;
  double xm3,xm4,q;
  double Etot,P1,P12,cos3,sin3,cos32,radn,den,P3,P31,E3,E31,E4,E41;
   
  //double xmli=7.01600;              //7Li mass
  //double xmbe=7.016928*u;           //7Be mass
  //double xmneu= 1.008664916*u;      //neutron mass
  //double xmprot=1.007825*u;         //proton mass
  //double xmli=7.01600*u;            //7Li mass
  //double q=-1.644;                  //q-value p(7Li,7Be)n
  
  if (peakcode=="hibe" || peakcode=="lowbe") 
    {
      q=-1.644;           //q-value p(7Li,7Be)n
      xm4=1.008644;       //neutron mass
      xm3=7.016929;       //beryllium mass
      Etot=beamEne+q;
    }
  else if (peakcode=="lip1" || peakcode=="lip2") 
    {   
      q=0.;               //q-value p(7Li,7Li)p
      xm4=1.007825;       //proton mass
      xm3=xmli;           //lithium mass
      
      Etot=beamEne+q;
    }
  else if (peakcode=="liC") 
    {
      q=0.;               //q-value 12C(7Li,7Li)12C
      xm3=xmli;           //lithium mass
      xm4=12.00;          //carbon mass
      
      Etot=beamEne+q;
    }
  else
    {
      std::cout << "Unknown code: " << peakcode << std::endl;
      return 0;
    }

    //Etot=beamEne+q;                              //total energy of outgoing particles

    //cout << "Etot: " << Etot << endl;

    P1=sqrt(2.*xmli*beamEne);                     //projectile impulse
    P12=P1*P1;

    cos3=cos(theta);
    sin3=sin(theta);
    cos32=cos3*cos3;
    radn=(4.*P12*cos32)-(4.*(1.+(xm4/xm3)))*(P12-(2.*xm4*Etot)); //xm3 ejectile mass, xm4 recoil mass
    if (radn<0) {
      cout << "Negative sqrt value!!" << endl;
      return 0;
    }
    den=2.*(1.+(xm4/xm3));
    
    P3=((2.*P1*cos3)+sqrt(radn))/den;          //ejectile impulse
    P31=((2.*P1*cos3)-sqrt(radn))/den;
    E3=(P3*P3)/(2.*xm3);                          //ejectile energy
    
    /* cout << "E3: " << E3 << endl; */
    
    E31=(P31*P31)/(2.*xm3);
    //cout << "E3: " << E3 << " " << "E31: " << E31 << endl;

    E4=Etot-E3;                                  //recoil energy
    E41=Etot-E31;
    //cout << "E4: " << E4 << " " << "E41: " << E41 << endl;
    
    if (peakcode=="lowbe" || peakcode=="lip2")
      return std::max(E3,E31);
    else if (peakcode == "hibe" || peakcode=="lip1")
      return std::min(E3,E31);
    else if (peakcode == "liC")
      return std::max(E3,E31);
    else 
      return 0;    
}
