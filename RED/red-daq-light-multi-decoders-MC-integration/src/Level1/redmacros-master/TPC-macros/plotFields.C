#include "larproperties.C"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "Math/SpecFunc.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"

const Int_t model = -1; //constant

Double_t GetLeakField(Double_t Edrift /* kV/cm*/, Double_t par)
{
  if (model < 0) //no correction
    return 0;
  //This is for a constant leak field
  else if (model == 0)
    return par/1000.;
  //This is for a Eleak = k/Edrift
  else if (model == 1)
    return par/(Edrift*1000.)/1000.;
  else if (model == 2) //Eleak = k(Eext-Edrift)
    return par*(3.8-Edrift);
}

void plotFields(){
  TGraphErrors* gr = new TGraphErrors("driftfield.dat");
  gr->SetMarkerStyle(20);

  TF1 *f1 = new TF1("f1", SpeedFunc, 0.05,1.5,1); 
  f1->SetParameter(0,88.);
  f1->SetLineColor(kBlue);

  TF1 *f2 = new TF1("f2", SpeedFunc, 0.05,1.5,1); 
  f2->SetParameter(0,85.);
  f2->SetLineColor(kGreen);
  
  TF1 *f3 = new TF1("f3", SpeedFunc, 0.05,1.5,1); 
  f3->SetParameter(0,91.);
  f3->SetLineColor(kRed);
  
  f1->GetXaxis()->SetTitle("Electric field (kV/cm)");
  f1->GetYaxis()->SetTitle("Electron drift velocity (mm/#mus)");


  TLegend* l1 = new TLegend(0.1,0.7,0.48,0.9);
  l1->AddEntry(f2,"85 K");
  l1->AddEntry(f1,"88 K");
  l1->AddEntry(f3,"91 K");

  f1->Draw();
  f2->Draw("same");
  f3->Draw("same");
  gr->Draw("ZP");
  l1->Draw("same");

  //Function to minimize
  auto chi2Function = [&](const Double_t *par)
    {
      //minimisation function computing the sum of squares of residuals
      // looping at the data points
      Double_t f=0;
      //Loop on data
      for (Int_t i=0;i<gr->GetN();i++) {
	Double_t x=0;
	Double_t y=0;
	gr->GetPoint(i,x,y);
	Double_t dy = gr->GetErrorY(i);

	Double_t offset = par[0];
	Double_t temperature = par[1];
	
	Double_t theory = LArProperties::GetLArDriftSpeed(temperature,
							  x+GetLeakField(x,offset));
	
	//cout << x << " " << offset << " " << theory << " " << y << endl;
	
	f += (y-theory)*(y-theory)/(dy*dy); //minimization!
      }   
      return f;
    };
  // wrap chi2 funciton in a function object for the fit
  // 2 is the number of fit parameters (size of array par)
  ROOT::Math::Functor fcn(chi2Function,2);
  ROOT::Fit::Fitter  fitter;
  double pStart[2] = {0.,88.};
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("LF"); //leak field
  fitter.Config().ParSettings(1).SetName("T"); //temperature
  if (model == 0)
    fitter.Config().ParSettings(0).SetLimits(-50.,50.);
  else if (model == 1)
    fitter.Config().ParSettings(0).SetLimits(-50.e3,50.e3); // (V/cm)^2
  else if (model == 2)
    fitter.Config().ParSettings(0).SetLimits(-0.05,0.05); // (V/cm)^2

  TGraph* chi2 = new TGraph();
  TGraph* field = new TGraph();
  Int_t counter = 0;
  ROOT::Fit::FitResult aResult; //best fit
  Double_t f0=DBL_MAX;
  for (Double_t T=84;T<94;T+=0.01)
    {
      fitter.Config().ParSettings(1).Release();
      fitter.Config().ParSettings(1).SetValue(T);
      fitter.Config().ParSettings(1).Fix();
      bool ok = fitter.FitFCN();
      if (!ok) {
	Error("line3Dfit","Line3D Fit failed");
      }   
      const ROOT::Fit::FitResult & result = fitter.Result();
      //result.Print(std::cout);
      
      Double_t LF = result.Parameter(0);
      Double_t Te = result.Parameter(1);
      Double_t f = chi2Function(result.GetParams());
      if (f<f0)
	{
	  f0=f;
	  aResult = result;
	}

      chi2->SetPoint(counter,T,f);
      field->SetPoint(counter,T,LF);
      
      counter++;
      //cout << LF << " " << Te << " " << f << endl;
    }
  cout << "Model for E_leak is " << ((model == 0) ? " constant field " : "k/Edrift") 
       << endl;
    
  cout << "Best fit is: " << aResult.Parameter(0) << " +/- " << aResult.Error(0) << ", T = " << 
    aResult.Parameter(1) << " " << 
    chi2Function(aResult.GetParams()) << endl;
  Double_t chi2min = chi2Function(aResult.GetParams());

  for (Int_t i=0;i<chi2->GetN()-1;i++)
    {
      Double_t x=0;
      Double_t y=0;
      Double_t x1=0;
      Double_t y1=0;
      chi2->GetPoint(i,x,y);
      chi2->GetPoint(i+1,x1,y1);
      if (y > (chi2min+1) && y1 < (chi2min+1))
	cout << "Left edge: " << (x+x1)/2. << endl;
      if (y < (chi2min+1) && y1 > (chi2min+1))
	cout << "Right edge: " << (x+x1)/2. << endl;
    }

  TCanvas* c2 = new TCanvas();
  c2->Divide(1,2);
  c2->cd(1);
  chi2->Draw("AZP");
  c2->cd(2);
  field->Draw("AZP");

  TGraphErrors* gr2 = new TGraphErrors();
  gr2->SetMarkerStyle(20);

  for (Int_t i=0;i<gr->GetN();i++) {
    Double_t x=0;
    Double_t y=0;
    gr->GetPoint(i,x,y);
    Double_t dy = gr->GetErrorY(i);
    
    Double_t offset =  aResult.Parameter(0);
    Double_t temperature = aResult.Parameter(1);
    
    Double_t theory = LArProperties::GetLArDriftSpeed(temperature,
						      x+GetLeakField(x,offset));
    
	//cout << x << " " << offset << " " << theory << " " << y << endl;
    gr2->SetPoint(i,x+GetLeakField(x,offset),y);
    cout << x << " " << GetLeakField(x,offset) << " " << theory << " " << y << endl;
    gr2->SetPointError(i,0.,dy);
  
  }   
  TLatex l;
  TCanvas* c3 = new TCanvas();
  c3->cd();
  TF1 *f4 = new TF1("f4", SpeedFunc, 0.05,1.5,1); 
  f4->SetParameter(0,aResult.Parameter(1));
  f4->SetLineColor(kRed);
  f4->GetXaxis()->SetTitle("Electric field (kV/cm)");
  f4->GetYaxis()->SetTitle("Electron drift velocity (mm/#mus)");
  f4->Draw();
  gr2->Draw("ZP");
  if (model == 0)
    l.DrawLatex(0.8,1,Form("Leak field: %5.1lf V/cm", aResult.Parameter(0)));
  else if (model == 1)
    l.DrawLatex(0.8,1,Form("k= %5.2e (V/cm)^2", aResult.Parameter(0)));
  else if (model == 2)
    l.DrawLatex(0.8,1,Form("k= %5.3lf ", aResult.Parameter(0)));
  l.DrawLatex(0.8,0.8,Form("Temperature: %5.1lf K", aResult.Parameter(1)));

  
  //Test full fit
  fitter.Config().ParSettings(1).Release();
  fitter.Config().ParSettings(1).SetLimits(80,100);
  //cout << fitter.Config().ParSettings(1).IsFixed() << endl;
  bool ok = fitter.FitFCN();
  if (!ok) {
    Error("line3Dfit","Line3D Fit failed");
  }   
  const ROOT::Fit::FitResult & resultF = fitter.Result();
  //resultF.Print(std::cout);
  
  Double_t LF = resultF.Parameter(0);
  Double_t Te = resultF.Parameter(1);
  Double_t f = chi2Function(resultF.GetParams());
   
  cout << "Full fit: (" << LF << " +/- " << resultF.Error(0) << "), (" <<
    Te << " +/- " << 
    resultF.Error(1) << "), chi2=" << f << endl;

  return;
}
