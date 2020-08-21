#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>

int plotVBias()
{
  gStyle->SetOptFit();
  vector<Int_t> runNumbers = {159, 161, 165, 171, 174, 175};
  vector<Double_t> Vtop = {32.,33.,34.,35.,36.,37.};
  vector<Double_t> Vbottom = {64.,66.,68.,70.,72.,74.};
  Int_t nlines=16;

  vector<TGraph*> par0(2*nlines,nullptr);
  vector<TGraph*> par1(2*nlines,nullptr);
  
  for (Int_t i=0;i<2*nlines;i++)
    {
      par0.at(i) = new TGraph();     
      par1.at(i) = new TGraph();
    }



  for (size_t ifile=0;ifile < runNumbers.size(); ifile++)
    {
      TString fname;
      fname.Form("ser_%d_b00.out",runNumbers[ifile]);
      cout << "Open file: " << fname << " vbias = " << Vtop[ifile] << endl;
      ifstream file(fname.Data());

      if (!file.good())
	{
	  cout << "file no buono" << endl;
	  return 0;
	}

      for (int iline=0;iline<nlines;iline++)
	{
	  Double_t aa,bb;
	  Int_t ichan;
	  file >> ichan >> aa >> bb;
	  if (bb > 0)
	    {
	      par0.at(ichan)->SetPoint(par0.at(ichan)->GetN(),Vtop.at(ifile),aa);
	      par1.at(ichan)->SetPoint(par1.at(ichan)->GetN(),Vtop.at(ifile),bb);
	    }
	  //cout << ichan << " " << aa << " " << bb << endl;
	}
      file.close();
      
      fname.Form("ser_%d_b01.out",runNumbers[ifile]);
      cout << "Open file: " << fname << " vbias = " << Vbottom[ifile] << endl;
      ifstream file2(fname.Data());

      if (!file2.good())
	{
	  cout << "file no buono" << endl;
	  return 0;
	}

      for (int iline=0;iline<nlines;iline++)
	{
	  Double_t aa,bb;
	  Int_t ichan;
	  file2 >> ichan >> aa >> bb;
	  if (bb > 0)
	    {
	      Int_t index = ichan+16;
	      if (index < 24)
		{
		  par0.at(index)->SetPoint(par0.at(index)->GetN(),Vtop.at(ifile),aa);
		  par1.at(index)->SetPoint(par1.at(index)->GetN(),Vtop.at(ifile),bb);
		}
	      else if (index < 28)
		{
		  if (!((index==24) && (ifile==5)))
		    {
		      par0.at(index)->SetPoint(par0.at(index)->GetN(),Vbottom.at(ifile),aa);
		      par1.at(index)->SetPoint(par1.at(index)->GetN(),Vbottom.at(ifile),bb);
		    }
		}
	    }
	  //cout << ichan << " " << aa << " " << bb << endl;
	}
      file2.close();


    }

  vector<TCanvas*> canvases;
  for (size_t i=0;i<7;i++)
    {
      canvases.push_back(new TCanvas());      
      canvases[i]->Divide(2,2);
    }
 
  TF1* fun = new TF1("fun","pol1");
  for (size_t i=0;i<28;i++)
     {
       Int_t canvasid=i/4;
       Int_t padId = i%4;
       canvases.at(canvasid)->cd(padId+1);
       par0.at(i)->GetXaxis()->SetTitle("V_{bias} (V)");
       par0.at(i)->GetYaxis()->SetTitle("Incercept");
       par1.at(i)->GetXaxis()->SetTitle("V_{bias} (V)");
       par1.at(i)->GetYaxis()->SetTitle("Slope");
       par1.at(i)->SetMarkerStyle(22);
       par1.at(i)->SetMarkerColor(kBlue);
       par0.at(i)->SetMarkerStyle(22);
       par0.at(i)->SetMarkerColor(kBlue);
       if (i<24)
	 par1.at(i)->SetTitle(Form("Chan %d (top)",i));
       else 
	 par1.at(i)->SetTitle(Form("Chan %d (bottom)",i));
       par1.at(i)->Fit(fun,"Q");
       par1.at(i)->Draw("AP");
       //cout << "Channel: " << i << " --> " << fun->GetParameter(0) << " + " <<
       // fun->GetParameter(1) << "*Vbias" << endl;
       cout << "Channel: " << i << " breakdown --> " << -1*fun->GetParameter(0)/fun->GetParameter(1) << endl;
     }
  
  //Save PDF
  TString pdffilename="bias.pdf";
  for (size_t l=0;l<canvases.size();l++)
    {
      if (!l)
	canvases[l]->Print(pdffilename+"(","pdf");
      else if (l == canvases.size()-1)
	canvases[l]->Print(pdffilename+")","pdf");
      else
	canvases[l]->Print(pdffilename,"pdf");
    }
  
  return 0;
}
