{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(42);
  //gStyle->SetTextFontSize(18);
  TLatex lat;
  TCanvas *c = new TCanvas("c","",0,0,800,600);
  c->Divide(2,1);  
  c_1->SetLogy();

  TFile g("run_715.root");
  TTree* t1 = (TTree*) g.Get("reco");

  TH1D* h1 = new TH1D("h1","Timing",150,-150,150);
  t1->Project("h1","(clusters.cdf_time-xmin[30])*2.",
	      "number_of_clusters>0 && clusters.charge > 40");
  h1->GetXaxis()->SetTitle("#Deltat (TPC-Si) (ns)");
  h1->GetYaxis()->SetTitle("Events");
  c->cd(1);
  h1->Draw();  

  
  TH1D* h2 = new TH1D("h2","F90",80,0.,0.8);
  t1->Project("h2","clusters.f90",
	      "number_of_clusters>0 && clusters.charge > 40 && clusters.cdf_time>2865 && clusters.cdf_time<2890");
  h2->GetXaxis()->SetTitle("PSD discriminator f_{90}");
  c->cd(2);
  h2->Draw();
  lat.DrawLatex(0.1,160,"#gamma-like");
  lat.DrawLatex(0.6,160,"n-like");
    
 
}
