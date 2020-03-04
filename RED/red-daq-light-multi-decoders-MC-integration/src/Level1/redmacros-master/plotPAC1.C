{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(42);
  //gStyle->SetTextFontSize(18);
  TLatex lat;
  TFile f("run_624.root");
  TCanvas *c = new TCanvas("c","",0,0,1200,600);
  //c->Divide(2,1);  
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.,0.,0.6,1.);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.6,0.,1.,1.);
  pad1->Draw();
  pad2->Draw();

  pad1->SetLogz();
  //c_1->SetLogz();
  //c_2->SetLogz();
  //c->SetLogz();
  TTree* t1 = (TTree*) f.Get("reco");
  
  TH2D* h2 = new TH2D("h2","Full Li plot",200,1000,11000,
		      200,1000,12000);

  t1->Project("h2","ymax[9]-baseline_mean[9]:ymax[11]-baseline_mean[11]");
  h2->GetXaxis()->SetTitle("E (a.u.)");
  h2->GetYaxis()->SetTitle("#DeltaE (a.u.)");
  
  TBox* box = new TBox(3500,4000,8200,11500);
  box->SetLineStyle(7);
  box->SetFillStyle(0);
  box->SetLineColor(kBlack);
  box->SetLineWidth(2);
  pad1->cd();
  h2->Draw("colz");
  box->Draw("same");
  lat.DrawLatex(4800,8500,"^{7}Be");
  lat.DrawLatex(7000,7000,"^{7}Be");
  lat.DrawLatex(1500,7000,"^{7}Li");
  lat.DrawLatex(5000,1500,"#alpha");
  // TLatex* lat1 = new TLatex(5000,9000,"^{7}Be");
  //lat1->Draw("same");

  TFile g("run_680.root");
  TTree* t2 = (TTree*) g.Get("reco");
  TH2D* h3 = new TH2D("h3","Empty box",200,3500,8200,
		      200,4000,11500);
  pad2->cd();
  h3->GetXaxis()->SetTitle("E (a.u.)");
  //h3->GetYaxis()->SetTitle("#DeltaE (a.u.)");

  h3->Draw();  
  Int_t n = 
    t2->Draw("ymax[9]-baseline_mean[9]:ymax[11]-baseline_mean[11]","","goff");
 
  TGraph *gr = new TGraph(n,t2->GetV2(),t2->GetV1()); 
  gr->SetMarkerColor(kRed);
  gr->Draw("zp"); 
  
  Int_t n2 = t2->Draw("ymax[9]-baseline_mean[9]:ymax[11]-baseline_mean[11]","number_of_clusters>0 && clusters[0].cdf_time>2800 && clusters[0].cdf_time<3200  && clusters[0].f90>0.3","goff");
  
  TGraph *gr2 = new TGraph(n2,t2->GetV2(),t2->GetV1()); 
  gr2->SetMarkerStyle(7);
  gr2->Draw("zp"); 

}
