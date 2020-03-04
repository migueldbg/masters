{
  TFile *_file0 = TFile::Open("run_74_int.root.light");
  TTree* reco = (TTree*) _file0->Get("t1");

  
  TH1D* light = new TH1D("light","h0",400,-20,20);
  reco->Project("light","tof","nchan==17 && chargeFar>25000 && chargeNear>12500 && f90Far<0.13 && f90Near<0.13");

  TFile *_file1 = TFile::Open("run_74_30_int.root");
  //TFile *_file0 = TFile::Open("run_74_int.root");
  TTree* reco2 = (TTree*) _file1->Get("reco");
  TH1D* orig = new TH1D("orig","h1",400,-20,20);
  reco2->Project("orig","(start_time[17]-start_time[0])*2.+16*(evheader.boardtimes[1]-evheader.boardtimes[0])/2.","charge[0]/7500*60>200. && charge[17]/7500.*60.>100. && f90[0]<0.13 && f90[17]<0.13 ");

  orig->SetLineColor(kRed);
  orig->DrawCopy();
  light->DrawCopy("same");

}
