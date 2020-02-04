//#include "Root_Plot.cc"

TTree *t;

TCut cutall;
TCut cutFits="";

void Load(int run, bool isMC=false, int ERorNR=0, bool ph2=true, bool s2Fit=false, bool speFit=false)
{
    //SetStyle();
    //--- Clustering tree
    TFile *f=new TFile(Form("run_%d%s%s%s.root",run, (isMC)?"_MC":"", (ERorNR == 1)?"ER":"", (ERorNR == 2)?"NR":""));
    f->GetObject("reco",t);
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.5);

    t->SetAlias("s1", "clusters[0].charge");
    t->SetAlias("s1_top", "clusters[0].tot_charge_top");
    t->SetAlias("s1_bot", "clusters[0].tot_charge_bottom");
    t->SetAlias("tba1", "(s1_top-s1_bot)/s1");
    t->SetAlias("ratio", "s2/s1");
    t->SetAlias("nc", "number_of_clusters");
    t->SetAlias("f1", "clusters[0].f90");
    t->SetAlias("f90", "f90_tot");
    t->SetAlias("r1", "clusters[0].rep");

    if (ph2)
    {
        t->SetAlias("s2", "clusters[1].charge");
        t->SetAlias("s2_top", "clusters[1].tot_charge_top");
        t->SetAlias("s2_bot", "clusters[1].tot_charge_bottom");
        t->SetAlias("tba2", "(s2_top-s2_bot)/s2");
        t->SetAlias("f2", "clusters[1].f90");
        t->SetAlias("r2", "clusters[1].rep");
        t->SetAlias("t1", "clusters[0].start_time");
        t->SetAlias("t2", "clusters[1].start_time");
        t->SetAlias("tdrift", "(t2-t1)*2.0/1000.0");
    }

    if (s2Fit)
    {
        t->SetAlias("nfits", "nfits");
        t->SetAlias("status", "fits.status");
        t->SetAlias("covstatus", "fits.covstatus");
        t->SetAlias("ndf", "fits.ndf");
        t->SetAlias("chi2", "fits.chi2");
        t->SetAlias("tau2", "fits.par[0][1]");
        t->SetAlias("p", "fits.par[0][2]");
        t->SetAlias("T", "fits.par[0][3]");
        t->SetAlias("sigma", "fits.par[0][4]");
        t->SetAlias("t0", "fits.par[0][6]");
        t->SetAlias("e_p", "fits.epar[0][2]");
        t->SetAlias("e_T", "fits.epar[0][3]");
        t->SetAlias("e_sigma", "fits.epar[0][4]");
        t->SetAlias("e_t0", "fits.epar[0][6]");
    }

    else if (speFit)
    {
        t->SetAlias("nfits", "nfits");
        t->SetAlias("status", "fits.status");
        t->SetAlias("covstatus", "fits.covstatus");
        t->SetAlias("ndf", "fits.ndf");
        t->SetAlias("chi2", "fits.chi2");
        t->SetAlias("tau", "fits.par[0][1]");
        t->SetAlias("p", "fits.par[0][2]");
        t->SetAlias("sigma", "fits.par[0][3]");
        t->SetAlias("e_tau", "fits.epar[0][1]");
        t->SetAlias("e_p", "fits.epar[0][2]");
        t->SetAlias("e_sigma", "fits.epar[0][4]");
    }

    if (ph2)
        cutall = "nc==2 && f1>0.2 && f2<0.2 && r1==1 && r2==1 && s1>200 && s1<1000";
    else
        cutall = "nc==1 && f1>0.1 && r1==1 && s1>200 && s1<1000";
    if (s2Fit)
        cutFits = "status==0 && covstatus==3 && ndf>0 && chi2/ndf<1.3";
    else if (speFit)
        cutFits = "status==0 && covstatus==3 && ndf>0 && chi2/ndf<1.3 && sigma>5 && sigma<15 && tau>450 && tau<750 && p>0.02 && p<0.1";
}
