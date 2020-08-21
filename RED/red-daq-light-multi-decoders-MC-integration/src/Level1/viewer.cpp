#include <iostream>
#include <cstring>
#include <map>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TTree.h>
#include <TF1.h>
#include <TApplication.h>
#include <TTimer.h>
#include <TSystem.h>

#include "getopt.h"

#include "Decoder.hh"
#include "MCDecoder.hh"
#include "Analyzer.hh"
#include "RDconfig.h"
#include "EvRaw0.hh"
#include "EvRec0.hh"
#include "RDconfig.h"
#include "RDBaseline.hh"
#include "RDFind_min.hh"
#include "RDIntegral.hh"
#include "RDFiltering.hh"
#include "RDClusterization.hh"
#include "RDPulseFitter.cc"

using namespace std;

void viewer(int run, string listfile, bool islaser, int sch, string customfile, string customMCfile, string serfile,
            bool isExtendedSER, string mappingfile, int customnsamples, bool generateTraces)
{
    TApplication* app = new TApplication("App",0,0);
    int ev = 0;
    std::string datadir = "";

    RDconfig::GetInstance()->SetIsExtendedSER(isExtendedSER);

    string configfile =(islaser) ?  "cfg/laser.cfg" :
                                    "cfg/global.cfg";
    if (customfile != "")
        configfile = customfile;

    // read configuration file
    RDconfig::GetInstance()->SetIsLaser(islaser);
    map<string, double> cfg  = RDconfig::GetInstance()->get_cfg(configfile);

    if (mappingfile != "")
        RDconfig::GetInstance()->SetConfigFile(mappingfile);
    if (serfile != "")
        RDconfig::GetInstance()->SetSERFile(serfile);

    string configmcfile = "cfg/MCconfig.cfg";
    if (customMCfile != "")
        configmcfile = customMCfile;

    VDecoder* theDecoder;
    if (generateTraces)
        theDecoder = new MCDecoder(configmcfile);
    else
        theDecoder = new Decoder(datadir);
    if (customnsamples>0)
        theDecoder->SetCustomNSamples(customnsamples);

    if (!(theDecoder->OpenFiles(listfile,run)))
    {
        fprintf (stderr,"Unable to open file. Exiting.\n");
        exit(1);
    }

    vector<double> ser = RDconfig::GetInstance()->GetSER();

    Analyzer* theAnalyzer = new Analyzer(islaser,configfile);

    if ( !isExtendedSER )
        for (int s = 0; s < ser.size(); s++)
            cout << "Ser " <<  s << " = " << ser[s] << endl;

    double scale = 0.75;
    TCanvas *c_top = new TCanvas("c_top", "c_top", 1000*scale, 1000*scale);
    TCanvas *c_bottom = new TCanvas("c_bottom", "c_bottom", 1000*scale, 1000*scale);
    TCanvas *c_sipm = new TCanvas("viewer_sipm", "viewer_sipm", 1000*scale, 620*scale);
    TCanvas *c_scint = new TCanvas("viewer_scint", "viewer_scint", 1000*scale, 1000*scale);

    c_top->Divide(6, 4); // final config
    c_bottom->Divide(2, 2);
    c_scint->Divide(4,4);

    TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);

    c_top->SetBatch(kFALSE);

    vector<TGraph*> g;
    vector<TGraph*> gf;
    vector<TGraph*> gp;
    vector<TGraph*> gw;
    vector<TGraph*> gc;

    vector<TF1*> fit;

    bool firsttime=true;

    RDBaseline* theBL = new RDBaseline;
    RDFind_min* theFM = new RDFind_min();
    Int_t evNumber = 0;
    if (sch >= 0)
        ev = sch;
    //cout << "Sono qui: "<< ev << " " << sch << endl;
    while (ev > -1) {
        // plot waveform with basic reconstruction

        EvRaw0* evRaw = theDecoder->ReadEvent();
        if (theDecoder->IsOver()) //File is over: rewind!
        {
            if (evNumber < ev)
            {
                cout << "The event #" << ev << " does not exist" << endl;
                exit(1);
            }
            theDecoder->CloseFiles(false);
            theDecoder->OpenFiles(listfile,run);
            continue;
        }
        evNumber = evRaw->GetEvHeader()->GetEventNum();
        if (evNumber < ev)
            continue;
        else if (evNumber > ev && !generateTraces) //rewind
        {
            delete theDecoder;

            VDecoder* theDecoder;
            theDecoder = new Decoder(datadir);
            theDecoder->OpenFiles(listfile,run);

            continue;
        }
        cout << "Display event #" << ev << endl;

        //--
        EvRec0 *theRecEvent = theAnalyzer->ProcessEvent(evRaw);

        //all waveforms here
        map<int,vector<double>* > wfs = evRaw->GetWFs();
        size_t n_samp = evRaw->GetSampleNumber(wfs.begin()->first);
        vector<double> wf_sipm(n_samp,0);

        Int_t n_channels = evRaw->GetChannelNumber();

        if (firsttime)
        {

            //if (n_channels <= 16)
            //  c->Divide(4, 4);
            //else
            //  c->Divide(4, 6);
            for (size_t k=0;k<ser.size();k++)
            {
                g.push_back(new TGraph);
                gf.push_back(new TGraph);
                gp.push_back(new TGraph);
                gw.push_back(new TGraph());
                gc.push_back(new TGraph());

                fit.push_back(new TF1(Form("fit%d", k),
                                      "[0]*x+[1]",
                                      cfg["baseline_start"],
                              cfg["baseline_stop"]));
                fit.at(k)->SetLineWidth(2);
                fit.at(k)->SetLineColor(kRed);
                g.at(k)->SetLineColor(kBlue);
                g.at(k)->SetLineWidth(2);
                gf.at(k)->SetLineColor(kCyan);
                gf.at(k)->SetLineWidth(2);
                gp.at(k)->SetMarkerStyle(23);
                gp.at(k)->SetMarkerColor(kGreen+1);
                gp.at(k)->SetMarkerSize(2);
                gw.at(k)->SetMarkerStyle(22);
                gw.at(k)->SetMarkerColor(kRed);
                gw.at(k)->SetMarkerSize(1.5);
                gc.at(k)->SetMarkerColor(kOrange);
                gc.at(k)->SetMarkerStyle(8);
                gc.at(k)->SetMarkerSize(1.5);
            }
            firsttime = false;
        }
        // single channels

        int ntop(0);
        int nbot(0);
        int nsci(0);

        //for (int k = 0; k < n_channels; k++) {
        for (map<int,vector<double>* >::iterator it=wfs.begin();
             it!=wfs.end(); ++it)
        {
            int k = it->first;
            vector<double>  *wf;
            if (wfs.count(k))
                wf = wfs.find(k)->second;
            else
                continue;

            //c->cd(k + 1);


            // baseline
            theBL->DoIt(wf,
                        (int) cfg["baseline_start"],
                    (int) cfg["baseline_stop"]);
            double bm = theBL->mean();

            fit.at(k)->SetParameters(theBL->m(), theBL->q());


            if (RDconfig::GetInstance()->IsSiPM(k) && ser[k]>0) {
                for (int i = 0; i < n_samp; i++) {
                    wf_sipm[i] += (wf->at(i) - bm)/ser[k];
                }
            }

            for (int i = 0; i < n_samp; i++)
                g.at(k)->SetPoint(i, i, wf->at(i));
            g.at(k)->SetTitle(Form("Event %d - Ch %d;samples;ADC", ev, k));
            // moving average
            RDFiltering *ff = new RDFiltering(wf, n_samp);
            vector<double> *wa = ff->mavg( cfg["moving_avg"]);
            vector<double>  waa(n_samp);
            vector<double>  wff(n_samp);

            for (int a = 0; a < n_samp; a++) {
                waa[a] = wa->at(a);
                wff[a] = wf->at(a);
            }

            delete ff;

            for (int i = 0; i < n_samp; i++)
                gf.at(k)->SetPoint(i, i, waa[i]);

            if (RDconfig::GetInstance()->IsTop(k)) {
                c_top->cd(1 + ntop);
                g.at(k)->SetTitle(Form("Ch%d - SiPM Top; Sa; ADC", k));
                g.at(k)->Draw("AL");
                gf.at(k)->Draw("Lsame");

                fit.at(k)->Draw("same");

                // find minimum
                theFM->DoIt(&wff, (int) cfg["charge_start"], (int) cfg["charge_stop"], (double) bm, (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);

                gp.at(k)->SetPoint(0, theFM->x_min(), theFM->y_min());
                gp.at(k)->Draw("Psame");
                gc.at(k)->SetPoint(0, theFM->s_time(), wff[theFM->s_time()]);
                gc.at(k)->Draw("Psame");


                // integration window
                gw.at(k)->SetPoint(0, (int) cfg["charge_start"], bm);
                gw.at(k)->SetPoint(1, (int) cfg["charge_stop"], bm);

                gw.at(k)->Draw("Psame");
                ntop++;
            }
            if (RDconfig::GetInstance()->IsBottom(k)) {
                c_bottom->cd(1 + nbot);
                g.at(k)->SetTitle(Form("Ch%d - SiPM Bottom; Sa; ADC", k));
                g.at(k)->Draw("AL");
                gf.at(k)->Draw("Lsame");

                fit.at(k)->Draw("same");

                // find minimum
                theFM->DoIt(&wff, (int) cfg["charge_start"], (int) cfg["charge_stop"],
                        (double) bm, (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);

                gp.at(k)->SetPoint(0, theFM->x_min(), theFM->y_min());
                gp.at(k)->Draw("Psame");
                gc.at(k)->SetPoint(0, theFM->s_time(), wff[theFM->s_time()]);
                gc.at(k)->Draw("Psame");

                // integration window
                gw.at(k)->SetPoint(0, (int) cfg["charge_start"], bm);
                gw.at(k)->SetPoint(1, (int) cfg["charge_stop"], bm);

                gw.at(k)->Draw("Psame");
                nbot++;
            }
            if (RDconfig::GetInstance()->IsLSci(k)) {
                c_scint->cd(1 + nsci);
                g.at(k)->SetTitle(Form("Ch%d - LSci; Sa; ADC", k));
                g.at(k)->Draw("AL");
                gf.at(k)->Draw("Lsame");

                fit.at(k)->Draw("same");

                // find minimum
                theFM->DoIt(&wff, (int) cfg["charge_start"], (int) cfg["charge_stop"],
                        (double) bm, (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);

                gp.at(k)->SetPoint(0, theFM->x_min(), theFM->y_min());
                gp.at(k)->Draw("Psame");
                gc.at(k)->SetPoint(0, theFM->s_time(), wff[theFM->s_time()]);
                gc.at(k)->Draw("Psame");

                // integration window
                gw.at(k)->SetPoint(0, (int) cfg["charge_start"], bm);
                gw.at(k)->SetPoint(1, (int) cfg["charge_stop"], bm);

                gw.at(k)->Draw("Psame");
                nsci++;
            }
            if (RDconfig::GetInstance()->IsSi(k)) {
                c_scint->cd(1 + nsci);
                g.at(k)->SetTitle(Form("Ch%d - Si; Sa; ADC", k));
                g.at(k)->Draw("AL");
                gf.at(k)->Draw("Lsame");

                fit.at(k)->Draw("same");

                // find minimum
                theFM->DoIt(&wff, (int) cfg["charge_start"], (int) cfg["charge_stop"],
                        (double) bm, (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);

                gp.at(k)->SetPoint(0, theFM->x_min(), theFM->y_min());
                gp.at(k)->Draw("Psame");
                gc.at(k)->SetPoint(0, theFM->s_time(), wff[theFM->s_time()]);
                gc.at(k)->Draw("Psame");

                // integration window
                gw.at(k)->SetPoint(0, (int) cfg["charge_start"], bm);
                gw.at(k)->SetPoint(1, (int) cfg["charge_stop"], bm);

                gw.at(k)->Draw("Psame");
                nsci++;
            }
        }
        c_top->Modified();
        c_top->Update();
        c_bottom->Modified();
        c_bottom->Update();
        c_scint->Modified();
        c_scint->Update();



        // all sipm sum
        c_sipm->cd();

        TGraph *g_sipm = new TGraph();
        g_sipm->SetTitle(Form("All SiPMs' Sum - Ev. number %d; samples; PE", ev));
        g_sipm->SetLineColor(kBlue);
        g_sipm->SetLineWidth(3);

        for (int i = 0; i < n_samp; i++)
            g_sipm->SetPoint(i, i, wf_sipm[i]);
        g_sipm->Draw("AL");


        // moving average
        TGraph *g_sipm_avg = new TGraph();
        g_sipm_avg->SetTitle(Form("All SiPMs' Sum - Ev. number %d; samples; PE", ev));
        g_sipm_avg->SetLineColor(kCyan);
        g_sipm_avg->SetLineWidth(3);

        RDFiltering *fff = new RDFiltering(&wf_sipm, n_samp);
        vector<double> *wa_sipm = fff->mavg(cfg["moving_avg"]);
        for (int i = 0; i < n_samp; i++)
            g_sipm_avg->SetPoint(i, i, wa_sipm->at(i));
        g_sipm_avg->Draw("Lsame");


        // baseline
        theBL->DoIt(&wf_sipm,
                    cfg["baseline_start"],
                cfg["baseline_stop"]);
        double bm = theBL->mean();

        // clusterization
        RDClusterization *cc = new
                RDClusterization(wa_sipm, (int) cfg["baseline_stop"], n_samp, cfg["cluster_thr"], (int) cfg["peak_sigma"],
                (int) cfg["cluster_start"],
                (int) cfg["cluster_stop"]);

        vector<int> start_time = cc->start_time();
        vector<int> c_start = cc->start();
        vector<int> c_stop = cc->stop();
        vector<int> c_rep = cc->rep();

        vector<double> min_x, min_y, s_time, f_time;
        int n_cl = start_time.size();

        if (n_cl == 0) {
            cout << "Cluster not found!" << endl;
            start_time.push_back((int) cfg["cluster_start"]);
            c_start.push_back(0);
            c_stop.push_back( (int) cfg["cluster_start"] + (int) cfg["cluster_stop"]);
            min_x.push_back(0);
            min_y.push_back(0);
            n_cl = 1;
        }
        else {
            TGraph *g_cl[n_cl];
            TGraph *gp_cl[n_cl];
            TGraph *gc_cl[n_cl];
            TGraph *gs_cl[n_cl];

            for (int n = 0; n < n_cl; n++) {
                theFM->DoIt(&wf_sipm, c_start[n], c_stop[n], (double) bm, (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);
                min_x.push_back(theFM->x_min());
                min_y.push_back(theFM->y_min());
                s_time.push_back(theFM->s_time());
                f_time.push_back(theFM->f_time());

                g_cl[n] = new TGraph();
                g_cl[n]->SetLineColor(kRed); g_cl[n]->SetLineWidth(2);
                g_cl[n]->SetPoint(0, c_start[n], bm);
                g_cl[n]->SetPoint(1, c_start[n], min_y[n]);
                g_cl[n]->SetPoint(2, c_stop[n], min_y[n]);
                g_cl[n]->SetPoint(3, c_stop[n], bm);
                g_cl[n]->Draw("Lsame");

                gp_cl[n] = new TGraph();
                gp_cl[n]->SetMarkerStyle(23);
                gp_cl[n]->SetMarkerColor(kGreen+1);
                gp_cl[n]->SetMarkerSize(2);
                gp_cl[n]->SetPoint(0, min_x[n], min_y[n]);
                gp_cl[n]->Draw("Psame");

                gc_cl[n] = new TGraph();
                gc_cl[n]->SetMarkerStyle(8);
                gc_cl[n]->SetMarkerColor(kOrange);
                gc_cl[n]->SetMarkerSize(2);
                gc_cl[n]->SetPoint(0, s_time[n], wf_sipm[(int) s_time[n]]);
                gc_cl[n]->Draw("Psame");

                gs_cl[n] = new TGraph();
                gs_cl[n]->SetMarkerStyle(33);
                gs_cl[n]->SetMarkerColor(kRed);
                gs_cl[n]->SetMarkerSize(2.5);
                gs_cl[n]->SetPoint(0, f_time[n], 0);
                gs_cl[n]->Draw("Psame");

            }
        }
        delete cc;

        // Draw S2 fits
        TGraph *gfitS2=new TGraph();
        TGraphErrors *gS2=new TGraphErrors();
        auto fFitter = new RDPulseFitter();

        for (int ifit=0;ifit<theRecEvent->GetFits().size();ifit++)
        {
            //if ( theRecEvent->GetFits(ifit)->type!=1 ) continue;
            int start = theRecEvent->GetFits(ifit)->start;
            int stop  = theRecEvent->GetFits(ifit)->end;
            fFitter->ReBin(&wf_sipm,start,stop,gS2,0.);

            double par[11];
            for (int ipar=0;ipar<11;ipar++) par[ipar]=theRecEvent->GetFits(ifit)->par[ipar];

            for (int sa=start;sa<stop;sa++)
            {
                double v[1]={double(sa)};
                double val=(islaser) ? fFitter->ImpulseFunc(v,par) : fFitter->yfit_sa(v,par);
                gfitS2->SetPoint(gfitS2->GetN(), (double) sa, val );
            }
        }
        gfitS2->SetLineWidth(2);

        if ( gS2->GetN()>0 ) { gS2->SetMarkerColor(kRed); gS2->SetMarkerSize(1.2); gS2->SetMarkerStyle(20); gS2->Draw("Psame"); }
        if ( gfitS2->GetN()>0 ) { gfitS2->SetLineWidth(3); gfitS2->SetLineColor(kBlack); gfitS2->Draw("LPsame"); }

        // cumulative sum

        TGraph *g_cum = new TGraph();
        g_cum->SetLineColor(kViolet); g_cum->SetLineWidth(2);

        vector<double> *cs = fff->cumsum();

        if (min_y[0] != 0) {
            for (int i = 0; i < n_samp; i++) {
                cs->at(i) = cs->at(i)/cs->at(n_samp - 1)*min_y[0];
                g_cum->SetPoint(i, i, cs->at(i));
            }
        }

        g_cum->Draw("same");

        // draw canvas
        c_sipm->Modified();
        c_sipm->Update();

        delete fff;

        timer->TurnOn(); // for intective modification of the canvas
        timer->Reset();

        // progress opt
        char cev[256];
        int rev;
        cout << "Insert event number, RETURN to continue or -1 to quit: ";
        cin.getline(cev, 256);
        rev = atoi(cev);
        if (rev < 0) {
            cout << "Goodbye!" << endl;
            ev = rev;
            // return;
        }
        else if (rev == 0) ev++;
        else ev = rev;

        timer->TurnOff();
    }

    //Cleanup
    for (size_t i=0;i<g.size();i++)
    {
        delete g.at(i);
        delete gf.at(i);
        delete gp.at(i);
        delete gw.at(i);
        delete fit.at(i);
    }
    delete theBL;
    delete theFM;
    //clean exit
    //theDecoder->CloseFiles(false);
    delete theDecoder;
    delete app;
    return;
}

int main(int argc, char *argv[]) {

    int c;
    std::string listfile = "";
    int nevent = -1;
    int runnr = 0;
    bool islaser = false;
    std::string configfile = "";
    std::string configmcfile = "";
    std::string serfile = "";
    std::string mappingfile = "";
    int customnsamples = -1;
    bool isExtendedSER = false;
    bool generateTraces = false;
    static struct option long_options[] =
    {
        {"run",required_argument,0,'r'},
        {"list",required_argument,0,'l'},
        {"nevent",required_argument,0,'n'},
        {"config",required_argument,0,'c'},
        {"configmc",required_argument,0,'i'},
        {"mapping",required_argument,0,'m'},
        {"ser",required_argument,0,'e'},
        {"extendedser",no_argument,0,'f'},
        {"laser",no_argument,0,'s'},
        {"nsamples",required_argument,0,'a'},
        {"generateMCtraces",no_argument,0,'g'},
        {"help",no_argument,0,'h'},
        {0, 0, 0, 0}
    };

    // getopt_long stores the option index here.
    int option_index = 0;

    // Parse options
    while ((c = getopt_long (argc, argv, "r:l:n:c:i:e:m:a:fsgh",long_options,
                             &option_index)) != -1) {
        switch (c)
        {
        case 'r':
            if (runnr!=0) {
                fprintf (stderr, "Error while processing option '-r'. Multiple runs specified.\n");
                exit(1);
            }
            if ( sscanf(optarg,"%d",&runnr) != 1 ) {
                fprintf (stderr, "Error while processing option '-r'. Wrong parameter '%s'.\n", optarg);
                exit(1);
            }
            if (runnr<0) {
                fprintf (stderr, "Error while processing option '-r'. Run number set to %d (must be >=0).\n", runnr);
                exit(1);
            }
            fprintf(stdout,"Merging files from run %d\n",runnr);
            break;
        case 'n':
            if ( sscanf(optarg,"%d",&nevent) != 1 ) {
                fprintf (stderr, "Error while processing option '-n'. Wrong parameter '%s'.\n", optarg);
                exit(1);
            }
            //nevent = atoi(optarg);
            fprintf(stdout,"Starting from event %d\n",nevent);
            break;
        case 'l':
            listfile = optarg;
            fprintf(stdout,"Data will be read from files listed in '%s'\n",listfile.c_str());
            break;
        case 'c':
            configfile = optarg;
            fprintf(stdout,"Using custom config file: %s \n",configfile.c_str());
            break;
        case 'i':
            configmcfile = optarg;
            fprintf(stdout,"Using custom MC config file: %s \n",configmcfile.c_str());
            break;
        case 'e':
            serfile = optarg;
            fprintf(stdout,"Using custom ser file: %s \n",serfile.c_str());
            break;
        case 'f':
            isExtendedSER = true;
            fprintf(stdout,"Using extended ser file.\n");
            break;
        case 's':
            islaser = true;
            fprintf(stdout,"This is a laser runs");
            break;
        case 'a':
            customnsamples = atoi(optarg);
            if (customnsamples<=0)
            {
                fprintf (stderr, "Error while processing option '--nsamples'. Wrong parameter '%s'.\n", optarg);
                exit(1);
            }
            fprintf(stdout,"Set custom number of samples %d\n",customnsamples);
            break;
        case 'm':
            mappingfile = optarg;
            fprintf(stdout,"Using channel mapping from file: %s \n",mappingfile.c_str());
            break;
        case 'g':
            fprintf(stdout,"Events will be generated from Monte Carlo simulation.\n");
            generateTraces=true;
            break;
        case 'h':
            fprintf(stdout,"\nviewer ([-r run_number]|[-l list_file]) [-n eventNb] [-c config file]  [--laser] [-h]\n\n");
            fprintf(stdout,"  -r: define run to process\n");
            fprintf(stdout,"  -l: define file with list of data files to process\n");
            fprintf(stdout,"      n.b. either -r or -l must be specified\n");
            fprintf(stdout,"  -n: event number to be displayed (default: %d)\n",nevent);
            fprintf(stdout,"  -m: define channel mapping from file (default: cfg/channelmapping.cfg)\n");
            fprintf(stdout,"  -c: define custom config file\n");
            fprintf(stdout,"  -i: define custom MC config file\n");
            fprintf(stdout,"  --ser (or -e): Read SERs from a custom file (default: cfg/ser.cfg)\n");
            fprintf(stdout,"  --extendedser (or -f): Read from a extended ser file (multi parameter ser file)\n");
            fprintf(stdout,"  --laser (or -s): laser run (default = normal run)\n");
            fprintf(stdout,"  --nsamples (or -a): override the number of samples from the event header (use with care!) \n");
            fprintf(stdout,"  -g: events will be generated from Monte Carlo simulation \n");
            fprintf(stdout,"  -h: show this help message and exit\n\n");
            exit(0);
        case '?':
            fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
            exit(1);
        default:
            abort();
        }
    }

    // Verify that some input was specified
    if (!generateTraces)
    {
        if ( listfile.compare("")==0 && runnr==0 ) {
            fprintf (stderr,"No run number and no list file were specified. Exiting.");
            exit(1);
        }
        if ( listfile.compare("")!=0 && runnr!=0 ) {
            fprintf (stderr,"Both run number and list file were specified. Exiting.");
            exit(1);
        }
    }
    else
    {
        if ( listfile.compare("")!=0 ) {
            fprintf (stderr,"Both generate MC traces and list file were specified. Exiting.");
            exit(1);
        }
        if ( runnr!=0 ) {
            fprintf (stderr,"Both generate MC traces and run number were specified. Exiting.");
            exit(1);
        }
        /* FIXME */
        if ( mappingfile.compare("cfg/channelmapping_Na.cfg")!=0 ) {
            fprintf (stderr,"For now to generate MC traces is only allowed the channel mapping of Naples:\n");
            mappingfile = "cfg/channelmapping_Na.cfg";
            fprintf(stdout,"- Using channel mapping from file: %s\n",mappingfile.c_str());
        }
    }

    //Main call
    viewer(runnr,listfile,islaser,nevent,configfile,configmcfile,serfile,isExtendedSER,mappingfile,customnsamples,generateTraces);

    return 0;
}
