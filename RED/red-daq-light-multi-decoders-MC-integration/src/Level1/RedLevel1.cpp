#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <string>

#include "getopt.h"

#include "ADCBoard.hh"
#include "Decoder.hh"
#include "MCDecoder.hh"
#include "Analyzer.hh"
#include "RunHeader.hh"
#include "RDconfig.h"

#include "RootIO.hh"

int main(int argc, char* argv[])
{
    int c;
    std::string listfile = "";
    std::string outfile = "reco/rawdata.root";
    std::string datadir = "";
    //
    int neventsperfile = 10000000; //very large number (= default)
    int runnr = 0;
    int verbose = 0;
    int skip = 0;
    int updatedb = 0;
    bool islaser = false;
    bool writewaveforms = false;
    std::string configfile = "";
    std::string configmcfile = "";
    std::string mappingfile = "";
    std::string serfile = "";
    bool isExtendedSER = false;
    bool customserfile=false;
    bool custommapfile=false;
    bool aligntimes = true;
    bool generateTraces = false;
    int customnsamples = -1;

    static struct option long_options[] =
    {
        {"datadir",required_argument, 0, 'd'},
        {"outdir",required_argument,0,'o'},
        {"run",required_argument,0,'r'},
        {"list",required_argument,0,'l'},
        {"verbose",required_argument,0,'v'},
        {"skip",required_argument,0,'k'},
        {"config",required_argument,0,'c'},
        {"configmc",required_argument,0,'i'},
        {"mapping",required_argument,0,'m'},
        {"ser",required_argument,0,'e'},
        {"extendedser",no_argument,0,'f'},
        {"neventsperfile",required_argument,0,'n'},
        {"update",no_argument,0,'u'},
        {"laser",no_argument,0,'s'},
        {"waveforms",no_argument,0,'w'},
        {"no-sync",no_argument,0,'x'},
        {"nsamples",required_argument,0,'a'},
        {"generateMCtraces",no_argument,0,'g'},
        {"help",no_argument,0,'h'},
        {0, 0, 0, 0}
    };

    // getopt_long stores the option index here.
    int option_index = 0;

    // Parse options
    while ((c = getopt_long (argc, argv, "d:r:l:o:v:k:n:c:i:m:e:a:fuwxgh",long_options,
                             &option_index)) != -1) 
    {
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
            if ( sscanf(optarg,"%d",&neventsperfile) != 1 ) {
                fprintf (stderr, "Error while processing option '-n'. Wrong parameter '%s'.\n", optarg);
                exit(1);
            }
            if (neventsperfile<0) {
                fprintf (stderr, "Error while processing option '-n'. Number of events must be >=0, found %d.\n", neventsperfile);
                exit(1);
            }
            fprintf(stdout,"Writing up to %d events per output file\n",neventsperfile);
            break;
        case 'd':
            datadir = optarg;
            fprintf(stdout,"Set input data directory to '%s'\n",datadir.c_str());
            break;
        case 'o':
            outfile = optarg;
            fprintf(stdout,"Set output data file to '%s'\n",outfile.c_str());
            break;
        case 'l':
            listfile = optarg;
            fprintf(stdout,"Data will be read from files listed in '%s'\n",listfile.c_str());
            break;
        case 'e':
            serfile = optarg;
            fprintf(stdout,"Using custom ser file: %s \n",serfile.c_str());
            customserfile = true;
            break;
        case 'f':
            isExtendedSER = true;
            fprintf(stdout,"Using extended ser file.\n");
            break;
        case 'c':
            configfile = optarg;
            fprintf(stdout,"Using custom config file: %s \n",configfile.c_str());
            break;
        case 'i':
            configmcfile = optarg;
            fprintf(stdout,"Using custom MC config file: %s \n",configmcfile.c_str());
            break;
        case 'm':
            mappingfile = optarg;
            fprintf(stdout,"Using channel mapping from file: %s \n",mappingfile.c_str());
            custommapfile = true;
            break;
        case 's':
            islaser = true;
            fprintf(stdout,"This is a laser run \n");
            break;
        case 'v':
            if ( sscanf(optarg,"%d",&verbose) != 1 ) {
                fprintf (stderr, "Error while processing option '-v'. Wrong parameter '%s'.\n", optarg);
                exit(1);
            }
            if (verbose<0) {
                fprintf (stderr, "Error while processing option '-v'. Verbose level set to %d (must be >=0).\n", verbose);
                exit(1);
            }
            fprintf(stdout,"Set verbose level to %d\n",verbose);
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
        case 'k':
            if ( sscanf(optarg,"%d",&skip) != 1 ) {
                fprintf (stderr, "Error while processing option '-k'. Wrong parameter '%s'.\n", optarg);
                exit(1);
            }
            if (skip<=0) {
                fprintf (stderr, "Error while processing option '-k'. Verbose level set to %d (must be >0).\n", verbose);
                exit(1);
            }
            fprintf(stdout,"Skip first %d events\n",skip);
            break;
        case 'u':
            fprintf(stdout,"Enabling DB update\n");
            updatedb = 1;
            break;
        case 'w':
            fprintf(stdout,"Writing waveforms on the output file\n");
            writewaveforms = true;
            break;
        case 'x':
            fprintf(stdout,"Events from the boards will *NOT* be aligned.\n");
            aligntimes=false;
            break;
        case 'g':
            fprintf(stdout,"Events will be generated from Monte Carlo simulation.\n");
            generateTraces=true;
            break;
        case 'h':
            fprintf(stdout,"\nRedLevel1 ([-r run_number]|[-l list_file]) [-d input files directory] [-o output root file] [-n events per file] [-c config file] [-m mappingfile] [-v verbosity] [--laser] [--skip] [-e serfile] [--extendedser] [-w] [-g] [-h]\n\n");
            fprintf(stdout,"  -r: define run to process\n");
            fprintf(stdout,"  -l: define file with list of data files to process\n");
            fprintf(stdout,"      n.b. either -r or -l must be specified\n");
            fprintf(stdout,"  -u: update DB with total number of events merged (default: no DB update)\n");
            fprintf(stdout,"  -d: define directory where input files are located (default: \"%s\")\n",datadir.c_str());
            fprintf(stdout,"  -o: define an output file in root format (default: \"%s\")\n",outfile.c_str());
            fprintf(stdout,"  -n: define max number of events per output file (0=no limit, default: %d)\n",neventsperfile);
            fprintf(stdout,"  -c: define custom config file\n");
            fprintf(stdout,"  -i: define custom MC config file\n");
            fprintf(stdout,"  -m: define channel mapping from file\n");
            fprintf(stdout,"  --ser (or -e): Read SERs from a custom file (default: cfg/ser.cfg)\n");
            fprintf(stdout,"  --extendedser (or -f): Read from a extended ser file (multi parameter ser file)\n");
            fprintf(stdout,"  --skip (or -k): Set the number of initial events to be skipped (default: 0)\n");
            fprintf(stdout,"  -v: define verbose level (default: %d)\n",verbose);
            fprintf(stdout,"  --laser (or -s): laser run (default = normal run)\n");
            fprintf(stdout,"  --no-sync (or -x): do not sync events from the board (default = sync)\n");
            fprintf(stdout,"  --waveforms (or -w): write waveforms on the output file (default = false)\n");
            fprintf(stdout,"  --nsamples (or -a): override the number of samples from the event header (use with care!)\n");
            fprintf(stdout,"  --generateMCtraces (or -g): events will be generated from Monte Carlo simulation\n");
            fprintf(stdout,"  -h: show this help message and exit\n\n");
            exit(0);
        case '?':
            if (optopt == 'v') {
                // verbose with no argument: just enable at minimal level
                verbose = 1;
                break;
            } else if (optopt == 'r' || optopt == 'l' || optopt == 'd' || optopt == 'o' || optopt == 'n')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint(optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
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
            fprintf (stderr,"No run number and no list file were specified. Exiting.\n");
            exit(1);
        }
        if ( listfile.compare("")!=0 && runnr!=0 ) {
            fprintf (stderr,"Both run number and list file were specified. Exiting.\n");
            exit(1);
        }        
    }
    else
    {
        if ( listfile.compare("")!=0 ) {
            fprintf (stderr,"Both generate MC traces and list file were specified. Exiting.\n");
            exit(1);
        }
        if ( runnr!=0 ) {
            fprintf (stderr,"Both generate MC traces and run number were specified. Exiting.\n");
            exit(1);
        }
        /* FIXME */
        if ( mappingfile.compare("cfg/channelmapping_Na.cfg")!=0 ) {
            fprintf (stderr,"For now to generate MC traces is only allowed the channel mapping of Naples.\n");
            custommapfile = true;
            mappingfile = "cfg/channelmapping_Na.cfg";           
        }
    }
    
    TFile* fout = new TFile(outfile.c_str(), "recreate");
    fout->cd();
    Int_t n_samp;
    Int_t n_channels;
    const int NMAXR = 48;
    Double_t sers[NMAXR];
    Int_t chanmap[NMAXR];
    std::vector<TString> chanID;
    TTree* metaevent = new TTree("metaevent", "Event meta info");
    metaevent->Branch("n_samples", &n_samp);
    metaevent->Branch("n_channels", &n_channels);
    metaevent->Branch("run_number", &runnr);
    metaevent->Branch("sers",sers,"sers[n_channels]/D");
    metaevent->Branch("chanmap",chanmap,"chanmap[n_channels]/I");
    metaevent->Branch("chanID",&chanID);
    for (int u=0;u<NMAXR;u++)
    {
        sers[u]=-1;
        chanmap[u]=5; //initialize to kUndefined
        chanID.push_back("NULL");
    }

    //Ok, now start real things
    RDconfig::GetInstance()->SetIsExtendedSER(isExtendedSER);
    
    if (custommapfile)
        RDconfig::GetInstance()->SetConfigFile(mappingfile);
    
    VDecoder* theDecoder;
    if (generateTraces)
        theDecoder = new MCDecoder(configmcfile);
    else
        theDecoder = new Decoder(datadir);
    
    theDecoder->SetVerbosity(verbose);
    if (aligntimes)
        theDecoder->SetAlignTimes(true);
    if (skip)
        theDecoder->SetEventsToSkip(skip);
    if (customnsamples>0)
        theDecoder->SetCustomNSamples(customnsamples);
    Analyzer* theAnalyzer = new Analyzer(islaser,configfile);
    theAnalyzer->SetVerbosity(verbose);

    TTree* raweventtree = new TTree("raw","Raw events");
    raweventtree->Branch("rawevent",theDecoder->GetRawEvent());

    TTree* recoeventtree = new TTree("reco","Reco events");
    recoeventtree->Branch("recoevent",theAnalyzer->GetRecEvent());

    if (!(theDecoder->OpenFiles(listfile,runnr)))
    {
        fprintf (stderr,"Unable to open file. Exiting.\n");
        exit(1);
    }
    RunHeader* theRunHeader = theDecoder->GetRunHeader();

    time_t timer;
    time(&timer);

    TTree* runheader = new TTree("runheader","Run Header");
    runheader->Branch("runheader",theRunHeader);

    // Loop over all events in files
    static bool firsttime = true;
    for(int iloop=0; !(theDecoder->IsOver());iloop++)
    {
        if (verbose>=1) printf("=== Processing event %8d ===\n",iloop);

        EvRaw0* theEvent = theDecoder->ReadEvent();

        if (theDecoder->IsOver())
            continue;

        //int eventnr = theDecoder->GetEventNumber();
        if (firsttime)
        {
            n_samp = theEvent->GetSampleNumber();
            n_channels = theEvent->GetChannelNumber();
            if (customserfile)
                RDconfig::GetInstance()->SetSERFile(serfile);
            vector<double> ser = RDconfig::GetInstance()->GetSER();
            if (ser.size() < n_channels)
            {
                cerr << "SER values are " << ser.size() << " and are hence insufficient for " <<
                        n_channels << " channels " << endl;
                cerr << "Please check " << endl;
                abort();
            }
            //Loop on the valid files
            map<int,vector<double>* > wfs = theEvent->GetWFs();
            for (map<int,vector<double>* >::iterator it=wfs.begin();
                 it!=wfs.end(); ++it)
            {
                sers[it->first] = ser.at(it->first);
                chanmap[it->first] = RDconfig::GetInstance()->GetChannelType(it->first);
                chanID.at(it->first) = RDconfig::GetInstance()->GetChannelName(it->first);
            }

            metaevent->Fill();
            metaevent->Write();
            string runtype = (islaser) ? "laser" : "physics";
            theRunHeader->SetRunType(runtype);
            firsttime = false;
        }
        time_t now; time(&now);
        if (!(iloop%500))
            printf("=== Processing event : iloop=%8d evtn=%8d    (%7.1f s) ===\n",iloop,theEvent->GetEvHeader()->GetEventNum(),difftime(now,timer));
        if ( theEvent->GetEvHeader()->GetEventNum() > neventsperfile ) break;

        EvRec0* theRecEvent = theAnalyzer->ProcessEvent(theEvent);
        if (writewaveforms)
            raweventtree->Fill();
        //if (theRecEvent->GetNClusters())
        recoeventtree->Fill();
    }
    theRunHeader->SetNumEvt(raweventtree->GetEntries());
    runheader->Fill();
    runheader->Write();
    if (raweventtree->GetEntries())
        raweventtree->Write();
    recoeventtree->Write();
    TTree* configtree = RDconfig::GetInstance()->GetConfigurationTree();
    if (configtree)
        configtree->Write();

    fout->Close();
    delete fout;

    int retcode = theDecoder->CloseFiles(updatedb);
    delete theDecoder;

    return retcode;
}
