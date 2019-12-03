#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "ADCBoard.hh"
//#include "DBService.hh"
#include "RootIO.hh"
#include "InitOutputRoot.hh"
#include "TApplication.h"
#include "TH2F.h"
#include <TCanvas.h>

#include "regex.h"

#include "baseline.h"
#include "peakfinder.h"
#include "integral.h"

using namespace std;


#define FILENAME "viewer_conf"
#define MAXBUF 1024
#define DELIM "="

struct config
{	
   int int_min, int_max;
};

struct config get_config(char *filename)
{
        struct config configstruct;
        FILE *file = fopen (filename, "r");

        if (file != NULL)
        {
                char line[MAXBUF];
                int i = 0;

                while(fgets(line, sizeof(line), file) != NULL)
                {	
                        char *cfline;
                        cfline = strstr((char *)line,DELIM);
                        cfline = cfline + strlen(DELIM);
   
                        if (i == 0){
                                configstruct.int_min = atoi(cfline);
                                //printf("%s",configstruct.imgserver);
                        } else if (i == 1){
				configstruct.int_max = atoi(cfline);
                                
                        }                       
                        i++;
                } // End while
                fclose(file);
        } // End if file
       
               
       
        return configstruct;

}

//MAIN___________________________________________________________________________
int main(int argc, char* argv[]){
int c;
char datadir[100] = "/home/MRescigno/rawdata";
//const char *dir = "ls ../PadmeDAQ/data -Art | tail -n 3 | head -n 1";//"ls ../RunControl/data/rawdata -Art | tail -n 3 | head -n 1";
const char *dir = "ls /home/MRescigno/rawdata -Art | tail -n 3 | head -n 1";
char temp[100];
char filePath[100];// = datadir+"/"+"data/rawdata/run_9_b00_2017_07_27_11_30_34";

char last_file[20], previous_file[20];

FILE *f_data;
//
//f_data = popen(dir, "r");
//fscanf(f_data, "%s", &last_file);
//pclose(f_data);

//printf("%s\n",last_file);

//std::string filePath = "";
//string datafile;
//int runnr = 9;
int b_id = 0;
while ((c = getopt (argc, argv, "f:")) != -1) {
    	switch (c){
		case 'f':
        		//filePath = datadir+"/"+optarg;
			strcpy(datadir, optarg);
        		fprintf(stdout,"Set input data directory to '%s'\n",datadir);
			strcpy(temp, "ls ");
			strcat(temp, datadir);
			strcat(temp, " -Art | tail -n 3 | head -n 1");
			dir = &temp[0];
        		break;
     		default:
			
        		break;
	}
}

strcpy(filePath, datadir);
strcat(filePath, "/");

ADCBoard* board;
std::vector<ADCBoard*> boards;
std::vector<ADCBoard*>::iterator it;
board = new ADCBoard(b_id);
double wf[NSAMPLE_MAX];
//int canale = 0;
int conto  = 0;

n_samp = 10000;

TH1F *charge_hist[8];
for(int i = 0; i< 8; i++){
  if (i>5)  charge_hist[i] = new TH1F (Form("charge_hist_ch%d", i),"charge",250,-1000,300000); //Lsci
  else if(i>3) charge_hist[i] = new TH1F (Form("charge_hist_ch%d", i),"charge",250,-1000,30000); //SenSL // laser
  else charge_hist[i] =    new TH1F (Form("charge_hist_ch%d", i),"charge",250,-1000,300000); //FBK

	//if(i>3) charge_hist[i] = new TH1F (Form("charge_hist_ch%d", i),"charge",250,-1000,600000); //SenSL // normal
	//else charge_hist[i] =    new TH1F (Form("charge_hist_ch%d", i),"charge",250,-1000,6000000); //FBK
}
// = new TH1F ("charge_hist","charge",500,0,50000);

TH1F *WF0 = new TH1F ("WF0","Waveform for ch 0",n_samp,0,n_samp-1);
TH1F *WF1 = new TH1F ("WF1","Waveform for ch 1",n_samp,0,n_samp-1);
TH1F *WF2 = new TH1F ("WF2","Waveform for ch 2",n_samp,0,n_samp-1);
TH1F *WF3 = new TH1F ("WF3","Waveform for ch 3",n_samp,0,n_samp-1);
TH1F *WF4 = new TH1F ("WF4","Waveform for ch 4",n_samp,0,n_samp-1);
TH1F *WF5 = new TH1F ("WF5","Waveform for ch 5",n_samp,0,n_samp-1);
TH1F *WF6 = new TH1F ("WF6","Waveform for ch 6",n_samp,0,n_samp-1);
TH1F *WF7 = new TH1F ("WF7","Waveform for ch 6",n_samp,0,n_samp-1);

TH2F *F90vsE = new TH2F("F90vsE", "F90vsE", 100,0, 50000, 100,0,0.3);

TApplication theApp("App",&argc, argv);

TCanvas *c1 = new TCanvas("c1","c1", 1200,700);
TCanvas *c2 = new TCanvas("c2","c2", 1200,700);

c1->Divide(4,2);
c2->Divide(4,2);


while(1 ){
	//GetEventCounter(); //number of events in the file
	board = new ADCBoard(b_id);
	board->Reset();

	strcpy(previous_file, last_file);

	f_data = popen(dir, "r");
	fscanf(f_data, "%s", &last_file);
	pclose(f_data);
        strcpy(filePath,datadir);
	strcat(filePath,"/");
	strcat(filePath,last_file);
	
	//printf("%c\n", filePath[0]);
	if(strstr(filePath,"_b0") != NULL){ 
	  //	   printf("OK\n");
	}

	if(strcmp(previous_file, last_file) == 0){
		usleep(10000);
	
	}
	else if(strstr(filePath,"_b0") != NULL) {
	
		c1->SetTitle(last_file);
		c2->SetTitle(last_file);
		// Add files to board
		//std::string filePath = datadir+"/"+"data/rawdata/run_9_b00_2017_07_27_11_30_34";
		printf("File %d - Reading from file %s\n",0,filePath);
		board->AddFile(filePath);
		// Add board to boards vector
		//boards.push_back(board);

		//Show list of known boards/files

		for (it = boards.begin(); it != boards.end(); ++it) {
			board = *it;
			printf("Board %d Files %d\n",board->GetBoardId(),board->GetNFiles());
			printf("File %d %s\n",0,board->GetFileName(0).c_str());
		}

		//_______________________________________________________________________________________________

	
		board->NextEvent();
		n_samp = board->Event()->GetAcceptedChannelMask() >> 16 ;
		n_channels=0;

		  for(int i=0; i<16; i++){
		    
		    if(board->Event()->GetAcceptedChannelMask() & (0x1<<i))
		       n_channels++;
		     
		  }

		//TH1F *charge_hist = new TH1F ("charge_hist","charge",500,0,100000);
		//charge_hist->Reset();
		WF0->Reset();
		WF0->SetBins(n_samp, 0, n_samp);
		WF1->Reset();
		WF1->SetBins(n_samp, 0, n_samp);
		WF2->Reset();
		WF2->SetBins(n_samp, 0, n_samp);
		WF3->Reset();
		WF3->SetBins(n_samp, 0, n_samp);
		WF4->Reset();
		WF4->SetBins(n_samp, 0, n_samp);
		WF5->Reset();
		WF5->SetBins(n_samp, 0, n_samp);
		WF6->Reset();
		WF6->SetBins(n_samp, 0, n_samp);
		WF7->Reset();
		WF7->SetBins(n_samp, 0, n_samp);
		//printf("%d\n",n_samp);

		for (int s=0; s < n_samp; s++){
			wf0[s]=board->Event()->GetADCChannelSample(0,s);
			wf1[s]=board->Event()->GetADCChannelSample(1,s);
			wf2[s]=board->Event()->GetADCChannelSample(2,s);
			wf3[s]=board->Event()->GetADCChannelSample(3,s);
			wf4[s]=board->Event()->GetADCChannelSample(4,s);
			wf5[s]=board->Event()->GetADCChannelSample(5,s);
			wf6[s]=board->Event()->GetADCChannelSample(6,s);
			wf7[s]=board->Event()->GetADCChannelSample(7,s);

			WF0->Fill(s,wf0[s]);
			WF1->Fill(s,wf1[s]);
			WF2->Fill(s,wf2[s]);
			WF3->Fill(s,wf3[s]);
			WF4->Fill(s,wf4[s]);
			WF5->Fill(s,wf5[s]);
			WF6->Fill(s,wf6[s]);
			WF7->Fill(s,wf7[s]);
		}

		//WF0->Draw("hist");
		  
		//int wf[NSAMPLE_MAX];
		//int canale = 0;
		conto  = 0;

		while(board->NextEvent()!= 0){
			//if(conto%2){
				charge_tot=0;
				for (int j=0; j < 8; j++) {
					memset(wf,0,sizeof(wf));
				    	for (int s=0; s < n_samp; s++){
						wf[s]=board->Event()->GetADCChannelSample(j,s);
				    	}
					
					baseline b(wf, 0, 400); // normal
					//baseline b(wf, 0,400, 1600, 2000); // laser

					baseline_mean[j] = b.mean();
					baseline_rms[j] = b.rms();
					// peakfinder
					//peakfinder pf(wf, 0, n_samp);
					//peak_x[j] = pf.x_peak();
					//peak_y[j] = pf.y_peak();
					// charge & f90
					integral *c;

					if(j>5) c = new integral(wf, b.m(), b.q(), 870, 1050, 910);
					else if(j>3) c = new integral(wf, b.m(), b.q(), 400, 3500, 200); // normal 
					else c = new integral(wf, b.m(), b.q(), 600, 3500, 200); 
	
					
					//if(j>3) c = new integral(wf, b.m(), b.q(), 400, 1600, 200); // laser
					//else c = new integral(wf, b.m(), b.q(), 400, 1600, 200);

					//charge[j] = c.charge()/ser[j];
					//f90[j] = c.f90();
   	   		if(j>6) F90vsE->Fill(c->charge(), c->f90()); 
					  //printf("f90: %f\n", c->f90());  
					charge_hist[j]->Fill(c->charge());


					//f90_tot += charge[j]*f90[j]; 
					//printf("%f   ",charge_tot);
                                       
					delete c;  
				}
			//}
			conto++;
		}

		//printf("%s\n\n",processName())

		c1->cd(1);
		WF0->Draw("hist");
		c1->cd(2);
		WF1->Draw("hist");
		c1->cd(3);
		WF2->Draw("hist");
		c1->cd(4);
		WF3->Draw("hist");
		c1->cd(5);
		WF4->Draw("hist");
		c1->cd(6);
		WF5->Draw("hist");
		c1->cd(7);
		WF6->SetAxisRange(800., 1100.,"X");
		WF6->Draw("hist");
		c1->cd(8);
		WF7->SetAxisRange(800., 1100.,"X");
		WF7->Draw("hist");
		gPad->Modified(); gPad->Update();
		
		for(int i = 0; i<8; i++){
			c2->cd(i+1);//->SetLogy();
			charge_hist[i]->Draw();
		}
				c2->cd(7);
				F90vsE->Draw();

		

		//theApp.Run();
		gPad->Modified(); gPad->Update();
		delete board;
		//delete WF0; delete WF1; delete WF2; delete WF3;
		usleep(500000); //useconds
	}
}
//theApp.Terminate();
return 0;

}

