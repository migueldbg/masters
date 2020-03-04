#include "root/rootstart.h"
void generateS1Corr(){
	SetMyStyle();
	gStyle->SetNumberContours(99);
	TFile *f = TFile::Open("allscan_Hist.root");
	TGraph2D *S1Map[10][2];
	double topMeanMax=0, botMeanMax=0;
	double topMeanMiddle, botMeanMiddle;
		TH3D * htop3D = (TH3D*)(f->Get("hS1Top_barXY_z4"));
		TProfile2D * ptop = (htop3D)->Project3DProfile("yx");
		topMeanMiddle = ptop->GetMean(3);
		TH3D * hbot3D = (TH3D*)(f->Get("hS1Bottom_barXY_z4"));
		TProfile2D * pbot = (hbot3D)->Project3DProfile("yx");
		botMeanMiddle = pbot->GetMean(3);
	double topMean[10], botMean[10];
	TCanvas *c = new TCanvas("c","c");
	c->Divide(3,2);
	int nbinsxy = 15;
	for(int iz=0; iz<10; iz++){
		S1Map[iz][0] = new TGraph2D();
		S1Map[iz][0]->SetNameTitle(Form("S1CorrMapTop_z%d",iz),Form("S1CorrMapTop_z%d",iz));
		S1Map[iz][1] = new TGraph2D();
		S1Map[iz][1]->SetNameTitle(Form("S1CorrMapBot_z%d",iz),Form("S1CorrMapBot_z%d",iz));
		//top
		TH3D * htop3D = (TH3D*)(f->Get(Form("hS1Top_barXY_z%d",iz)));
		htop3D->FitSlicesZ(0,1,0,1,0,10,"QNLG1S");
		TH2D * htopMeanFit = (TH2D*)gDirectory->Get(Form("hS1Top_barXY_z%d_1",iz));
		TH2D * htopConstFit = (TH2D*)gDirectory->Get(Form("hS1Top_barXY_z%d_0",iz));
		TProfile2D * ptop = (htop3D)->Project3DProfile("yx");
		TH2D * htop = (TH2D*)((htop3D)->Project3D("yx"));
		topMean[iz] = ptop->GetMean(3);
		cout<<"top mean"<<topMean[iz]<<" ";
		//ptop->Scale(1./topMean[iz]);
		ptop->Scale(1./topMeanMiddle);
		ptop->SetMaximum(1.5); ptop->SetMinimum(0.5);
		c->cd(1); ptop->Draw("COLZ");
		c->cd(2); htop->Draw("COLZtext");
		if(topMean[iz]>topMeanMax)topMeanMax = topMean[iz];
		for(int ix=0; ix<nbinsxy; ix++){
			for(int iy=0; iy<nbinsxy; iy++){
				if(htop->GetBinContent(htop->GetBin(ix+1,iy+1)) > 15){
					double smooth = 
						ptop->GetBinContent(htop->GetBin(ix+1,iy+1)) * htop->GetBinContent(htop->GetBin(ix+1,iy+1))+
						ptop->GetBinContent(htop->GetBin(ix+0,iy+1)) * htop->GetBinContent(htop->GetBin(ix+0,iy+1))+
						ptop->GetBinContent(htop->GetBin(ix+2,iy+1)) * htop->GetBinContent(htop->GetBin(ix+2,iy+1))+
						ptop->GetBinContent(htop->GetBin(ix+1,iy+0)) * htop->GetBinContent(htop->GetBin(ix+1,iy+0))+
						ptop->GetBinContent(htop->GetBin(ix+1,iy+2)) * htop->GetBinContent(htop->GetBin(ix+1,iy+2));
					smooth /= htop->GetBinContent(htop->GetBin(ix+1,iy+1))+
                              htop->GetBinContent(htop->GetBin(ix+0,iy+1))+
                              htop->GetBinContent(htop->GetBin(ix+2,iy+1))+
                              htop->GetBinContent(htop->GetBin(ix+1,iy+0))+
                              htop->GetBinContent(htop->GetBin(ix+1,iy+2));
					S1Map[iz][0]->SetPoint(S1Map[iz][0]->GetN(), ix*0.13333+0.05+1.5, iy*0.13333+0.05+1.5, 1./smooth);//1./ptop->GetBinContent(htop->GetBin(ix+1,iy+1)));
				}
			}
		}
		if(S1Map[iz][0]->GetN()>3){
		c->cd(3);
		S1Map[iz][0]->SetMaximum(1.5);
		S1Map[iz][0]->SetMinimum(0.5);
		S1Map[iz][0]->SetMarginBinsContent(topMeanMiddle/topMean[iz]);
		S1Map[iz][0]->Draw("colz");
		}
		//bot
		TH3D * hbot3D = (TH3D*)(f->Get(Form("hS1Bottom_barXY_z%d",iz)));
		hbot3D->FitSlicesZ(0,1,0,1,0,10,"QNLG1S");
		TH2D * hbotMeanFit = (TH2D*)gDirectory->Get(Form("hS1Bottom_barXY_z%d_1",iz));
		TH2D * hbotConstFit = (TH2D*)gDirectory->Get(Form("hS1Bottom_barXY_z%d_0",iz));
		TProfile2D * pbot = (hbot3D)->Project3DProfile("yx");
		TH2D * hbot = (TH2D*)((hbot3D)->Project3D("yx"));
		botMean[iz] = pbot->GetMean(3);
		cout<<"bottom mean"<<botMean[iz]<<endl;
		//pbot->Scale(1./botMean[iz]);
		pbot->Scale(1./botMeanMiddle);
		pbot->SetMaximum(1.5); pbot->SetMinimum(0.5);
		c->cd(4); pbot->Draw("COLZ");
		c->cd(5); hbot->Draw("COLZtext");
		if(botMean[iz]>botMeanMax)botMeanMax = botMean[iz];
		for(int ix=0; ix<nbinsxy; ix++){
			for(int iy=0; iy<nbinsxy; iy++){
				if(hbot->GetBinContent(hbot->GetBin(ix+1,iy+1)) > 15){
					double smooth = 
						pbot->GetBinContent(hbot->GetBin(ix+1,iy+1)) * hbot->GetBinContent(hbot->GetBin(ix+1,iy+1))+
						pbot->GetBinContent(hbot->GetBin(ix+0,iy+1)) * hbot->GetBinContent(hbot->GetBin(ix+0,iy+1))+
						pbot->GetBinContent(hbot->GetBin(ix+2,iy+1)) * hbot->GetBinContent(hbot->GetBin(ix+2,iy+1))+
						pbot->GetBinContent(hbot->GetBin(ix+1,iy+0)) * hbot->GetBinContent(hbot->GetBin(ix+1,iy+0))+
						pbot->GetBinContent(hbot->GetBin(ix+1,iy+2)) * hbot->GetBinContent(hbot->GetBin(ix+1,iy+2));
					smooth /= hbot->GetBinContent(hbot->GetBin(ix+1,iy+1))+
                              hbot->GetBinContent(hbot->GetBin(ix+0,iy+1))+
                              hbot->GetBinContent(hbot->GetBin(ix+2,iy+1))+
                              hbot->GetBinContent(hbot->GetBin(ix+1,iy+0))+
                              hbot->GetBinContent(hbot->GetBin(ix+1,iy+2));
					S1Map[iz][1]->SetPoint(S1Map[iz][1]->GetN(), ix*0.13333+0.05+1.5, iy*0.13333+0.05+1.5, 1./smooth);//1./pbot->GetBinContent(hbot->GetBin(ix+1,iy+1)));
				}
			}
		}
		c->cd(6);
		if(S1Map[iz][1]->GetN()>3){
		S1Map[iz][1]->SetMaximum(1.5);
		S1Map[iz][1]->SetMinimum(0.5);
		S1Map[iz][1]->SetMarginBinsContent(botMeanMiddle/botMean[iz]);
		S1Map[iz][1]->Draw("colz");
		}
		c->Update();
		if(iz==0){c->Print("S1Map_fieldScan.pdf(");}
		else if(iz==9){ c->Print("S1Map_fieldScan.pdf)");}
		else{  c->Print("S1Map_fieldScan.pdf");}
	}
	TFile *fout = new TFile("S1Map_fieldScan.root","RECREATE");
	fout->cd();
	for(int iz=0; iz<10; iz++){
		S1Map[iz][0]->Write();
		S1Map[iz][1]->Write();
	}
	fout->Close();
	////normalize in z to maximum.
	//TCanvas *c2 = new TCanvas("c2");
	//c2->Divide(5,4);
	//for(int iz=0; iz<10; iz++){
	//	for(int ip=0; ip<S1Map[iz][0]->GetN(); ip++){
	//		S1Map[iz][0]->SetPoint(ip, S1Map[iz][0]->GetX()[ip], S1Map[iz][0]->GetY()[ip], S1Map[iz][0]->GetZ()[ip]*topMeanMax/topMean[iz]);
	//	}
	//	for(int ip=0; ip<S1Map[iz][1]->GetN(); ip++){
	//		S1Map[iz][1]->SetPoint(ip, S1Map[iz][1]->GetX()[ip], S1Map[iz][1]->GetY()[ip], S1Map[iz][1]->GetZ()[ip]*botMeanMax/botMean[iz]);
	//	}
	//	c2->cd(iz*2+1);
	//	S1Map[iz][0]->Draw("colz");
	//	c2->cd(iz*2+2);
	//	S1Map[iz][1]->Draw("colz");
	//}
	//c2->Print("S1Map_fieldScan.pdf)");
}
