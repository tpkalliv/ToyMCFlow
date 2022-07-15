#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include <vector>
/* #include "../include/toyflowinputs.h" */
#include "/home/tkallio/softwares/root_macros/toyMCFlowcode/ToyMCFlow2022/include/rootcommon.h"

using namespace std;

const Int_t NH = 5; 

Int_t gMarkers[]= {20,24,21,25,22,26,23,27,32,28};
Int_t gColors[]={kRed+1, kOrange+2, kCyan+2, kSpring-6, kRed-7, kOrange+1,kCyan-6,kGreen+7,kRed-9,kOrange-9,kAzure+6,kGreen-9};
Int_t gStyles[]={1,2,3,4,5,6,7,8,9,10};

TH1D *hPhiPsi[NH];
TH1D *hEventPlane[NH];

TF1 *fFit[NH];
Double_t vn[NH];
Double_t vnError[NH];
Double_t v_psi[NH];
Double_t fvns[NH];


void LoadData(TString);
void FitDrawPhi();
void SaveVns();

//-----------Main Function------------------
void FitSingle(TString infile="outputs/out_tytk.root")
{
	LoadData(infile);
	FitDrawPhi();
	SaveVns();
}

//-------Member Functions-----------------
void LoadData(TString inputname)
{
	TFile *fIn = TFile::Open(inputname,"read");
	for(int ih=1; ih<=NH; ih++){
		hPhiPsi[ih-1]=(TH1D*)fIn->Get(Form("hPhiPsiH%02d", ih));
	}

	
	for(int ih = 0; ih < NH; ih++){
		hEventPlane[ih]=(TH1D*)fIn->Get(Form("hEventPlane%02d", ih+1));
		fvns[ih] = hEventPlane[ih]->GetMean();
		cout << fvns[ih] << endl;
	}
	
}

void FitDrawPhi()
{	
	for (Int_t ih=1; ih<=NH; ih++){
		TString formula = Form("[0]*(1 + 2*[1]*TMath::Cos(%d*x))",ih);
		fFit[ih] = new TF1(Form("fFitH%02d",ih), formula,0, 2.0*TMath::Pi());
		fFit[ih]->SetParameter(0,1000);
		fFit[ih]->SetParameter(ih, fvns[ih-1]);
	}
	
	for (Int_t ih=1; ih<=NH; ih++) hPhiPsi[ih-1]->Fit(fFit[ih]->GetName());

	for (Int_t ih=1; ih<=NH; ih++){
		vn[ih]=fFit[ih]->GetParameter(1);
		vnError[ih]=fFit[ih]->GetParError(1);
	}


	gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas("C2","canvas",1024,740);
	
	can->SetFillStyle(4000);
	can->SetLeftMargin(0.15);
   	can->SetBottomMargin(0.15);
   	
   	//For editing canvas #include "include/rootcommon.h"
   	double lowx = 0.,highx=2*TMath::Pi();
  	//double ly=hDeltaPhiSum[ic]->GetMinimum()*0.99,hy=hDeltaPhiSum[ic]->GetMaximum()*1.01;
  	Double_t ly = 8500, hy=16500;
  	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
  	hset( *hfr, "#Delta#phi=#phi-#psi_{n}", "dN/d#Delta#phi", 1.1, 1.1, 0.05,0.05, 0.01,0.01, 0.03,0.03, 510, 505);//settings of the upper pad: x-axis, y-axis
  	hfr->Draw();
  	TLegend *legend = new TLegend(0.5,0.6,0.8,0.85,"","brNDC");
    legend->SetTextSize(0.02);legend->SetBorderSize(0);legend->SetFillStyle(0);//legend settings;
  	/* legend->AddEntry((TObjArray*)NULL,Form("Centrality %s",strCentrality[ic].Data())," "); */


  	for(int ih=2; ih<=NH; ih++) {

  		cout << ih << endl;
  		//hPhiPsi[ih][ic]->SetTitle(Form("Cent%02d",ic));
  		hPhiPsi[ih-1]->SetMarkerStyle(20);
		hPhiPsi[ih-1]->SetMarkerColor(gColors[ih]);
		hPhiPsi[ih-1]->Draw("psame");		
		fFit[ih]->SetLineColor(gColors[ih]);
		fFit[ih]->Draw("same");
		legend->AddEntry(fFit[ih],Form("v_{%d} = %.3f #pm %.4f ",ih, vn[ih], vnError[ih]),"l");

  	}
 

	fFit[2]->Draw();

	for (int i = 3; i < NH+1; i++) {
		fFit[i]->Draw("same");
	} 	
	

  	legend->Draw();
	gPad->GetCanvas()->SaveAs(Form("SingleParticle.pdf"));



}



void SaveVns()
{
	TFile *output = new TFile("out_VnFitSingle.root", "recreate");
	for (int i = 0; i < NH; i++) { hPhiPsi[i]->Write(); } 
	output->Write();
	output->Close();
}