#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TFile.h"
#include "TRatioPlot.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include <sstream>
#include <string>

double DeltaPhi(double phi1, double phi2);

void FlowCoefficients() {
	const Int_t NH = 5; //Number of harmonics
	const Int_t NP = 1000; //Number of particles
	const Int_t NE = 100; //Number of events

	cout << "Number of harmonics: " << NH << endl;
	cout << "Number of particles: " << NP << endl;
	cout << "Number of events: " << NE << "\n" << endl;

	Double_t v_n[NH] = {0.0,0.12, 0.06, 0.03, 0.013}; //




	Double_t v_psi[NH] = {0.0}; // Analytical EP Method: Values of Cos(n*(psi_i-Psi_n))
	Double_t v_obs[NH] = {0.0}; // Q-vector EP Method: Values of Cos(n*(psi_i-Psi_EP))
	Double_t v_deltaphi[NH] = {0.0}; // Two-Particle Method: Values of Cos(n*(psi_i-psi_j))
	Double_t Psi_n[NH] = {0.0}; // Symmetry plane values (in radians from Reaction Plane) created by sampling Uniform distribution

	// Mean values
	Double_t MeanEventPlane[NH];
	Double_t MeanResolution[NH];
	Double_t MeanTwoParticle[NH];
	Double_t MeanEPObs[NH];

	// Standard deviation (sigma)
	Double_t ErrorEP[NH];
	Double_t ErrorResolution[NH];
	Double_t ErrorTP[NH];
	Double_t ErrorEPObs[NH];

	// Histogram initializations
	TH1D *hEventPlane[NH];
	TH1D *hEventPlaneEP[NH];
	TH1D *hPhiPsi[NH];
	TH1D *hPhiPsiQ[NH];
	TH1D *hEPObs[NH];
	TH1D *hResolution[NH];
	TH1D *hTwoParticle[NH];

	TH1D *hDeltaPhiSum = new TH1D ("hDeltaPhiSum", "hDeltaPhiSum", 200, 0.0, 2*TMath::Pi());

	// Histograms for summing up calculated Phi-Psi -angle for every method
	TH1D *hPhiPsi_Sum = new TH1D ("hPhiPsi_Sum", "hPhiPsi_Sum", 200, 0.0, 2*TMath::Pi());
	TH1D *hPhiPsiQ_Sum = new TH1D ("hPhiPsiQ_Sum", "hPhiPsiQ_Sum", 200, 0.0, 2*TMath::Pi());
	TH1D *hPhiPsi2P_Sum = new TH1D ("hPhiPsi2P_Sum", "hPhiPsi2P_Sum", 200, 0.0, 2*TMath::Pi());


	// Produces a string of Fourier decomposition about the anisotropic flow as azimuthal distribution
	string cosine = "[0]*(1";
	for (int i=1; i<NH+1; i++) {
		ostringstream app;
		app << "+2*[" << i << "]*TMath::Cos(" << i << "*(x-[" << i+NH << "]))";
		string append = app.str();
		cosine = cosine + append;
	}
	cosine = cosine + ")";
	const char* cos = cosine.c_str();
	cout << cosine << endl;

	TFile *output  = TFile::Open("outputs/out_tytk.root", "recreate");

	TF1 *fourier = new TF1("Fourier", cos, 0.0, 2.0*TMath::Pi());

	TF1 *uniform[NH];


	TStopwatch *timer = new TStopwatch();
	Double_t time;
	timer->Start();

	// Creating histogram for every observable harmonic
	for (int i=0; i<NH; i++) {
		uniform[i]= new TF1(Form("uniform%02d", i+1), "1.0", -1.0*TMath::Pi()/(i+1), 1.0*TMath::Pi()/(i+1));
	

		hEventPlane[i]=new TH1D(Form("hEventPlane%02d",i+1),Form("hEventPlane%02d",i+1),200, -1.0, 1.0);
		hEPObs[i]=new TH1D(Form("hEpobs%02d",i+1),Form("hEpobs%02d",i+1),200,-1.0,1.0);
		hEventPlaneEP[i]=new TH1D(Form("hEventPlaneEP%02d",i+1),Form("hEventPlaneEP%02d",i+1),200,-1.0, 1.0);
		hPhiPsiQ[i]=new TH1D(Form("hPhiPsiQ%02d",i+1),Form("hPhiPsiQ%02d",i+1),200,0.0, 2.0*TMath::Pi());
		hTwoParticle[i]=new TH1D(Form("hTwoParticle%02d",i+1),Form("hTwoParticle%02d",i+1),200, -1.0, 1.0);
		hPhiPsi[i] = new TH1D(Form("hPhiPsiH%02d",i+1),Form("hPhiPsiH%02d",i+1),200,0.0, 2.0*TMath::Pi());
		hResolution[i] = new TH1D(Form("hResolution%02d",i+1),Form("hResolution%02d",i+1),200,-100, 100);
	}

	// Filling up Fourier decomposition with input values Vn and Psi_n
	fourier->SetParameter(0,NP);
	

	

	// Going through all the events (event: nucleus-nucleus collision)
	for (int i=0; i<NE;i++) {
		Double_t Qn_x[NH] = {0.0};
		Double_t Qn_y[NH] = {0.0};
		Double_t phis[NP] = {0.0};
		Double_t Resolution[NH] = {0.0};
		Double_t Psi_EP[NH] = {0.0};


		for (int n=0; n<NH; n++) {
			Psi_n[n] = uniform[n]->GetRandom();
			fourier->SetParameter(n+1,v_n[n]);
			fourier->SetParameter(n+1+NH, Psi_n[n]);
		}


		// Going through 1000 particles (from one event nucleus-nucleus collision)
		for (int j=0;j<NP;j++) {
			phis[j] = fourier->GetRandom();

			
			// Going through every harmonic for each particle
			for (int k=0; k<NH; k++) {

				// Analytical Event Plane Method:

				hPhiPsi[k]->Fill(DeltaPhi(phis[j],Psi_n[k]));

				v_psi[k]= TMath::Cos((k+1)*DeltaPhi(phis[j],Psi_n[k]));

				hEventPlane[k]->Fill(v_psi[k]);
				/* 
				fourier_analyticalEP->SetParameter(k+1, v_psi[k]);
				fourier_analyticalEP->SetParameter(k+NH+1, Psi_n[k]);
				*/

				Qn_x[k] += TMath::Cos((k+1)*phis[j]);
				Qn_y[k] += TMath::Sin((k+1)*phis[j]);
			}

		}


		// Calculating Event Plane with Q-vector average
		for (Int_t n=0; n<NH;n++) {
			Psi_EP[n] = TMath::ATan2(Qn_y[n],Qn_x[n])/(n+1); 
		}

		// Going through 1000 particles for an Event
		for (int i=0; i<NP; i++) {

			// Event Plane Method: With Q-vectors
			for (int n=0; n<NH; n++) {

				/*
				hPhiPsiQ_Sum->Fill(DeltaPhi(phis[i],Psi_EP[n]));
				*/

				v_obs[n] = TMath::Cos((n+1)*DeltaPhi(phis[i],Psi_EP[n]));
				hEPObs[n]->Fill(v_obs[n]);

				/*
				fourier_obs->SetParameter(n+1, v_obs[n]);
				fourier_obs->SetParameter(n+NH+1, Psi_n[n]);
				*/

			}

			// Two Particle Method
			for (int j=0; j<NP; j++) {
				if (j==i) continue; // Excluding auto-correlations (double counting)
				for (int n=0; n<NH; n++) {
					
					hDeltaPhiSum->Fill(DeltaPhi(phis[i],phis[j]));

					v_deltaphi[n] = TMath::Cos((n+1)*DeltaPhi(phis[i],phis[j]));
					hTwoParticle[n]->Fill(v_deltaphi[n]); 

					/*
					fourier_2p->SetParameter(n+1, v_deltaphi[n]);
					fourier_2p->SetParameter(n+NH+1, Psi_n[n]);
					*/
				}
			}
		}

		// Correction 
		for (int n=0; n<NH; n++) {
			Resolution[n] = TMath::Cos((n+1)*(DeltaPhi(Psi_n[n],Psi_EP[n])));
			hResolution[n]->Fill(Resolution[n]);
		}

	}


	for (int i = 0; i < NH; i++){
		hPhiPsi[i]->Write();
		hEventPlane[i]->Write();
	}


	
	TCanvas *hDeltaPhiSum_Canvas = new TCanvas ("hDeltaPhiSum", "hDeltaPhiSum", 900, 700); 
	hDeltaPhiSum->Write();
	hDeltaPhiSum->Draw();
	

	output->cd();
	output->Write();
	output->Close();


	time = timer->CpuTime();
	cout << "Time: " << time << "s" << endl;
}


// Calculating angle difference using ATan2-function

double DeltaPhi(double phi1, double phi2) {
	// dphi
	double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
	return res>0 ? res : 2.*TMath::Pi()+res ;

}

/*

Number of harmonics: 5
Number of particles: 1000
Number of events: 100

v_EP1: 0.00313163
v_EPObs1: 0.337481
v_phi1: 0.00981959

v_EP2: 0.11792
v_EPObs2: 0.121513
v_phi2: 0.117196

v_EP3: 0.0621664
v_EPObs3: 0.0702543
v_phi3: 0.0617024

v_EP4: 0.0297577
v_EPObs4: 0.0560131
v_phi4: 0.0283041

v_EP5: 0.0307163
v_EPObs5: 0.0610436
v_phi5: 0.0316459

Time: 55.1s

*/