#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TH2.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TText.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom.h"
#include "/home/doornen/ROOT/Macros/Fitter/Fitter.h"

// Actual fitting procedure ---------------------------------------
void EfficiencyFitterAngleCut(int caseNumber)  {
  //c2->Divide(3,5);
  //c2->cd();
  TCanvas *fcanvas;  
  char dummytext[1000];
  char root_file[200];
  TH1F *h[200];
 
  sprintf(dummytext,"./output/FittingResults/Case%iAngleCut.txt",caseNumber);

  fFitOutput = fopen(dummytext,"w");
  fprintf(fFitOutput,"BeamEnergy_GammaEnergy_FittedEnergy_sigma_Resolution_sigma_Counts_sigma \n");
  
  sprintf(dummytext,"Canvas Fitting Case %i",caseNumber);
  fcanvas = new TCanvas(dummytext,dummytext,1400,1000);
  fcanvas->Divide(5,12);
  
  sprintf(root_file,"../../../SimulationResults/Reconstructor/Shogun/Efficiency/Case%i_100MeV1000keV.root",caseNumber);
  TFile *infile = new TFile(root_file,"READ");
  
  for(int k=0;k<18;k++) {
    sprintf(dummytext,"doppler_angle_cut[%i];1",k);
    infile->GetObject(dummytext,h[k]);
    cout<<" I am here"<<endl;
    Int_t minvalue = 0.9*1000;
    Int_t maxvalue = 1.08*1000+10;
    
    Int_t dummy3 = 0.98*1000;
    Int_t dummy4 = 1.02*1000;
    
    double range_min = 0.82*1000;
    double range_max = 1.18*1000;
    fcanvas->cd(k+1);
    h[k]->Draw();
    h[k]->GetXaxis()->SetRangeUser(range_min,range_max);
    fcanvas->Update();
    Double_t CPos[5]={minvalue,maxvalue,1000,dummy3,dummy4};
    Double_t FRes[100];
    RT_FitManPar(h[k],CPos,FRes,0,1,1,0,100,1000);
  } 
  fclose(fFitOutput); 
}
