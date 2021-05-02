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
void EfficiencyFitter(int caseNumber)  {
  //c2->Divide(3,5);
  //c2->cd();
  TCanvas *fcanvas;  
  char dummytext[1000];
  char root_file[200];
  TH1F *h[200];
  //TH1F *h2[12][14];
  //TH1F *h2;
  sprintf(dummytext,"./output/FittingResults/Case%ic.txt",caseNumber);

  fFitOutput = fopen(dummytext,"w");
  fprintf(fFitOutput,"BeamEnergy_GammaEnergy_FittedEnergy_sigma_Resolution_sigma_Counts_sigma \n");
  int k = 0;

  sprintf(dummytext,"Canvas Fitting Case %i",caseNumber);
  fcanvas = new TCanvas(dummytext,dummytext,1400,1000);
  fcanvas->Divide(5,12);

  for(int i=100;i<=200;i+=50)  {
    for(int j=250;j<=3000;j+=250)  {
   
      cout<<"Reading Root-tree"<<i<<endl;

      if(j<1000)
        sprintf(root_file,"../../../SimulationResults/Reconstructor/Shogun/Efficiency/Case%ic_%iMeV%ikeV.root",caseNumber,i,j);
      else
        sprintf(root_file,"../../../SimulationResults/Reconstructor/Shogun/Efficiency/Case%ic_%iMeV%ikeV.root",caseNumber,i,j);
      TFile *infile = new TFile(root_file,"READ");
      sprintf(dummytext,"doppler_total;1");
      infile->GetObject(dummytext,h[k]);
      
      //c2->Update(); 
      // Values of Fitting Range 
      Int_t minvalue = 0.9*j;
      Int_t maxvalue = 1.08*j+10;
       
      Int_t dummy3 = 0.99*j;
      Int_t dummy4 = 1.01*j;
      
      double range_min = 0.82*j;
      double range_max = 1.18*j;
      fcanvas->cd(k+1);
      h[k]->Draw();
      h[k]->GetXaxis()->SetRangeUser(range_min,range_max);
      fcanvas->Update();
      Double_t CPos[5]={minvalue,maxvalue,j,dummy3,dummy4};
      Double_t FRes[100];
      RT_FitManPar(h[k],CPos,FRes,0,1,1,0,i,j);
      //cout<<"Fine till here 5 "<<j<<endl;
      k++;
    }
  }
  fclose(fFitOutput); 
}
