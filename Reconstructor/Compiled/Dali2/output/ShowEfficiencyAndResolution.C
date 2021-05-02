#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TH1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TText.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"

FILE *fDataInput;

void ShowEfficiencyAndResolution(int eventSimulated,int binning)
{  
  char dummyText[1000];
  char inputFile[1000];
  char outputFile[1000];
  float detector;
  float energyBeam[100] = {0.0};
  float energyGamma[100] = {0};
  float energyGammaFitted[100] = {0.0};
  float energyGammaFittedError[100] = {0.0};
  float energResolutionFitted[100],energyResolutionFittedError[100],countsFitted[100],countsFittedError[100];
    
  float energyResolutionPercent[100];
  float efficiencyPercent[100];
  sprintf(inputFile,"./FittingResults/BeamPipe1mmPb6pRes.txt");
  sprintf(outputFile,"./FittingResults/EfficiencyBeamPipe1mmPb6pRes.txt");
  fDataInput = fopen(inputFile,"r");
  int distinguisher[10] = {0.};
 
  int counterBeamEnergyGammasFitted[6]={0};
  int i = 0;
  fscanf(fDataInput,"%s",dummyText);
  
  while(!feof(fDataInput))  {
    fscanf(fDataInput,"%f %f %f %f %f %f %f %f %f",&energyBeam[i],&energyGamma[i],&detector,&energyGammaFitted[i],
           &energyGammaFittedError[i],&energResolutionFitted[i],&energyResolutionFittedError[i],
           &countsFitted[i],&countsFittedError[i]);
    
    cout<<"I am working on line: "<<i+2<<endl;
    energyResolutionPercent[i] = 100. * 2.35 *  energResolutionFitted[i]/energyGammaFitted[i];
    efficiencyPercent[i] = 100. * countsFitted[i]/ binning / eventSimulated;

    for(int j =0;j<=4;j++)  {
      int dummy = 100 + 50 * i;
      if(energyBeam[i]==dummy) counterBeamEnergyGammasFitted[j]++;
    }
    i++;
  }
  fclose(fDataInput);

  
  TCanvas *fCanvas;  
  sprintf(dummyText,"Canvas");
  fCanvas = new TCanvas(dummyText,dummyText,500,800);
  fCanvas->cd();
  fCanvas->SetFillColor(0);
  //fCanvas->SetLogy();
  //fCanvas->SetLogx();
  fCanvas->SetBorderSize(0);
  fCanvas->SetBorderMode(0);
  fCanvas->SetFrameBorderMode(0);
  fCanvas->SetFrameFillColor(0);

  // ------------>Primitives in pad: c_up
  TPad *c_up = new TPad("c_up", "c_up",0.0,0.5,1.0,1.0);
  c_up->Draw();
  c_up->cd();
  c_up->SetFillColor(0);
  c_up->SetBorderSize(0);
  c_up->SetRightMargin(0.03);
  c_up->SetTopMargin(0.15);
  c_up->SetBottomMargin(0.01); 

  int style;
  // The Efficiency
  TGraph *grEfficiency[5];
  for(int k=0;k<5;k++)  {
    //int dummy = 100 + 50 * i;
    //sprintf(dummytext,"");
    style = 21;
    grEfficiency[k] = new TGraph(counterBeamEnergyGammasFitted[k]);
    grEfficiency[k]->SetMarkerStyle(style+k);
    //grEfficiency[k]->SetMarkerColor(1+k);
  }

// The Resolution
  TGraph *grResolution[5];
  for(int k=0;k<5;k++)  {
    //int dummy = 100 + 50 * i;
    //sprintf(dummytext,"");
    style = 21;
    grResolution[k] = new TGraph(counterBeamEnergyGammasFitted[k]);
    grResolution[k]->SetMarkerStyle(style+k);
    //grResolution[k]->SetMarkerColor(1+k);
  }
  int a[5]={0};
  for(int j=0;j<i;j++)  {
    for(int k =0;k<=3;k++)  {
      int dummy = 100 + 50 * k;
      if(energyBeam[j]==dummy) {
        grEfficiency[k]->SetPoint(a[k],energyGamma[j],efficiencyPercent[j]); 
        grResolution[k]->SetPoint(a[k],energyGamma[j],energyResolutionPercent[j]); 
        a[k]++;
      }
    }
  }

  TH1F *GraphTotal = new TH1F("Efficiency","",2000,0.0,3250.0);
  GraphTotal->SetStats(0);
  GraphTotal->SetLabelFont(42);
  GraphTotal->GetXaxis()->SetTitle("");
  GraphTotal->GetXaxis()->SetNdivisions(205);
  GraphTotal->GetXaxis()->SetRangeUser(100,3250);
  GraphTotal->GetXaxis()->SetLabelFont(42);
  GraphTotal->GetXaxis()->SetTitleFont(42);
  GraphTotal->GetYaxis()->SetRangeUser(5,60);
  GraphTotal->GetYaxis()->SetNdivisions(605);
  GraphTotal->GetYaxis()->SetTitle("FEP Efficiency / %");
  GraphTotal->GetYaxis()->SetTitleOffset(1.0);
  GraphTotal->GetYaxis()->SetLabelFont(42);
  GraphTotal->GetYaxis()->SetTitleFont(42);
 
  for(int k =0;k<=3;k++)  grEfficiency[k]->SetHistogram(GraphTotal);

  grEfficiency[0]->Draw("ACP");
  for(int k =1;k<=3;k++)  grEfficiency[k]->Draw("CP");
  
  TLegend *leg = new TLegend(0.75,0.52,0.97,0.85,NULL,"brNDC");
  leg->SetTextAlign(32);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  
   for(int k =0;k<=3;k++) {
     int dummy = 100 + 50 * k;
     sprintf(dummyText,"%i MeV/u",dummy);
     TLegendEntry *entry=leg->AddEntry("NULL",dummyText,"cpl");
     style = 21;
     entry->SetMarkerStyle(style+k);
     //entry->SetMarkerColor(1+k);
     entry->SetTextColor(1);
     entry->SetLineColor(1);
     entry->SetTextFont(42);
     entry->SetLineWidth(1.);
   }
  leg->Draw();

  c_up->Modified();
  //__________________________________________________________________________________________________
  fCanvas->cd();
  // ------------>Primitives in pad: c_down
  TPad *c_down = new TPad("c_down", "c_down",0.0,0.0,1.0,0.5);
  c_down->Draw();
  c_down->cd();
  c_down->SetFillColor(0);
  c_down->SetBorderSize(0);
  c_down->SetRightMargin(0.03);
  c_down->SetTopMargin(0.01);
  c_down->SetBottomMargin(0.15); 

  TH1F *GraphTotalRes = new TH1F("Energy Resolution","",2000,0.0,3250.0);
  GraphTotalRes->SetStats(0);
  GraphTotalRes->SetLabelFont(42);
  GraphTotalRes->GetXaxis()->SetNdivisions(205);
  GraphTotalRes->GetXaxis()->SetTitle("Energy (keV)");
  GraphTotalRes->GetXaxis()->SetRangeUser(100,3250);
  GraphTotalRes->GetXaxis()->SetLabelFont(42);
  GraphTotalRes->GetXaxis()->SetTitleFont(42);
  GraphTotalRes->GetYaxis()->SetRangeUser(5,15.5);
  GraphTotalRes->GetYaxis()->SetNdivisions(605);
  GraphTotalRes->GetYaxis()->SetTitle("Energy Resolution / %");  
  GraphTotalRes->GetYaxis()->SetTitleOffset(1.00);
  GraphTotalRes->GetYaxis()->SetLabelFont(42);
  GraphTotalRes->GetYaxis()->SetTitleFont(42);

  for(int k =0;k<=3;k++)  grResolution[k]->SetHistogram(GraphTotalRes);
  
  grResolution[0]->Draw("ACP");
  for(int k =1;k<=3;k++)  grResolution[k]->Draw("CP");

  TLegend *leg2 = new TLegend(0.75,0.66,0.97,0.99,NULL,"brNDC");
  leg2->SetTextAlign(32);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  leg2->SetLineColor(1);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  
   for(int k =0;k<=3;k++) {
     int dummy = 100 + 50 * k;
     sprintf(dummyText,"%i MeV/u",dummy);
     TLegendEntry *entry=leg2->AddEntry("NULL",dummyText,"cpl");
     style = 21;
     entry->SetMarkerStyle(style+k);
     //entry->SetMarkerColor(1+k);
     entry->SetTextColor(1);
     entry->SetLineColor(1);
     entry->SetTextFont(42);
     entry->SetLineWidth(1.);
   }
  leg2->Draw();

  
  //c_down->Modified();

  // Writing the resolution to file:
  FILE* fEfficiencyOut = fopen(outputFile,"w");
  fprintf(fEfficiencyOut,"BeamEnergy GammaEnergy Resolution Efficiency \n");
  i=0;
  for(i=0;i<100;i++)  {
    if(energyGamma[i]>0)
      fprintf(fEfficiencyOut,"%f %f %f %f \n",energyBeam[i],energyGamma[i],energyResolutionPercent[i],efficiencyPercent[i]);
  }
  fclose(fEfficiencyOut);  
}
