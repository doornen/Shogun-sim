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

void ShowEfficiencyAndResolutionStopped(int caseNumber,int eventSimulated,int binning)  {  
  char dummyText[1000];
  float energyBeam[100],energyGamma[100],energyGammaFitted[100],energyGammaFittedError[100];
  float energResolutionFitted[100],energyResolutionFittedError[100],countsFitted[100],countsFittedError[100];
    
  float energyResolutionPercent[100];
  float efficiencyPercent[100];
  sprintf(dummyText,"./FittingResults/Case%iStopped.txt",caseNumber);
  fDataInput = fopen(dummyText,"r");
  int distinguisher[10] = {0.};

  int counterBeamEnergyGammasFitted[6]={0};
  int i = 0;
  fscanf(fDataInput,"%s",dummyText);
  int dummyInput;
  while(!feof(fDataInput)&&i<6)  {
    fscanf(fDataInput,"%i %f %f %f %f %f %f %f %f",&dummyInput,&energyBeam[i],&energyGamma[i],&energyGammaFitted[i],
           &energyGammaFittedError[i],&energResolutionFitted[i],&energyResolutionFittedError[i],
           &countsFitted[i],&countsFittedError[i]);
    
    cout<<"I am working on line: "<<i+2<<endl;
    energyResolutionPercent[i] = 100. * 2.35 *  energResolutionFitted[i]/energyGammaFitted[i];
    efficiencyPercent[i] = 100. * countsFitted[i]/ binning / eventSimulated;
    cout<<"Resolution, Efficiency: "<<energyResolutionPercent[i]<<"; "<<efficiencyPercent[i]<<endl;


    for(int j =0;j<=4;j++)  {
      int dummy = 0;
      if(energyBeam[i]==dummy) counterBeamEnergyGammasFitted[j]++;
    }
    i++;
  }
  fclose(fDataInput);

  TCanvas *fCanvas;  
  sprintf(dummyText,"Canvas");
  fCanvas = new TCanvas(dummyText,dummyText,600,400);
  fCanvas->cd();
  fCanvas->SetFillColor(0);
  //fCanvas->SetLogy();
  //fCanvas->SetLogx();
  fCanvas->SetBorderSize(0);
  fCanvas->SetBorderMode(0);
  fCanvas->SetFrameBorderMode(0);
  fCanvas->SetFrameFillColor(0);
  
  int style;
  // The Efficiency
  TGraph *grEfficiency[1];
  for(int k=0;k<1;k++)  {
    //int dummy = 100 + 50 * i;
    //sprintf(dummytext,"");
    style = 21;
    grEfficiency[k] = new TGraph(counterBeamEnergyGammasFitted[k]);
    grEfficiency[k]->SetMarkerStyle(style+k);
    //grEfficiency[k]->SetMarkerColor(1+k);
  }
  
  // The Efficiency
  TGraph *grResolution[1];
  for(int k=0;k<1;k++)  {
    //int dummy = 100 + 50 * i;
    //sprintf(dummytext,"");
    style = 21;
    grResolution[k] = new TGraph(counterBeamEnergyGammasFitted[k]);
    grResolution[k]->SetMarkerStyle(style+k);
    //grResolution[k]->SetMarkerColor(1+k);
  }

  int a[5]={0};
  for(int j=0;j<i;j++)  {
    for(int k =0;k<=0;k++)  {
      int dummy = 0;
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
  GraphTotal->GetXaxis()->SetTitle("Energy [keV]");
  GraphTotal->GetYaxis()->SetRangeUser(10,100);
  GraphTotal->GetXaxis()->SetRangeUser(100,3250);
  GraphTotal->GetYaxis()->SetTitle("FEP Efficiency / %");
  GraphTotal->GetYaxis()->SetTitleOffset(1.0);
  GraphTotal->GetXaxis()->SetLabelFont(42);
  GraphTotal->GetXaxis()->SetTitleFont(42);
  GraphTotal->GetYaxis()->SetLabelFont(42);
  GraphTotal->GetYaxis()->SetTitleFont(42);
 
  for(int k =0;k<=0;k++)  grEfficiency[k]->SetHistogram(GraphTotal);

  grEfficiency[0]->Draw("ACP");
  //for(int k =1;k<=0;k++)  grEfficiency[k]->Draw("CP");
  
  TLegend *leg = new TLegend(0.70,0.54,0.88,0.88,NULL,"brNDC");
  leg->SetTextAlign(32);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  
  //for(int k =0;k<=2;k++) {
  //   int dummy = 100 + 50 * k;
  //   sprintf(dummyText,"%i MeV/u",dummy);
  //   TLegendEntry *entry=leg->AddEntry("NULL",dummyText,"cpl");
  //   style = 21;
  //   entry->SetMarkerStyle(style+k);
  //   //entry->SetMarkerColor(1+k);
  //   entry->SetTextColor(1);
  //   entry->SetLineColor(1);
  //   entry->SetTextFont(42);
  //   entry->SetLineWidth(1.);
  // }
  // leg->Draw();
   
   fCanvas->Modified();
   //__________________________________________________________________________________________________
   TCanvas *fCanvas2;  
   sprintf(dummyText,"Canvas2");
   fCanvas2 = new TCanvas(dummyText,dummyText,600,400);
   fCanvas2->cd();
   fCanvas2->SetFillColor(0);
   //fCanvas->SetLogy();
   //fCanvas->SetLogx();
   fCanvas2->SetBorderSize(0);
   fCanvas2->SetBorderMode(0);
   fCanvas2->SetFrameBorderMode(0);
  fCanvas2->SetFrameFillColor(0);

  TH1F *GraphTotalRes = new TH1F("Energy Resolution","",2000,0.0,3250.0);
  GraphTotalRes->SetStats(0);
  GraphTotalRes->SetLabelFont(42);
  GraphTotalRes->GetXaxis()->SetTitle("Energy [keV]");
  GraphTotalRes->GetYaxis()->SetRangeUser(2,12.);
  GraphTotalRes->GetXaxis()->SetRangeUser(100,3250);
  GraphTotalRes->GetYaxis()->SetTitle("Energy Resolution (FWHM) / %");
  GraphTotalRes->GetYaxis()->SetTitleOffset(1.);
  GraphTotalRes->GetXaxis()->SetLabelFont(42);
  GraphTotalRes->GetXaxis()->SetTitleFont(42);
  GraphTotalRes->GetYaxis()->SetLabelFont(42);
  GraphTotalRes->GetYaxis()->SetTitleFont(42);
  
  for(int k =0;k<=0;k++)  grResolution[k]->SetHistogram(GraphTotalRes);
  grResolution[0]->Draw("ACP");
  //for(int k =1;k<=2;k++)  grResolution[k]->Draw("CP");

  //leg->Draw();
  fCanvas2->Modified();
}
