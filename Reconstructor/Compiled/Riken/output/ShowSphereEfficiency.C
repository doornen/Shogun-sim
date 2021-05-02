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

void ShowSphereEfficiency()  {  
  char dummyText[1000];
  int thickness[100];
  float energyBeam[100],energyGamma[100];
  float counts[100];
  
  float efficiencyPercent[100];
  sprintf(dummyText,"./FittingResults/SphereLaBr3.txt");
  fDataInput = fopen(dummyText,"r");
  
  
  fscanf(fDataInput,"%s",dummyText);
  int dummyInput;
  int i;
  while(!feof(fDataInput)&&i<72)  {
    fscanf(fDataInput,"%i %f %f %f",&thickness[i],&energyBeam[i],&energyGamma[i],&counts[i]);
    
    cout<<"I am working on line: "<<i+2<<endl;
    efficiencyPercent[i] = counts[i]/200000*100;
    cout<<"Efficiency: "<<efficiencyPercent[i]<<endl;

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
  TGraph *grEfficiency[6];
  for(int k=0;k<6;k++)  {
    //int dummy = 100 + 50 * i;
    //sprintf(dummytext,"");
    style = 21;
    grEfficiency[k] = new TGraph(12);
    grEfficiency[k]->SetMarkerStyle(style+k);
    //grEfficiency[k]->SetMarkerColor(1+k);
  }
  
  int a[6]={0};
  for(int j=0;j<i;j++)  {
    for(int k =0;k<3;k++)  {
      int dummy = 100 + 50 * k;
      for(int l =0;l<2;l++)  {
        int dummy2 = 8 + 2*l; 
        if(energyBeam[j]==dummy &&  thickness[j]==dummy2) {
          grEfficiency[k+3*l]->SetPoint(a[k+3*l],energyGamma[j],efficiencyPercent[j]); 
          a[k+3*l]++;
        }
      }
    }
  }
  TH1F *GraphTotal = new TH1F("Efficiency","",2000,0.0,3250.0);
  GraphTotal->SetStats(0);
  GraphTotal->SetLabelFont(42);
  GraphTotal->GetXaxis()->SetTitle("Energy [keV]");
  GraphTotal->GetYaxis()->SetRangeUser(30,100);
  GraphTotal->GetXaxis()->SetRangeUser(100,3250);
  GraphTotal->GetYaxis()->SetTitle("FEP Efficiency / %");
  GraphTotal->GetYaxis()->SetTitleOffset(1.0);
  GraphTotal->GetXaxis()->SetLabelFont(42);
  GraphTotal->GetXaxis()->SetTitleFont(42);
  GraphTotal->GetYaxis()->SetLabelFont(42);
  GraphTotal->GetYaxis()->SetTitleFont(42);
 
  for(int k =0;k<=5;k++)  grEfficiency[k]->SetHistogram(GraphTotal);

  grEfficiency[0]->Draw("ACP");
  for(int k =1;k<=5;k++)  grEfficiency[k]->Draw("CP");
  
  TLegend *leg = new TLegend(0.60,0.54,0.88,0.88,NULL,"brNDC");
  leg->SetTextAlign(32);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  
 
  for(int k =0;k<=2;k++) {
    for(int j=0;j<=1;j++) {
      int dummy = 100 + 50 * k;
      int dummy2 = 8 + 2*j;
      sprintf(dummyText,"%i MeV/u, %i cm",dummy,dummy2);
      TLegendEntry *entry=leg->AddEntry("NULL",dummyText,"cpl");
      style = 21;
      entry->SetMarkerStyle(style+k+3*j);
      //entry->SetMarkerColor(1+k);
      entry->SetTextColor(1);
      entry->SetLineColor(1);
      entry->SetTextFont(42);
      entry->SetLineWidth(1.);
    }
  }
  leg->Draw();

  fCanvas->Modified();
}
