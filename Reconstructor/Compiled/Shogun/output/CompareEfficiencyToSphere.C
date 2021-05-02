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
FILE *fDataInputSphere;

void CompareEfficiencyToSphere(int caseNumber,int eventSimulated,int binning)  {  
  char dummyText[1000];
  float energyBeam[100],energyGamma[100],energyGammaFitted[100],energyGammaFittedError[100];
  float energResolutionFitted[100],energyResolutionFittedError[100],countsFitted[100],countsFittedError[100];
    
  float energyResolutionPercent[100];
  float efficiencyPercent[100];
  sprintf(dummyText,"./FittingResults/Case%i.txt",caseNumber);
  fDataInput = fopen(dummyText,"r");
  int distinguisher[10] = {0.};

  int counterBeamEnergyGammasFitted[6]={0};
  int i = 0;
  fscanf(fDataInput,"%s",dummyText);
  int dummyInput;
  while(!feof(fDataInput)&&i<37)  {
    fscanf(fDataInput,"%i %f %f %f %f %f %f %f %f",&dummyInput,&energyBeam[i],&energyGamma[i],&energyGammaFitted[i],
           &energyGammaFittedError[i],&energResolutionFitted[i],&energyResolutionFittedError[i],
           &countsFitted[i],&countsFittedError[i]);
    
    cout<<"I am working on line: "<<i+2<<endl;
    energyResolutionPercent[i] = 100. * 2.35 *  energResolutionFitted[i]/energyGammaFitted[i];
    efficiencyPercent[i] = 100. * countsFitted[i]/ binning / eventSimulated;
    cout<<"Resolution, Efficiency: "<<energyResolutionPercent[i]<<"; "<<efficiencyPercent[i]<<endl;


    for(int j =0;j<=4;j++)  {
      int dummy = 100 + 50 * i;
      if(energyBeam[i]==dummy) counterBeamEnergyGammasFitted[j]++;
    }
    i++;
  }
  fclose(fDataInput);

  //----------------------------------------------------------------
  //Reading the sphere data:
  int spherethickness[100];
  float sphereenergybeam[100],sphereenergygamma[100];
  float spherecounts[100];
  float sphereefficiencypercent[100];

  sprintf(dummyText,"../../Riken/output/FittingResults/SphereLaBr3.txt");
  fDataInputSphere = fopen(dummyText,"r");
  fscanf(fDataInputSphere,"%s",dummyText);
 
  int m=0;
  while(!feof(fDataInputSphere)&&m<72)  {
    fscanf(fDataInputSphere,"%i %f %f %f",&spherethickness[m],&sphereenergybeam[m],&sphereenergygamma[m],&spherecounts[m]);
    
    cout<<"I am working on line: "<<m+2<<endl;
    sphereefficiencypercent[m] = spherecounts[m]/200000*100;
    cout<<"Sphere Efficiency: "<<sphereefficiencypercent[m]<<endl;
    m++;
  }
  fclose(fDataInputSphere);

  //____________________________________________________________________


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
  TGraph *grEfficiency[3];
  for(int k=0;k<3;k++)  {
    //int dummy = 100 + 50 * i;
    //sprintf(dummytext,"");
    style = 21;
    grEfficiency[k] = new TGraph(counterBeamEnergyGammasFitted[k]);
    grEfficiency[k]->SetMarkerStyle(style+k);
    //grEfficiency[k]->SetMarkerColor(1+k);
  }
  
  int a[5]={0};
  for(int j=0;j<i;j++)  {
    for(int k =0;k<=2;k++)  {
      int dummy = 100 + 50 * k;
      if(energyBeam[j]==dummy) {
        
        
        for(int o=0;o<m;o++)  {
          for(int n =0;n<3;n++)  {
            int dummy3 = 100 + 50 * n;
            if(sphereenergybeam[o]== energyBeam[j]&& 
               spherethickness[o]==8  &&
               sphereenergygamma[o] == energyGamma[j])
              {
                grEfficiency[k]->SetPoint(a[k],energyGamma[j],(efficiencyPercent[j]/
                                                               sphereefficiencypercent[o])); 
                a[k]++;      
                cout<<"Relative efficiency: "<<efficiencyPercent[j]/sphereefficiencypercent[o]<<endl;
              
            }
          }
        }   
      }
    }
  }
  
  TH1F *GraphTotal = new TH1F("Efficiency","",2000,0.0,3250.0);
  GraphTotal->SetStats(0);
  GraphTotal->SetLabelFont(42);
  GraphTotal->GetXaxis()->SetTitle("Energy [keV]");
  GraphTotal->GetYaxis()->SetRangeUser(0.4,0.8);
  GraphTotal->GetXaxis()->SetRangeUser(100,3250);
  GraphTotal->GetYaxis()->SetTitle("SHOGUN / SPHERE Efficiency");
  GraphTotal->GetYaxis()->SetTitleOffset(1.0);
  GraphTotal->GetXaxis()->SetLabelFont(42);
  GraphTotal->GetXaxis()->SetTitleFont(42);
  GraphTotal->GetYaxis()->SetLabelFont(42);
  GraphTotal->GetYaxis()->SetTitleFont(42);
 
  for(int k =0;k<=2;k++)  grEfficiency[k]->SetHistogram(GraphTotal);
  grEfficiency[0]->Draw("ACP");
  for(int k =1;k<=2;k++)  grEfficiency[k]->Draw("CP");
  
  TLegend *leg = new TLegend(0.70,0.54,0.88,0.88,NULL,"brNDC");
  leg->SetTextAlign(32);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  
   for(int k =0;k<=2;k++) {
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

  fCanvas->Modified();
}
