//
//#include "TFile.h"
//#include "TH1F.h"
//#include "TMath.h"
//#include "TCanvas.h"
//#include "TF1.h"
//#include "TGraphErrors.h"
//#include "TGraph.h"
//#include "TPad.h"
//#include "TLegend.h"
//#include <TLatex.h>
//#include <iostream.h>
//#include <test.h>

//void test()
{
  //gStyle->Reset();
  gStyle->SetOptStat(kFALSE);
  
  gROOT->ProcessLine( ".L fit.h" );
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  
  //****************************************************************************
  // The simulated peaks
  TFile *peak1 = new TFile("./134sn_be_128pd_925kev.root");
 
  //The experimental data:
  char dummyText[1000];
  sprintf(dummyText,"128pd_spectrum.csv");
  //xsprintf(dummyText,"126pd.dat");

  FILE *fDataInput= fopen(dummyText,"r");
  
  TH1 *h1 = new TH1F("128pd","",600,0.0,3000.0);
  int i=0;

  float x[600];
  float y[600];

  while(!feof(fDataInput)&&i<600) {
    fscanf(fDataInput,"%f %f",&x[i],&y[i]);
    h1->SetBinContent(i,(int)y[i]);
    i++;
  }
  //for(i=1;i<600;i++)
  //  cout<<"i :"<<i<<" y[i]: "<<y[i]<<endl;

  fclose(fDataInput);
  h1->Rebin(10);

  //  TFile *OutFile = new TFile( "whole_doppler.root", "RECREATE" );
  TCanvas *fCanvas=new TCanvas("Canvas","Responsive function",500,350);
  fCanvas->SetBorderSize(0);
  fCanvas->SetBorderMode(0);
  fCanvas->SetFrameBorderMode(0);
  fCanvas->SetFrameFillColor(0);
  fCanvas->SetBottomMargin(0.15);
  fCanvas->SetLeftMargin(0.15);
  fCanvas->cd();
  
  //****************************************************************************
  TH1F *hpeak1=new TH1F("hpeak1","temp1",160,0,4000);
  
  hpeak1 = (TH1F*)peak1->Get("doppler;1");
  hpeak1->Rebin(2);

  h1->GetXaxis()->SetRangeUser(300,2470);
  h1->GetXaxis()->SetNdivisions(205);
  h1->GetYaxis()->SetTitle("Counts / 50 keV");
  h1->GetYaxis()->SetRangeUser(1,15);
  h1->GetYaxis()->SetNdivisions(304);
  
  h1->GetXaxis()->SetTitleOffset(0.9);  
  h1->GetYaxis()->SetTitleOffset(0.7);
  
  h1->GetXaxis()->SetTitleFont(132);
  h1->GetYaxis()->SetTitleFont(132);
   
  h1->GetXaxis()->SetTitleSize(0.08);
  h1->GetYaxis()->SetTitleSize(0.08);
  
  h1->GetXaxis()->SetLabelSize(0.08);
  h1->GetYaxis()->SetLabelSize(0.08);
  
  //How to get error bar for each bin
  h1->SetDefaultSumw2(kTRUE);
  
  h1->SetTitle("");
  hpeak1->Draw("");

  //*****************************************************************************
  TGraph *peak1g = new TGraph(hpeak1);
  
  //******************Function Definition****************************************
  const Double_t fitmin=200.;
  const Double_t fitmax=1500.;
  TF1 *whole = new TF1("whole",ex_respf,fitmin,fitmax,3);
  whole->SetParameters(0.0012,3.32-2.0,-0.00125);
  whole->SetParLimits(0,0.0012,0.0012);
  whole->SetParLimits(1,3.32-2.0,3.32-2.0);
  whole->SetParLimits(2,-0.00125,-0.00125);
  whole->SetLineColor(1);
  whole->SetLineWidth(2);
  whole->SetNpx(200);
  h1->Fit(whole,"R");
  whole->Draw("same");
  
  //*****************************************************************************
  TF1 * peak1f= new TF1( "peak1f", resp1,fitmin,fitmax,1);
  peak1f->SetParameter(0,whole->GetParameter(0));
  peak1f->SetLineColor(2);
  peak1f->SetLineWidth(2);
  peak1f->SetLineStyle(2);
  peak1f->Draw("same");
  
  TF1 * expon= new TF1( "expon",expf ,fitmin,fitmax,2);
  expon->SetParameters(whole->GetParameter(1),whole->GetParameter(2));
  expon->SetLineColor(4);
  expon->SetLineWidth( 2 );
  expon->SetLineStyle(2);
  expon->Draw("same");
  TLegend *legend = new TLegend(0.5,0.6,0.9,0.9);
  legend->SetFillColor(0);
  legend->SetTextFont(132);
  legend->SetTextSize(0.06);
  legend->AddEntry(h1,"Experiment","lpe");
  legend->AddEntry(peak1f,"Simulation 925 keV","l");
  //legend->AddEntry(peak2f,"777 keV","l");
  //legend->AddEntry(peak3f,"504keV","l");
  legend->AddEntry(expon,"Background Fit","l");
  legend->AddEntry(whole,"Sim. + Background","l");
  legend->Draw();
  h1->Draw("e SAME");
  fCanvas->SetLogy(0);
  
  TLatex *   tex = new TLatex(600,12,"^{128}Pd");
  tex->SetTextFont(132);
  tex->SetTextSize(0.08);
  tex->SetLineWidth(2);
  tex->Draw();
  /*
  TLatex *   tex = new TLatex(400,53,"a)");
  tex->SetTextFont(132);
  tex->SetTextSize(0.07);
  tex->SetLineWidth(2);
  tex->Draw();
  */

  //peak1->Close();
}

