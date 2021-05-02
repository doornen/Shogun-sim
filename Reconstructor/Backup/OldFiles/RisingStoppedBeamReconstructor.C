#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TH1.h"
#include "TH2F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TText.h"
#include "TMath.h"
#include "TF1.h"

float GetTOFGammaNs(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                    float target_x,float target_y,float target_z)
{
  float vx1 = gamma_det_x - target_x;                         //cout<<"vx1: "<<vx1<<endl;
  float vy1 = gamma_det_y - target_y;                         //cout<<"vy1: "<<vy1<<endl;
  float vz1 = gamma_det_z - target_z;                         //cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  return  l1/29.97925;
} 

void GetGammaEnergiesSorted(float gammaEnergy[])
{
  float dummyGammaEnergy[10];
  int rankOfGamma[10]={0};
  for(int i=0;i<10;i++) dummyGammaEnergy[i] = gammaEnergy[i];
  for(int i=0;i<10;i++)
  {
    if(gammaEnergy[i]>0)
    for(int j=0;j<10;j++){if(gammaEnergy[i]>dummyGammaEnergy[j]) rankOfGamma[i]++;}
  }
  for(int i=0;i<10;i++) gammaEnergy[i] = dummyGammaEnergy[rankOfGamma[i]];
}
void RisingStoppedBeamReconstructor()
{ 
  //float eventNumberCheck = 0.;
  const int numberOfClusterDetectors = 15;
  float eventNumberPrevious = 0.;
  float gammaEnergy[20];
  char root_in[200];
  char root_out[200];
  char temp[200];
  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
  float beta_average    = 0.4295;        //The average value at the moment of decay for the correction
  float beta_mean_tof   = 0.4295;//The mean beta from Tof
  float decay_z         = 0.0;  //The lifetime of the excited state moves the average decay point along the beam axis
  //------------------------------------------------------------------------------------------------------
  //Reading the input file, which can change the values
  FILE *fin = fopen("RisingStoppedBeamReconstructor.in","r");
  while(!feof(fin))
  {
    fscanf(fin,"%s ",temp); 
    if(strcmp(temp,"INPUTFILE")==0)
    {
      fscanf(fin,"%s",&root_in);
      printf("%s %s \n",temp,root_in);
    }
    else if(strcmp(temp,"OUTPUTFILE")==0)  
    {
      fscanf(fin,"%s ",&root_out); 
      printf("%s %s \n",temp,root_out);
    }
    else if(strcmp(temp,"SPECTRABINANDRANGE")==0)
    {
      fscanf(fin,"%i %f %f",&numBin,&firstBin,&lastBin);
      printf("%s %i %f %f\n",temp, numBin,firstBin,lastBin);
    }
    else if(strcmp(temp,"BETADOPPLERAVERAGE")==0)
    {
      fscanf(fin,"%f",&beta_average); 
      printf("%s %f \n",temp,beta_average);
    }
    else if(strcmp(temp,"BETATOFAVERAGE")==0)
    {
      fscanf(fin,"%f",&beta_mean_tof); 
      printf("%s %f \n",temp,beta_mean_tof);
    }
    else if(strcmp(temp,"DECAYPOSITIONZ")==0)
    {
      fscanf(fin,"%f",&decay_z); 
      printf("%s %f \n",temp,decay_z);
    }
    else if(strcmp(temp,"STATISTICSREDUCTIONFACTOR")==0)
    {
      fscanf(fin,"%f",&reduction_factor); 
      printf("%s %f \n",temp,reduction_factor);
      if(reduction_factor<1) return;
    }
    else if(strcmp(temp,"END")==0) break;
    else 
    {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      abort();
    }
  }

//-----------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------
// These variables are the same as for the eventbuilder 
float theta_low, theta_high; //Angle covered
int mass,z,charge, mass_f, z_f, charge_f;
float thickness; 
float evnum;
float decayIDevnum;
float x_p0,y_p0,z_p0,x_pvb,y_pvb,z_pvb,x_pva,y_pva,z_pva;
float energy_p;
float beta_b,beta_r,beta_a,halflife,decay_time_after_interaction;
float x_g0,y_g0,z_g0,x_gv,y_gv,z_gv;
float e_rest,e_doppler;
float theta_gamma_rest, theta_gamma_lab;
float energy_vertex_stored;
int kind_target;
//***********************************************
float total_event_num;//energy_notcor;
int gamma_det_type;  // 1 for dali, 2 for grape, 3 for sgt 
float pos_det_x,pos_det_y,pos_det_z;     //reconstructed HI position from detector after target
float ver_rec_x, ver_rec_y, ver_rec_z;
float beta_rec;
//-----------------------------------------------
//For the Cluster array
int cluster_flag[numberOfClusterDetectors][7];
float cluster_x[numberOfClusterDetectors][7],cluster_y[numberOfClusterDetectors][7],cluster_z[numberOfClusterDetectors][7];
float cluster_energy_not_cor[numberOfClusterDetectors][7];
float cluster_time[numberOfClusterDetectors][7] = {{0.}};
//Input of the resolutions:
int cluster_en_res_opt;
float cluster_en_res[2];
float cluster_time_res[2];
float pos_det_at_target_res;
float pos_det_after_target_res; // Pos Resolution in mm FWHM!
float beta_res; //Resolution of beta in FWHM!
//End, same variables
//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------

  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();
  
  TH1F *h_number_of_observed_decays = new TH1F("h_number_of_observed_decays","h_number_of_observed_decays",100,-0.5,99.5);
  TH2F *h_decayId_vs_energy_sum = new TH2F("h_decayId_vs_energy_sum","h_decayId_vs_energy_sum",10,0,10,1000,0,5000);
  
  TH1F *h_cluster_crystal_mult;
  TH1F *h_cluster_single;
  TH1F *h_cluster_energy_sum;
  TH2F *h_cluster_gamma_gamma[10];

  // Crystal multiplicity,ok
  h_cluster_crystal_mult = new TH1F("cluster_Crystal_Mult","cluster_Crystal_mult",100,0,100);
  //Single hits, OK
  h_cluster_single = new TH1F("cluster_single","cluster_single",numBin,firstBin,lastBin);
  // All, OK, This includes addback events.
  h_cluster_energy_sum = new TH1F("h_cluster_energy_sum","h_cluster_energy_sum",numBin,firstBin,lastBin);
  //---------------------------------------
  for(int i=0;i<10;i++)
  {
    sprintf(temp,"h_cluster_gamma_gamma[%i]",i);
    h_cluster_gamma_gamma[i]= new TH2F(temp,temp,numBin,firstBin,lastBin,numBin,firstBin,lastBin);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile *infile = new TFile(root_in,"READ");
  TTree *t = (TTree*)infile->Get("ObservedEvents");

  t->SetBranchAddress("EventNumber",&evnum);
  t->SetBranchAddress("DecayIDOfEventNumber",&decayIDevnum);
  t->SetBranchAddress("X_vertex",&x_p0);    //X Position at the fragmentation point
  t->SetBranchAddress("Y_vertex",&y_p0);    //Y Position at the fragmentation point
  t->SetBranchAddress("Z_vertex",&z_p0);    //Z Position at the fragmentation point
  t->SetBranchAddress("XV_projectile_before_target",&x_pvb);  // Normalized Vector of beam before the target
  t->SetBranchAddress("YV_projectile_before_target",&y_pvb);
  t->SetBranchAddress("ZV_projectile_before_target",&z_pvb);
  t->SetBranchAddress("XV_projectile_after_target",&x_pva);   // Normalized Vector of beam after the target
  t->SetBranchAddress("YV_projectile_after_target",&y_pva);
  t->SetBranchAddress("ZV_projectile_after_target",&z_pva);
  t->SetBranchAddress("Energy_projectile",&energy_p); // energy of beam before the target in MeV/u
  t->SetBranchAddress("Beta_before_target",&beta_b);    // Beta before the target
  t->SetBranchAddress("Beta_real",&beta_r);             // Beta during deexcitation	
  t->SetBranchAddress("Beta_after_target",&beta_a);     // Beta After Target
  t->SetBranchAddress("Halflife",&halflife);          // Halflife
  t->SetBranchAddress("Decay_Time_after_interaction",&decay_time_after_interaction);
  t->SetBranchAddress("X_vertex_gamma",&x_g0);            //X Position at the gamma emmittance point
  t->SetBranchAddress("Y_vertex_gamma",&y_g0);            //Y Position at the gamma emmittance point
  t->SetBranchAddress("Z_vertex_gamma",&z_g0);            //Z Position at the gamma emmittance point
  t->SetBranchAddress("XV_gamma",&x_gv);		         // Gamma vector
  t->SetBranchAddress("YV_gamma",&y_gv);
  t->SetBranchAddress("ZV_gamma",&z_gv);
  t->SetBranchAddress("E_gamma_rest",&e_rest);		 //Energy at rest
  t->SetBranchAddress("E_gamma_Doppler",&e_doppler);//Theta of doppler boosted gamma
  t->SetBranchAddress("Theta_gamma_rest",&theta_gamma_rest);
  t->SetBranchAddress("Theta_gamma_lab",&theta_gamma_lab);
  //Energy of fragment at fragmentation
  t->SetBranchAddress("Energy_Vertex",&energy_vertex_stored); 
  //So far identical to the ROOT-Tree of step one
  //The new stuff:
  //Setting the type of the gamma detector
  t->SetBranchAddress("Gamma_det_type",&gamma_det_type);
  t->SetBranchAddress("Cluster_flag",cluster_flag);//ID of which det recorded a gamma in the event
  t->SetBranchAddress("Cluster_energy_not_cor",cluster_energy_not_cor);
  t->SetBranchAddress("Cluster_time",cluster_time);
  
  //-------------------------------
  //Position detector after targetpos_det_after_target_res
  t->SetBranchAddress("pos_det_rec_X",&pos_det_x); // reconstructed position from a det after the secondary target
  t->SetBranchAddress("pos_det_rec_Y",&pos_det_y); //Including the resolution of the detectors
  t->SetBranchAddress("pos_det_rec_Z",&pos_det_z);
  t->SetBranchAddress("Vertex_Reconstructed_X",&ver_rec_x); // reconstructed vertex postion, 
  t->SetBranchAddress("Vertex_Reconstructed_Y",&ver_rec_y); //including resolution
  t->SetBranchAddress("Vertex_Reconstructed_Z",&ver_rec_z);    
  //reconstructed beta 
  t->SetBranchAddress("Beta_Reconstructed",&beta_rec); // reconstructed beta, including resolution

  //------------------------------------------------------------------------------------------------------
  //Going to the Header file, which has the information of constant values.
  infile->cd();
  TTree *tHeader = (TTree*)infile->Get("Header");
  tHeader->SetBranchAddress("Mass",&mass);                      // Beam mass
  tHeader->SetBranchAddress("Z",&z);                               // Beam z	
  tHeader->SetBranchAddress("Charge",&charge);                // Beam charge	
  tHeader->SetBranchAddress("Mass_Fragment",&mass_f);         // Fragment mass
  tHeader->SetBranchAddress("Z_Fragment",&z_f);	                 // Fragment z
  tHeader->SetBranchAddress("Charge_Fragment",&charge_f);   // Fragment charge
  tHeader->SetBranchAddress("TargetKind",&kind_target);
  tHeader->SetBranchAddress("TargetThicknessCM",&thickness);
  tHeader->SetBranchAddress("Total_Event_Number",&total_event_num);
  tHeader->SetBranchAddress("ThetaLow",&theta_low);
  tHeader->SetBranchAddress("ThetaHigh",&theta_high);
  tHeader->SetBranchAddress("cluster_en_res_opt",&cluster_en_res_opt);
  tHeader->SetBranchAddress("cluster_en_res",cluster_en_res);
  tHeader->SetBranchAddress("cluster_time_res",cluster_time_res);
  tHeader->SetBranchAddress("cluster_x",cluster_x);
  tHeader->SetBranchAddress("cluster_y",cluster_y);
  tHeader->SetBranchAddress("cluster_z",cluster_z);
  //beta, and position resolution 
  tHeader->SetBranchAddress("Beta_Resolution",&beta_res); 
  tHeader->SetBranchAddress("Pos_Det_at_Target_res",&pos_det_at_target_res); 
  tHeader->SetBranchAddress("Pos_Det_After_Target_Res",&pos_det_after_target_res);
  tHeader->GetEntry(0);
  //-----------------------------------------------------------------------------------------------------------
  //Finished reading the Root file

  Int_t nentries = (Int_t)(t->GetEntries()/reduction_factor);
  cout <<"Entries: "<< nentries << endl;
  //Variables that are determined event by event:
  int cluster_crystal_mult=0;
  float cluster_energy_sum = 0.;
  float cluster_energy_max = -999.9;
  
  float cluster_x_thatsit;
  float cluster_y_thatsit;
  float cluster_z_thatsit;
  float gamma_tof;

  int observedDecayNumber = 0;
  //---------------------------------------------------------------------------------------------------------
  //Gates for the analysis:
  float energy_gate_lower[20]={120.,240.,315.,695.};
  float energy_gate_upper[20]={140.,265.,340.,750.};

  cout<<"Starting the analysis"<<endl;
  float dummyEnergy = 0.; 
  for(int iii=0;iii<nentries;iii++)
  {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    t->GetEntry(iii);

    //-------------------------------------------------------------
    // Starting with the cluster analysis:
    //---------------------------------------------------------------------------------------------------------
    for(int nnn=0;nnn<numberOfClusterDetectors;nnn++)
    {
      for(int j=0;j<7;j++)
      {
        if(cluster_flag[nnn][j]==1.0) 
        {
          cluster_crystal_mult++;
          //summed energy
          cluster_energy_sum = cluster_energy_sum+cluster_energy_not_cor[nnn][j];
          // getting the MI
          if(cluster_energy_not_cor[nnn][j]>cluster_energy_max)
          {
            cluster_energy_max=cluster_energy_not_cor[nnn][j];
            cluster_x_thatsit=cluster_x[nnn][j];
            cluster_y_thatsit=cluster_y[nnn][j];
            cluster_z_thatsit=cluster_z[nnn][j];
          }
        }
      }
    }
    h_cluster_crystal_mult->Fill(cluster_crystal_mult);
    if(cluster_crystal_mult>0)
    {
      //Performing the Doppler-correction from the angles
      //of the detector that registered the highest gamma-ray energy.
      dummyEnergy = cluster_energy_sum;
      if(evnum == eventNumberPrevious) gammaEnergy[observedDecayNumber+1] = dummyEnergy;
      
      h_cluster_energy_sum->Fill(dummyEnergy);
      h_decayId_vs_energy_sum->Fill(decayIDevnum,dummyEnergy);
      if(cluster_crystal_mult==1) h_cluster_single->Fill(dummyEnergy);
    }
    // The variables have to be reset.
    cluster_crystal_mult=0;
    cluster_energy_sum=0.0;
    cluster_energy_max = -999.9;
    observedDecayNumber++;
    if(evnum != eventNumberPrevious)
    { 
      h_number_of_observed_decays->Fill(observedDecayNumber);
      //Filling the gamma-gamma spectra:
      //Sorting in the right order
      for(int i =1;i<observedDecayNumber;i++)
      {
        if(gammaEnergy[0]>=gammaEnergy[i])
        {
	   h_cluster_gamma_gamma[i]->Fill(gammaEnergy[0],gammaEnergy[i]);
           h_cluster_gamma_gamma[0]->Fill(gammaEnergy[0],gammaEnergy[i]);
        }
        else
        {
          h_cluster_gamma_gamma[i]->Fill(gammaEnergy[i],gammaEnergy[0]);
          h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[0]);
        }
        if(i>1)
        {
          if(gammaEnergy[1]>=gammaEnergy[i])
            h_cluster_gamma_gamma[0]->Fill(gammaEnergy[1],gammaEnergy[i]);
          else h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[1]);
        }
        if(i>2)
        { 
          if(gammaEnergy[2]>=gammaEnergy[i])
            h_cluster_gamma_gamma[0]->Fill(gammaEnergy[2],gammaEnergy[i]);
          else h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[2]);
        }
        if(i>3)
          if(gammaEnergy[3]>=gammaEnergy[i])
            h_cluster_gamma_gamma[0]->Fill(gammaEnergy[3],gammaEnergy[i]);
          else h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[3]);
      } 
      observedDecayNumber = 0;
      for(int i =0;i<20;i++)
      {
        gammaEnergy[i]=0.;
      }
      gammaEnergy[observedDecayNumber] = dummyEnergy;
    }
    eventNumberPrevious = evnum;
  }
  //---------------------------------------------------------------
  //Writing into the root file:
  //---------------------------------------------------------------
  //cluster-----------------------------------------------------------
  //---------------------------------------------------------------
  h_number_of_observed_decays->Write();
  h_cluster_crystal_mult->Write();
  h_cluster_single->Write();
  h_cluster_energy_sum->Write();
  for(int i=0;i<10;i++)h_cluster_gamma_gamma[i]->Write();
}
