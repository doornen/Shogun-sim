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

TH1F *h_cluster_doppler_beta_var[20];
TH1F *h_cluster_doppler_outer_ring_beta_var[20];
TH1F *h_cluster_doppler_ring_beta_var[3][20];
TH1F *h_miniball_doppler_beta_var[20];
TH1F *h_miniball_doppler_ring_beta_var[2][20];
TFile *rootfile;
TFile *infile;

//_________________________________________________________________________________________________
float GetDopplerCorrectedEnergy(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                                float target_x,float target_y,float target_z,
                                float pos_det_x, float pos_det_y, float pos_det_z, 
                                float gamma_energy, float beta_rec,float beta_mean_tof,float beta_average,float decay_z)
{
  float vx1 = gamma_det_x - target_x;                                      //cout<<"vx1: "<<vx1<<endl;
  float vy1 = gamma_det_y - target_y;                                      //cout<<"vy1: "<<vy1<<endl;
  float vz1 = gamma_det_z - target_z - decay_z;                            //cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  float vx2 = pos_det_x - target_x;                                        //cout<<"vx2: "<<vx2<<endl;
  float vy2 = pos_det_y - target_y;                                        //cout<<"vy2: "<<vy2<<endl;
  float vz2 = pos_det_z - target_z - decay_z;                              //cout<<"vy2: "<<vz2<<endl;
 
  double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
          
  // The theta angle:
  double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);             //cout<<"Theta lab: "<<theta_lab<<endl;
  // the beta value:
  float beta                     = beta_average + (beta_rec-beta_mean_tof); //cout<<"beta: "<<beta<<endl;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

//_________________________________________________________________________________________________
// Starting the main program://////////////////////////////////////////////////////////////////////
//_________________________________________________________________________________________________
void RisingFastBeamReconstructorBetaVar36K()
{ 
  //float eventNumberCheck = 0.;
  const int numberOfClusterDetectors  = 15;
  const int numberOfMiniballDetectors = 8;
  int clusterInclude                  = 1;
  int miniballInclude                 = 1;
  float eventNumberPrevious           = 0.;
  char root_in[200];
  char root_out[200];
  char temp[200];
  int numBin             = 1000;
  float firstBin         = -2.0;
  float lastBin          = 3998.0;
  float reduction_factor = 1.;
  float beta_mean_tof    = 0.56311;        // The mean beta from Tof
  float decayZ[20]          = 
{0.0,0.03068,0.05838,0.08597,0.1134,0.1409,0.1683,0.1957,0.2231,0.2505,0.2779,0.3053,0.3326,0.360,0.3874,0.4147,0.4421,0.4695,0.4968,0.5242};
  float betaVar[20] = {0.5459,0.5434,0.5415,0.5400,0.5387,0.5376,0.5368,0.5361,0.5355,0.5350,0.5346,0.5342,0.5338,0.5336,0.5333,0.5331,0.5329,0.5327,0.5325,0.5323};
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
float total_event_num;
int gamma_det_type;                      // 1 for cluster, 2 for miniball, 3 for hector 
float pos_det_x,pos_det_y,pos_det_z;     // reconstructed HI position from detector after target
float ver_rec_x, ver_rec_y, ver_rec_z;
float beta_rec;
//-----------------------------------------------
// For the Cluster array
int cluster_flag[numberOfClusterDetectors][7];
float cluster_x[numberOfClusterDetectors][7],cluster_y[numberOfClusterDetectors][7],cluster_z[numberOfClusterDetectors][7];
float cluster_energy_not_cor[numberOfClusterDetectors][7] = {{0.}};
float cluster_time[numberOfClusterDetectors][7] = {{0.}};
// For the Miniball array
int miniball_flag[numberOfMiniballDetectors][3][7];
float miniball_x[numberOfMiniballDetectors][3][7],miniball_y[numberOfMiniballDetectors][3][7],miniball_z[numberOfMiniballDetectors][3][7];
float miniball_energy_not_cor[numberOfMiniballDetectors][3][7] = {{{0.}}};
float miniball_time[numberOfMiniballDetectors][3][7] = {{{0.}}};
//-----------------------------------------------
// Input of the resolutions:
int cluster_en_res_opt;
float cluster_en_res[2];
float cluster_time_res[2];
int miniball_en_res_opt;
float miniball_en_res[2];
float miniball_time_res[2];
//-----------------------------------------------
float pos_det_at_target_res;
float pos_det_after_target_res; // Pos Resolution in mm FWHM!
float beta_res; //Resolution of beta in FWHM!

// End, same variables
//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------

  for(int i = 798; i<818;i++)
  {
    for(int j = 0; j<20;j++)
    {
      sprintf(root_out,"/home/doornen/GEANT4/g4.9.0.p01/SimulationResults/Reconstructor/RisingFast/36k/out_events_built_37ca_36k_195mev_%ikev_%ips.root",i,j);

      rootfile = new TFile(root_out,"RECREATE");
      rootfile->cd();
    
      if(clusterInclude==1)
      {
        for(int k=0;k<20;k++)
        {
          sprintf(temp,"h_cluster_doppler_beta_var[%i]",k);
          h_cluster_doppler_beta_var[k] = new TH1F(temp,"",numBin,firstBin,lastBin);
        }
        for(int k=0;k<20;k++)
        {
          sprintf(temp,"h_cluster_doppler_outer_ring_beta_var[%i]",k);
          h_cluster_doppler_outer_ring_beta_var[k] = new TH1F(temp,"",numBin,firstBin,lastBin);
        }
        for(int k=0;k<3;k++)
        {
          for(int l=0;l<20;l++)
          {
            sprintf(temp,"h_cluster_doppler_ring_beta_var[%i][%i]",k,l);
            h_cluster_doppler_ring_beta_var[k][l] = new TH1F(temp,"",numBin,firstBin,lastBin);
          }
        }
      }
      //_______________________________________________________________________________________________ 

      if(miniballInclude==1)
      {
        for(int k=0;k<20;k++)
        {
          sprintf(temp,"h_miniball_doppler_beta_var[%i]",k);
          h_miniball_doppler_beta_var[k] = new TH1F(temp,"",numBin,firstBin,lastBin);
        }
        for(int k=0;k<2;k++)
        {
          for(int l=0;l<20;l++)
          {
            sprintf(temp,"h_miniball_doppler_ring_beta_var[%i][%i]",k,l);
            h_miniball_doppler_ring_beta_var[k][l] = new TH1F(temp,"",numBin,firstBin,lastBin);
          }
        }
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////
      sprintf(root_in,"/media/BOOTCAMP/SimulationResults/Builder/36k/events_built_37ca_36k_195mev_%ikev_%ips.root",i,j);

      infile = new TFile(root_in,"READ");
      TTree *t = (TTree*)infile->Get("ObservedEvents");

      t->SetBranchAddress("EventNumber",&evnum);
      t->SetBranchAddress("DecayIDOfEventNumber",&decayIDevnum);
      t->SetBranchAddress("X_vertex",&x_p0);                      // Position at the fragmentation point
      t->SetBranchAddress("Y_vertex",&y_p0);    
      t->SetBranchAddress("Z_vertex",&z_p0);                      
      t->SetBranchAddress("XV_projectile_before_target",&x_pvb);  // Normalized Vector of beam before the target
      t->SetBranchAddress("YV_projectile_before_target",&y_pvb);
      t->SetBranchAddress("ZV_projectile_before_target",&z_pvb);
      t->SetBranchAddress("XV_projectile_after_target",&x_pva);   // Normalized Vector of beam after the target
      t->SetBranchAddress("YV_projectile_after_target",&y_pva);
      t->SetBranchAddress("ZV_projectile_after_target",&z_pva);
      t->SetBranchAddress("Energy_projectile",&energy_p);         // Energy of beam before the target in MeV/u
      t->SetBranchAddress("Beta_before_target",&beta_b);          // Beta of beam before the target
      t->SetBranchAddress("Beta_real",&beta_r);                   // Beta of fragment during deexcitation	
      t->SetBranchAddress("Beta_after_target",&beta_a);           // Beta of fragment After Target
      t->SetBranchAddress("Halflife",&halflife);                  // Halflife
      t->SetBranchAddress("Decay_Time_after_interaction",&decay_time_after_interaction);
      t->SetBranchAddress("X_vertex_gamma",&x_g0);                // Position at the gamma emmittance point
      t->SetBranchAddress("Y_vertex_gamma",&y_g0);            
      t->SetBranchAddress("Z_vertex_gamma",&z_g0);            
      t->SetBranchAddress("XV_gamma",&x_gv);		          // Gamma vector
      t->SetBranchAddress("YV_gamma",&y_gv);
      t->SetBranchAddress("ZV_gamma",&z_gv);
      t->SetBranchAddress("E_gamma_rest",&e_rest);	          // Energy at rest
      t->SetBranchAddress("E_gamma_Doppler",&e_doppler);          // Theta of doppler boosted gamma
      t->SetBranchAddress("Theta_gamma_rest",&theta_gamma_rest);
      t->SetBranchAddress("Theta_gamma_lab",&theta_gamma_lab);
      t->SetBranchAddress("Energy_Vertex",&energy_vertex_stored); // Energy of fragment at fragmentation 
      t->SetBranchAddress("Gamma_det_type",&gamma_det_type);
      if(clusterInclude==1)
      {
        t->SetBranchAddress("Cluster_flag",cluster_flag);         // The Cluster detectors
        t->SetBranchAddress("Cluster_energy_not_cor",cluster_energy_not_cor);
        t->SetBranchAddress("Cluster_time",cluster_time);
      }
      if(miniballInclude==1)
      {
        t->SetBranchAddress("miniball_flag",miniball_flag);       // The Miniball detectors
        t->SetBranchAddress("miniball_energy_not_cor",miniball_energy_not_cor);
        t->SetBranchAddress("miniball_time",miniball_time);
      }
      t->SetBranchAddress("pos_det_rec_X",&pos_det_x);            // Reconstructed position from a det after the secondary target
      t->SetBranchAddress("pos_det_rec_Y",&pos_det_y);            // including the resolution of the detectors
      t->SetBranchAddress("pos_det_rec_Z",&pos_det_z);
      t->SetBranchAddress("Vertex_Reconstructed_X",&ver_rec_x);   // reconstructed vertex postion, 
      t->SetBranchAddress("Vertex_Reconstructed_Y",&ver_rec_y);   // including resolution
      t->SetBranchAddress("Vertex_Reconstructed_Z",&ver_rec_z);    
      t->SetBranchAddress("Beta_Reconstructed",&beta_rec);        // reconstructed beta, including resolution

      //------------------------------------------------------------------------------------------------------
      //Going to the Header file, which has the information of constant values.
      infile->cd();
      TTree *tHeader = (TTree*)infile->Get("Header");
      tHeader->SetBranchAddress("Mass",&mass);                    // Beam mass
      tHeader->SetBranchAddress("Z",&z);                          // Beam z	
      tHeader->SetBranchAddress("Charge",&charge);                // Beam charge	
      tHeader->SetBranchAddress("Mass_Fragment",&mass_f);         // Fragment mass
      tHeader->SetBranchAddress("Z_Fragment",&z_f);	              // Fragment z
      tHeader->SetBranchAddress("Charge_Fragment",&charge_f);     // Fragment charge
      tHeader->SetBranchAddress("TargetKind",&kind_target);
      tHeader->SetBranchAddress("TargetThicknessCM",&thickness);
      tHeader->SetBranchAddress("Total_Event_Number",&total_event_num);
      tHeader->SetBranchAddress("ThetaLow",&theta_low);
      tHeader->SetBranchAddress("ThetaHigh",&theta_high);
      if(clusterInclude==1)
      {
        tHeader->SetBranchAddress("cluster_en_res_opt",&cluster_en_res_opt);
        tHeader->SetBranchAddress("cluster_en_res",cluster_en_res);
        tHeader->SetBranchAddress("cluster_time_res",cluster_time_res);
        tHeader->SetBranchAddress("cluster_x",cluster_x);
        tHeader->SetBranchAddress("cluster_y",cluster_y);
        tHeader->SetBranchAddress("cluster_z",cluster_z);
      }
      if(miniballInclude==1)
      {
        tHeader->SetBranchAddress("miniball_en_res_opt",&miniball_en_res_opt);
        tHeader->SetBranchAddress("miniball_en_res",miniball_en_res);
        tHeader->SetBranchAddress("miniball_time_res",miniball_time_res);
        tHeader->SetBranchAddress("miniball_x",miniball_x);
        tHeader->SetBranchAddress("miniball_y",miniball_y);
        tHeader->SetBranchAddress("miniball_z",miniball_z);
      }
      tHeader->SetBranchAddress("Beta_Resolution",&beta_res); 
      tHeader->SetBranchAddress("Pos_Det_at_Target_res",&pos_det_at_target_res); 
      tHeader->SetBranchAddress("Pos_Det_After_Target_Res",&pos_det_after_target_res);
      tHeader->GetEntry(0);
      //-----------------------------------------------------------------------------------------------------------
      // Finished reading the Root file

      Int_t nentries = (Int_t)(t->GetEntries()/reduction_factor);
      cout <<"Entries: "<< nentries << endl;
      // Variables that are determined event by event:
      int cluster_crystal_mult  = 0;
      int cluster_detector_mult = 0;
      float cluster_energy_sum  = 0.;
      float cluster_energy_max  = -999.9;
      bool cluster_fired[numberOfClusterDetectors] = {false};
  
      int miniball_segment_mult  = 0;
      int miniball_crystal_mult  = 0;
      int miniball_detector_mult = 0; 
      float miniball_energy_sum  = 0.;
      float miniball_energy_max  = -999.9;
      bool miniball_fired[numberOfMiniballDetectors] = {false};

      float cluster_x_thatsit  = 0.;
      float cluster_y_thatsit  = 0.;
      float cluster_z_thatsit  = 0.;
      float miniball_x_thatsit = 0.;
      float miniball_y_thatsit = 0.;
      float miniball_z_thatsit = 0.;

      cout<<"Starting the analysis"<<endl;

      float dummyEnergy = 0.; 
      for(int iii=0;iii<nentries;iii++)
      {
        if((iii+1)%10000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
        t->GetEntry(iii);

        //-------------------------------------------------------------
        // Starting with the cluster analysis:-------------------------
        //-------------------------------------------------------------
        int dummy = -1;
        if(clusterInclude==1)
        {
          for(int nnn=0;nnn<numberOfClusterDetectors;nnn++)
          {
            for(int jjj=0;jjj<7;jjj++)
            {
              if(cluster_flag[nnn][jjj]==1.0) 
              {
                if(dummy!=nnn) cluster_detector_mult++;
                dummy  = nnn;
                cluster_crystal_mult++;
                cluster_fired[nnn] = true;
                //summed energy
                cluster_energy_sum = cluster_energy_sum+cluster_energy_not_cor[nnn][jjj];
                // getting the MI
                if(cluster_energy_not_cor[nnn][jjj]>cluster_energy_max)
                {
                  cluster_energy_max = cluster_energy_not_cor[nnn][jjj];
                  cluster_x_thatsit  = cluster_x[nnn][jjj];
                  cluster_y_thatsit  = cluster_y[nnn][jjj];
                  cluster_z_thatsit  = cluster_z[nnn][jjj];
                }
              }
            }
          }
          //____________________________________________________________________________________________
          if(cluster_crystal_mult==1 && cluster_detector_mult==1) // Checking if the gamma-ray energy was released within the same cluster!!
          {
            for(int kkk=0;kkk<20;kkk++)
            {
              dummyEnergy = GetDopplerCorrectedEnergy(cluster_x_thatsit,cluster_y_thatsit,cluster_z_thatsit,
              ver_rec_x,ver_rec_y,ver_rec_z,pos_det_x,pos_det_y,pos_det_z,cluster_energy_sum,beta_rec,beta_mean_tof,betaVar[kkk],decayZ[kkk]);

              h_cluster_doppler_beta_var[kkk]->Fill(dummyEnergy);
              for(int lll=0;lll<15;lll++)
              {
                if(lll<5 && cluster_fired[lll] == true)
                {
                  h_cluster_doppler_ring_beta_var[0][kkk]->Fill(dummyEnergy);
                }
                if(lll>4 && lll<10 && cluster_fired[lll] == true)
                {
                  h_cluster_doppler_ring_beta_var[1][kkk]->Fill(dummyEnergy);
                  h_cluster_doppler_outer_ring_beta_var[kkk]->Fill(dummyEnergy);
                }
                if(lll>10 && cluster_fired[lll] == true)
                {
                  h_cluster_doppler_ring_beta_var[2][kkk]->Fill(dummyEnergy);
                  h_cluster_doppler_outer_ring_beta_var[kkk]->Fill(dummyEnergy);
                }
              }
            }
          }
        }
        //____________________________________________________________________________________________
        // Clusters Finished. Analyzing Miniball
        dummy = -1;
        if(miniballInclude==1)
        {
          for(int nnn=0;nnn<numberOfMiniballDetectors;nnn++)
          {
            for(int jjj=0;jjj<3;jjj++)
            {
              for(int kkk=0;kkk<7;kkk++)
              {
                if(miniball_flag[nnn][jjj][kkk]==1.0 && kkk==0) 
                {
                  if(dummy!=nnn) miniball_detector_mult++;
                  miniball_crystal_mult++;
                  miniball_fired[nnn] = true;
                  //summed energy
                  miniball_energy_sum = miniball_energy_sum+miniball_energy_not_cor[nnn][jjj][0];
                  dummy  = nnn;
                }
                if(miniball_flag[nnn][jjj][kkk]==1.0 && kkk!=0)
                {  
                  miniball_segment_mult++;
                  // getting the MI
                  if(miniball_energy_not_cor[nnn][jjj][kkk]>miniball_energy_max)
                  {
                    miniball_energy_max = miniball_energy_not_cor[nnn][jjj][kkk];
                    miniball_x_thatsit  = miniball_x[nnn][jjj][kkk];
                    miniball_y_thatsit  = miniball_y[nnn][jjj][kkk];
                    miniball_z_thatsit  = miniball_z[nnn][jjj][kkk];
                  }
                }
              }
            }
          }
          if(miniball_crystal_mult==1 && miniball_detector_mult==1) // Checking if the gamma-ray energy was released within the same cluster!!
          {
          // Performing the Doppler-correction from the angles
          // of the detector that registered the highest gamma-ray energy.
            for(int kkk=0;kkk<20;kkk++)
            {
              dummyEnergy = GetDopplerCorrectedEnergy(miniball_x_thatsit,miniball_y_thatsit,miniball_z_thatsit,
              ver_rec_x,ver_rec_y,ver_rec_z,pos_det_x,pos_det_y,pos_det_z,miniball_energy_sum,beta_rec,beta_mean_tof,betaVar[kkk],decayZ[kkk]);

              h_miniball_doppler_beta_var[kkk]->Fill(dummyEnergy);
           
              if(miniball_fired[0]==true || miniball_fired[2]==true || miniball_fired[5]==true || miniball_fired[7]==true)
              {
                h_miniball_doppler_ring_beta_var[0][kkk]->Fill(dummyEnergy);
              }
              if(miniball_fired[1]==true || miniball_fired[3]==true || miniball_fired[4]==true || miniball_fired[6]==true)
              {
                h_miniball_doppler_ring_beta_var[1][kkk]->Fill(dummyEnergy);
              }
            }
          }
        }     
        //_____________________________________________________________________________________________
        // The variables have to be reset.
        cluster_crystal_mult  = 0;
        cluster_detector_mult = 0;
        cluster_energy_sum    = 0.0;
        cluster_energy_max    = -999.9;
        for(int nnn=0;nnn<numberOfClusterDetectors;nnn++)
        {
          cluster_fired[nnn]  = false;
        }

        miniball_segment_mult  = 0;
        miniball_crystal_mult  = 0;
        miniball_detector_mult = 0;
        miniball_energy_sum    = 0.0;
        miniball_energy_max    = -999.9;
        for(int nnn=0;nnn<numberOfMiniballDetectors;nnn++)
        {
          miniball_fired[nnn] = false;
        }
        eventNumberPrevious = evnum;
      }
      //--------------------------------------------------------------------------
      // Writing into the root file:----------------------------------------------
      //--------------------------------------------------------------------------
      rootfile->cd();
      // Cluster------------------------------------------------------------------
      if(clusterInclude==1)
      {
        for(int m=0;m<20;m++)
        {
          h_cluster_doppler_beta_var[m]                                 ->Write();
          h_cluster_doppler_outer_ring_beta_var[m]                      ->Write();
          for(int nn=0;nn<3;nn++)h_cluster_doppler_ring_beta_var[nn][m] ->Write();
        }
      }
      // Miniball-----------------------------------------------------------------
      if(miniballInclude==1)
      {
        for(int m=0;m<20;m++)
        {
          h_miniball_doppler_beta_var[m]                                ->Write();
          for(int nn=0;nn<2;nn++)h_miniball_doppler_ring_beta_var[nn][m]->Write();
        }
      }
      rootfile->Close("");
      infile->Close("");
    }
  }
}
