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
#include "TTree.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TText.h"
#include "TMath.h"
#include "TF1.h"


void LifetimeAnalyser()
{ 

//char root_in[200]="/home/doornen/GEANT4/SimulationResults/Builder/36mg/events_built_36mg_150mev_0005mgBe_660_kev_50000fs_borrelbindenergy.root";
char root_in[200]="/media/BOOTCAMP/SimulationResults/Builder/lifetime/events_built_100mev_1000kev_90ps.root";
//char root_in[200]="/media/BOOTCAMP/SimulationResults/Builder/34cl/lifetime/events_built_34cl_195mev_461kev_19ps.root";
//char root_in[200]="/media/BOOTCAMP/SimulationResults/Builder/36k/lifetime/events_built_36k_195mev_810kev_18ps.root";
//char root_out[200]="/home/doornen/GEANT4/SimulationResults/Reconstructor/36mg/out_events_built_36mg_150mev_0005mgBe_660_kev_50000fs_borrelbindenergy.root";
char root_out[200]="/media/BOOTCAMP/SimulationResults/Reconstructor/lifetime/out_events_built_100mev_1000kev_90ps.root";
//char root_out[200]="/media/BOOTCAMP/SimulationResults/Reconstructor/34cl/lifetime/out_events_built_34cl_195mev_461kev_19ps.root";
//char root_out[200]="/media/BOOTCAMP/SimulationResults/Reconstructor/36k/lifetime/out_events_built_36k_195mev_810kev_18ps.root";
  char name[100];
  int  numBin = 1000;
  float firstBin   = -2;
  float lastBin   = 3998;
  //Beta for 34cl at 195 meV in steps of 1 ps
  //double beta_var[21]={0.5459,0.544,0.5424,0.541,0.54,0.5391,0.5384,0.5378,0.5373,0.5369,
  //0.5365,0.5362,0.536,0.5357,0.5355,0.5353,0.5351,0.535,0.5349,0.5347,0.5346}; 
  //Beta for 36k at 195mev in steps of 1 ps

// double beta_var[41]={0.546443,0.544162,0.542243,0.540674,0.53938,
//  0.538325,0.537461,0.536745,0.53615,0.535643,0.535207,0.534831,0.5345,0.534208,0.533947,0.533715,0.533508,0.53332,0.53315,0.533,
//0.5328,0.5327,0.5326,0.5325,0.5324,0.5323,0.5322,0.5321,0.532,0.532,0.5319,0.5318,0.5318,0.5317,0.5317,0.5316,0.5316,0.5315,0.5315,0.5315,0.5314
//};
  
  //double beta_mean_tof = 0.563082;//The mean beta from Tof
  //The lifetime of the excited state moves the average deacay point in zet
  //for 34cl in steps of 1ps
  //double decay_move[21] ={0.0,0.02787,0.056,0.08413,0.1123,0.1404,0.1685,0.1967,0.2248,
  //0.2529,0.281,0.3092,0.3373,0.3654,0.3936,0.4217,0.4498,0.4779,0.5061,0.5342,0.5623};
  //for 36k in steps of 1ps

  //double decay_move[41]={0.0,0.2605,0.5431,0.8256,1.108,1.391,1.673,1.956,2.238,
  //2.521,2.804,3.086,3.369,3.651,3.934,4.216,4.499,4.782,5.064,5.347,
  //5.774,6.054,6.287,6.574,6.86,7.146,7.432,7.719,8.005,8.291,8.577,8.864,9.152,9.436,
   //9.722,10.01,10.29,10.58,10.87,11.15,11.42};

  // Beta for 100 MeV 
  double beta_var[21]; 
  for(int i=0;i<21;i++)
  {  
     //100 MeV
     beta_var[i]=0.4295;
     //150 MeV
     //beta_var[i]=0.5081;
     //200 MeV
     //beta_var[i]=0.5677;
     //250 MeV
     //beta_var[i]=0.6152;
  }  
  double beta_mean_tof = 0.4295;//The mean beta from Tof
  double decay_move[21]; 
  for(int i=0;i<21;i++)
  {
     // 100 MeV
     decay_move[i]=i*0.1042;
     // 150 MeV
     //decay_move[i]=i*0.1284;
     // 200 MeV
     //decay_move[i]=i*0.1497;
     // 250 MeV
     //decay_move[i]=i*0.169;
  } 
  
  float evnum, kind;
  float x_p0, y_p0,z_p0; //Position at fragmentation point
  float x_pvb,y_pvb,z_pvb; //Normalized beam vector before the target
  float x_pva,y_pva,z_pva; //Normalized beam vector after the target
  float energy_p,beta_b,beta_r,beta_a,halflife;

  float x_g0,y_g0,z_g0;  //Position at gamma emmittance point
  float x_gv,y_gv,z_gv;  //Gamma-vector
  float e_rest,e_doppler, theta_gamma_rest, theta_gamma_lab;
  int z, mass, charge;   //Beam specifications
  float  energy_vertex_stored; //energy at fragmentation point
  int z_f, mass_f, charge_f; //fragment specifications
  float thickness,kind_target_f; //target
  float total_event_num;   //Number of total events
  float energy_notcor;  //Observed energy in detector (any)

  float x_gamma_det,y_gamma_det,z_gamma_det; //Position of gamma detector
  float gamma_det_type;  //Type of gamma detector
  //Cluster
  float cluster_res;
  float id;
  float id_cluster;
  float id_crystal_in_cluster;
 //miniball 
  float miniball_res;
  float miniball_cluster_thatsit;
  float id_miniball;
  float id_cluster_in_miniball;
  float id_crystal_in_miniball;
  float id_segment_in_miniball;
  //cate
  float cate_pos_res;
  float x_cate,y_cate,z_cate;  //position at Cate including resolution;
  float x_ver_rec, y_ver_rec, z_ver_rec; //reconstruction of decay at target point
  float cate_si_res;//resolution of si
  float delta_e_cate;//measured delta E
  float cate_csi_res;//resolution of csi
  float e_cate;//measured E
  //Tof
  float beta_res;//resolution tof
  float beta_rec;//reconstructed beta, including resolution
  float total_gamma_num;

  float theta_low,theta_high,solidangle_factor;
  TFile *rootfile = new TFile(root_out,"RECREATE");
  //Creating cluster spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //Single hit
  rootfile->cd();
  char tempname[100];
  TH1F *h_cluster_dopplercorrected_single_beta_var[21];
  TH1F *h_cluster_dopplercorrected_single_inner_ring_beta_var[21];
  TH1F *h_cluster_dopplercorrected_single_middle_ring_beta_var[21];
  TH1F *h_cluster_dopplercorrected_single_outer_ring_beta_var[21];
  TH1F *h_cluster_x_position_of_gamma_at_moment_of_decay = new TH1F("Cluster Gamma X position",
                                "Cluster Gamma X position",500,-5.0,20.0); 
  for(int i=0;i<21;i++)
  {
    rootfile->cd();
    sprintf(tempname,"Cluster_DopplerCorrected_SingleHit_Beta_Var_%i",i);
    h_cluster_dopplercorrected_single_beta_var[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
    rootfile->cd();
    sprintf(tempname,"Cluster_DopplerCorrected_SingleHit_Inner_Ring_Beta_Var_%i",i);
    h_cluster_dopplercorrected_single_inner_ring_beta_var[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
    rootfile->cd();
    sprintf(tempname,"Cluster_DopplerCorrected_SingleHit_Middle_Ring_Beta_Var_%i",i);
    h_cluster_dopplercorrected_single_middle_ring_beta_var[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
    rootfile->cd();
    sprintf(tempname,"Cluster_DopplerCorrected_SingleHit_Outer_Ring_Beta_Var_%i",i);
    h_cluster_dopplercorrected_single_outer_ring_beta_var[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }   
  
  
  //Not doppler corrected
  rootfile->cd();
  TH1F *h_cluster_notdopplercorrected = new TH1F("Cluster_NotDopplerCorrected",
                                "Cluster_NotDopplerCorrected",numBin,firstBin,lastBin); 
  //////////////////////////////////////////////
  //Creating Miniball spectra:

  ////////////////////////////////////////////
   rootfile->cd(); 
  TH1F *h_miniball_dopplercorrected_single_beta_var[21];
  TH1F *h_miniball_dopplercorrected_single_inner_ring_beta_var[21];
  TH1F *h_miniball_dopplercorrected_single_outer_ring_beta_var[21];
  TH1F *h_miniball_x_position_of_gamma_at_moment_of_decay = new TH1F("Miniball Gamma X position",
                                "Miniball Gamma X position",500,-5.0,20);
  
  for(int i=0;i<21;i++)
  {
    rootfile->cd();
    sprintf(tempname,"Miniball_DopplerCorrected_SingleHit_Beta_Var_%i",i);
    h_miniball_dopplercorrected_single_beta_var[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
    rootfile->cd();
    sprintf(tempname,"Miniball_DopplerCorrected_SingleHit_Inner_Ring_Beta_Var_%i",i);
    h_miniball_dopplercorrected_single_inner_ring_beta_var[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
    rootfile->cd();
    sprintf(tempname,"Miniball_DopplerCorrected_SingleHit_Outer_Ring_Beta_Var_%i",i);
    h_miniball_dopplercorrected_single_outer_ring_beta_var[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }   
  //Not doppler corrected
  rootfile->cd();
  TH1F *h_miniball_notdopplercorrected = new TH1F("Miniball_NotDopplerCorrected",
                                "Miniball_NotDopplerCorrected",numBin,firstBin,lastBin); 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile *infile = new TFile(root_in,"READ");
  infile->cd();
  TTree *t = (TTree*)infile->Get("ObservedEvents");
  t->SetBranchAddress("EventNumber",&evnum);
  t->SetBranchAddress("EventTag",&kind);
  t->SetBranchAddress("X_vertex",&x_p0);
  t->SetBranchAddress("Y_vertex",&y_p0); 
  t->SetBranchAddress("Z_vertex",&z_p0);
  t->SetBranchAddress("XV_projectile_before_target",&x_pvb);
  t->SetBranchAddress("YV_projectile_before_target",&y_pvb);
  t->SetBranchAddress("ZV_projectile_before_target",&z_pvb);
  t->SetBranchAddress("XV_projectile_after_target",&x_pva);
  t->SetBranchAddress("YV_projectile_after_target",&y_pva);
  t->SetBranchAddress("ZV_projectile_after_target",&z_pva);
  t->SetBranchAddress("Energy_projectile",&energy_p);
  t->SetBranchAddress("Beta_before_target",&beta_b);
  t->SetBranchAddress("Beta_Reconstructed",&beta_rec);  //This is the recontructed beta which inlcudes the tof resolution!!!!
  t->SetBranchAddress("Beta_real",&beta_r);
  t->SetBranchAddress("Beta_after_target",&beta_a);
  t->SetBranchAddress("Halflife",&halflife);
  t->SetBranchAddress("X_vertex_gamma",&x_g0);
  t->SetBranchAddress("Y_vertex_gamma",&y_g0);
  t->SetBranchAddress("Z_vertex_gamma",&z_g0);
  t->SetBranchAddress("XV_gamma",&x_gv);
  t->SetBranchAddress("YV_gamma",&y_gv);
  t->SetBranchAddress("ZV_gamma",&z_gv);
  t->SetBranchAddress("E_gamma_rest",&e_rest);
  t->SetBranchAddress("E_gamma_Doppler",&e_doppler);
  t->SetBranchAddress("Theta_gamma_rest",&theta_gamma_rest);
  t->SetBranchAddress("Theta_gamma_lab",&theta_gamma_lab);
  t->SetBranchAddress("Mass",&mass);
  t->SetBranchAddress("Z",&z);
  t->SetBranchAddress("Charge",&charge);
  t->SetBranchAddress("Mass_Fragment",&mass_f);
  t->SetBranchAddress("Z_Fragment",&z_f);
  t->SetBranchAddress("Charge_Fragment",&charge_f);
  t->SetBranchAddress("Energy_Vertex",&energy_vertex_stored);
  // energy released and position fo any of the gamma detectors
  t->SetBranchAddress("Observed_Energy",&energy_notcor);
  t->SetBranchAddress("X_Gamma_det",&x_gamma_det);
  t->SetBranchAddress("Y_Gamma_det",&y_gamma_det);
  t->SetBranchAddress("Z_Gamma_det",&z_gamma_det);
  t->SetBranchAddress("Gamma_det_type",&gamma_det_type);
  //clusters
  t->SetBranchAddress("Cluster_res",&cluster_res);
  t->SetBranchAddress("ID_Ge",&id);
  t->SetBranchAddress("ID_Cluster",&id_cluster);//Which Cluster recorded 
  t->SetBranchAddress("ID_Crystal_In_Cluster",&id_crystal_in_cluster);
  //miniball
  t->SetBranchAddress("Miniball_res",&miniball_res);
  t->SetBranchAddress("ID_Miniball",&id_miniball);
  t->SetBranchAddress("ID_Cluster_In_Miniball",&id_cluster_in_miniball);
  t->SetBranchAddress("ID_Crystal_In_Miniball",&id_crystal_in_miniball);
  t->SetBranchAddress("ID_Segment_In_Miniball",&id_segment_in_miniball);
  //Cate
  t->SetBranchAddress("Cate_Pos_Res",&cate_pos_res); //Position at Cate reconstructed
  t->SetBranchAddress("X_Cate_Reconstructed",&x_cate); //Position at Cate reconstructed
  t->SetBranchAddress("Y_Cate_Reconstructed",&y_cate); //Including the resolution of the detectors
  t->SetBranchAddress("Z_Cate_Reconstructed",&z_cate);
  t->SetBranchAddress("X_Vertex_Reconstructed",&x_ver_rec); // reconstructed vertex, 
  t->SetBranchAddress("Y_Vertex_Reconstructed",&y_ver_rec); //including resolution
  t->SetBranchAddress("Z_Vertex_Reconstructed",&z_ver_rec);    
  t->SetBranchAddress("Cate_Si_Resolution",&cate_si_res);
  t->SetBranchAddress("Delta_E_CATE",&delta_e_cate);//Energy loss in Cate,
  t->SetBranchAddress("Cate_CsI_Resolution",&cate_csi_res);
  t->SetBranchAddress("E_CATE",&e_cate);

  infile->cd();
  TTree *theader = (TTree*)infile->Get("Header");
  theader->SetBranchAddress("TargetKind",&kind_target_f);
  theader->SetBranchAddress("TargetThicknessCM",&thickness);
  theader->SetBranchAddress("Total_Event_Number",&total_event_num);
  theader->SetBranchAddress("Total_Gamma_Number",&total_gamma_num);

  theader->GetEntry(0);
 
  Int_t nentries = (Int_t)t->GetEntries();
  cout << nentries << endl;
  
  int evnum_previous = -1;
  int mul_crystal_cluster=0;
  int mul_cluster_cluster=0;
  int flag_cluster_cluster[15];

  int mul_segment_miniball=0;
  int mul_crystal_miniball=0;
  int flag_crystal_miniball[8][3];

  float energy_sum_cluster;
  float energy_sum_miniball;
  float betaCorrection;

  float energy_max_cluster = -999.9;
  float energy_max_miniball = -999.9;
  float x_thatsit_cluster;
  float y_thatsit_cluster;
  float z_thatsit_cluster;
  float x_thatsit_miniball;
  float y_thatsit_miniball;
  float z_thatsit_miniball;
  float x_cate_thatsit;
  float y_cate_thatsit;
  float z_cate_thatsit;
  
  float x_target_thatsit;
  float y_target_thatsit;
  float z_target_thatsit;
  
  int id_thatsit;

  for(int iii=0;iii<nentries;iii++)
  {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    t->GetEntry(iii);

    if(evnum!=evnum_previous)
    {
      if(mul_crystal_cluster!=0)
      {
        h_cluster_notdopplercorrected->Fill(energy_sum_cluster);
        double vx1 = x_thatsit_cluster - x_target_thatsit;
        double vy1 = y_thatsit_cluster - y_target_thatsit;    
        double vx2 = x_cate_thatsit - x_target_thatsit;
        double vy2 = y_cate_thatsit - y_target_thatsit;
 	  
        if(mul_crystal_cluster==1) 
	{ 
          //Varying beta and decay position for lifetime
          for(int i=0;i<21;i++)
	  {
	  
            double vz1 = z_thatsit_cluster - z_target_thatsit - decay_move[i];
            double vz2 = z_cate_thatsit - z_target_thatsit - decay_move[i] ;
            double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);
            double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
            //cout<<"I am here"<<j<<endl;
            double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);
            double betaCorrectionVar = beta_var[i] + (beta_rec-beta_mean_tof);
            double energy_dopplercorrectedVar = energy_sum_cluster*(1-betaCorrectionVar*cos(theta_lab))
                                         /sqrt(1.0-betaCorrectionVar*betaCorrectionVar);
            h_cluster_dopplercorrected_single_beta_var[i]->Fill(energy_dopplercorrectedVar);
	    if(i==9 && energy_dopplercorrectedVar>600) 
	    h_cluster_x_position_of_gamma_at_moment_of_decay->Fill(z_g0);
            if(id_cluster==1 || id_cluster==2 || id_cluster==7 || id_cluster==8 || id_cluster==9)
            { 
              h_cluster_dopplercorrected_single_inner_ring_beta_var[i]->Fill(energy_dopplercorrectedVar);
            }
            if(id_cluster==3 || id_cluster==5 || id_cluster==11 || id_cluster==13 || id_cluster==15 )
            {
              h_cluster_dopplercorrected_single_middle_ring_beta_var[i]->Fill(energy_dopplercorrectedVar);
            }
            if(id_cluster==4 || id_cluster==6 || id_cluster==10 || id_cluster==12 || id_cluster==14 )
            {
              h_cluster_dopplercorrected_single_outer_ring_beta_var[i]->Fill(energy_dopplercorrectedVar);
            }
          } 
	}
      }
      if(mul_segment_miniball>=1)
      {
        //cout<<"Entering Miniball"<<endl;
        for(int j=0;j<8;j++)
        {
          for(int k=0;k<3;k++)
          {
            mul_crystal_miniball = mul_crystal_miniball + flag_crystal_miniball[j][k];
          }
	}
        
        h_miniball_notdopplercorrected->Fill(energy_sum_miniball);

        double vx1 = x_thatsit_miniball - x_target_thatsit;
        double vy1 = y_thatsit_miniball - y_target_thatsit;
       
        double vx2 = x_cate_thatsit - x_target_thatsit;
        double vy2 = y_cate_thatsit - y_target_thatsit;
        
	if(mul_crystal_miniball==1)
	{  
          //Varying beta and decay position for lifetime
          for(int i=0;i<21;i++)
	  {
            double vz1 = z_thatsit_miniball - z_target_thatsit - decay_move[i];
            double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);
            double vz2 = z_cate_thatsit - z_target_thatsit - decay_move[i];
            double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);

            double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);
            double betaCorrectionVar = beta_var[i] + (beta_rec-beta_mean_tof);
            double energy_dopplercorrectedVar = energy_sum_miniball*(1-betaCorrectionVar*cos(theta_lab))
                                         /sqrt(1.0-betaCorrectionVar*betaCorrectionVar);
          	
            h_miniball_dopplercorrected_single_beta_var[i]->Fill(energy_dopplercorrectedVar);
	    if(i==9 && energy_dopplercorrectedVar>600) 
	    h_miniball_x_position_of_gamma_at_moment_of_decay->Fill(z_g0);
	    
            if(miniball_cluster_thatsit==0||miniball_cluster_thatsit==2||miniball_cluster_thatsit==5||miniball_cluster_thatsit==7)
            {
              h_miniball_dopplercorrected_single_inner_ring_beta_var[i]->Fill(energy_dopplercorrectedVar);
            }
            if(miniball_cluster_thatsit==1||miniball_cluster_thatsit==3||miniball_cluster_thatsit==4||miniball_cluster_thatsit==6)
            {
              h_miniball_dopplercorrected_single_outer_ring_beta_var[i]->Fill(energy_dopplercorrectedVar);
            }
	  }  
        }
      }
      evnum_previous = evnum;

      mul_crystal_cluster = 0;
      energy_sum_cluster = 0.0;
      energy_max_cluster = 0.0;

      for(int j=0;j<8;j++)
      {
        for(int k=0;k<3;k++)
        {
          flag_crystal_miniball[j][k]=0;
        }
      }
      mul_segment_miniball = 0;
      mul_crystal_miniball = 0;
      energy_sum_miniball = 0.0;
      energy_max_miniball = 0.0;

    }
    if(energy_notcor>0.0)
    {
      if(gamma_det_type==1)
      {
        energy_sum_cluster = energy_sum_cluster + energy_notcor;
       mul_crystal_cluster = mul_crystal_cluster+1;
        if(energy_notcor > energy_max_cluster) 
        {
          energy_max_cluster = energy_notcor;

          x_thatsit_cluster = x_gamma_det;
          y_thatsit_cluster = y_gamma_det;
          z_thatsit_cluster = z_gamma_det;

          x_target_thatsit = x_ver_rec;
          y_target_thatsit = y_ver_rec;
          z_target_thatsit = z_ver_rec;
 
          x_cate_thatsit = x_cate;
          y_cate_thatsit = y_cate;
          z_cate_thatsit = z_cate;
          id_thatsit = id;
	}
      }
      if(gamma_det_type==2)
      {
        //cout<<"gamma_det_type==2"<<endl;
        int id_cl = (int)id_cluster_in_miniball;
        //cout<<"id_cl = "<<id_cl<<endl;
        int id_cr = (int)id_crystal_in_miniball;
        flag_crystal_miniball[id_cl][id_cr] = 1;

        mul_segment_miniball = mul_segment_miniball+1;
        //cout<<"mul_segment_miniball = "<<mul_segment_miniball<<endl;
        energy_sum_miniball = energy_sum_miniball + energy_notcor;
      
        if(energy_notcor > energy_max_miniball) 
        {
          energy_max_miniball = energy_notcor;
          miniball_cluster_thatsit = id_cl;
          x_thatsit_miniball = x_gamma_det;
          y_thatsit_miniball = y_gamma_det;
          z_thatsit_miniball = z_gamma_det;

          x_target_thatsit = x_ver_rec;
          y_target_thatsit = y_ver_rec;
          z_target_thatsit = z_ver_rec;
 
          x_cate_thatsit = x_cate;
          y_cate_thatsit = y_cate;
          z_cate_thatsit = z_cate;
	}
      }
    }
  }
  //Writing into the root file:
  /////////////////////////////
  //Cluster!!!!!!!!!!!!!!!!!!!!
  /////////////////////////////
  rootfile->cd();
  h_cluster_notdopplercorrected->Write();
  rootfile->cd();
  h_cluster_x_position_of_gamma_at_moment_of_decay->Write();
  for(int i=0;i<21;i++)
  {
    rootfile->cd();
    h_cluster_dopplercorrected_single_beta_var[i]->Write();
    rootfile->cd();
    h_cluster_dopplercorrected_single_inner_ring_beta_var[i]->Write();
    rootfile->cd();
    h_cluster_dopplercorrected_single_middle_ring_beta_var[i]->Write();
    rootfile->cd();
    h_cluster_dopplercorrected_single_outer_ring_beta_var[i]->Write();
  }
  /////////////////////////////
  //Miniball!!!!!!!!!!!!!!!!!!!!
  /////////////////////////////
  rootfile->cd();
  h_miniball_notdopplercorrected->Write();
  rootfile->cd();
  h_miniball_x_position_of_gamma_at_moment_of_decay->Write();
  for(int i=0;i<21;i++)
  {
    rootfile->cd();
    h_miniball_dopplercorrected_single_beta_var[i]->Write();
    rootfile->cd();
    h_miniball_dopplercorrected_single_inner_ring_beta_var[i]->Write();
    rootfile->cd();
    h_miniball_dopplercorrected_single_outer_ring_beta_var[i]->Write();
  }
}
