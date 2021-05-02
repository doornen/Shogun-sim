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

void RisingFastBeamReconstructor()
{ 
//char root_in[200]="/d/rising03/GEANT4/Builder/efficiency/events_built_efficiency_200mev_2500kev.root";
//char root_out[200]="/d/rising03/GEANT4/Reconstructor/efficiency/out_events_built_efficiency_200mev_2500kev.root";
//char root_in[200]="/d/rising03/GEANT4/Builder/36ca/events_built_36ca_195mev_0700mgBe_3015kev_100fs_borrelbindenergy.root";
//char root_out[200]="/d/rising03/GEANT4/Reconstructor/36ca/out_events_built_36ca_195mev_0700mgBe_3015kev_100fs_borrelbindenergy.root";
//char root_in[200]="/d/rising03/GEANT4/Builder/32s/events_built_32s_195mev_700mgBe_2230kev_00168fs_borrelbindenergy_3mm_pos.root";
//char root_out[200]="/d/rising03/GEANT4/Reconstructor/32s/out_events_built_32s_195mev_700mgBe_2230kev_00168fs_borrelbindenergy_3mm_pos.root";


//char root_in[200]="/d/rising03/GEANT4/Builder/36k/new/events_built_36k_195mev_0700mgBe_810kev_7500fs_borrelbindenergy.root";
//char root_out[200]="/d/rising03/GEANT4/Reconstructor/36k/new/out_events_built_36k_195mev_0700mgBe_810kev_7500fs_borrelbindenergy.root";
//char root_in[200]="/d/rising03/GEANT4/Builder/34ar/events_built_34ar_195mev_0700mgBe_2090kev_100fs_borrelbindenergy.root";
//char root_out[200]="/d/rising03/GEANT4/Reconstructor/34ar/out_events_built_34ar_195mev_0700mgBe_2090kev_100fs_borrelbindenergy.root";
//char root_in[200]="/d/rising03/GEANT4/Builder/34cl/events_built_34cl_195mev_0700mgBe_461kev_10000fs_borrelbindenergy.root";
//char root_out[200]="/d/rising03/GEANT4/Reconstructor/34cl/out_events_built_34cl_195mev_0700mgBe_461kev_10000fs_borrelbindenergy.root";

char root_in[200]="/media/BOOTCAMP/events_built_dummy.root";
char root_out[200]="/media/BOOTCAMP/out_events_built_dummy.root";


//char root_out[200]="./newCalculations/34ar/out_events_built_34ar_2090kev_320fs.root";
  char name[100];
  int  numBin = 500;
  float firstBin   = 0.0 - 2000/(2*numBin);
  float lastBin   = 2000.0 - 2000/(2*numBin);
  double beta = 0.533;        //The average value for the correction
  double beta_mean_tof = 0.5631;//The mean beta from Tof
  //double beta = 0.5455;        //The average value for the correction
  //double beta_mean_tof = 0.5631;//The mean beta from Tof
  double decay_move_target = 1.0;  //The lifetime of the excited state moves the average decay point
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
  //hector
  float hector_res;
  float id_hector;
  float hector_thatsit;
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
  rootfile->cd();

  //Creating CATE spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TH2F *h_beta_before_beta_after = new TH2F("h_beta_before_beta_after","h_beta_before_beta_after",200,0.50,0.57,200,0.50,0.57);
  TH2F *h_beta_real_beta_before = new TH2F("h_beta_real_beta_before","h_beta_real_beta_before",200,0.50,0.57,200,0.50,0.57);
  TH2F *h_beta_real_gamma_dopp = new TH2F("h_beta_real_e_gamma","h_beta_real_e_gamma",200,0.50,0.57,1250,0,5000);
  TH2F *h_beta_real_beta_after = new TH2F("h_beta_real_beta_after","h_beta_real_beta_after",200,0.50,0.57,200,0.50,0.57);
  TH2F *h_cate_pos = new TH2F("CatePos","CatePos",200,-20.0,20.0,200,-20.0,20.0);
  TH2F *h_e_de = new TH2F("Cate_e_de","Cate_e_de",400,4000.0,8000.0,400,0.0,200.0);
  
  //Creating cluster spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rootfile->cd();
  TH1F *h_cluster_crystal_mult = new TH1F("Cluster_Crystal_Mult",
                                        "Cluster_Crystal_mult",15,0,15);
  TH1F *h_cluster_cluster_mult = new TH1F("Cluster_Cluster_Mult",
                                        "Cluster_Cluster_mult",15,0,15);
  //Single hits
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_single = new TH1F("Cluster_DopplerCorrected_SingleHit",
                                        "Cluster_DopplerCorrected_SingleHit",numBin,firstBin,lastBin);
  TH1F *h_cluster_dopplercorrected_mult2 = new TH1F("Cluster_DopplerCorrected_mul2",
                                        "Cluster_DopplerCorrected_mult2",numBin,firstBin,lastBin);					
  TH1F *h_cluster_dopplercorrected_mult3 = new TH1F("Cluster_DopplerCorrected_mul3",
                                        "Cluster_DopplerCorrected_mult3",numBin,firstBin,lastBin);
  //Single hit ringwise
  char tempname[100];
  TH1F *h_cluster_dopplercorrected_single_ring[3];
  for(int i=0;i<3;i++)
  {
    rootfile->cd();
    sprintf(tempname,"Cluster_DopplerCorrected_SingleHit_Ring_%i",i+1);
    h_cluster_dopplercorrected_single_ring[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }  
  //Outer rings
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_single_outerrings = new TH1F("Cluster_DopplerCorrected_SingleHit_OuterRings",
                                              "Cluster_DopplerCorrected_SingleHit_OuterRings",numBin,firstBin,lastBin);

  //Making no tracking but adjusting beta from tof
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_single_no_tracking = new TH1F("Cluster_DopplerCorrected_SingleHit_NoTracking",
                                        "Cluster_DopplerCorrected_SingleHit_NoTracking",numBin,firstBin,lastBin);
  
//Making tracking but no adjustment of beta from tof
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_single_beta_const = new TH1F("Cluster_DopplerCorrected_SingleHit_BetaConstant",
                                        "Cluster_DopplerCorrected_SingleHit_BetaConstant",numBin,firstBin,lastBin);

 //Making no tracking AND no adjustment of beta from tof
 //The most primitive way:
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_single_no_tracking_beta_const = new TH1F("Cluster_DopplerCorrected_SingleHit_NoTracking_BetaConstant",
                                        "Cluster_DopplerCorrected_SingleHit_NoTracking_BetaConstant",numBin,firstBin,lastBin);

  //Not doppler corrected
  rootfile->cd();
  TH1F *h_cluster_notdopplercorrected_total = new TH1F("Cluster_NotDopplerCorrected_Total",
                                "Cluster_NotDopplerCorrected_Total",numBin,firstBin,lastBin);
  TH1F *h_cluster_notdopplercorrected_total_ring[3];
  for(int i=0;i<3;i++)
  {
    rootfile->cd();
    sprintf(tempname,"Cluster_NotDopplerCorrected_Total_Ring_%i",i+1);
    h_cluster_notdopplercorrected_total_ring[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }  				
				
  // single hit AND addback
  
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_total = new TH1F("Cluster_DopplerCorrected_Total",
                                            "Cluster_DopplerCorrected_Total",numBin,firstBin,lastBin);
  TH1F *h_cluster_dopplercorrected_total_ring[3];
  for(int i=0;i<3;i++)
  {
    rootfile->cd();
    sprintf(tempname,"Cluster_DopplerCorrected_Total_Ring_%i",i+1);
    h_cluster_dopplercorrected_total_ring[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }  		

  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_total_outerrings = new TH1F("Cluster_DopplerCorrected_Total_OuterRings",
                                              "Cluster_DopplerCorrected_Total_OuterRings",numBin,firstBin,lastBin);
   
  //Only addback
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_addback = new TH1F("Cluster_DopplerCorrected_Addback",
                                      "Cluster_DopplerCorrected_Addback",numBin,firstBin,lastBin);
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_addback_65 = new TH1F("Cluster_DopplerCorrected_Addback_65",
                                         "Cluster_DopplerCorrected_Only_Addback_65",numBin,firstBin,lastBin);
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_addback_70 = new TH1F("Cluster_DopplerCorrected_Addback_70",
                                          "Cluster_DopplerCorrected_Addback_70",numBin,firstBin,lastBin);
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_addback_75 = new TH1F("Cluster_DopplerCorrected_Addback_75",
                                          "Cluster_DopplerCorrected_Only_Addback_75",numBin,firstBin,lastBin);
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_addback_80 = new TH1F("Cluster_DopplerCorrected_Addback_80",
                                          "Cluster_DopplerCorrected_Only_Addback_80",numBin,firstBin,lastBin);
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_addback_85 = new TH1F("Cluster_DopplerCorrected_Addback_85",
                                         "Cluster_DopplerCorrected_Only_Addback_85",numBin,firstBin,lastBin);
  rootfile->cd();
  TH1F *h_cluster_dopplercorrected_addback_90 = new TH1F("Cluster_DopplerCorrected_Addback_90",
                                         "Cluster_DopplerCorrected_Addback_90",numBin,firstBin,lastBin);
  
  ///////////////////////////////
  //Miniball!!!!!!!!!!!!!!!!!!!!!
  //////////////////////////////
  TH1F *h_miniball_segment_mult = new TH1F("Miniball_Segment_Mult",
                                        "Miniball_Segment_mult",15,0,15);
  TH1F *h_miniball_crystal_mult = new TH1F("Miniball_Crystal_Mult",
                                        "Miniball_Crystal_mult",15,0,15);
  TH1F *h_miniball_cluster_mult = new TH1F("Miniball_Cluster_Mult",
                                        "Miniball_Cluster_mult",15,0,15);

  //Single hits
  TH1F *h_miniball_dopplercorrected_single = new TH1F("Miniball_DopplerCorrected_Single",
                                        "Miniball_DopplerCorrected_Single",numBin,firstBin,lastBin);
  
  TH1F *h_miniball_dopplercorrected_single_innerring = new TH1F("Miniball_DopplerCorrected_Single_InnerRing",
                                        "Miniball_DopplerCorrected_Single_InnerRing",numBin,firstBin,lastBin);
  
  TH1F *h_miniball_dopplercorrected_single_outerring = new TH1F("Miniball_DopplerCorrected_Single_OuterRing",
                                        "Miniball_DopplerCorrected_Single_OuterRing",numBin,firstBin,lastBin);
  //Segment Mult=1 
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_single_mulseg1 = new TH1F("Miniball_DopplerCorrected_SingleHit_MulSeg1",
                                        "Miniball_DopplerCorrected_SingleHit_MulSeg1",numBin,firstBin,lastBin);
  //Mult==2
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_single_mulseg2 = new TH1F("Miniball_DopplerCorrected_SingleHit_MulSeg2",
                                        "Miniball_DopplerCorrected_SingleHit_MulSeg2",numBin,firstBin,lastBin);
  TH1F *h_miniball_dopplercorrected_single_mulseg3 = new TH1F("Miniball_DopplerCorrected_SingleHit_MulSeg3",
                                        "Miniball_DopplerCorrected_SingleHit_MulSeg3",numBin,firstBin,lastBin);				
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_single_innerring_mulseg2 = new TH1F("Miniball_DopplerCorrected_SingleHit_InnerRing_MulSeg2",
                                              "Miniball_DopplerCorrected_SingleHit_InnerRing_MulSeg2",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_single_outerring_mulseg2 = new TH1F("Miniball_DopplerCorrected_SingleHit_OuterRing_MulSeg2",
                                              "Miniball_DopplerCorrected_SingleHit_OuterRing_MulSeg2",numBin,firstBin,lastBin);

  //Taking both, mulseg1 and mulseg2 and restrictions on the tracking and beta determination:
  //No Tracking
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_single_no_tracking = new TH1F("Miniball_DopplerCorrected_SingleHit_NoTracking",
                                        "Miniball_DopplerCorrected_SingleHit_NoTracking",numBin,firstBin,lastBin);
  //Beta Const
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_single_beta_const = new TH1F("Miniball_DopplerCorrected_SingleHit_BetaConstant",
                                        "Miniball_DopplerCorrected_SingleHit_BetaConstant",numBin,firstBin,lastBin);
  //No tracking, Beta Const
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_single_no_tracking_beta_const = new TH1F("Miniball_DopplerCorrected_SingleHit_NoTracking_BetaConstant",
                                        "Miniball_DopplerCorrected_SingleHit_NoTracking_BetaConstant",numBin,firstBin,lastBin);

  //Not doppler corrected
  rootfile->cd();
  TH1F *h_miniball_notdopplercorrected = new TH1F("Miniball_NotDopplerCorrected",
                                "Miniball_NotDopplerCorrected",numBin,firstBin,lastBin);
  
  TH1F *h_miniball_notdopplercorrected_mulseg2_ring[2];
  TH1F *h_miniball_notdopplercorrected_ring[2];

  for(int i=0;i<2;i++)
  {
    rootfile->cd();
    sprintf(name,"h_miniball_notdopplercorrected_mulseg2_ring[%i]",i);
    h_miniball_notdopplercorrected_mulseg2_ring[i] = new TH1F(name,name,numBin,firstBin,lastBin);

  }
  for(int i=0;i<2;i++)
  {
    rootfile->cd();
    sprintf(name,"h_miniball_notdopplercorrected_ring[%i]",i);
    h_miniball_notdopplercorrected_ring[i] = new TH1F(name,name,numBin,firstBin,lastBin);

  }
  // Single hit AND addback
  //No Cut on segments:
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_total = new TH1F("Miniball_DopplerCorrected_Total",
                                            "Minball_DopplerCorrected_Total",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_total_innerring = new TH1F("Miniball_DopplerCorrected_Total_InnerRing",
                                               "Miniball_DopplerCorrected_Total_InnerRing",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_total_outerring = new TH1F("Miniball_DopplerCorrected_Total_OuterRing",
                                               "Miniball_DopplerCorrected_Total_OuterRing",numBin,firstBin,lastBin); 				  
  //Mulseg==2
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_total_mulseg2 = new TH1F("Miniball_DopplerCorrected_Total_MulSeg2",
                                            "Minball_DopplerCorrected_Total_MulSeg2",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_total_innerring_mulseg2 = new TH1F("Miniball_DopplerCorrected_Total_InnerRing_MulSeg2",
                                               "Miniball_DopplerCorrected_Total_InnerRing_MulSeg2",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_total_outerring_mulseg2 = new TH1F("Miniball_DopplerCorrected_Total_OuterRing_MulSeg2",
                                               "Miniball_DopplerCorrected_Total_OuterRing_MulSeg2",numBin,firstBin,lastBin);   
 
 //Only addback
 
  TH1F *h_miniball_dopplercorrected_addback = new TH1F("Miniball_DopplerCorrected_Addback",
                                      "Miniball_DopplerCorrected_Only_Addback",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_addback_mulseg2 = new TH1F("Miniball_DopplerCorrected_Addback_MulSeg2",
                                      "Miniball_DopplerCorrected_Only_Addback_MulSeg2",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_addback_innerring_mulseg2 = new TH1F("Miniball_DopplerCorrected_Addback_InnerRing_MulSeg2",
                                               "Miniball_DopplerCorrected_AddBack_InnerRing_MulSeg2",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_miniball_dopplercorrected_addback_outerring_mulseg2 = new TH1F("Miniball_DopplerCorrected_Addback_OuterRing_MulSeg2",
                                               "Miniball_DopplerCorrected_Addback_OuterRing_MulSeg2",numBin,firstBin,lastBin);
					       		      
  ///////////////////////////////
  //Hector!!!!!!!!!!!!!!!!!!!!!
  ///////////////////////////////
    
  TH1F *h_hector_crystal_mult = new TH1F("Hector_Crystal_Mult",
                                        "Hector_Crystal_mult",15,0,15);
  //Single hits
  rootfile->cd();
  TH1F *h_hector_dopplercorrected = new TH1F("Hector_DopplerCorrected",
                                        "Hector_DopplerCorrected",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_hector_dopplercorrected_90 = new TH1F("Hector_DopplerCorrected_90",
                                              "Hector_DopplerCorrected_90",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_hector_dopplercorrected_142 = new TH1F("Hector_DopplerCorrected_142",
                                              "Hector_DopplerCorrected_142",numBin,firstBin,lastBin);
  //No tracking
  rootfile->cd();
  TH1F *h_hector_dopplercorrected_no_tracking = new TH1F("Hector_DopplerCorrected_NoTracking",
                                        "Hector_DopplerCorrected_NoTracking",numBin,firstBin,lastBin);
  //Beta constant
  rootfile->cd();
  TH1F *h_hector_dopplercorrected_beta_const = new TH1F("Hector_DopplerCorrected_BetaConstant",
                                        "Hector_DopplerCorrected_BetaConstant",numBin,firstBin,lastBin);
  //NO tracking, Beta constant
  rootfile->cd();
  TH1F *h_hector_dopplercorrected_no_tracking_beta_const = new TH1F("Hector_DopplerCorrected_NoTracking_BetaConstant",
                                        "Hector_DopplerCorrected_NoTracking_BetaConstant",numBin,firstBin,lastBin);

  //Not doppler corrected
  rootfile->cd();
  TH1F *h_hector_notdopplercorrected = new TH1F("Hector_NotDopplerCorrected",
                                "Hector_NotDopplerCorrected",numBin,firstBin,lastBin);
  rootfile->cd();
  
  TH1F *h_hector_notdopplercorrected_90 = new TH1F("Hector_NotDopplerCorrected_90",
                                              "Hector_NotDopplerCorrected_90",numBin,firstBin,lastBin);
  
  rootfile->cd();
  TH1F *h_hector_notdopplercorrected_142 = new TH1F("Hector_NotDopplerCorrected_142",
                                              "Hector_NotDopplerCorrected_142",numBin,firstBin,lastBin);
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
  t->SetBranchAddress("Beta_Reconstructed",&beta_rec);  //This is the reconstructed beta which inlcudes the tof resolution!!!!
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
  //Hector:
  t->SetBranchAddress("Hector_res",&hector_res);
  t->SetBranchAddress("ID_Hector",&id_hector);
  //Miniball
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
  cout <<"Entries: "<< nentries << endl;
  
  int evnum_previous = -1;
  int mul_crystal_cluster=0;
  int mul_cluster_cluster=0;
  int flag_cluster_cluster[15];

  int mul_segment_miniball=0;
  int mul_crystal_miniball=0;
  int mul_cluster_miniball=0;
  int flag_cluster_miniball[8];
  int flag_crystal_miniball[8][3];

  int mul_crystal_hector=0;

  float energy_sum_cluster;
  float energy_sum_miniball;
  float energy_sum_hector;
  float energy_0_cluster;
  float energy_1_cluster;
  float energy_2_cluster;
  float energy_0_miniball;
  float energy_1_miniball;
  float energy_2_miniball;
  float betaCorrection;

  float energy_max_cluster = -999.9;
  float energy_max_miniball = -999.9;
  float energy_max_hector = -999.9;
  float x_thatsit_cluster;
  float y_thatsit_cluster;
  float z_thatsit_cluster;
  float x_thatsit_miniball;
  float y_thatsit_miniball;
  float z_thatsit_miniball;
  float x_thatsit_hector;
  float y_thatsit_hector;
  float z_thatsit_hector;

  float x_cate_thatsit;
  float y_cate_thatsit;
  float z_cate_thatsit;
  
  float x_target_thatsit;
  float y_target_thatsit;
  float z_target_thatsit;
  
  float beta_dummy;
  int id_thatsit;
  //int id_miniball_thatsit;
  //int id_hector_thatsit;
  //int kind_thatsit;

  //double gamma_rec;
  //double energy_rec;
  //double energy_cor;
  //double gamma_cor;

  for(int iii=0;iii<nentries;iii++)
  {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    t->GetEntry(iii);

    if(evnum!=evnum_previous)
    {
     h_beta_real_gamma_dopp->Fill(beta_r,e_doppler);
    //checking Beta
     h_beta_before_beta_after->Fill(beta_b,beta_a);
    h_beta_real_beta_after->Fill(beta_r,beta_a);
    h_beta_real_beta_before->Fill(beta_r,beta_b);
     //Filling cate spectra
    h_cate_pos->Fill(x_cate,y_cate);
    h_e_de->Fill(e_cate,delta_e_cate);
      if(mul_crystal_cluster!=0)
      {
        //cout<<"Entering Cluster"<<endl;
        for(int j=0;j<15;j++)
        {
          mul_cluster_cluster = mul_cluster_cluster + flag_cluster_cluster[j];
	}
        h_cluster_notdopplercorrected_total->Fill(energy_sum_cluster);
        h_cluster_crystal_mult->Fill(mul_crystal_cluster);
        h_cluster_cluster_mult->Fill(mul_cluster_cluster);

        if(mul_cluster_cluster==1)
	{
          if(id_cluster==1 || id_cluster==2 || id_cluster==7 || id_cluster==8 || id_cluster==9 )
          { 
            h_cluster_notdopplercorrected_total_ring[0]->Fill(energy_sum_cluster);
          }
          if(id_cluster==3 || id_cluster==5 || id_cluster==11 || id_cluster==13 || id_cluster==15 )
          {
            h_cluster_notdopplercorrected_total_ring[1]->Fill(energy_sum_cluster);
          }
          if(id_cluster==4 || id_cluster==6 || id_cluster==10 || id_cluster==12 || id_cluster==14 )
          {
            h_cluster_notdopplercorrected_total_ring[2]->Fill(energy_sum_cluster);
          }
	}  
        double vx1 = x_thatsit_cluster - x_target_thatsit;
        double vy1 = y_thatsit_cluster - y_target_thatsit;
        double vz1 = z_thatsit_cluster - z_target_thatsit;
        double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);
        double vx2 = x_cate_thatsit - x_target_thatsit;
        double vy2 = y_cate_thatsit - y_target_thatsit;
        double vz2 = z_cate_thatsit - z_target_thatsit;
        double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
          
        //With normal tracking:
        double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);
        //Without any tracking:
        double l1_no_tr = sqrt(x_thatsit_cluster*x_thatsit_cluster+y_thatsit_cluster*y_thatsit_cluster+z_thatsit_cluster*z_thatsit_cluster);
      
      double theta_lab_no_tracking = acos((z_thatsit_cluster/l1_no_tr));
        betaCorrection = beta + (beta_dummy-beta_mean_tof);
	//betaCorrection = beta_dummy;
        //dopplerCorrection including tracking and precise beta determination
        double energy_dopplercorrected = energy_sum_cluster*(1-betaCorrection*cos(theta_lab))
                                         /sqrt(1.0-betaCorrection*betaCorrection);
        //No tracking
        double energy_dopplercorrected_no_tracking =  energy_sum_cluster*(1-betaCorrection*cos(theta_lab_no_tracking))
                                         /sqrt(1.0-betaCorrection*betaCorrection);
        //Beta constant
        double energy_dopplercorrected_beta_const = energy_sum_cluster*(1-beta*cos(theta_lab))
                                         /sqrt(1.0-beta*beta);
        //No tracking and beta constant
        double energy_dopplercorrected_no_tracking_beta_const = energy_sum_cluster*(1-beta*cos(theta_lab_no_tracking))
	                                 /sqrt(1.0-beta*beta);
        if(mul_crystal_cluster==1) 
	{
          //  cout << id_thatsit << "   " << x_ge_thatsit << " " << y_ge_thatsit << " " << z_ge_thatsit << endl;
          h_cluster_dopplercorrected_single->Fill(energy_dopplercorrected);
          h_cluster_dopplercorrected_total->Fill(energy_dopplercorrected);
          h_cluster_dopplercorrected_single_no_tracking->Fill(energy_dopplercorrected_no_tracking);
          h_cluster_dopplercorrected_single_beta_const->Fill(energy_dopplercorrected_beta_const);
          h_cluster_dopplercorrected_single_no_tracking_beta_const->Fill(energy_dopplercorrected_no_tracking_beta_const);

          if(id_cluster==1 || id_cluster==2 || id_cluster==7 || id_cluster==8 || id_cluster==9 )
          {
            h_cluster_dopplercorrected_single_ring[0]->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_total_ring[0]->Fill(energy_dopplercorrected);
          }
          if(id_cluster==3 || id_cluster==5 || id_cluster==11 || id_cluster==13 || id_cluster==15 )
          {
            h_cluster_dopplercorrected_single_ring[1]->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_total_ring[1]->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_single_outerrings->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_total_outerrings->Fill(energy_dopplercorrected);
          }
          if(id_cluster==4 || id_cluster==6 || id_cluster==10 || id_cluster==12 || id_cluster==14 )
          {
            h_cluster_dopplercorrected_single_ring[2]->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_total_ring[2]->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_single_outerrings->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_total_outerrings->Fill(energy_dopplercorrected);
          }
	}
        if(mul_crystal_cluster==2 && mul_cluster_cluster==1) 
        {
	  h_cluster_dopplercorrected_mult2->Fill(energy_dopplercorrected);
          h_cluster_dopplercorrected_total->Fill(energy_dopplercorrected);
	  h_cluster_dopplercorrected_addback->Fill(energy_dopplercorrected);
          if(id_cluster==1 || id_cluster==2 || id_cluster==7 || id_cluster==8 || id_cluster==9 )
          {
            h_cluster_dopplercorrected_total_ring[0]->Fill(energy_dopplercorrected);
          }
          if(id_cluster==3 || id_cluster==5 || id_cluster==11 || id_cluster==13 || id_cluster==15 )
          {
            h_cluster_dopplercorrected_total_ring[1]->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_total_outerrings->Fill(energy_dopplercorrected);
          }
          if(id_cluster==4 || id_cluster==6 || id_cluster==10 || id_cluster==12 || id_cluster==14 )
          {
            h_cluster_dopplercorrected_total_ring[2]->Fill(energy_dopplercorrected);
            h_cluster_dopplercorrected_total_outerrings->Fill(energy_dopplercorrected);
          }
          if(mul_crystal_cluster==2)
          {
            float energy_tmp = energy_0_cluster + energy_1_cluster;
            if(energy_0_cluster>=0.65*(energy_tmp) || energy_1_cluster>=0.65*(energy_tmp))
              h_cluster_dopplercorrected_addback_65->Fill(energy_dopplercorrected);
            if(energy_0_cluster>=0.70*(energy_tmp) || energy_1_cluster>=0.70*(energy_tmp))
              h_cluster_dopplercorrected_addback_70->Fill(energy_dopplercorrected);
            if(energy_0_cluster>=0.75*(energy_tmp) || energy_1_cluster>=0.75*(energy_tmp))
              h_cluster_dopplercorrected_addback_75->Fill(energy_dopplercorrected);    
            if(energy_0_cluster>=0.80*(energy_tmp) || energy_1_cluster>=0.80*(energy_tmp))
              h_cluster_dopplercorrected_addback_80->Fill(energy_dopplercorrected);
            if(energy_0_cluster>=0.85*(energy_tmp) || energy_1_cluster>=0.85*(energy_tmp))
              h_cluster_dopplercorrected_addback_85->Fill(energy_dopplercorrected);
            if(energy_0_cluster>=0.90*(energy_tmp) || energy_1_cluster>=0.90*(energy_tmp))
              h_cluster_dopplercorrected_addback_90->Fill(energy_dopplercorrected);
          }	   
        }
	if(mul_crystal_cluster==3 && mul_cluster_cluster==1) 
        {
	  h_cluster_dopplercorrected_mult3->Fill(energy_dopplercorrected);
	}  
      }
      if(mul_segment_miniball>=1)
      {
        //cout<<"Entering Miniball"<<endl;
        for(int j=0;j<8;j++)
        {
          mul_cluster_miniball = mul_cluster_miniball + flag_cluster_miniball[j];
          for(int k=0;k<3;k++)
          {
            mul_crystal_miniball = mul_crystal_miniball + flag_crystal_miniball[j][k];
          }
	}
        h_miniball_segment_mult->Fill(mul_segment_miniball);
        h_miniball_crystal_mult->Fill(mul_crystal_miniball);
        h_miniball_cluster_mult->Fill(mul_cluster_miniball);

        h_miniball_notdopplercorrected->Fill(energy_sum_miniball);
	if(mul_cluster_miniball==1)
	{
          if(miniball_cluster_thatsit==0||miniball_cluster_thatsit==2||miniball_cluster_thatsit==5||miniball_cluster_thatsit==7)
          {
            h_miniball_notdopplercorrected_ring[0]->Fill(energy_sum_miniball);
	  }
          if(miniball_cluster_thatsit==1||miniball_cluster_thatsit==3||miniball_cluster_thatsit==4||miniball_cluster_thatsit==6)
          {
            h_miniball_notdopplercorrected_ring[1]->Fill(energy_sum_miniball);
	  }
	}   
        double vx1 = x_thatsit_miniball - x_target_thatsit;
        double vy1 = y_thatsit_miniball - y_target_thatsit;
        double vz1 = z_thatsit_miniball - z_target_thatsit;
        double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);
        double vx2 = x_cate_thatsit - x_target_thatsit;
        double vy2 = y_cate_thatsit - y_target_thatsit;
        double vz2 = z_cate_thatsit - z_target_thatsit;
        double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
	    
        double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);
        //Without any tracking:
        double l1_no_tr=  sqrt(x_thatsit_miniball*x_thatsit_miniball+
                               y_thatsit_miniball*y_thatsit_miniball+
			       z_thatsit_miniball*z_thatsit_miniball);
        double theta_lab_no_tracking = acos((z_thatsit_miniball/l1_no_tr));

	betaCorrection = beta + (beta_dummy-beta_mean_tof);
        double energy_dopplercorrected = energy_sum_miniball*(1-betaCorrection*cos(theta_lab))
                                         /sqrt(1.0-betaCorrection*betaCorrection);	
        //No tracking
        double energy_dopplercorrected_no_tracking =  energy_sum_miniball*(1-betaCorrection*cos(theta_lab_no_tracking))
                                         /sqrt(1.0-betaCorrection*betaCorrection);
        //Beta constant
        double energy_dopplercorrected_beta_const = energy_sum_miniball*(1-beta*cos(theta_lab))
                                         /sqrt(1.0-beta*beta);
        //No tracking and beta constant
        double energy_dopplercorrected_no_tracking_beta_const = energy_sum_miniball*(1-beta*cos(theta_lab_no_tracking))
                                         /sqrt(1.0-beta*beta);
	if(mul_crystal_miniball==1)
	{			 
          h_miniball_dopplercorrected_single->Fill(energy_dopplercorrected);
          h_miniball_dopplercorrected_total->Fill(energy_dopplercorrected);
          h_miniball_dopplercorrected_single_no_tracking->Fill(energy_dopplercorrected_no_tracking);
          h_miniball_dopplercorrected_single_beta_const->Fill(energy_dopplercorrected_beta_const);
          h_miniball_dopplercorrected_single_no_tracking_beta_const->Fill(energy_dopplercorrected_no_tracking_beta_const);
          if(miniball_cluster_thatsit==0||miniball_cluster_thatsit==2||miniball_cluster_thatsit==5||miniball_cluster_thatsit==7)
          {
            h_miniball_dopplercorrected_single_innerring->Fill(energy_dopplercorrected);
            h_miniball_dopplercorrected_total_innerring->Fill(energy_dopplercorrected);
          }
          if(miniball_cluster_thatsit==1||miniball_cluster_thatsit==3||miniball_cluster_thatsit==4||miniball_cluster_thatsit==6)
          {
            h_miniball_dopplercorrected_single_outerring->Fill(energy_dopplercorrected);
            h_miniball_dopplercorrected_total_outerring->Fill(energy_dopplercorrected);
          }
	 }
        // MultSeg==1
	if(mul_segment_miniball==1 && mul_crystal_miniball==1 && mul_cluster_miniball==1)
	{ 
          h_miniball_dopplercorrected_single_mulseg1->Fill(energy_dopplercorrected);
        }
        //MultSeg==2
        if(mul_segment_miniball==2 && mul_crystal_miniball==1 && mul_cluster_miniball==1) 
	{
          h_miniball_dopplercorrected_single_mulseg2->Fill(energy_dopplercorrected);
          h_miniball_dopplercorrected_total_mulseg2->Fill(energy_dopplercorrected);
          if(miniball_cluster_thatsit==0||miniball_cluster_thatsit==2||miniball_cluster_thatsit==5||miniball_cluster_thatsit==7)
          {
            h_miniball_notdopplercorrected_mulseg2_ring[0]->Fill(energy_sum_miniball);
            h_miniball_dopplercorrected_single_innerring_mulseg2->Fill(energy_dopplercorrected);
            h_miniball_dopplercorrected_total_innerring_mulseg2->Fill(energy_dopplercorrected);
          }
          if(miniball_cluster_thatsit==1||miniball_cluster_thatsit==3||miniball_cluster_thatsit==4||miniball_cluster_thatsit==6)
          {
            h_miniball_notdopplercorrected_mulseg2_ring[1]->Fill(energy_sum_miniball);
            h_miniball_dopplercorrected_single_outerring_mulseg2->Fill(energy_dopplercorrected);
            h_miniball_dopplercorrected_total_outerring_mulseg2->Fill(energy_dopplercorrected);
          }
        }
	//MultSeg==3
	if(mul_segment_miniball==3 && mul_crystal_miniball==1 && mul_cluster_miniball==1) 
	{
	  h_miniball_dopplercorrected_single_mulseg3->Fill(energy_dopplercorrected);
	}
        //The addback condition:
        if(mul_segment_miniball>=2 && mul_crystal_miniball==2 && mul_cluster_miniball==1) 
        {
          h_miniball_dopplercorrected_total->Fill(energy_dopplercorrected);
          h_miniball_dopplercorrected_addback->Fill(energy_dopplercorrected);
          if(miniball_cluster_thatsit==0||miniball_cluster_thatsit==2||miniball_cluster_thatsit==5||miniball_cluster_thatsit==7)
          {
            h_miniball_dopplercorrected_total_innerring->Fill(energy_dopplercorrected);
          }
          if(miniball_cluster_thatsit==1||miniball_cluster_thatsit==3||miniball_cluster_thatsit==4||miniball_cluster_thatsit==6)
          {
            h_miniball_dopplercorrected_total_outerring->Fill(energy_dopplercorrected);
          }
          if(mul_segment_miniball==2)
	  {
            h_miniball_dopplercorrected_total_mulseg2->Fill(energy_dopplercorrected);
            h_miniball_dopplercorrected_addback_mulseg2->Fill(energy_dopplercorrected);
            if(miniball_cluster_thatsit==0||miniball_cluster_thatsit==2||miniball_cluster_thatsit==5||miniball_cluster_thatsit==7)
            {
              h_miniball_dopplercorrected_addback_innerring_mulseg2->Fill(energy_dopplercorrected);
              h_miniball_dopplercorrected_total_innerring_mulseg2->Fill(energy_dopplercorrected);
            }
            if(miniball_cluster_thatsit==1||miniball_cluster_thatsit==3||miniball_cluster_thatsit==4||miniball_cluster_thatsit==6)
            {
              h_miniball_dopplercorrected_addback_outerring_mulseg2->Fill(energy_dopplercorrected);
              h_miniball_dopplercorrected_total_outerring_mulseg2->Fill(energy_dopplercorrected);
            }
          }
	}
      }
      if(mul_crystal_hector!=0)
      {
        //cout<<"Entering hector"<<endl;
        h_hector_notdopplercorrected->Fill(energy_sum_hector);
        h_hector_crystal_mult->Fill(mul_crystal_hector);
        double vx1 = x_thatsit_hector - x_target_thatsit;
        double vy1 = y_thatsit_hector - y_target_thatsit;
        double vz1 = z_thatsit_hector - z_target_thatsit;
        double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);
        double vx2 = x_cate_thatsit - x_target_thatsit;
        double vy2 = y_cate_thatsit - y_target_thatsit;
        double vz2 = z_cate_thatsit - z_target_thatsit;
        double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
	    
        double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);
        //Without any tracking:
        double l1_no_tr = sqrt(x_thatsit_hector*x_thatsit_hector+
                               y_thatsit_hector*y_thatsit_hector+
		               z_thatsit_hector*z_thatsit_hector);
        double theta_lab_no_tracking = acos((z_thatsit_hector/l1_no_tr));
        betaCorrection = beta + (beta_dummy-beta_mean_tof);
        double energy_dopplercorrected = energy_sum_hector*(1-betaCorrection*cos(theta_lab))
                                         /sqrt(1.0-betaCorrection*betaCorrection);
        //No tracking
        double energy_dopplercorrected_no_tracking =  energy_sum_cluster*(1-betaCorrection*cos(theta_lab_no_tracking))
                                         /sqrt(1.0-betaCorrection*betaCorrection);
        //Beta constant
        double energy_dopplercorrected_beta_const = energy_sum_cluster*(1-beta*cos(theta_lab))
                                         /sqrt(1.0-beta*beta);
        //No tracking and beta constant
        double energy_dopplercorrected_no_tracking_beta_const = energy_sum_cluster*(1-beta*cos(theta_lab_no_tracking))
                                         /sqrt(1.0-beta*beta);
	  
        if(mul_crystal_hector==1) 
	{
          h_hector_dopplercorrected->Fill(energy_dopplercorrected);
          h_hector_dopplercorrected_no_tracking->Fill(energy_dopplercorrected_no_tracking);
          h_hector_dopplercorrected_beta_const->Fill(energy_dopplercorrected_beta_const);
          h_hector_dopplercorrected_no_tracking_beta_const->Fill(energy_dopplercorrected_no_tracking_beta_const);
          if(hector_thatsit==5 || hector_thatsit==6)
	  { 
            h_hector_dopplercorrected_90->Fill(energy_dopplercorrected);
            h_hector_notdopplercorrected_90->Fill(energy_sum_hector);
          }
          else
	  {
            h_hector_dopplercorrected_142->Fill(energy_dopplercorrected);
            h_hector_notdopplercorrected_142->Fill(energy_sum_hector);
          }
        }
      }
      evnum_previous = evnum;
      for(int j=0;j<15;j++)
      {
        flag_cluster_cluster[j]=0;
      }
      mul_crystal_cluster = 0;
      mul_cluster_cluster = 0;
      energy_sum_cluster = 0.0;
      energy_max_cluster = 0.0;
      energy_0_cluster = 0.0;
      energy_1_cluster = 0.0;
      energy_2_cluster = 0.0;


      for(int j=0;j<8;j++)
      {
        for(int k=0;k<3;k++)
        {
          flag_crystal_miniball[j][k]=0;
        }
        flag_cluster_miniball[j]=0;
      }
      mul_segment_miniball = 0;
      mul_crystal_miniball = 0;
      mul_cluster_miniball = 0;
      energy_sum_miniball = 0.0;
      energy_max_miniball = 0.0;
      energy_0_miniball = 0.0;
      energy_1_miniball = 0.0;
      energy_2_miniball = 0.0;

      mul_crystal_hector=0;
      energy_max_hector = 0.0;
      energy_sum_hector=0.0;
    }
    if(energy_notcor>0.0)
    {
      if(gamma_det_type==1)
      {
        //cout<<"gamma_det_type==1"<<endl;
        int id_cl = (int)id_cluster-1;
        flag_cluster_cluster[id_cl] = 1;
        mul_crystal_cluster = mul_crystal_cluster+1;
        energy_sum_cluster = energy_sum_cluster + energy_notcor;
      
        if(mul_crystal_cluster==1) energy_0_cluster = energy_notcor;
        if(mul_crystal_cluster==2) energy_1_cluster = energy_notcor;
        if(mul_crystal_cluster==3) energy_2_cluster = energy_notcor;

        if(energy_notcor > energy_max_cluster) 
        {
          energy_max_cluster = energy_notcor;

          x_thatsit_cluster = x_gamma_det;
          y_thatsit_cluster = y_gamma_det;
          z_thatsit_cluster = z_gamma_det;

          x_target_thatsit = x_ver_rec;
          y_target_thatsit = y_ver_rec;
          z_target_thatsit = z_ver_rec+decay_move_target;
 
          x_cate_thatsit = x_cate;
          y_cate_thatsit = y_cate;
          z_cate_thatsit = z_cate;
	    
          id_thatsit = id;
	  beta_dummy=beta_rec;
          //kind_thatsit = kind;

          //gamma_rec = 1.0/sqrt(1.0-beta_rec*beta_rec);
          //energy_rec = (gamma_rec-1.0)*931.5;
	}
      }
      if(gamma_det_type==2)
      {
        //cout<<"gamma_det_type==2"<<endl;
        int id_cl = (int)id_cluster_in_miniball;
        //cout<<"id_cl = "<<id_cl<<endl;
        flag_cluster_miniball[id_cl] = 1;
        int id_cr = (int)id_crystal_in_miniball;
        //cout<<"id_cr = "<<id_cr<<endl;
        flag_crystal_miniball[id_cl][id_cr] = 1;

        mul_segment_miniball = mul_segment_miniball+1;
        //cout<<"mul_segment_miniball = "<<mul_segment_miniball<<endl;
        energy_sum_miniball = energy_sum_miniball + energy_notcor;
      
        if(mul_segment_miniball==1) energy_0_miniball = energy_notcor;
        if(mul_segment_miniball==2) energy_1_miniball = energy_notcor;
        if(mul_segment_miniball==3) energy_2_miniball = energy_notcor;

        if(energy_notcor > energy_max_miniball) 
        {
          energy_max_miniball = energy_notcor;
          miniball_cluster_thatsit = id_cl;
          x_thatsit_miniball = x_gamma_det;
          y_thatsit_miniball = y_gamma_det;
          z_thatsit_miniball = z_gamma_det;

          x_target_thatsit = x_ver_rec;
          y_target_thatsit = y_ver_rec;
          z_target_thatsit = z_ver_rec+decay_move_target;
 
          x_cate_thatsit = x_cate;
          y_cate_thatsit = y_cate;
          z_cate_thatsit = z_cate;
	  
	  beta_dummy=beta_rec;
	}
      }
      if(gamma_det_type==3)
      {
        //cout<<"gamma_det_type==3"<<endl;
        mul_crystal_hector = mul_crystal_hector+1;
        //cout<<"mul_crystal_hector =  "<<mul_crystal_hector<<endl;
        //cout<<" id_hector = "<<id_hector<<endl;
        energy_sum_hector = energy_sum_hector + energy_notcor;
      
        if(energy_notcor > energy_max_hector) 
        {
          energy_max_hector = energy_notcor;
          hector_thatsit = id_hector;

          x_thatsit_hector = x_gamma_det; //cout<<" x_thatsit_hector = "<< x_thatsit_hector<<endl;
          y_thatsit_hector = y_gamma_det; //cout<<" y_thatsit_hector = "<< y_thatsit_hector<<endl;
          z_thatsit_hector = z_gamma_det; //cout<<" z_thatsit_hector = "<< z_thatsit_hector<<endl;

          x_target_thatsit = x_ver_rec;
          y_target_thatsit = y_ver_rec;
          z_target_thatsit = z_ver_rec + decay_move_target;
 
          x_cate_thatsit = x_cate;
          y_cate_thatsit = y_cate;
          z_cate_thatsit = z_cate;
	  
	  beta_dummy=beta_rec;
	}
      }
    }
  }
  //Writing into the root file:
  //BETA
  rootfile->cd();
  h_beta_real_gamma_dopp->Write();
  h_beta_before_beta_after->Write();
  h_beta_real_beta_before->Write();
  h_beta_real_beta_after->Write();
  //CATE
    rootfile->cd();
  h_cate_pos->Write();
  rootfile->cd();
  h_e_de->Write();
  /////////////////////////////
  //Cluster!!!!!!!!!!!!!!!!!!!!
  /////////////////////////////
  rootfile->cd();
  h_cluster_crystal_mult->Write();
  rootfile->cd();
  h_cluster_cluster_mult->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_single->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_mult2->Write();
   rootfile->cd();
  h_cluster_dopplercorrected_mult3->Write();
  for(int i=0;i<3;i++)
  {
    rootfile->cd();
    h_cluster_dopplercorrected_single_ring[i]->Write();
  }
  rootfile->cd();
  h_cluster_dopplercorrected_single_outerrings->Write();

  rootfile->cd();
  h_cluster_dopplercorrected_single_no_tracking->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_single_beta_const->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_single_no_tracking_beta_const->Write();

  rootfile->cd();
  h_cluster_notdopplercorrected_total->Write();
  for(int i=0;i<3;i++)
  {
    rootfile->cd();
    h_cluster_notdopplercorrected_total_ring[i]->Write();
  }

  rootfile->cd();
  h_cluster_dopplercorrected_total->Write();
  for(int i=0;i<3;i++)
  {
    rootfile->cd();
    h_cluster_dopplercorrected_total_ring[i]->Write();
  }
  rootfile->cd();
  h_cluster_dopplercorrected_total_outerrings->Write();


  rootfile->cd();
  h_cluster_dopplercorrected_addback->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_addback->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_addback_65->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_addback_70->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_addback_75->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_addback_80->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_addback_85->Write();
  rootfile->cd();
  h_cluster_dopplercorrected_addback_90->Write();
  ////////////////////
  //Miniball!!!!!!!!!!
  ////////////////////
  rootfile->cd();
  h_miniball_segment_mult->Write();
  rootfile->cd();
  h_miniball_crystal_mult->Write();
  rootfile->cd();
  h_miniball_cluster_mult->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_innerring->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_outerring->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_mulseg1->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_mulseg2->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_mulseg3->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_innerring_mulseg2->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_outerring_mulseg2->Write();

  rootfile->cd();
  h_miniball_dopplercorrected_single_no_tracking->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_beta_const->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_single_no_tracking_beta_const->Write();

  rootfile->cd();
  h_miniball_notdopplercorrected->Write();
  for(int i=0;i<2;i++)
  {
    h_miniball_notdopplercorrected_mulseg2_ring[i]->Write();

  }
  for(int i=0;i<2;i++)
  {
    h_miniball_notdopplercorrected_ring[i]->Write();

  }
  rootfile->cd();
  h_miniball_dopplercorrected_total->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_total_innerring->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_total_outerring->Write();
  
  rootfile->cd();
  h_miniball_dopplercorrected_total_mulseg2->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_total_innerring_mulseg2->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_total_outerring_mulseg2->Write();

  rootfile->cd();
  h_miniball_dopplercorrected_addback->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_addback_mulseg2->Write();
  h_miniball_dopplercorrected_addback_innerring_mulseg2->Write();
  rootfile->cd();
  h_miniball_dopplercorrected_addback_outerring_mulseg2->Write();

  ///////////////////////////////
  //Hector!!!!!!!!!!!!!!!!!!!!!
  ///////////////////////////////
  rootfile->cd();
  h_hector_crystal_mult->Write();

  rootfile->cd();
  h_hector_dopplercorrected->Write();
  rootfile->cd();
  h_hector_dopplercorrected_90->Write();
  rootfile->cd();
  h_hector_dopplercorrected_142->Write();
 
  rootfile->cd();
  h_hector_dopplercorrected_no_tracking->Write();
  rootfile->cd();
  h_hector_dopplercorrected_beta_const->Write();
  rootfile->cd();
  h_hector_dopplercorrected_no_tracking_beta_const->Write();

  rootfile->cd();
  h_hector_notdopplercorrected->Write();
  rootfile->cd();
  h_hector_notdopplercorrected_90->Write();
  rootfile->cd();
  h_hector_notdopplercorrected_142->Write(); 
}
