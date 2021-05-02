#include "Globals.hh"
#include "MyNamespace.hh"

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

using namespace std;
using namespace MyNamespace;

//_________________________________________________________________________________________________
float GetDopplerCorrectedEnergy(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                                float target_x,float target_y,float target_z,
                                float pos_det_x, float pos_det_y, float pos_det_z, 
                                float gamma_energy, float beta_rec,float beta_mean_tof,float beta_ave,float decay_z)  {
  float vx1 = gamma_det_x - target_x;                                           // cout<<"vx1: "<<vx1<<endl;
  float vy1 = gamma_det_y - target_y;                                           // cout<<"vy1: "<<vy1<<endl;
  float vz1 = gamma_det_z - target_z - decay_z;                                 // cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  float vx2 = pos_det_x - target_x;                                             // cout<<"vx2: "<<vx2<<endl;
  float vy2 = pos_det_y - target_y;                                             // cout<<"vy2: "<<vy2<<endl;
  float vz2 = pos_det_z - target_z;                                   // cout<<"vy2: "<<vz2<<endl;
 
  double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
          
  // The theta angle:
  double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);                 // cout<<"Theta lab: "<<theta_lab<<endl;
  // the beta value:
  float beta                     = beta_ave + (beta_rec-beta_mean_tof);     // cout<<"beta: "<<beta<<endl;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

//_________________________________________________________________________________________________
float GetTOFGammaNs(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                    float target_x,float target_y,float target_z)  {
  float vx1 = gamma_det_x - target_x;                                           // cout<<"vx1: "<<vx1<<endl;
  float vy1 = gamma_det_y - target_y;                                           // cout<<"vy1: "<<vy1<<endl;
  float vz1 = gamma_det_z - target_z;                                           // cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  return  l1/29.97925;
} 

//_________________________________________________________________________________________________
void GetGammaEnergiesSorted(float gammaEnergy[])  {
  float dummyGammaEnergy[10];
  int rankOfGamma[10]={0};
  for(int i=0;i<10;i++) dummyGammaEnergy[i] = gammaEnergy[i];
  for(int i=0;i<10;i++)  {
    if(gammaEnergy[i]>0)
      for(int j=0;j<10;j++){if(gammaEnergy[i]>dummyGammaEnergy[j]) rankOfGamma[i]++;}
  }
  for(int i=0;i<10;i++) gammaEnergy[i] = dummyGammaEnergy[rankOfGamma[i]];
}


void GetAverageInteractionPointMiniball();
void GetAverageInteractionPointCluster();

double cluster_fi_ave[3][NUMBEROFCLUSTERDETECTORS][7] = {{{0.0}}};
int    cluster_fi_interactions[NUMBEROFCLUSTERDETECTORS][7] = {{0}};
float  cluster_pos_ave[3][NUMBEROFCLUSTERDETECTORS][7];

double miniball_fi_ave[3][NUMBEROFMINIBALLDETECTORS][3][7] = {{{{0.0}}}};
int    miniball_fi_interactions[NUMBEROFMINIBALLDETECTORS][3][7] = {{{0}}};
float  miniball_pos_ave[3][NUMBEROFMINIBALLDETECTORS][3][7];
//_________________________________________________________________________________________________
// Starting the main program://////////////////////////////////////////////////////////////////////
//_________________________________________________________________________________________________
int RisingReconstructor()  { 
  //int main(int argc,char** argv) { 
  //float eventNumberCheck = 0.;
  int clusterInclude                  = 0;
  int miniballInclude                 = 0;
  int hectorInclude                   = 0;
  float eventNumberPrevious           = 0.;
  float gammaEnergy[20];
  
  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
 
  float beta_ave         = 0.4295;        // The average value at the moment of decay for the correction
  float beta_mean_tof    = 0.4295;        // The mean beta from Tof
  float beta_ave_at      = 0.4295;        // The average beta after the target
  float decay_z          = 0.0;           // The lifetime of the excited state moves the average decay point along the beam axis
  int cluster_fi_option  = 0;
  int miniball_fi_option = 0;

  //------------------------------------------------------------------------------------------------------
  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/RisingReconstructor.in","r");
  while(!feof(fin)) {
    fscanf(fin,"%s ",temp); 
    if(strcmp(temp,"INPUTFILE")==0)  {
      fscanf(fin,"%s",&root_in);
      printf("%s %s \n",temp,root_in);
    }
    else if(strcmp(temp,"OUTPUTFILE")==0)  {
      fscanf(fin,"%s ",&root_out); 
      printf("%s %s \n",temp,root_out);
    }
    else if(strcmp(temp,"SPECTRABINANDRANGE")==0)  {
      fscanf(fin,"%i %f %f",&numBin,&firstBin,&lastBin);
      printf("%s %i %f %f\n",temp, numBin,firstBin,lastBin);
    }
    else if(strcmp(temp,"BETADOPPLERAVERAGE")==0)  {
      fscanf(fin,"%f",&beta_ave); 
      printf("%s %f \n",temp,beta_ave);
    }
    else if(strcmp(temp,"BETATOFAVERAGE")==0)  {
      fscanf(fin,"%f",&beta_mean_tof); 
      printf("%s %f \n",temp,beta_mean_tof);
    }
    else if(strcmp(temp,"BETAAFTERTARGETAVERAGE")==0)  {
      fscanf(fin,"%f",&beta_ave_at); 
      printf("%s %f \n",temp,beta_ave_at);
    }
    else if(strcmp(temp,"DECAYPOSITIONZ")==0)  {
      fscanf(fin,"%f",&decay_z); 
      printf("%s %f \n",temp,decay_z);
    }
    else if(strcmp(temp,"STATISTICSREDUCTIONFACTOR")==0)  {
      fscanf(fin,"%f",&reduction_factor); 
      printf("%s %f \n",temp,reduction_factor);
      if(reduction_factor<1) return 0;
    }
    else if(strcmp(temp,"CLUSTERINCLUDE")==0)  {
      fscanf(fin,"%i",&clusterInclude); 
      printf("%s %i \n",temp,clusterInclude);
    }
    else if(strcmp(temp,"MINIBALLINCLUDE")==0)  {
      fscanf(fin,"%i",&miniballInclude); 
      printf("%s %i \n",temp,miniballInclude);
    }
    else if(strcmp(temp,"HECTORINCLUDE")==0)  {
      fscanf(fin,"%i",&hectorInclude); 
      printf("%s %i \n",temp,hectorInclude);
    }
    else if(strcmp(temp,"CLUSTERFIFIND")==0)  {
      fscanf(fin,"%i",&cluster_fi_option); 
      printf("%s %i \n",temp,cluster_fi_option);
    }
    else if(strcmp(temp,"MINIBALLFIFIND")==0)  {
      fscanf(fin,"%i",&miniball_fi_option); 
      printf("%s %i \n",temp,miniball_fi_option);
    }
    else if(strcmp(temp,"END")==0) break;
    else  {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      return 0;
    }
  } 
 
  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();
  
  TH1F *h_number_of_observed_decays = new TH1F("number_of_observed_decays","",100,-0.5,99.5);
  TH2F *h_decayId_vs_energy_sum = new TH2F("decayId_vs_energy_sum","",10,0,10,1000,0,5000);
  
  TH1F *h_cluster_crystal_mult;
  TH1F *h_cluster_detector_mult;
  TH1F *h_cluster_single;
  TH1F *h_cluster_energy_sum;
  TH1F *h_cluster_ring[3];
  TH2F *h_cluster_gamma_gamma[10];

  TH1F *h_cluster_doppler_single;
  TH1F *h_cluster_doppler_total;
  TH1F *h_cluster_doppler_total_ring[3];

  TH1F *h_cluster_doppler_single_ave_z;
  TH1F *h_cluster_doppler_single_ave_z_ave_fi;
  TH1F *h_cluster_doppler_single_fi;
  TH1F *h_cluster_doppler_single_fi_ave_z;

  TH1F *h_cluster_doppler_single_at;
  TH1F *h_cluster_doppler_single_at_fep;
  TH1F *h_cluster_doppler_single_fi_at;

  TH1F *h_cluster_doppler_single_fep;  // Take only events with FEP.
  TH1F *h_cluster_doppler_single_ave_z_fep;  // Take only events with FEP.
  TH1F *h_cluster_doppler_single_ave_z_ave_fi_fep;  // Take only events with FEP.

  if(clusterInclude==1)  {
    // Crystal multiplicity,ok
    h_cluster_crystal_mult = new TH1F("cluster_crystal_mult","",100,0,100);
    h_cluster_detector_mult = new TH1F("cluster_detector_mult","",100,0,100);
    //Single hits, OK
    h_cluster_single = new TH1F("cluster_single","",numBin,firstBin,lastBin);
    // All, OK, This includes addback events.
    h_cluster_energy_sum = new TH1F("cluster_energy_sum","",numBin,firstBin,lastBin);

    h_cluster_doppler_single = new TH1F("cluster_doppler_single","",numBin,firstBin,lastBin);
    h_cluster_doppler_single_ave_z = new TH1F("cluster_doppler_single_ave_z","",numBin,firstBin,lastBin);
    h_cluster_doppler_single_at = new TH1F("cluster_doppler_single_at","",numBin,firstBin,lastBin);
    h_cluster_doppler_single_at_fep = new TH1F("cluster_doppler_single_at_fep","",numBin,firstBin,lastBin);

    h_cluster_doppler_single_fep = new TH1F("cluster_doppler_single_fep","",numBin,firstBin,lastBin);   
    h_cluster_doppler_single_ave_z_fep = new TH1F("cluster_doppler_single_ave_z_fep","",numBin,firstBin,lastBin);   
    h_cluster_doppler_single_ave_z_ave_fi_fep = new TH1F("cluster_doppler_single_ave_z_ave_fi_fep","",numBin,firstBin,lastBin);   

    //______________________________________________
    if(cluster_fi_option==2)  {
      FILE *fInAve  = fopen("./output/ClusterAverageInteractionPoint.out","r");
      int counter = 0;
      float x,y,z;
      int counts;
      int crystalNumber = 0;
      while(!feof(fInAve) && counter<NUMBEROFCLUSTERDETECTORS && crystalNumber<(NUMBEROFCLUSTERDETECTORS*7)){
        fscanf(fInAve,"%i %i %f %f %f",&counter,&counts,&x,&y,&z); 
        //counter++;                                            
        cluster_pos_ave[0][counter][crystalNumber % 7] = x;
        cluster_pos_ave[1][counter][crystalNumber % 7] = y;
        cluster_pos_ave[2][counter][crystalNumber % 7] = z;
        //cout<<counter<<" "<<cluster_pos_ave[0][counter][crystalNumber % 7]
        //    <<" "<<cluster_pos_ave[1][counter][crystalNumber % 7]
        //    <<" "<<cluster_pos_ave[2][counter][crystalNumber % 7] <<endl;
        crystalNumber++;
        //cout<<crystalNumber<<endl;
      }
      h_cluster_doppler_single_ave_z_ave_fi = new TH1F("cluster_doppler_single_ave_z_ave_fi","",numBin,firstBin,lastBin);
      h_cluster_doppler_single_fi = new TH1F("cluster_doppler_single_fi","",numBin,firstBin,lastBin);
      h_cluster_doppler_single_fi_ave_z = new TH1F("cluster_doppler_single_fi_ave_z","",numBin,firstBin,lastBin);
      h_cluster_doppler_single_fi_at = new TH1F("cluster_doppler_single_fi_at","",numBin,firstBin,lastBin);
    }
    h_cluster_doppler_total = new TH1F("cluster_doppler_total","",numBin,firstBin,lastBin);

    for(int i=0;i<3;i++)  {
      sprintf(temp,"cluster_doppler_total_ring[%i]",i);
      h_cluster_doppler_total_ring[i] = new TH1F(temp,"",numBin,firstBin,lastBin);
    }
    for(int i=0;i<3;i++)  {
      sprintf(temp,"cluster_ring[%i]",i);
      h_cluster_ring[i] = new TH1F(temp,"",numBin,firstBin,lastBin);
    }
    //---------------------------------------
    for(int i=0;i<10;i++)  {
      sprintf(temp,"cluster_gamma_gamma[%i]",i);
      h_cluster_gamma_gamma[i]= new TH2F(temp,temp,numBin,firstBin,lastBin,numBin,firstBin,lastBin);
    }
  }

  //_______________________________________________________________________________________________ 
  TH1F *h_miniball_segment_mult;
  TH1F *h_miniball_crystal_mult;
  TH1F *h_miniball_detector_mult;

  TH1F *h_miniball_doppler_single;              // i.e. single crystal
  //TH1F *h_miniball_doppler_single_innerring;
  //TH1F *h_miniball_doppler_single_outerring;

  //TH1F *h_miniball_doppler_single_mulseg1;
  //TH1F *h_miniball_doppler_single_mulseg2;
  //TH1F *h_miniball_doppler_single_mulseg3;

  //TH1F *h_miniball_doppler_single_innerring_mulseg2;
  //TH1F *h_miniball_doppler_single_outerring_mulseg2;

  TH1F *h_miniball;
  //TH1F *h_miniball_mulseg2_ring[2];
  TH1F *h_miniball_ring[2];

  TH1F *h_miniball_doppler_total;
  TH1F *h_miniball_doppler_total_ring[2];
  //TH1F *h_miniball_doppler_total_mulseg2;
  //TH1F *h_miniball_doppler_addback;

  TH1F *h_miniball_doppler_single_ave_z;
  TH1F *h_miniball_doppler_single_ave_z_ave_fi;
  TH1F *h_miniball_doppler_single_fi;
  TH1F *h_miniball_doppler_single_fi_ave_z;

  TH1F *h_miniball_doppler_single_at;       // Using the beta after the taret (at)
  TH1F *h_miniball_doppler_single_fi_at;
  TH1F *h_miniball_doppler_single_at_fep;

  TH1F *h_miniball_doppler_single_fep;  // Take only events with FEP.
  TH1F *h_miniball_doppler_single_ave_z_fep;  // Take only events with FEP
  TH1F *h_miniball_doppler_single_ave_z_ave_fi_fep;  // Take only events with FEP

  if(miniballInclude==1)  {
    h_miniball_segment_mult  = new TH1F("miniball_segment_mult" ,"",15,0,15);
    h_miniball_crystal_mult  = new TH1F("miniball_crystal_mult" ,"",15,0,15);
    h_miniball_detector_mult = new TH1F("miniball_detector_mult","",15,0,15);

    //Single hits, i.e. only one crystal, several segments possible
    h_miniball_doppler_single = new TH1F("miniball_doppler_single","",numBin,firstBin,lastBin);
    h_miniball_doppler_single_ave_z = new TH1F("miniball_doppler_single_ave_z","",numBin,firstBin,lastBin);
    h_miniball_doppler_single_fep = new TH1F("miniball_doppler_single_fep","",numBin,firstBin,lastBin);
    h_miniball_doppler_single_ave_z_fep = new TH1F("miniball_doppler_single_ave_z_fep","",numBin,firstBin,lastBin);
    h_miniball_doppler_single_ave_z_ave_fi_fep = new TH1F("miniball_doppler_single_ave_z_ave_fi_fep","",numBin,firstBin,lastBin);
 
    h_miniball_doppler_single_at = new TH1F("miniball_doppler_single_at","",numBin,firstBin,lastBin);
    h_miniball_doppler_single_at_fep = new TH1F("miniball_doppler_single_at_fep","",numBin,firstBin,lastBin);
    //______________________________________________
    if(miniball_fi_option==2)  {
      FILE *fInAve  = fopen("./output/MiniballAverageInteractionPoint.out","r");
      int counter = 0;
      float x,y,z;
      int counts;
      int numberOfSegments = 0;
      while(!feof(fInAve) && counter<NUMBEROFMINIBALLDETECTORS && numberOfSegments<(NUMBEROFMINIBALLDETECTORS*3*7)){
        fscanf(fInAve,"%i %i %f %f %f",&counter,&counts,&x,&y,&z); 
        //counter++;                                            
        miniball_pos_ave[0][counter][((int)numberOfSegments/7) % 3][numberOfSegments % 7] = x;
        miniball_pos_ave[1][counter][((int)numberOfSegments/7) % 3][numberOfSegments % 7] = y;
        miniball_pos_ave[2][counter][((int)numberOfSegments/7) % 3][numberOfSegments % 7] = z;
        cout<<counter<<" "<<miniball_pos_ave[0][counter][((int)numberOfSegments/7) %3][numberOfSegments % 7]
            <<" "<<miniball_pos_ave[1][counter][((int)numberOfSegments/7) %3][numberOfSegments % 7]
            <<" "<<miniball_pos_ave[2][counter][((int)numberOfSegments/7) %3][numberOfSegments % 7] <<endl;
        numberOfSegments++;
      }
      h_miniball_doppler_single_ave_z_ave_fi = new TH1F("miniball_doppler_single_ave_z_ave_fi","",numBin,firstBin,lastBin);
      h_miniball_doppler_single_fi = new TH1F("miniball_doppler_single_fi","",numBin,firstBin,lastBin);
      h_miniball_doppler_single_fi_ave_z = new TH1F("miniball_doppler_single_fi_ave_z","",numBin,firstBin,lastBin);

      h_miniball_doppler_single_fi_at = new TH1F("miniball_doppler_single_fi_at","",numBin,firstBin,lastBin);
    }
    //h_miniball_doppler_single_innerring = new TH1F("Miniball_Doppler_Single_InnerRing","",numBin,firstBin,lastBin);
    //h_miniball_doppler_single_outerring = new TH1F("Miniball_Doppler_Single_OuterRing","",numBin,firstBin,lastBin);
 
    //h_miniball_doppler_single_mulseg1 = new TH1F("Miniball_Doppler_SingleHit_MulSeg1","",numBin,firstBin,lastBin);
    //h_miniball_doppler_single_mulseg2 = new TH1F("Miniball_Doppler_SingleHit_MulSeg2","",numBin,firstBin,lastBin);
    //h_miniball_doppler_single_mulseg3 = new TH1F("Miniball_Doppler_SingleHit_MulSeg3","",numBin,firstBin,lastBin);				
  
    //h_miniball_doppler_single_innerring_mulseg2 = new TH1F("Miniball_Doppler_SingleHit_InnerRing_MulSeg2","",numBin,firstBin,lastBin);
    //h_miniball_doppler_single_outerring_mulseg2 = new TH1F("Miniball_Doppler_SingleHit_OuterRing_MulSeg2","",numBin,firstBin,lastBin);

    //Not doppler corrected
    h_miniball = new TH1F("miniball","",numBin,firstBin,lastBin);
    /*
      for(int i=0;i<2;i++)
      {
      sprintf(name,"h_miniball_mulseg2_ring[%i]",i);
      h_miniball_mulseg2_ring[i] = new TH1F(name,"",numBin,firstBin,lastBin);
      }*/
    for(int i=0;i<2;i++)  {
      sprintf(temp,"miniball_ring[%i]",i);
      h_miniball_ring[i] = new TH1F(temp,"",numBin,firstBin,lastBin);
    }
    // Single hit AND addback
    //No Cut on segments
    h_miniball_doppler_total = new TH1F("miniball_doppler_total","",numBin,firstBin,lastBin);
    for(int i=0;i<2;i++)  {
      sprintf(temp,"miniball_doppler_total_ring[%i]",i);
      h_miniball_doppler_total_ring[i] = new TH1F(temp,"",numBin,firstBin,lastBin);
    }
    //h_miniball_doppler_total_mulseg2 = new TH1F("Miniball_Doppler_Total_MulSeg2"  ,"",numBin,firstBin,lastBin);
  
    //h_miniball_doppler_addback = new TH1F("Miniball_Doppler_Addback","",numBin,firstBin,lastBin);
  }
					       		      
  ///////////////////////////////
  //Hector!!!!!!!!!!!!!!!!!!!!!
  ///////////////////////////////
  TH1F *h_hector_crystal_mult;
  TH1F *h_hector_doppler;
  TH1F *h_hector_doppler_ring[2];
  
  TH1F *h_hector;
  TH1F *h_hector_ring[2];

  if(hectorInclude==1)  {
    h_hector_crystal_mult = new TH1F("hector_crystal_mult","",15,0,15);
    
    h_hector_doppler = new TH1F("hector_doppler","",numBin,firstBin,lastBin);
    for(int i=0;i<2;i++)  {
      sprintf(temp,"h_hector_doppler_ring[%i]",i);
      h_hector_doppler_ring[i] = new TH1F(temp,"",numBin,firstBin,lastBin);
    }
    h_hector = new TH1F("hector","",numBin,firstBin,lastBin);
    for(int i=0;i<2;i++)  {
      sprintf(temp,"hector_ring[%i]",i);
      h_hector_ring[i] = new TH1F(temp,"",numBin,firstBin,lastBin);
    }
  }
   
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile *infile = new TFile(root_in,"READ");
  TTree *t = (TTree*)infile->Get("ObservedEvents");

  t->SetBranchAddress("EventNumber",&evnum);
  t->SetBranchAddress("DecayIDOfEventNumber",&decayIDevnum);
  t->SetBranchAddress("ProjectileVertex",p0);                                   // Position at the fragmentation point
  t->SetBranchAddress("VProjectileBeforeTarget",pvb);                           // Normalized Vector of beam before the target
  t->SetBranchAddress("VProjectileAfterTarget",pva);                            // Normalized Vector of beam after the target
  t->SetBranchAddress("EnergyProjectile",&energy_p);                            // Energy of beam before the target in MeV/u
  t->SetBranchAddress("BetaBeforeTarget",&beta_b);                              // Beta of beam before the target
  t->SetBranchAddress("BetaReal",&beta_r);                                      // Beta of fragment during deexcitation	
  t->SetBranchAddress("BetaAfterTarget",&beta_a);                               // Beta of fragment After Target
  t->SetBranchAddress("Halflife",&halflife);                                    // Halflife
  t->SetBranchAddress("DecayTimeAfterInteraction",&decay_time_after_interaction);
  t->SetBranchAddress("VertexGamma",g0);                                        // Position at the gamma emmittance point
  t->SetBranchAddress("VGamma",gv);		                                // Gamma vector
  t->SetBranchAddress("EGammaRest",&e_rest);		                        // Energy at rest
  t->SetBranchAddress("EGammaDoppler",&e_doppler);                              // Theta of doppler boosted gamma
  t->SetBranchAddress("ThetaGammaRest",&theta_gamma_rest);
  t->SetBranchAddress("ThetaGammaLab",&theta_gamma_lab);
  t->SetBranchAddress("EnergyVertex",&energy_vertex_stored);                    // Energy of fragment at fragmentation 
  t->SetBranchAddress("GammaDetType",&gamma_det_type);
  if(clusterInclude==1)  {
    t->SetBranchAddress("ClusterFlag",cluster_flag);                            // The Cluster detectors
    t->SetBranchAddress("ClusterEnergyNotCor",cluster_energy_not_cor);
    t->SetBranchAddress("ClusterTime",cluster_time);
    if(cluster_fi_option >= 1)  {
      t->SetBranchAddress("ClusterFI",cluster_fi);
    }
  }
  if(miniballInclude==1)  {
    t->SetBranchAddress("MiniballFlag",miniball_flag);                          // The Miniball detectors
    t->SetBranchAddress("MiniballEnergyNotCor",miniball_energy_not_cor);
    t->SetBranchAddress("MiniballTime",miniball_time);
    if(miniball_fi_option >= 1)  {
      t->SetBranchAddress("MiniballFI",miniball_fi);
    }
  }
  if(hectorInclude==1)  {
    t->SetBranchAddress("HectorFlag",hector_flag);                              // The Hector detectors
    t->SetBranchAddress("HectorEnergyNotCor",hector_energy_not_cor);
    t->SetBranchAddress("HectorTime",hector_time);  
  }
  t->SetBranchAddress("PosDet",pos_det);                                        // Reconstructed position 
  t->SetBranchAddress("VertexReconstructed",ver_rec);                           // reconstructed vertex postion, 
  t->SetBranchAddress("BetaReconstructed",&beta_rec);                           // reconstructed beta, including resolution

  //------------------------------------------------------------------------------------------------------
  //Going to the Header file, which has the information of constant values.
  infile->cd();
  TTree *tHeader = (TTree*)infile->Get("Header");
  tHeader->SetBranchAddress("Mass",&mass);                                      // Beam mass
  tHeader->SetBranchAddress("Z",&z);                                            // Beam z	
  tHeader->SetBranchAddress("Charge",&charge);                                  // Beam charge	
  tHeader->SetBranchAddress("MassFragment",&mass_f);                            // Fragment mass
  tHeader->SetBranchAddress("ZFragment",&z_f);	                                // Fragment z
  tHeader->SetBranchAddress("ChargeFragment",&charge_f);                        // Fragment charge
  tHeader->SetBranchAddress("TargetKind",&kind_target);
  tHeader->SetBranchAddress("TargetThicknessCM",&thickness);
  tHeader->SetBranchAddress("TotalEventNumber",&total_event_num);
  tHeader->SetBranchAddress("ThetaLow",&theta_low);
  tHeader->SetBranchAddress("ThetaHigh",&theta_high);
  if(clusterInclude==1)  {
    tHeader->SetBranchAddress("ClusterEnResOpt",&cluster_en_res_opt);
    tHeader->SetBranchAddress("ClusterEnRes",cluster_en_res);
    tHeader->SetBranchAddress("ClusterTimeRes",cluster_time_res);
    tHeader->SetBranchAddress("ClusterPos",cluster_pos);
  }
  if(miniballInclude==1)  {
    tHeader->SetBranchAddress("MiniballEnResOpt",&miniball_en_res_opt);
    tHeader->SetBranchAddress("MiniballEnRes",miniball_en_res);
    tHeader->SetBranchAddress("MiniballTimeRes",miniball_time_res);
    tHeader->SetBranchAddress("MiniballPos",miniball_pos);
  }
  if(hectorInclude==1)
    {
      tHeader->SetBranchAddress("HectorEnResOpt",&hector_en_res_opt);
      tHeader->SetBranchAddress("HectorEnRes",hector_en_res);
      tHeader->SetBranchAddress("HectorTimeRes",hector_time_res);
      tHeader->SetBranchAddress("HectorPos",hector_pos);
    }
  tHeader->SetBranchAddress("BetaResolution",&beta_res); 
  tHeader->SetBranchAddress("PosDetAtTargetRes",&pos_det_at_target_res); 
  tHeader->SetBranchAddress("PosDetAfterTargetRes",&pos_det_after_target_res);
  tHeader->GetEntry(0);
  //-----------------------------------------------------------------------------------------------------------
  // Finished reading the Root file

  nentries = (Int_t)(t->GetEntries()/reduction_factor);
  cout <<"Entries: "<< nentries << endl;
  // Variables that are determined event by event:
  int cluster_crystal_mult  = 0;
  int cluster_detector_mult = 0;
  float cluster_energy_sum  = 0.;
  float cluster_energy_max  = -999.9;
  bool cluster_fired[NUMBEROFCLUSTERDETECTORS] = {false};
  
  int miniball_segment_mult  = 0;
  int miniball_crystal_mult  = 0;
  int miniball_detector_mult = 0; 
  float miniball_energy_sum  = 0.;
  float miniball_energy_max  = -999.9;
  bool miniball_fired[NUMBEROFMINIBALLDETECTORS] = {false};

  int hector_crystal_mult  = 0;
  float hector_energy_sum  = 0.;
  float hector_energy_max  = -999.9;
  bool hector_fired[NUMBEROFHECTORDETECTORS] = {false};

  float cluster_x_thatsit  = 0.;
  float cluster_y_thatsit  = 0.;
  float cluster_z_thatsit  = 0.;
  int   cluster_detector_HE= 0;  // The detector that recorded the highest energy in the Event;
  int   cluster_crystal_HE = 0;  // The according crystal
  float miniball_x_thatsit = 0.;
  float miniball_y_thatsit = 0.;
  float miniball_z_thatsit = 0.;
  int   miniball_detector_HE= 0;  // The detector that recorded the highest energy in the Event;
  int   miniball_crystal_HE = 0;  // The according crystal
  int   miniball_segment_HE = 0;  // The according segment
  float hector_x_thatsit   = 0.;
  float hector_y_thatsit   = 0.;
  float hector_z_thatsit   = 0.;

  float gamma_tof = 0.;

  int observedDecayNumber = 0;
  //---------------------------------------------------------------------------------------------------------
  //Gates for the analysis:
  float energy_gate_lower[20] = {120.,240.,315.,695.};
  float energy_gate_upper[20] = {140.,265.,340.,750.};

  cout<<"Starting the analysis"<<endl;

  float doppler_energy              = 0.; 
  float doppler_energy_ave_z        = 0.;
  float doppler_energy_ave_z_ave_fi = 0.;
  float doppler_energy_fi           = 0.;
  float doppler_energy_fi_ave_z     = 0.;
  float doppler_energy_fi_fe        = 0.;    // fe = full energy
  float doppler_energy_at           = 0.;
  float doppler_energy_fi_at        = 0.;
  
  for(int iii=0;iii<nentries;iii++) {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    t->GetEntry(iii);

    //-------------------------------------------------------------
    // Starting with the cluster analysis:-------------------------
    //-------------------------------------------------------------
    int dummy = -1;
    if(clusterInclude==1)  {
      for(int nnn=0;nnn<NUMBEROFCLUSTERDETECTORS;nnn++)  {
        for(int j=0;j<7;j++)  {
          if(cluster_flag[nnn][j]==1.0)  {
            if(dummy!=nnn) cluster_detector_mult++;
            dummy  = nnn;
            cluster_crystal_mult++;
            cluster_fired[nnn] = true;
            //summed energy
            cluster_energy_sum = cluster_energy_sum+cluster_energy_not_cor[nnn][j];
            // getting the MI
            if(cluster_energy_not_cor[nnn][j]>cluster_energy_max)  {
              cluster_energy_max = cluster_energy_not_cor[nnn][j];
              cluster_x_thatsit  = cluster_pos[0][nnn][j];
              cluster_y_thatsit  = cluster_pos[1][nnn][j];
              cluster_z_thatsit  = cluster_pos[2][nnn][j];
              cluster_crystal_HE = j;
              cluster_detector_HE = nnn;
            }
          }
        }
      }
      h_cluster_detector_mult->Fill(cluster_detector_mult);
      h_cluster_crystal_mult ->Fill(cluster_crystal_mult);

      if(cluster_crystal_mult>0 && cluster_detector_mult==1)  {
        // Performing the Doppler-correction from the angles
        // of the detector that registered the highest gamma-ray energy.
        //if(cluster_fi_option !=2)
        doppler_energy = GetDopplerCorrectedEnergy(cluster_x_thatsit,cluster_y_thatsit,cluster_z_thatsit,
                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2],
                                                   cluster_energy_sum,beta_rec,beta_mean_tof,beta_ave,0);
        
        doppler_energy_ave_z = GetDopplerCorrectedEnergy(cluster_x_thatsit,cluster_y_thatsit,cluster_z_thatsit,
                                                         ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2],
                                                         cluster_energy_sum,beta_rec,beta_mean_tof,beta_ave,decay_z);
        doppler_energy_at = GetDopplerCorrectedEnergy(cluster_x_thatsit,cluster_y_thatsit,cluster_z_thatsit,
                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2],
                                                   cluster_energy_sum,beta_rec,beta_mean_tof,beta_ave_at,0);

        //else 
        if(cluster_fi_option ==2) {
          doppler_energy_ave_z_ave_fi = GetDopplerCorrectedEnergy(cluster_pos_ave[0][cluster_detector_HE][cluster_crystal_HE],
                                                                  cluster_pos_ave[1][cluster_detector_HE][cluster_crystal_HE],
                                                                  cluster_pos_ave[2][cluster_detector_HE][cluster_crystal_HE],
                                                                  0,0,ver_rec[2],0,0,pos_det[2],
                                                                  cluster_energy_sum,beta_rec,beta_mean_tof,beta_ave,decay_z);
          doppler_energy_fi = GetDopplerCorrectedEnergy(cluster_fi[0][cluster_detector_HE][cluster_crystal_HE],
                                                        cluster_fi[1][cluster_detector_HE][cluster_crystal_HE],
                                                        cluster_fi[2][cluster_detector_HE][cluster_crystal_HE],
                                                        0,0,ver_rec[2],0,0,pos_det[2],
                                                        cluster_energy_sum,beta_rec,beta_mean_tof,beta_ave,0);
          doppler_energy_fi_ave_z = GetDopplerCorrectedEnergy(cluster_fi[0][cluster_detector_HE][cluster_crystal_HE],
                                                              cluster_fi[1][cluster_detector_HE][cluster_crystal_HE],
                                                              cluster_fi[2][cluster_detector_HE][cluster_crystal_HE],
                                                              0,0,ver_rec[2],0,0,pos_det[2],
                                                              cluster_energy_sum,beta_rec,beta_mean_tof,beta_ave,decay_z);
          doppler_energy_fi_fe = GetDopplerCorrectedEnergy(cluster_fi[0][cluster_detector_HE][cluster_crystal_HE],
                                                           cluster_fi[1][cluster_detector_HE][cluster_crystal_HE],
                                                           cluster_fi[2][cluster_detector_HE][cluster_crystal_HE],
                                                           0,0,ver_rec[2],0,0,pos_det[2],
                                                           e_doppler,beta_rec,beta_mean_tof,beta_ave,0);
          doppler_energy_fi_at = GetDopplerCorrectedEnergy(cluster_fi[0][cluster_detector_HE][cluster_crystal_HE],
                                                        cluster_fi[1][cluster_detector_HE][cluster_crystal_HE],
                                                        cluster_fi[2][cluster_detector_HE][cluster_crystal_HE],
                                                        0,0,ver_rec[2],0,0,pos_det[2],
                                                        cluster_energy_sum,beta_rec,beta_mean_tof,beta_ave_at,0);
        }
        
        h_cluster_doppler_total->Fill(doppler_energy);
        if(evnum == eventNumberPrevious) gammaEnergy[observedDecayNumber+1] = cluster_energy_sum;
      
        h_cluster_energy_sum->Fill(cluster_energy_sum);
        h_decayId_vs_energy_sum->Fill(decayIDevnum,cluster_energy_sum);
        if(cluster_crystal_mult==1)  {
          h_cluster_single->Fill(cluster_energy_sum);

          h_cluster_doppler_single->Fill(doppler_energy);
          h_cluster_doppler_single_ave_z->Fill(doppler_energy_ave_z);
          h_cluster_doppler_single_at->Fill(doppler_energy_at);
          //_________________________________________________
          // Checking the effect of getting the average position:
          if(cluster_fi_option ==2)  {
            h_cluster_doppler_single_ave_z_ave_fi->Fill(doppler_energy_ave_z_ave_fi);
            h_cluster_doppler_single_fi->Fill(doppler_energy_fi);
            h_cluster_doppler_single_fi_ave_z->Fill(doppler_energy_fi_ave_z);
            h_cluster_doppler_single_fi_at->Fill(doppler_energy_fi_at);

            if(doppler_energy_fi_fe-5< doppler_energy_fi){
              h_cluster_doppler_single_fep->Fill(doppler_energy);
              h_cluster_doppler_single_ave_z_fep->Fill(doppler_energy_ave_z);
              h_cluster_doppler_single_ave_z_ave_fi_fep->Fill(doppler_energy_ave_z_ave_fi);
              h_cluster_doppler_single_at_fep->Fill(doppler_energy_at);
            }
          }
          //________________________________________________
          // Only if almost all the energy was released within one crystal the average position is determined:
          if(cluster_fi_option >= 1 && doppler_energy>= e_rest/1.1)  {
            for(int jjj=0;jjj<3;jjj++)  {
              cluster_fi_ave[jjj][cluster_detector_HE][cluster_crystal_HE] = 
                cluster_fi_ave[jjj][cluster_detector_HE][cluster_crystal_HE] + cluster_fi[jjj][cluster_detector_HE][cluster_crystal_HE]; 
              if(jjj==0) cluster_fi_interactions[cluster_detector_HE][cluster_crystal_HE]++;
            }
          }
        }
        
        for(int j=0;j<15;j++)  {
          if(j<5 && cluster_fired[j] == true)  {
            h_cluster_ring[0]              ->Fill(cluster_energy_sum);
            h_cluster_doppler_total_ring[0]->Fill(doppler_energy);
          } 
          if(j>4 && j<10 && cluster_fired[j] == true)  {
            h_cluster_ring[1]              ->Fill(cluster_energy_sum);
            h_cluster_doppler_total_ring[1]->Fill(doppler_energy);
          }
          if(j>10 && cluster_fired[j] == true)  {
            h_cluster_ring[2]              ->Fill(cluster_energy_sum);
            h_cluster_doppler_total_ring[2]->Fill(doppler_energy);
          }
        }
      }
    }
    //____________________________________________________________________________________________
    // Clusters Finished. Analyzing Miniball
    dummy = -1;
    if(miniballInclude==1)  {
      for(int nnn=0;nnn<NUMBEROFMINIBALLDETECTORS;nnn++)  {
        for(int j=0;j<3;j++)  {
          for(int k=0;k<7;k++)  {
            if(miniball_flag[nnn][j][k]==1.0 && k==0)  {
              if(dummy!=nnn) miniball_detector_mult++;
              miniball_crystal_mult++;
              miniball_fired[nnn] = true;
              //summed energy
              miniball_energy_sum = miniball_energy_sum+miniball_energy_not_cor[nnn][j][0];
              dummy  = nnn;
            }
            if(miniball_flag[nnn][j][k]==1.0 && k!=0)  {  
              miniball_segment_mult++;
              // getting the MI
              if(miniball_energy_not_cor[nnn][j][k]>miniball_energy_max)  {
                miniball_energy_max = miniball_energy_not_cor[nnn][j][k];
                miniball_x_thatsit  = miniball_pos[0][nnn][j][k];
                miniball_y_thatsit  = miniball_pos[1][nnn][j][k];
                miniball_z_thatsit  = miniball_pos[2][nnn][j][k];
                miniball_detector_HE= nnn; //cout<<"nnn: "<<nnn<<endl;
                miniball_crystal_HE = j;   //cout<<"j: "<<j<<endl;
                miniball_segment_HE = k;   //cout<<"l: "<<k<<endl;
              }
             }
          }
        }
      }
      h_miniball_detector_mult->Fill(miniball_detector_mult);
      h_miniball_crystal_mult ->Fill(miniball_crystal_mult);
      h_miniball_segment_mult ->Fill(miniball_segment_mult);

      if(miniball_crystal_mult>0 && miniball_detector_mult==1) {
        // Performing the Doppler-correction from the angles
        // of the detector that registered the highest gamma-ray energy.
        doppler_energy = GetDopplerCorrectedEnergy(miniball_x_thatsit,miniball_y_thatsit,miniball_z_thatsit,
                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2],
                                                   miniball_energy_sum,beta_rec,beta_mean_tof,beta_ave,0);
        
        doppler_energy_ave_z = GetDopplerCorrectedEnergy(miniball_x_thatsit,miniball_y_thatsit,miniball_z_thatsit,
                                                         ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2],
                                                         miniball_energy_sum,beta_rec,beta_mean_tof,beta_ave,decay_z);
        doppler_energy_at = GetDopplerCorrectedEnergy(miniball_x_thatsit,miniball_y_thatsit,miniball_z_thatsit,
                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2],
                                                   miniball_energy_sum,beta_rec,beta_mean_tof,beta_ave_at,0);
        //else 
        if(miniball_fi_option ==2) {
          doppler_energy_ave_z_ave_fi=GetDopplerCorrectedEnergy(miniball_pos_ave[0][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                                miniball_pos_ave[1][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                                miniball_pos_ave[2][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                                0,0,ver_rec[2],0,0,pos_det[2],
                                                                miniball_energy_sum,beta_rec,beta_mean_tof,beta_ave,decay_z);
          doppler_energy_fi = GetDopplerCorrectedEnergy(miniball_fi[0][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                        miniball_fi[1][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                        miniball_fi[2][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                        0,0,ver_rec[2],0,0,pos_det[2],
                                                        miniball_energy_sum,beta_rec,beta_mean_tof,beta_ave,0);
          doppler_energy_fi_ave_z = GetDopplerCorrectedEnergy(miniball_fi[0][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                              miniball_fi[1][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                              miniball_fi[2][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                              0,0,ver_rec[2],0,0,pos_det[2],
                                                              miniball_energy_sum,beta_rec,beta_mean_tof,beta_ave,decay_z);
          doppler_energy_fi_fe = GetDopplerCorrectedEnergy(miniball_fi[0][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                           miniball_fi[1][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                           miniball_fi[2][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                           0,0,ver_rec[2],0,0,pos_det[2],
                                                           e_doppler,beta_rec,beta_mean_tof,beta_ave,0);
          doppler_energy_fi_at = GetDopplerCorrectedEnergy(miniball_fi[0][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                        miniball_fi[1][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                        miniball_fi[2][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE],
                                                        0,0,ver_rec[2],0,0,pos_det[2],
                                                        miniball_energy_sum,beta_rec,beta_mean_tof,beta_ave_at,0);
        }

        h_miniball              ->Fill(miniball_energy_sum);
        h_miniball_doppler_total->Fill(doppler_energy);

        if(miniball_fired[0]==true || miniball_fired[2]==true || miniball_fired[5]==true || miniball_fired[7]==true)  {
          h_miniball_ring[0]              ->Fill(miniball_energy_sum);
          h_miniball_doppler_total_ring[0]->Fill(doppler_energy);
        }
        if(miniball_fired[1]==true || miniball_fired[3]==true || miniball_fired[4]==true || miniball_fired[6]==true)  {
          h_miniball_ring[1]              ->Fill(miniball_energy_sum);
          h_miniball_doppler_total_ring[1]->Fill(doppler_energy);
        }
        if(miniball_crystal_mult==1 && miniball_segment_mult==1 
           //&&( miniball_fired[0]==true || miniball_fired[1]==true ||  miniball_fired[6]==true || miniball_fired[7]==true) 
           )  {
          h_miniball_doppler_single->Fill(doppler_energy);
          h_miniball_doppler_single_at->Fill(doppler_energy_at);
          h_miniball_doppler_single_ave_z->Fill(doppler_energy_ave_z);
          //_________________________________________________
          // Checking the effect of getting the average position:
          if(miniball_fi_option ==2)  {
            h_miniball_doppler_single_ave_z_ave_fi->Fill(doppler_energy_ave_z_ave_fi);
            h_miniball_doppler_single_fi->Fill(doppler_energy_fi);
            h_miniball_doppler_single_fi_at->Fill(doppler_energy_fi_at);
            h_miniball_doppler_single_fi_ave_z->Fill(doppler_energy_fi_ave_z);
            if(doppler_energy_fi_fe-5<doppler_energy_fi)  {
              h_miniball_doppler_single_fep->Fill(doppler_energy);
              h_miniball_doppler_single_ave_z_fep->Fill(doppler_energy_ave_z);
              h_miniball_doppler_single_ave_z_ave_fi_fep->Fill(doppler_energy_ave_z_ave_fi);
              h_miniball_doppler_single_at_fep->Fill(doppler_energy_at);
            } 
          }
          //________________________________________________
          // Only if almost all the energy was released within one crystal the average position is determined:
          if(miniball_fi_option >= 1 && doppler_energy>= e_rest/1.1)  {
            for(int jjj=0;jjj<3;jjj++)  {
              miniball_fi_ave[jjj][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE] = 
                miniball_fi_ave[jjj][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE] + 
                miniball_fi[jjj][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE]; 
              //cout<<"miniball_fi_ave["<<jjj<<"]["<<miniball_detector_HE<<"]["<<miniball_crystal_HE<<"]["
              //    <<miniball_segment_HE<<"] = "
              //    <<miniball_fi_ave[jjj][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE]<<endl;
              //cout<<miniball_fi[jjj][miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE]<<endl;
              if(jjj==0) miniball_fi_interactions[miniball_detector_HE][miniball_crystal_HE][miniball_segment_HE]++;
            }
          }
        }
      }
    }
    //____________________________________________________________________________________________
    // Miniball Finished. Analyzing Hector
    dummy = -1;
    if(hectorInclude==1)  {
      for(int nnn=0;nnn<NUMBEROFHECTORDETECTORS;nnn++)  {
        if(hector_flag[nnn]==1.0) {
          hector_crystal_mult++;
          hector_fired[nnn] = true;
          //summed energy
          hector_energy_sum = hector_energy_sum + hector_energy_not_cor[nnn];

          if(hector_energy_not_cor[nnn] > hector_energy_max)  {
            hector_energy_max = hector_energy_not_cor[nnn];
            hector_x_thatsit  = hector_pos[0][nnn];
            hector_y_thatsit  = hector_pos[1][nnn];
            hector_z_thatsit  = hector_pos[2][nnn];
          }
        }
      }
      h_hector_crystal_mult->Fill(hector_crystal_mult);
      doppler_energy = GetDopplerCorrectedEnergy(hector_x_thatsit,hector_y_thatsit,hector_z_thatsit,
                                              ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2],
                                              hector_energy_sum,beta_rec,beta_mean_tof,beta_ave,decay_z);
      
      if(hector_crystal_mult==1)  {
        h_hector        ->Fill(hector_energy_sum);
        h_hector_doppler->Fill(doppler_energy);
        
        if(hector_fired[5] == true || hector_fired[6] == true)  {
          h_hector_ring[0]        ->Fill(hector_energy_sum);
          h_hector_doppler_ring[0]->Fill(doppler_energy);
        }
        else  {
          h_hector_ring[1]        ->Fill(hector_energy_sum);
          h_hector_doppler_ring[1]->Fill(doppler_energy);
        }
      }
    }
    //_____________________________________________________________________________________________
    // The variables have to be reset.
    cluster_crystal_mult  = 0;
    cluster_detector_mult = 0;
    cluster_energy_sum    = 0.0;
    cluster_energy_max    = -999.9;
    for(int nnn=0;nnn<NUMBEROFCLUSTERDETECTORS;nnn++)  {
      cluster_fired[nnn]  = false;
    }
   
    miniball_segment_mult  = 0;
    miniball_crystal_mult  = 0;
    miniball_detector_mult = 0;
    miniball_energy_sum    = 0.0;
    miniball_energy_max    = -999.9;
    for(int nnn=0;nnn<NUMBEROFMINIBALLDETECTORS;nnn++)  {
      miniball_fired[nnn] = false;
    }

    hector_crystal_mult  = 0;
    hector_energy_sum    = 0.0;
    hector_energy_max    = -999.9;
    for(int nnn=0;nnn<NUMBEROFHECTORDETECTORS;nnn++)  {
      hector_fired[nnn] = false;
    }

    observedDecayNumber++;
    //_____________________________________________________________________________________________
    // Stopped beam stuff. In particular, the gamma gamma matrices.
    if(evnum != eventNumberPrevious)  { 
      h_number_of_observed_decays->Fill(observedDecayNumber);
      //Filling the gamma-gamma spectra:
      //Sorting in the right order
      for(int i =1;i<observedDecayNumber;i++)  {
        if(gammaEnergy[0]>=gammaEnergy[i])  {
          h_cluster_gamma_gamma[i]->Fill(gammaEnergy[0],gammaEnergy[i]);
          h_cluster_gamma_gamma[0]->Fill(gammaEnergy[0],gammaEnergy[i]);
        }
        else  {
          h_cluster_gamma_gamma[i]->Fill(gammaEnergy[i],gammaEnergy[0]);
          h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[0]);
        }
        if(i>1)  {
          if(gammaEnergy[1]>=gammaEnergy[i])
            h_cluster_gamma_gamma[0]->Fill(gammaEnergy[1],gammaEnergy[i]);
          else h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[1]);
        }
        if(i>2)  { 
          if(gammaEnergy[2]>=gammaEnergy[i])
            h_cluster_gamma_gamma[0]->Fill(gammaEnergy[2],gammaEnergy[i]);
          else h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[2]);
        }
        if(i>3) {
          if(gammaEnergy[3]>=gammaEnergy[i])
            h_cluster_gamma_gamma[0]->Fill(gammaEnergy[3],gammaEnergy[i]);
          else h_cluster_gamma_gamma[0]->Fill(gammaEnergy[i],gammaEnergy[3]);
h_miniball_doppler_single_ave_z_ave_fi            ->Write();        }
      } 
      observedDecayNumber = 0;
      for(int i =0;i<20;i++)
        {
          gammaEnergy[i]=0.;
        }
      gammaEnergy[observedDecayNumber] = doppler_energy;
    }
    eventNumberPrevious = evnum;
  }
  //_______________________________________________________________
  // Getting the average point:
  if(cluster_fi_option >= 1) GetAverageInteractionPointCluster();
  if(miniball_fi_option >= 1) GetAverageInteractionPointMiniball();

  //--------------------------------------------------------------
  //Writing into the root file:-----------------------------------
  //--------------------------------------------------------------
  rootfile->cd();
  h_decayId_vs_energy_sum                               ->Write();
  h_number_of_observed_decays                           ->Write();
  // Cluster------------------------------------------------------
  if(clusterInclude==1)  {
    h_cluster_crystal_mult                              ->Write();
    h_cluster_detector_mult                             ->Write();
    h_cluster_single                                    ->Write();
    h_cluster_energy_sum                                ->Write();
    for(int i=0;i<3;i++)h_cluster_ring[i]               ->Write();
    for(int i=0;i<10;i++)h_cluster_gamma_gamma[i]       ->Write();
    h_cluster_doppler_single                            ->Write();
    h_cluster_doppler_single_at                         ->Write();
    h_cluster_doppler_single_ave_z                      ->Write();
    h_cluster_doppler_single_fep                        ->Write();
    h_cluster_doppler_single_at_fep                     ->Write();
    h_cluster_doppler_single_ave_z_fep                  ->Write();
    if(cluster_fi_option==2)  {  
      h_cluster_doppler_single_ave_z_ave_fi             ->Write();
      h_cluster_doppler_single_fi                       ->Write();
      h_cluster_doppler_single_fi_at                    ->Write();
      h_cluster_doppler_single_fi_ave_z                 ->Write();
      h_cluster_doppler_single_ave_z_ave_fi_fep         ->Write();
    }
    h_cluster_doppler_total                             ->Write();
    for(int i=0;i<3;i++)h_cluster_doppler_total_ring[i] ->Write();
  }
  // Miniball-----------------------------------------------------
  if(miniballInclude==1)  {
    h_miniball_segment_mult                             ->Write();
    h_miniball_crystal_mult                             ->Write();
    h_miniball_detector_mult                            ->Write();
    h_miniball_doppler_single                           ->Write();
    h_miniball_doppler_single_at                        ->Write(); 
    h_miniball_doppler_single_ave_z                     ->Write();
    h_miniball_doppler_single_fep                       ->Write();
    h_miniball_doppler_single_at_fep                    ->Write();
    h_miniball_doppler_single_ave_z_fep                 ->Write();
    if(miniball_fi_option==2)  {  
      h_miniball_doppler_single_ave_z_ave_fi            ->Write();
      h_miniball_doppler_single_fi                      ->Write();
      h_miniball_doppler_single_fi_ave_z                ->Write();
      h_miniball_doppler_single_ave_z_ave_fi_fep        ->Write();
      h_miniball_doppler_single_fi_at                   ->Write();
    }
    h_miniball                                          ->Write();
    for(int i=0;i<2;i++)h_miniball_ring[i]              ->Write();
    h_miniball_doppler_total                            ->Write();
    for(int i=0;i<2;i++)h_miniball_doppler_total_ring[i]->Write();
  }
  // Hector-------------------------------------------------------
  if(hectorInclude==1)  {
    h_hector_crystal_mult                               ->Write();
    h_hector_doppler                                    ->Write();
    for(int i=0;i<2;i++)h_hector_doppler_ring[i]        ->Write();
    h_hector                                            ->Write();
    for(int i=0;i<2;i++)h_hector_ring[i]                ->Write();
  }
  return 0;
}

//_______________________________________________________________
void GetAverageInteractionPointMiniball(){
  cout<<" Getting the average FI point for all the Miniball crystals..."<<endl;

  FILE *fOutAve  = fopen("./output/MiniballAverageInteractionPoint.out","w");
  FILE *fOutDiff = fopen("./output/MiniballInteractionDifferenceToCenterOfGravity.out","w");
  float difference[3];
  float distance1, distance2;
  float thetaDifference;
  for( int i = 0;i<NUMBEROFMINIBALLDETECTORS;i++)  {
    for( int j = 0;j<3;j++)  {
      for( int k = 0;k<7;k++)  {
        for(int jjj=0;jjj<3;jjj++)  {
          if(k!=0)miniball_fi_ave[jjj][i][j][k] = miniball_fi_ave[jjj][i][j][k]/miniball_fi_interactions[i][j][k];
          difference[jjj] = miniball_fi_ave[jjj][i][j][k] - miniball_pos[jjj][i][j][k];
        }
        // The average interaction point:
        fprintf(fOutAve,"%4i %4i %10.2f %10.2f %10.2f \n", i, miniball_fi_interactions[i][j][k], 
                miniball_fi_ave[0][i][j][k], miniball_fi_ave[1][i][j][k], miniball_fi_ave[2][i][j][k]); 
        //cout<<i<<" "<<fi_interactions[i]<<" "<<fi_ave[0][i]<<" "<<fi_ave[1][i]<<" "<<fi_ave[2][i];
    
        distance1 = sqrt(miniball_pos[0][i][j][k]*miniball_pos[0][i][j][k] + 
                         miniball_pos[1][i][j][k]*miniball_pos[1][i][j][k] + 
                         miniball_pos[2][i][j][k]*miniball_pos[2][i][j][k]);
        distance2 = sqrt(miniball_fi_ave[0][i][j][k]*miniball_fi_ave[0][i][j][k] + 
                         miniball_fi_ave[1][i][j][k]*miniball_fi_ave[1][i][j][k] + 
                         miniball_fi_ave[2][i][j][k]*miniball_fi_ave[2][i][j][k]);

        //   thetaDifference = acos(dali2_pos[0][i]*fi_ave[0][i] + dali2_pos[1][i]*fi_ave[1][i] + dali2_pos[2][i]*fi_ave[2][i]
        //                      /distance1/distance2)

        // Comparison to the center of gravity points:
        fprintf(fOutDiff,"%4i %10.2f %10.2f %10.2f %10.2f\n", i, difference[0],difference[1],difference[2],thetaDifference); 
        cout<<i<<" "<<difference[0]<<" "<<difference[1]<<" "<<difference[2]<<endl;
      }
    }
  }
  fclose(fOutAve); 
  fclose(fOutDiff);
}

//__________________________________________________________
void GetAverageInteractionPointCluster(){
  cout<<" Getting the average FI point for all the Cluster crystals..."<<endl;

  FILE *fOutAve  = fopen("./output/ClusterAverageInteractionPoint.out","w");
  FILE *fOutDiff = fopen("./output/ClusterInteractionDifferenceToCenterOfGravity.out","w");
  float difference[3];
  float distance1, distance2;
  float thetaDifference;
  for( int i = 0;i<NUMBEROFCLUSTERDETECTORS;i++)  {
    for( int j = 0;j<7;j++)  {
      for(int jjj=0;jjj<3;jjj++)  {
        cluster_fi_ave[jjj][i][j] = cluster_fi_ave[jjj][i][j]/cluster_fi_interactions[i][j];
        difference[jjj] = cluster_fi_ave[jjj][i][j] - cluster_pos[jjj][i][j];
      }
      // The average interaction point:
      fprintf(fOutAve,"%4i %4i %10.2f %10.2f %10.2f \n", i, cluster_fi_interactions[i][j], 
              cluster_fi_ave[0][i][j], cluster_fi_ave[1][i][j], cluster_fi_ave[2][i][j]); 
            //cout<<i<<" "<<fi_interactions[i]<<" "<<fi_ave[0][i]<<" "<<fi_ave[1][i]<<" "<<fi_ave[2][i];
    
            distance1 = sqrt(cluster_pos[0][i][j]*cluster_pos[0][i][j] + 
                             cluster_pos[1][i][j]*cluster_pos[1][i][j] + 
                             cluster_pos[2][i][j]*cluster_pos[2][i][j]);
            distance2 = sqrt(cluster_fi_ave[0][i][j]*cluster_fi_ave[0][i][j] + 
                             cluster_fi_ave[1][i][j]*cluster_fi_ave[1][i][j] + 
                             cluster_fi_ave[2][i][j]*cluster_fi_ave[2][i][j]);
            
            //   thetaDifference = acos(dali2_pos[0][i]*fi_ave[0][i] + dali2_pos[1][i]*fi_ave[1][i] + dali2_pos[2][i]*fi_ave[2][i]
            //                      /distance1/distance2);
            
            // Comparison to the center of gravity points:
            fprintf(fOutDiff,"%4i %10.2f %10.2f %10.2f %10.2f\n", i, difference[0],difference[1],difference[2],thetaDifference); 
            cout<<i<<" "<<difference[0]<<" "<<difference[1]<<" "<<difference[2]<<endl;
    }
  }
  fclose(fOutAve); 
  fclose(fOutDiff);
}












