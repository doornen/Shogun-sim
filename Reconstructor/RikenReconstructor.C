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

#define MAXFOLD 20
class GammaDetector;
//________________________________________________
class GammaDetector {
 public:
  GammaDetector(){
  };
  ~GammaDetector(){
  };
  
  float event;              // The number of the event
  int id[2][MAXFOLD];      // Detector number, first index: 0=dali2, 1=labr3
  int fold[2];             // Gamma ray fold given separately
  
  float e[2][MAXFOLD];     // The detected energy in the lab system
  float d[2][MAXFOLD];     // Doppler corrected energy
  float reste[2][MAXFOLD]; // The energy of the transition at rest
  
  float x[2][MAXFOLD];     // Position of the detector
  float y[2][MAXFOLD];
  float z[2][MAXFOLD];

  ClassDef(GammaDetector,1)
};
ClassImp(GammaDetector)

GammaDetector *gammaDet;
void ResetGammaDetectorValues();


//_________________________________________________________________________________________________
float GetDopplerCorrectedEnergy(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                                float target_x,float target_y,float target_z,
                                float pos_det_x, float pos_det_y, float pos_det_z, 
                                float gamma_energy, float beta_rec,float beta_mean_tof,float beta_average,float decay_z){ 
  double vx1 = gamma_det_x - target_x;                         //cout<<"vx1: "<<vx1<<endl;
  double vy1 = gamma_det_y - target_y;                         //cout<<"vy1: "<<vy1<<endl;
  double vz1 = gamma_det_z - target_z - decay_z;               //cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  double vx2 = pos_det_x - target_x;                           //cout<<"vx2: "<<vx2<<endl;
  double vy2 = pos_det_y - target_y;                           //cout<<"vy2: "<<vy2<<endl;
  double vz2 = pos_det_z - target_z;                           //cout<<"vy2: "<<vz2<<endl;
 
  double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
          
  // The theta angle:
  double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);//cout<<"Theta lab: "<<theta_lab<<endl;
  // the beta value:
  double beta = beta_average + (beta_rec-beta_mean_tof);         //cout<<"beta: "<<beta<<endl;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

int RikenReconstructor()  {
  // int main(int argc,char** argv)  { 
  
  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
  
  float beta_average     = 0.0;        // The average value at the moment of decay for the correction
  float beta_mean_tof    = 0.0;        // The mean beta from Tof
  float decay_z          = 0.0;        // The lifetime moves the average decay point along the beam axis
  int dali2_fi_option    = 0;
 
  float dali2_energy[NUMBEROFDALI2CRYSTALS]                  = {0.0};  // Unsorted energy, filled according to the loop from 0 to the last crystal
  float dali2_e_rest[NUMBEROFDALI2CRYSTALS]                  = {0.0};  // Energy in the rest system of detected energy in crystal
  float dali2_time[NUMBEROFDALI2CRYSTALS]                    = {0.0};  // same for the time
  float dali2_dopplerEnergy[NUMBEROFDALI2CRYSTALS]           = {0.0};  // Unsorted dopplerEnergy
  int  dali2_crystalFired[NUMBEROFDALI2CRYSTALS]             = {0};    // The crystal that fired the unsorted energy
  
  float labr3_energy[NUMBEROFLABR3ARRAYCRYSTALS]             = {0.0};  // Unsorted energy, filled according to the loop from 0 to the last crystal
  float labr3_e_rest[NUMBEROFLABR3ARRAYCRYSTALS]              = {0.0};  // Energy in the rest system of detected energy in crystal
  float labr3_time[NUMBEROFLABR3ARRAYCRYSTALS]                = {0.0};  // same for the time
  float labr3_dopplerEnergy[NUMBEROFLABR3ARRAYCRYSTALS]       = {0.0};  // Unsorted dopplerEnergy
  int  labr3_crystalFired[NUMBEROFLABR3ARRAYCRYSTALS]         = {0};    // The crystal that fired the unsorted energy

  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/RikenReconstructor.in","r");
  while(!feof(fin))  {
    fscanf(fin,"%s ",temp);  
    if(strcmp(temp,"INPUTFILE")==0)  {
      fscanf(fin,"%s",&root_in);
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
      fscanf(fin,"%f",&beta_average); 
      printf("%s %f \n",temp,beta_average);
    }
    else if(strcmp(temp,"BETATOFAVERAGE")==0)  {
      fscanf(fin,"%f",&beta_mean_tof); 
      printf("%s %f \n",temp,beta_mean_tof);
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
    else if(strcmp(temp,"DALI2FIFIND")==0)  {
      fscanf(fin,"%i",&dali2_fi_option); 
      printf("%s %i \n",temp,dali2_fi_option);
    }
    else if(strcmp(temp,"DALI2INCLUDE")==0) {
      fscanf(fin,"%i ",&dali2_opt); 
      printf("%s %i \n",temp,dali2_opt);
    }
    else if(strcmp(temp,"SHOGUNINCLUDE")==0) {
      fscanf(fin,"%i ",&Shogun_opt); 
      printf("%s %i \n",temp,Shogun_opt);
    }
    else if(strcmp(temp,"GRAPEINCLUDE")==0) {
      fscanf(fin,"%i ",&grape_opt); 
      printf("%s %i \n",temp,grape_opt);
    }
    else if(strcmp(temp,"SGTINCLUDE")==0) {
      fscanf(fin,"%i ",&sgt_opt); 
      printf("%s %i \n",temp,sgt_opt);
    }
    else if(strcmp(temp,"SPHEREINCLUDE")==0)  {
      fscanf(fin,"%i ",&sphere_opt); 
      printf("%s %i \n",temp,sphere_opt);
    }
    else if(strcmp(temp,"LABR3INCLUDE")==0)  {
      fscanf(fin,"%i ",&LaBr3Opt); 
      printf("%s %i \n",temp,LaBr3Opt);
    }
    else if(strcmp(temp,"END")==0) break;
    else {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      return 0;
    }
  }
  //--------------------------------------------------------------------------------------------------
  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();

  TH1F *h_dali2_crystal_mult;
  TH1F *h_dali2_doppler;
  TH1F *h_dali2_doppler_total; // Taking the sum energy
  TH1F *h_dali2_energy;

  TH2F *h_dali2_crystal_fired_doppler;
  TH2F *h_dali2_crystal_fired_energy;
  

  if(dali2_opt==1)  {
    //Creating dali2 spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Crystal multiplicity,ok
    h_dali2_crystal_mult = new TH1F("h_dali2_crystal_Mult","h_dali2_crystal_mult",NUMBEROFDALI2CRYSTALS,0,NUMBEROFDALI2CRYSTALS);  
    // All, OK
    h_dali2_doppler = new TH1F("h_dali2_doppler","h_dali2_doppler",numBin,firstBin,lastBin);
    
    // All, taking the sum energy....
    h_dali2_doppler_total = new TH1F("h_dali2_doppler_total","h_dali2_doppler_total",numBin,firstBin,lastBin);
    
    //Not doppler corrected, OK
    h_dali2_energy = new TH1F("h_dali2_energy","h_dali2_energy",numBin,firstBin,lastBin);

    //crystal fired vs. doppler
    h_dali2_crystal_fired_doppler = new TH2F("h_dali2_crystal_fired_doppler","",NUMBEROFDALI2CRYSTALS,0,NUMBEROFDALI2CRYSTALS,
                                             numBin,firstBin,lastBin);

    //crystal fired vs. energy
    h_dali2_crystal_fired_energy = new TH2F("h_dali2_crystal_fired_energy","",NUMBEROFDALI2CRYSTALS,0,NUMBEROFDALI2CRYSTALS,
                                             numBin,firstBin,lastBin);    
  }

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // SHOGUN:
  TH1F *h_shogun_crystal_mult;
  TH1F *h_shogun_doppler_single;
  TH1F *h_shogun_doppler;
  TH1F *h_shogun_doppler_total;
  TH1F *h_shogun_energy_total;
  if(Shogun_opt==1)  {
    //Creating shogun spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Crystal multiplicity,ok
    h_shogun_crystal_mult = new TH1F("h_shogun_crystal_Mult","h_shogun_crystal_mult",100,0,100);
    
    h_shogun_doppler = new TH1F("h_shogun_doppler","h_shogun_doppler",numBin,firstBin,lastBin);
        //Single hits, OK
    h_shogun_doppler_single = new TH1F("h_shogun_doppler_single","h_shogun_doppler_single",numBin,firstBin,lastBin);
    // All, OK, This includes addback events.
    h_shogun_energy_total = new TH1F("h_shogun_energy_total","h_shogun_energy_total",numBin,firstBin,lastBin);   
    h_shogun_doppler_total = new TH1F("h_shogun_doppler_total","h_shogun_doppler_total",numBin,firstBin,lastBin);
  }
  //---------------------------------------------------------------
  ///////////////////////////////
  //Grape!!!!!!!!!!!!!!!!!!!!!
  //////////////////////////////
  TH1F *h_grape_segment_mult;
  TH1F *h_grape_crystal_mult;
  TH1F *h_grape_detector_mult;
  TH1F *h_grape_doppler_cor;
  TH1F *h_grape_doppler_cor_single_segment;
  TH1F *h_grape_doppler_cor_single_crystal;
  TH1F *h_grape_doppler_cor_single_detector;
  TH1F *h_grape_doppler_cor_two_segments;
  TH1F *h_grape_doppler_cor_two_crystals;
  TH1F *h_grape_doppler_cor_det[NUMBEROFGRAPEDETECTORS];
  TH1F *h_grape_energy;
  TH1F *h_grape_energy_det[NUMBEROFGRAPEDETECTORS];

  if(grape_opt==1)  {
    h_grape_segment_mult = new TH1F("h_grape_segment_mult",
                                    "h_grape_segment_mult",15,0,15);
    h_grape_crystal_mult = new TH1F("h_grape_crystal_mult",
                                    "h_grape_crystal_mult",15,0,15);
    h_grape_detector_mult = new TH1F("h_grape_detector_mult",
                                     "h_grape_detector_mult",15,0,15);
    // Hit in only one segment, OK
    h_grape_doppler_cor = new TH1F("h_grape_doppler_cor",

                                   "h_grape_doppler_cor",numBin,firstBin,lastBin); 
    // Hit in only one segment, OK
    h_grape_doppler_cor_single_segment = new TH1F("h_grape_doppler_cor_single_segment",
                                                  "h_grape_doppler_cor_single_segment",numBin,firstBin,lastBin); 
    // Hit in only one crystal, OK
    h_grape_doppler_cor_single_crystal = new TH1F("h_grape_doppler_cor_single_crystal",
                                                  "h_grape_doppler_cor_single_crystal",numBin,firstBin,lastBin);
    // Hit in only one detector
    h_grape_doppler_cor_single_detector = new TH1F("h_grape_doppler_cor_single_detector",
                                                   "h_grape_doppler_cor_single_detector",numBin,firstBin,lastBin); 
    // Hit in two segments in the same detector
    h_grape_doppler_cor_two_segments = new TH1F("h_grape_doppler_cor_two_segments",
                                                "h_grape_doppler_cor_two_segments",numBin,firstBin,lastBin);
    //Energy added up from two crystals in the same detector
    h_grape_doppler_cor_two_crystals = new TH1F("h_grape_doppler_cor_two_crystals",
                                                "h_grape_doppler_cor_two_crystals",numBin,firstBin,lastBin);
    //All the detectors individually
    //
    for(int i=0;i<NUMBEROFGRAPEDETECTORS;i++)  {
      sprintf(temp,"h_grape_doppler_cor_det[%i]",i);
      h_grape_doppler_cor_det[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
    }
    //Not Doppler corrected, OK
    h_grape_energy = new TH1F("h_grape_energy","h_grape_energy",numBin,firstBin,lastBin);
 
    for(int i=0;i<NUMBEROFGRAPEDETECTORS;i++)  {
      sprintf(temp,"h_grape_energy_det[%i]",i);
      h_grape_energy_det[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
    }
  }
  //---------------------------------------------------------------	      
  ///////////////////////////////
  //SGT!!!!!!!!!!!!!!!!!!!!!!!!!!
  ///////////////////////////////
  // multiplicity of the coaxial crystals within the array
  TH1F *h_sgt_coax_mult;
  TH1F *h_sgt_strip_mult[NUMBEROFSGTDETECTORS];
  TH1F *h_sgt_strip_fan[NUMBEROFSGTDETECTORS];
  TH1F *h_sgt_doppler_cor_det[NUMBEROFSGTDETECTORS];
  TH1F *h_sgt_doppler_cor;
  TH1F *h_sgt_coax_energy_det[NUMBEROFSGTDETECTORS];

  if(sgt_opt==1)  {
    cout<<"Creating the SGT spectra"<<endl;
    h_sgt_coax_mult = new TH1F("h_sgt_coax_mult","h_sgt_coax_mult",NUMBEROFSGTDETECTORS,0,NUMBEROFSGTDETECTORS);
    
    // multiplicity and fan spectra of strips of the detectors
    for(int i=0;i<NUMBEROFSGTDETECTORS;i++)  {
      sprintf(temp,"h_sgt_strip_mult[%i]",i);
      h_sgt_strip_mult[i] = new TH1F(temp,temp,25,0,25);

      sprintf(temp,"h_sgt_strip_fan[%i]",i);
      h_sgt_strip_fan[i] = new TH1F(temp,temp,26,0,26);
    }

    // Doppler corrected, detectorwise
    for(int i=0;i<NUMBEROFSGTDETECTORS;i++)  {
      sprintf(temp,"h_sgt_doppler_cor_det[%i]",i);
      h_sgt_doppler_cor_det[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
    }
    //cout<<"Doppler corrected, total array"<<endl;
    sprintf(temp,"h_sgt_doppler_cor");
    h_sgt_doppler_cor = new TH1F(temp,temp,numBin,firstBin,lastBin);
    
    //Not Doppler-corrected
    for(int i=0;i<NUMBEROFSGTDETECTORS;i++)  {
      sprintf(temp,"h_sgt_coax_energy_det[%i]",i);
      h_sgt_coax_energy_det[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
    }
  }
  //-----------------------------------------------------------------------------
  // The Sphere:
  TH1F *h_sphere_doppler_cor;
  TH1F *h_sphere_doppler_cor_angle_cut[18];
  if(sphere_opt==1)  {
    sprintf(temp,"h_sphere_doppler_cor");
    h_sphere_doppler_cor = new TH1F(temp,temp,numBin,firstBin,lastBin);
    for(int i=0;i<18;i++)  {
      sprintf(temp,"h_sphere_doppler_cor_angle_cut[%i]",i);
      h_sphere_doppler_cor_angle_cut[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
    }
  } 

  //-----------------------------------------------------------------------------
  // The LaBr3 array:
  TH1F *h_labr3_doppler_cor;
  if(LaBr3Opt==1)  {
    sprintf(temp,"h_labr3_doppler_cor");
    h_labr3_doppler_cor = new TH1F(temp,temp,numBin,firstBin,lastBin);
  } 

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  infile = new TFile(root_in,"READ");
  t = (TTree*)infile->Get("ObservedEvents");
  TH1F *h_evnum;
  h_evnum = new TH1F("h_evnum","h_evnum",1000,0,t->GetEntries());

  t->SetBranchAddress("EventNumber",&evnum);
  t->SetBranchAddress("DecayIDOfEventNumber",&decayIDevnum);
  t->SetBranchAddress("ProjectileVertex",p0);                                   // Position at the fragmentation point
  t->SetBranchAddress("VProjectileBeforeTarget",pvb);                           // Normalized Vector of beam before the target
  t->SetBranchAddress("VProjectileAfterTarget",pva);                            // Normalized Vector of beam after the target
  t->SetBranchAddress("EnergyProjectile",&energy_p);                            // energy of beam before the target in MeV/u
  t->SetBranchAddress("BetaBeforeTarget",&beta_b);                              // Beta before the target
  t->SetBranchAddress("BetaReal",&beta_r);                                      // Beta during deexcitation	
  t->SetBranchAddress("BetaAfterTarget",&beta_a);                               // Beta After Target
  t->SetBranchAddress("Halflife",&halflife);                                    // Halflife
  t->SetBranchAddress("DecayTimeAfterInteraction",&decay_time_after_interaction);
  t->SetBranchAddress("VertexGamma",g0);                                        // Position at the gamma emmittance point
  t->SetBranchAddress("VGamma",gv);		                                // Gamma vector
  t->SetBranchAddress("EGammaRest",&e_rest);		                        // Energy at rest
  t->SetBranchAddress("EGammaDoppler",&e_doppler);                              // Theta of doppler boosted gamma
  t->SetBranchAddress("ThetaGammaRest",&theta_gamma_rest);
  t->SetBranchAddress("ThetaGammaLab",&theta_gamma_lab);
  
  //Energy of fragment at fragmentation
  t->SetBranchAddress("EnergyVertex",&energy_vertex_stored); 
  t->SetBranchAddress("BGGamma",&bgGamma);                                      // To know if the gamma is from bg or a good event
  //So far identical to the ROOT-Tree of step one
  //The new stuff:
  //Setting the type of the gamma detector
  t->SetBranchAddress("GammaDetType",&gamma_det_type);
  //-------------------------------
  //DALI2:
  if(dali2_opt==1)  {
    t->SetBranchAddress("DALI2Flag",dali2_flag);                                // ID of which det recorded a gamma in the event
    t->SetBranchAddress("DALI2EnergyNotCor",dali2_energy_not_cor);
    t->SetBranchAddress("DALI2Time",dali2_time);
    if(dali2_fi_option >= 1)  {
      t->SetBranchAddress("DALI2FI",dali2_fi);
    }
  }
  //-------------------------------
  //Shogun:
  if(Shogun_opt==1)  {
    t->SetBranchAddress("ShogunFlag",Shogun_flag);                              //ID of which det recorded a gamma in the event
    t->SetBranchAddress("ShogunEnergyNotCor",Shogun_energy_not_cor);
    t->SetBranchAddress("ShogunTime",Shogun_time);
  }
  //-----------------------------------------------------------
  //GRAPE:
  if(grape_opt==1)  {
    t->SetBranchAddress("GRAPEFlag",grape_flag);
    t->SetBranchAddress("GRAPECrystalFlag",grape_crystal_flag);
    t->SetBranchAddress("GRAPEEnergyNotCor",grape_energy_not_cor);
    t->SetBranchAddress("GRAPETime",grape_time);
  }
  //-----------------------------------------------------------
  //SGT
  if(sgt_opt==1)  {
    t->SetBranchAddress("SGT_flag",sgt_flag);
    t->SetBranchAddress("SGT_coax_crystal_flag",sgt_coax_crystal_flag);
    t->SetBranchAddress("SGT_planar_strip_flag",sgt_planar_strip_flag); //central plus 25 strips
    t->SetBranchAddress("SGT_energy_not_cor",sgt_energy_not_cor);
    t->SetBranchAddress("SGT_time",sgt_time);
  }
  if(sphere_opt==1)  {
    //fi: First interaction
    t->SetBranchAddress("sphere_fi",sphere_fi); // reconstructed position from a det after the secondary target
    t->SetBranchAddress("sphereEnergyNotCor",&sphere_energy_not_cor);
  }
  //-------------------------------
  //Labr3:
  if(LaBr3Opt==1)  {
    t->SetBranchAddress("LaBr3Flag",LaBr3Flag);
    t->SetBranchAddress("LaBr3EnergyNotCor",LaBr3EnergyNotCor);
    t->SetBranchAddress("LaBr3Time",LaBr3Time);
  }

  //-------------------------------
  //Position detector after targetpos_det_after_target_res
  t->SetBranchAddress("PosDet",pos_det);                                        // Reconstructed position after the secondary target
  t->SetBranchAddress("VertexReconstructed",ver_rec);                           // Reconstructed vertex postion, 
 
  //reconstructed beta 
  t->SetBranchAddress("BetaReconstructed",&beta_rec);                           // Reconstructed beta, including resolution

  //------------------------------------------------------------------------------------------------------
  //Going to the Header file, which has the information of constant values.
  infile->cd();
  tHeader = (TTree*)infile->Get("Header");
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
  //DALI2:
  if(dali2_opt==1)  {
    tHeader->SetBranchAddress("DALI2EnResOpt",&dali2_en_res_opt);
    tHeader->SetBranchAddress("DALI2EnRes",dali2_en_res);
    tHeader->SetBranchAddress("DALI2TimeRes",dali2_time_res);
    tHeader->SetBranchAddress("DALI2Pos",dali2_pos);
    tHeader->SetBranchAddress("BetaResolution",&beta_res);  
    tHeader->SetBranchAddress("PosDetAtTargetRes",&pos_det_at_target_res); 
    tHeader->SetBranchAddress("PosDetAfterTargetRes",&pos_det_after_target_res);
  }
  //Shogun
  if(Shogun_opt==1)  {
    tHeader->SetBranchAddress("ShogunEnResOpt",&Shogun_en_res_opt);
    tHeader->SetBranchAddress("ShogunEnRes",Shogun_en_res);
    tHeader->SetBranchAddress("ShogunTimeRes",Shogun_time_res);
    tHeader->SetBranchAddress("ShogunPos",Shogun_pos);
  }
  //GRAPE:
  if(grape_opt==1)  {
    tHeader->SetBranchAddress("GRAPEEnResOpt",&grape_en_res_opt);
    tHeader->SetBranchAddress("GRAPEEnRes",grape_en_res);
    tHeader->SetBranchAddress("GRAPETimeRes",grape_time_res);
    tHeader->SetBranchAddress("GRAPEPos",grape_pos);
  }
  //SGT
  if(sgt_opt==1)  {
    tHeader->SetBranchAddress("SGT_en_res_opt",&sgt_en_res_opt);
    tHeader->SetBranchAddress("SGT_en_res",sgt_en_res);
    tHeader->SetBranchAddress("SGT_time_res",sgt_time_res);
    tHeader->SetBranchAddress("SGTPOS",sgt_pos);   //first the coax, then the fictive planar central, then the strips
  }
  //SGT
  if(LaBr3Opt==1)  {
    tHeader->SetBranchAddress("LaBr3EnResOpt",&LaBr3EnResOpt);
    tHeader->SetBranchAddress("LaBr3EnRes",LaBr3EnRes);
    tHeader->SetBranchAddress("LaBr3TimeRes",LaBr3TimeRes);
    tHeader->SetBranchAddress("LaBr3Pos",LaBr3Pos);
  }

  tHeader->GetEntry(0);
  //------------------------------------------------------
  //Finished reading the Root file

  nentries = (Int_t)(t->GetEntries()/reduction_factor);
  cout <<"Entries: "<< nentries << endl;

  //______________________________________________
  TFile  *outfile2 = new TFile("OutputTree.root","RECREATE");
  outfile2->cd();
  
  gammaDet = new GammaDetector();
  
  t2 = new TTree("ReconstructedEvents","ReconstructedEvents");
  t2->Branch("gammaDet",&gammaDet);
  //______________________________________________

  //Dividing old dali2 into 16 rings of equal z-position:
  //int dali2_ring_low[16]={0,6,12,20,32,45,56,70,80,90,104,116,128,140,148,154};
  //int dali2_ring_high[16]={5,11,19,31,44,55,69,79,89,103,115,127,139,147,153,159};
  //Dividing new dali2 into 17 rings of equal theta angle:
  int dali2_ring_low[17] ={0, 6,13,20,30,40,52,66,80,94,108,122,136,148,156,164,170};
  int dali2_ring_high[17]=   {5,12,19,29,39,51,65,79,93,107,121,135,147,155,163,169,181};

  int dali2_ring_energy_max = -1;
  float dali2DopplerCorrectedEnergy[NUMBEROFDALI2CRYSTALS] = {0.0};
 
  int dali2_crystal_mult  = 0;
  int Shogun_crystal_mult = 0;
  int grape_segment_mult  = 0;
  int grape_crystal_mult  = 0;
  int grape_detector_mult = 0;
  int sgt_strip_mult      = 0;
  int sgt_coax_mult       = 0;
  int sgt_detector_mult   = 0;
  int labr3_crystal_mult  = 0;

  int sgt_strip_mult_det[10] = {0.};

  float dali2_energy_sum  = 0.0;
  float Shogun_energy_sum = 0.0;
  float grape_energy_sum  = 0.0;
  float sgt_energy_sum    = 0.0;
  float labr3_energy_sum  = 0.0;
  float grape_energy_sum_det[NUMBEROFGRAPEDETECTORS] = {0.};
  float sgt_energy_sum_det[NUMBEROFSGTDETECTORS]     = {0.};
  
  float dali2_energy_max  = -999.9;
  float Shogun_energy_max = -999.9;
  float grape_energy_max  = -999.9;
  float sgt_energy_max    = -999.9;
  float grape_energy_max_det[NUMBEROFGRAPEDETECTORS] = {0.0}; 
  float sgt_energy_max_det[NUMBEROFSGTDETECTORS]     = {0.};
 
  float dali2_pos_thatsit[3]  = {0.0};
  float Shogun_pos_thatsit[3] = {0.0};
  float grape_pos_thatsit[3]  = {0.0};
  float sgt_pos_thatsit[3]    = {0.0};
  float grape_pos_thatsit_det[3][NUMBEROFGRAPEDETECTORS] = {{0.}}; 
  float sgt_pos_thatsit_det[3][NUMBEROFSGTDETECTORS]     = {{0.}};
  int eventNumberNext = 0; 


  cout<<"Starting the analysis"<<endl;
  for(int iii=0;iii<nentries;iii++) {        
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    
    if(iii<nentries-1)  {
      t->GetEntry(iii+1);
      
      eventNumberNext = evnum;
    }

    t->GetEntry(iii);
    

    // Reset all the values for the root tree filling:
    ResetGammaDetectorValues();    
    //-------------------------------------------------------------
    // Starting with the gamma ray detector spectra:
    if(dali2_opt==1)  { 
      for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++) {
        if(dali2_flag[nnn]==1)   {
          dali2_energy[nnn] += dali2_energy_not_cor[nnn];
          dali2_e_rest[nnn] += e_rest;
          dali2_energy_sum += dali2_energy_not_cor[nnn];
          dali2_crystalFired[dali2_crystal_mult]=nnn;
          //Filling the energy spectra without Doppler-correction        
          h_dali2_energy->Fill(dali2_energy_not_cor[nnn]);
          h_dali2_crystal_fired_energy->Fill(nnn,dali2_energy_not_cor[nnn]);
          
          //Performing the Doppler correction,
          //This one is detectorwise!
          float corrected_energy_crystal = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                                     ver_rec[0],ver_rec[1],ver_rec[2], 
                                                                     pos_det[0],pos_det[1],pos_det[2], 
                                                                     dali2_energy_not_cor[nnn],beta_rec,beta_mean_tof,
                                                                     beta_average,decay_z);
          h_dali2_doppler->Fill(corrected_energy_crystal);
          h_dali2_crystal_fired_doppler->Fill(nnn,corrected_energy_crystal);
          
          dali2_crystal_mult++;
          
          dali2DopplerCorrectedEnergy[dali2_crystal_mult-1] = corrected_energy_crystal;
                    
          //See how addback procedures work
          if(dali2_energy_not_cor[nnn]>dali2_energy_max){
            dali2_energy_max=dali2_energy_not_cor[nnn];
            for(int a=0;a<3;a++)
              dali2_pos_thatsit[a]=dali2_pos[a][nnn];
            //for(int ppp=0;ppp<15;ppp++) {
            //  if(nnn>=dali2_ring_low[ppp] && nnn<=dali2_ring_high[ppp]) 
            //    dali2_ring_energy_max = ppp;
            //}
          }
        }
      }
      h_dali2_crystal_mult->Fill(dali2_crystal_mult);
      gammaDet->fold[0] = dali2_crystal_mult;
      if(dali2_crystal_mult>0)  {
        //Performing the Doppler-correction from the angles
        //of the detector that registered the highest gamma-ray energy.
        float corrected_energy = GetDopplerCorrectedEnergy(dali2_pos_thatsit[0],dali2_pos_thatsit[1],dali2_pos_thatsit[2],
                                                           ver_rec[0],ver_rec[1],ver_rec[2], 
                                                           pos_det[0],pos_det[1],pos_det[2], 
                                                           dali2_energy_sum,beta_rec,beta_mean_tof,beta_average,decay_z);
        h_dali2_doppler_total->Fill(corrected_energy);
                
      }
    }
    // Ok, need to go a second time through the loop when all gamma detectors were present:
    if(eventNumberNext!=evnum)  {  
      gammaDet->event = evnum;
      for(int i=0;i<dali2_crystal_mult&&i<20;i++) {
        

        int nnn = dali2_crystalFired[i];
        // Have to calculate the Doppler energy:
        float corrected_energy_crystal = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                                     ver_rec[0],ver_rec[1],ver_rec[2], 
                                                                     pos_det[0],pos_det[1],pos_det[2], 
                                                                     dali2_energy[nnn],beta_rec,beta_mean_tof,
                                                                     beta_average,decay_z);

          //____________
          //if((iii+1)%1000 == 0) cout<<"Dali2, evnum: "<<evnum<<endl;
          h_evnum->Fill(evnum);

          gammaDet->id[0][i] = nnn;
          gammaDet->e[0][i] = dali2_energy[nnn];
          gammaDet->d[0][i] = corrected_energy_crystal;
          gammaDet->reste[0][i] = dali2_e_rest[nnn];
          
          gammaDet->x[0][i] = dali2_pos[0][nnn];
          gammaDet->y[0][i] = dali2_pos[1][nnn];
          gammaDet->z[0][i] = dali2_pos[2][nnn];
          //___________
        }
    }
         
    //End Dali2
    //---------------------------------------------------------------------------------------------
    //Begin SHOGUN
    if(Shogun_opt==1)  { 
      for(int nnn=0;nnn<NUMBEROFSHOGUNCRYSTALS;nnn++)  {
        if(Shogun_flag[nnn]==1) {
          Shogun_crystal_mult++;
          Shogun_energy_sum += Shogun_energy_not_cor[nnn];
          if(Shogun_energy_not_cor[nnn]>Shogun_energy_max)  {
            Shogun_energy_max=Shogun_energy_not_cor[nnn];
            for(int a=0;a<3;a++)  
              Shogun_pos_thatsit[a]=Shogun_pos[a][nnn];
          }
          float corrected_energy = GetDopplerCorrectedEnergy(Shogun_pos_thatsit[0],Shogun_pos_thatsit[1],Shogun_pos_thatsit[2], 
                                                           ver_rec[0],ver_rec[1],ver_rec[2], 
                                                           pos_det[0],pos_det[1],pos_det[2], 
                                                           Shogun_energy_not_cor[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
          h_shogun_doppler->Fill(corrected_energy);
        }
      }
      h_shogun_crystal_mult->Fill(Shogun_crystal_mult);
      if(Shogun_crystal_mult>0)  {
        //Performing the Doppler-correction from the angles
        //of the detector that registered the highest gamma-ray energy.
        float corrected_energy = GetDopplerCorrectedEnergy(Shogun_pos_thatsit[0],Shogun_pos_thatsit[1],Shogun_pos_thatsit[2], 
                                                           ver_rec[0],ver_rec[1],ver_rec[2], 
                                                           pos_det[0],pos_det[1],pos_det[2], 
                                                           Shogun_energy_sum,beta_rec,beta_mean_tof,beta_average,decay_z);
        h_shogun_doppler_total->Fill(corrected_energy);      
        h_shogun_energy_total->Fill(Shogun_energy_sum);
        if(Shogun_crystal_mult==1) h_shogun_doppler_single->Fill(corrected_energy);
      }
    }
    //End Shogun
    //-------------------------------------------------------------
    // Begin Grape
    if(grape_opt==1) {
      for(int nnn=0;nnn<NUMBEROFGRAPEDETECTORS;nnn++)  {
        if(grape_flag[nnn]==1.0) grape_detector_mult++;
        for(int ppp=0;ppp<2;ppp++)  {
          if(grape_crystal_flag[nnn][ppp][0]==1)  {
            //cout<<"nnn :"<<nnn<<endl;
            grape_crystal_mult++;
            grape_energy_sum = grape_energy_sum+ grape_energy_not_cor[nnn][ppp][0];
            grape_energy_sum_det[nnn] = grape_energy_sum_det[nnn] + grape_energy_not_cor[nnn][ppp][0];
          }
          for(int qqq=0;qqq<9;qqq++)  {
            if(grape_crystal_flag[nnn][ppp][qqq+1]==1.0) grape_segment_mult++;
            // The highest energy release in the segment of the total array
            if(grape_energy_not_cor[nnn][ppp][qqq+1]>grape_energy_max)  {
              grape_energy_max=grape_energy_not_cor[nnn][ppp][qqq+1];
              for(int a=0;a<3;a++)  
                grape_pos_thatsit[a]=grape_pos[a][nnn][ppp][qqq+1];
            }
            // The highest energy release in the segment of the detector [nnn]
            if(grape_energy_not_cor[nnn][ppp][qqq+1]>grape_energy_max_det[nnn])  {
              grape_energy_max_det[nnn]=grape_energy_not_cor[nnn][ppp][qqq+1];
              for(int a=0;a<3;a++)  
                grape_pos_thatsit_det[a][nnn]=grape_pos[a][nnn][ppp][qqq+1];
            }
          }
        }
        // The uncorrected energy per detector
        if(grape_flag[nnn]==1)  {
          h_grape_energy_det[nnn]->Fill(grape_energy_sum_det[nnn]);
          // Doppler correction per detector
          float corrected_energy_det = GetDopplerCorrectedEnergy(grape_pos_thatsit_det[0][nnn],
                                                                 grape_pos_thatsit_det[1][nnn],
                                                                 grape_pos_thatsit_det[2][nnn],
                                                                 ver_rec[0],ver_rec[1],ver_rec[2], 
                                                                 pos_det[0],pos_det[1],pos_det[2], 
                                                                 grape_energy_sum_det[nnn],beta_rec,beta_mean_tof,
                                                                 beta_average,decay_z); 
          h_grape_doppler_cor_det[nnn]->Fill(corrected_energy_det);
        }
      }
      h_grape_segment_mult->Fill(grape_segment_mult);
      h_grape_crystal_mult->Fill(grape_crystal_mult);
      h_grape_detector_mult->Fill(grape_detector_mult);
      h_grape_energy->Fill(grape_energy_sum);
      float corrected_energy = GetDopplerCorrectedEnergy(grape_pos_thatsit[0],grape_pos_thatsit[1],grape_pos_thatsit[2],
                                                         ver_rec[0],ver_rec[1],ver_rec[2], 
                                                         pos_det[0],pos_det[1],pos_det[2], 
                                                         grape_energy_sum,beta_rec,beta_mean_tof,beta_average,decay_z);
      
      h_grape_doppler_cor->Fill(corrected_energy);
      // Checking resolution, efficiency with condition on gamma-ray multiplicity
      if(grape_segment_mult==1) h_grape_doppler_cor_single_segment->Fill(corrected_energy);
      if(grape_segment_mult==2) h_grape_doppler_cor_two_segments->Fill(corrected_energy);
      if(grape_crystal_mult==1) h_grape_doppler_cor_single_crystal->Fill(corrected_energy);
      if(grape_crystal_mult==2) h_grape_doppler_cor_two_crystals->Fill(corrected_energy);
      if(grape_detector_mult==1) h_grape_doppler_cor_single_detector->Fill(corrected_energy);
    } 
    //End Grape
    //-------------------------------------------------------------
    // Begin SGT
    // Reminder of the spectra to be filled:
    //h_sgt_coax_mult; OK
    //h_sgt_strip_mult[1]; OK
    //h_sgt_strip_fan[1]; OK
    //h_sgt_doppler_cor_det[1]; OK
    //h_sgt_doppler_cor; OK
    // Remingder of the variables of the ROOT-tree:
    // sgt_flag[1], sgt_coax_crystal_flag[1], sgt_planar_strip_flag[1][26]; //fictive central contact plus the 25 strips
    // sgt_x[1][27],sgt_y[1][27],sgt_z[1][27]; //first the coax, then the fictive central planar contact, then the strips
    // sgt_energy_not_cor[1][27]; //first the coax, then the planar (sum of the strips), then the strips
    if(sgt_opt==1)  {
      for(int nnn=0;nnn<NUMBEROFSGTDETECTORS;nnn++)  {
        if(sgt_flag[nnn]==1) sgt_detector_mult++;
        if(sgt_coax_crystal_flag[nnn]==1)  {
          sgt_coax_mult++;
          sgt_energy_sum = sgt_energy_sum + sgt_energy_not_cor[nnn][0];
          sgt_energy_sum_det[nnn] = sgt_energy_sum_det[nnn] + sgt_energy_not_cor[nnn][0];
          h_sgt_coax_energy_det[nnn]->Fill(sgt_energy_not_cor[nnn][0]);
        }
        for(int ppp=1;ppp<26;ppp++)  {
          if(sgt_planar_strip_flag[nnn][ppp]==1)  {
            sgt_strip_mult++;
            sgt_strip_mult_det[nnn]++;
            //cout<<"StripFlag true. sgt_x["<<nnn<<"]["<<ppp+1<<"]: "<<sgt_x[nnn][ppp+1]<<endl;
            h_sgt_strip_fan[nnn]->Fill(ppp);
            // If there is only one gammaray simulated and only one detector, the next line is fine:
            sgt_energy_sum = sgt_energy_sum + sgt_energy_not_cor[nnn][ppp+1];
            sgt_energy_sum_det[nnn] = sgt_energy_sum_det[nnn] + sgt_energy_not_cor[nnn][ppp+1];
            if(sgt_energy_not_cor[nnn][ppp+1]>= sgt_energy_max)  {
              sgt_energy_max=sgt_energy_not_cor[nnn][ppp+1];
              for(int a=0;a<3;a++) 
                sgt_pos_thatsit[a]=sgt_pos[a][nnn][ppp+1];
            }
            // For individuall detectors:
            if(sgt_energy_not_cor[nnn][ppp+1]> sgt_energy_max_det[nnn])  {
              sgt_energy_max_det[nnn]=sgt_energy_not_cor[nnn][ppp+1];
              for(int a=0;a<3;a++)  
                sgt_pos_thatsit_det[a][nnn]=sgt_pos[a][nnn][ppp+1];
            } 
          }
        }
        if(sgt_strip_mult_det[nnn]>=1)  {
          float corrected_energy = GetDopplerCorrectedEnergy(sgt_pos_thatsit_det[0][nnn],
                                                             sgt_pos_thatsit_det[1][nnn],sgt_pos_thatsit_det[2][nnn],
                                                             ver_rec[0],ver_rec[1],ver_rec[2], 
                                                             pos_det[0],pos_det[1],pos_det[2], 
                                                             sgt_energy_sum_det[nnn],beta_rec,beta_mean_tof,
                                                             beta_average,decay_z);
          h_sgt_doppler_cor_det[nnn]->Fill(corrected_energy);
          h_sgt_doppler_cor->Fill(corrected_energy);
        }
        h_sgt_strip_mult[nnn]->Fill(sgt_strip_mult);
        sgt_strip_mult=0;
      }
      h_sgt_coax_mult->Fill(sgt_coax_mult);
    }
    // SGT finished
    //-------------------------------------------------------------
    // Trying a perfect sphere
    if(sphere_opt==1)  {
      //float corrected_energy = GetDopplerCorrectedEnergy(sphere_fi_x,sphere_fi_y,sphere_fi_z,
      float corrected_energy = GetDopplerCorrectedEnergy(gv[0]*1000,gv[1]*1000,gv[2]*1000,
                                                         ver_rec[0],ver_rec[1],ver_rec[2], 
                                                         pos_det[0],pos_det[1],pos_det[2], 
                                                         sphere_energy_not_cor,beta_rec,beta_mean_tof,beta_average,decay_z);
      h_sphere_doppler_cor->Fill(corrected_energy);

      for(int i = 0;i<18;i++){
        if(theta_gamma_lab>(i*10*3.14159/180) && theta_gamma_lab<=((i+1)*10*3.14159/180)) 
          h_sphere_doppler_cor_angle_cut[i]->Fill(corrected_energy);
      }
    }
    // Sphere finished
    //-------------------------------------------------------------
    // The LaBr3 array
    if(LaBr3Opt==1)  {
      for(int nnn=0;nnn<NUMBEROFLABR3ARRAYCRYSTALS;nnn++) {
        if(LaBr3Flag[nnn]==1)   {
          labr3_energy[nnn] += LaBr3EnergyNotCor[nnn];
          labr3_e_rest[nnn] += e_rest;

          labr3_energy_sum += LaBr3EnergyNotCor[nnn];
          labr3_crystalFired[labr3_crystal_mult]=nnn;
          //Performing the Doppler correction,
          //This one is detectorwise!
          float corrected_energy = GetDopplerCorrectedEnergy(LaBr3Pos[0][nnn],LaBr3Pos[1][nnn],LaBr3Pos[2][nnn],
                                                             ver_rec[0],ver_rec[1],ver_rec[2], 
                                                             pos_det[0],pos_det[1],pos_det[2], 
                                                             LaBr3EnergyNotCor[nnn],beta_rec,beta_mean_tof,
                                                             beta_average,decay_z);
          h_labr3_doppler_cor->Fill(corrected_energy);

          labr3_crystal_mult++;
        }
      }
      gammaDet->fold[1] = labr3_crystal_mult;
    }
    // Ok, need to go a second time through the loop when all gamma detectors were present:
    if(eventNumberNext!=evnum)  {  
      gammaDet->event = evnum;
      for(int i=0;i<labr3_crystal_mult&&i<20;i++) {
      
        int nnn = labr3_crystalFired[i];
        // Have to calculate the Doppler energy:
        float corrected_energy_crystal = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                                   ver_rec[0],ver_rec[1],ver_rec[2], 
                                                                   pos_det[0],pos_det[1],pos_det[2], 
                                                                   labr3_energy[nnn],beta_rec,beta_mean_tof,
                                                                   beta_average,decay_z);
        
        //____________
        //if((iii+1)%1000 == 0) cout<<"LaBr3, evnum: "<<evnum<<endl;
        h_evnum->Fill(evnum);
        
        gammaDet->id[1][i] = nnn;
        gammaDet->e[1][i] = labr3_energy[nnn];
        gammaDet->d[1][i] = corrected_energy_crystal;
        gammaDet->reste[1][i] = labr3_e_rest[nnn];
        
        gammaDet->x[1][i] = LaBr3Pos[0][nnn];
        gammaDet->y[1][i] = LaBr3Pos[1][nnn];
        gammaDet->z[1][i] = LaBr3Pos[2][nnn];
        //___________
      }
      
    }
    // LaBr3 finished
    //-------------------------------------------------------------
    // Filling the tree
    if((dali2_crystal_mult>0 || labr3_crystal_mult > 0) && eventNumberNext!=evnum)  {
      t2->Fill();
      //if((iii+1)%1000 == 0) cout<<"Filling evnum: "<<evnum<<endl;
    }

    // The variables have to be reset.
    dali2_ring_energy_max = -1;
    for(int i =0;i<NUMBEROFSGTDETECTORS;i++) 
      dali2DopplerCorrectedEnergy[i] = 0.;

    dali2_crystal_mult  = 0;
    Shogun_crystal_mult = 0;
    grape_segment_mult  = 0;
    grape_crystal_mult  = 0;
    grape_detector_mult = 0;
    sgt_strip_mult = 0;
    for(int i =0;i<10;i++)
      sgt_strip_mult_det[i]=0;
    sgt_coax_mult = 0;
    sgt_detector_mult = 0;
    labr3_crystal_mult=0;

    dali2_energy_sum  = 0.0;
    Shogun_energy_sum = 0.0;
    grape_energy_sum  = 0.0;
    sgt_energy_sum    = 0.0;
    labr3_energy_sum  = 0.;

    dali2_energy_max = -999.9;
    Shogun_energy_max = -999.9;
    grape_energy_max = -999.9;
    for(int i =0;i<NUMBEROFGRAPEDETECTORS;i++)  {
      grape_energy_max_det[i] = -999;
      grape_energy_sum_det[i] = 0.0;
    }
    sgt_energy_max = -999.9;
    for(int i =0;i<10;i++)  {
      sgt_energy_max_det[i] = -999;
      sgt_energy_sum_det[i] = 0.0;
    }

    for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
      dali2_energy[i]         = 0.0; 
      dali2_e_rest[i]         = 0.0;
      dali2_time[i]           = 0.0; 
      dali2_dopplerEnergy[i]  = 0.0; 
      dali2_crystalFired[i]   = -1;  
    }  

    for(int i=0;i<NUMBEROFLABR3ARRAYCRYSTALS;i++)  {
      labr3_energy[i]        = 0.0; 
      labr3_e_rest[i]        = 0.0;
      labr3_time[i]          = 0.0; 
      labr3_dopplerEnergy[i] = 0.0; 
      labr3_crystalFired[i]  = -1;   
    }

  }//End of loop
  outfile2->cd();  
  outfile2->Close();

  //---------------------------------------------------------------
  //Writing into the root file:
  rootfile->cd();
  h_evnum->Write();
  //---------------------------------------------------------------
  //Dali2-----------------------------------------------------------
  //---------------------------------------------------------------
  if(dali2_opt==1)  {
    h_dali2_crystal_mult->Write();     
    h_dali2_doppler->Write();
    h_dali2_doppler_total->Write(); 
    h_dali2_energy->Write();             
                            
    h_dali2_crystal_fired_doppler->Write();
    h_dali2_crystal_fired_energy->Write();
  }    
  //---------------------------------------------------------------
  //SHOGUN---------------------------------------------------------
  //---------------------------------------------------------------
  if(Shogun_opt==1)  {
    h_shogun_crystal_mult->Write();
    h_shogun_doppler->Write();
    h_shogun_doppler_single->Write();
    h_shogun_doppler_total->Write();
    h_shogun_energy_total->Write();
  }   
  //---------------------------------------------------------------
  ///////////////////////////////
  //Grape!!!!!!!!!!!!!!!!!!!!!
  //////////////////////////////
  if(grape_opt==1)  {
    h_grape_segment_mult->Write();
    h_grape_crystal_mult->Write();
    h_grape_detector_mult->Write();
    h_grape_doppler_cor->Write();
    h_grape_doppler_cor_single_segment->Write();
    h_grape_doppler_cor_single_crystal->Write();
    h_grape_doppler_cor_single_detector->Write();
    h_grape_doppler_cor_two_segments->Write();
    h_grape_doppler_cor_two_crystals->Write();
    for(int i=0;i<NUMBEROFGRAPEDETECTORS;i++)
      h_grape_doppler_cor_det[i]->Write();
    h_grape_energy->Write();
    for(int i=0;i<NUMBEROFGRAPEDETECTORS;i++)
      h_grape_energy_det[i]->Write();
  }
  //---------------------------------------------------------------	      
  ///////////////////////////////
  //SGT!!!!!!!!!!!!!!!!!!!!!!!!!!
  ///////////////////////////////
  if(sgt_opt==1)  {
    h_sgt_coax_mult->Write();
    for(int i=0;i<NUMBEROFSGTDETECTORS;i++)  {
      h_sgt_strip_mult[i]->Write();
      h_sgt_strip_fan[i]->Write();
      h_sgt_doppler_cor_det[i]->Write();
      h_sgt_coax_energy_det[i]->Write();
    }
    h_sgt_doppler_cor->Write();
  }
  //---------------------------------------------------------------	      
  ///////////////////////////////
  //Sphere///////////////////////
  ///////////////////////////////
  if(sphere_opt==1)  {
    h_sphere_doppler_cor->Write();
    for(int i = 0;i<18;i++)
      h_sphere_doppler_cor_angle_cut[i]->Write();
  }
  //---------------------------------------------------------------	      
  ///////////////////////////////
  //LaBr3////////////////////////
  ///////////////////////////////
  if(LaBr3Opt==1)  {
    h_labr3_doppler_cor->Write();
  }
}

//_________________________________________________________________________________________________
void ResetGammaDetectorValues()  {
  for(int i=0;i<2;i++)  {
    for(int j=0;j<MAXFOLD;j++)  {
      gammaDet->id[i][j] = -1;
      gammaDet->fold[i] = 0;
      
      gammaDet->e[i][j] = 0;
      gammaDet->d[i][j] = 0;
      gammaDet->reste[i][j] = 0;
      
      gammaDet->x[i][j] = -999.0;
      gammaDet->y[i][j] = -999.0;
      gammaDet->z[i][j] = -999.0;
    }
  }
}

