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

double thetaLast;

//_________________________________________________________________________________________________
float GetDopplerCorrectedEnergy(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                                float target_x,float target_y,float target_z,
                                float pos_det_x, float pos_det_y, float pos_det_z, 
                                float gamma_energy, float beta_rec,float beta_mean_tof,float beta_average,float decay_z){ 
  double vx1 = gamma_det_x - target_x;                         
  double vy1 = gamma_det_y - target_y;                         
  double vz1 = gamma_det_z - target_z - decay_z;               

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  double vx2 = pos_det_x - target_x;                           
  double vy2 = pos_det_y - target_y;                           
  double vz2 = pos_det_z - target_z;                           
 
  double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
          
  // The theta angle:
  double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);
  thetaLast = theta_lab * 180./3.14159;                        //the last calculated theta value in degrees;
  // the beta value:
  double beta = beta_average + (beta_rec-beta_mean_tof);         //cout<<"beta: "<<beta<<endl;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

const  int maxDecays = 5;
double fi_average[3][NUMBEROFSHOGUNCRYSTALS] = {{0.0}};
int    fi_interactions[NUMBEROFSHOGUNCRYSTALS] = {0};
float  Shogun_pos_ave[3][NUMBEROFSHOGUNCRYSTALS];

int ShogunReconstructorSimple(){ 

  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
  float beta_average     = 0.0;        // The average value at the moment of decay for the correction
  float beta_mean_tof    = 0.0;        // The mean beta from Tof
  float decay_z          = 0.0;        // The lifetime moves the average decay point along the beam axis
  int fi_option          = 0;
 
  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/ShogunReconstructor.in","r");
  while(!feof(fin)){
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
    else if(strcmp(temp,"FIFIND")==0)  {
      fscanf(fin,"%i",&fi_option); 
      printf("%s %i \n",temp,fi_option);
    }

    else if(strcmp(temp,"END")==0) break;
    else {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      return 0;
    }
  }

  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();
  
  TH1F *h_crystal_mult;
  TH1F *h_energy;

  TH1F *h_doppler;  

  TH1F *h_doppler_total;
  TH1F *h_doppler_total_individual_decays; 

  TH2F *h_doppler_gamma_gamma;

  //Creating SHOGUN spectra
  // Crystal multiplicity,ok
  h_crystal_mult = new TH1F("crystal_Mult","crystal_mult",NUMBEROFSHOGUNCRYSTALS,0,NUMBEROFSHOGUNCRYSTALS);

  // Not doppler corrected,
  h_energy = new TH1F("energy","energy",numBin,firstBin,lastBin);
  
  // All crystals, doppler corrected 
  h_doppler = new TH1F("doppler","doppler",numBin,firstBin,lastBin);

  // Taking the sum energy and making the Doppler correction with the crystal that has the highest energy:
  h_doppler_total = new TH1F("doppler_total","doppler_total",numBin,firstBin,lastBin);
  h_doppler_total_individual_decays = new TH1F("doppler_total_individual_decays","doppler_total_individual_decays",numBin,firstBin,lastBin);
  
  //For gamma-gamma analysis:
  h_doppler_gamma_gamma = new TH2F("doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  infile = new TFile(root_in,"READ");
  t = (TTree*)infile->Get("ObservedEvents");

  t->SetBranchAddress("EventNumber",&evnum);
  t->SetBranchAddress("DecayIDOfEventNumber",&decayIDevnum);                    // Labeling of the different decays
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
  //So far identical to the ROOT-Tree of step one
  //The new stuff:
  //Setting the type of the gamma detector
  t->SetBranchAddress("GammaDetType",&gamma_det_type);
  //-------------------------------
  //Shogun:
  t->SetBranchAddress("ShogunFlag",Shogun_flag);                                  // ID of which det recorded a gamma in the event
  t->SetBranchAddress("ShogunEnergyNotCor",Shogun_energy_not_cor);
  t->SetBranchAddress("ShogunTime",Shogun_time);
  if(fi_option >= 1)  {
    t->SetBranchAddress("ShogunFI",Shogun_fi);
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
  //SHOGUN:
  tHeader->SetBranchAddress("ShogunEnResOpt",&Shogun_en_res_opt);
  tHeader->SetBranchAddress("ShogunEnRes",Shogun_en_res);
  tHeader->SetBranchAddress("ShogunTimeRes",Shogun_time_res);
  tHeader->SetBranchAddress("ShogunPos",Shogun_pos);
  tHeader->SetBranchAddress("BetaResolution",&beta_res);  
  tHeader->SetBranchAddress("PosDetAtTargetRes",&pos_det_at_target_res); 
  tHeader->SetBranchAddress("PosDetAfterTargetRes",&pos_det_after_target_res);
  tHeader->GetEntry(0);
  //-----------------------------------------------------------------------------------------------------------
  //Finished reading the Root file

  nentries = (Int_t)t->GetEntries()/reduction_factor;
  cout <<"Entries: "<< nentries << endl;

  int crystalMult = 0;
  
  float energy[NUMBEROFSHOGUNCRYSTALS]                 = {0.0};
  float time[NUMBEROFSHOGUNCRYSTALS]                   = {0.0};
  float dopplerEnergy[NUMBEROFSHOGUNCRYSTALS]          = {0.0};  // Doppler corrected energy for the crystals
  float dopplerEnergyIndividualDecays[maxDecays]       = {0.0};  // Every individual gamma decay branch
 
  float energySum = 0.0;
  float energySumIndividualDecays[maxDecays] = {0.0};
  float energyMax = -999.9;
  float energyMaxIndividualDecays[maxDecays] = {-999.};

  float posThatsIt[3] = {0.0};
  float posThatsItIndividualDecays[maxDecays][3] = {{0.0}};
  int eventNumberNext = 0;

  cout<<"Starting the analysis"<<endl;
  
  for(int iii=0;iii<nentries;iii++)  {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    if(iii<nentries-1)  {
      t->GetEntry(iii+1);
      eventNumberNext = evnum;
    }

    t->GetEntry(iii);
      
    //-------------------------------------------------------------
    for(int nnn=0;nnn<NUMBEROFSHOGUNCRYSTALS;nnn++)  {
      if(Shogun_flag[nnn]==1)   {
        energy[nnn] = energy[nnn] + Shogun_energy_not_cor[nnn];
        energySum   = energySum +  Shogun_energy_not_cor[nnn];
        if(maxDecays>decayIDevnum) {
          energySumIndividualDecays[(int)decayIDevnum] = energySumIndividualDecays[(int)decayIDevnum] + Shogun_energy_not_cor[nnn];
          // Getting the right position for the Doppler correction of individual decays:
          if(Shogun_energy_not_cor[nnn]>energyMaxIndividualDecays[(int)decayIDevnum])  {
	    energyMaxIndividualDecays[(int)decayIDevnum] = Shogun_energy_not_cor[nnn];
	    posThatsItIndividualDecays[(int)decayIDevnum][0]  = Shogun_pos[0][nnn];
            posThatsItIndividualDecays[(int)decayIDevnum][1]  = Shogun_pos[1][nnn];
            posThatsItIndividualDecays[(int)decayIDevnum][2]  = Shogun_pos[2][nnn];
	  }
        }
      }
    }
    if(eventNumberNext!=evnum)  {                 // The last observed decay we can start to fill the spectra
      for(int nnn=0;nnn<NUMBEROFSHOGUNCRYSTALS;nnn++)  {
        if(energy[nnn]>0)  {
       	  // Filling the energy spectrum without Doppler-correction        
	  h_energy->Fill(energy[nnn]);
          
   	  //Performing the Doppler correction,
	  //This one is detectorwise
          float correctedEnergyCrystal;
          float correctedEnergySimple;
          correctedEnergyCrystal = GetDopplerCorrectedEnergy(Shogun_pos[0][nnn],Shogun_pos[1][nnn],Shogun_pos[2][nnn],
                                                             ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                             energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
          
          dopplerEnergy[crystalMult] = correctedEnergyCrystal;
          h_doppler->Fill(correctedEnergyCrystal);
          
	  dopplerEnergy[crystalMult] = correctedEnergyCrystal ;
          //Set the threshold
          //if(energy[nnn]>150)
          crystalMult++;
	  
	  if(energy[nnn]>energyMax)  {
	    energyMax = energy[nnn];
	    posThatsIt[0]  = Shogun_pos[0][nnn];
            posThatsIt[1]  = Shogun_pos[1][nnn];
            posThatsIt[2]  = Shogun_pos[2][nnn];
	  }
	}
      }
      h_crystal_mult->Fill(crystalMult);
      
      //____________________________________________________
      // Making a Doppler correction with the Sum of the Energy detected in the array
      // The Position is taken from the crystal that recorded the highest energy
      float correctedEnergyTotal = GetDopplerCorrectedEnergy(posThatsIt[0],posThatsIt[1],posThatsIt[2],
                                                             ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                             energySum,beta_rec,beta_mean_tof,beta_average,decay_z);
      h_doppler_total->Fill(correctedEnergyTotal);
      
                
      int counter = 0;
      // Treating the decays individually:
      while(energySumIndividualDecays[counter]>0 && counter<maxDecays)  {
        dopplerEnergyIndividualDecays[counter] =  GetDopplerCorrectedEnergy(posThatsItIndividualDecays[counter][0],
                                                                            posThatsItIndividualDecays[counter][1],
                                                                            posThatsItIndividualDecays[counter][2],
                                                                            ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                                            energySumIndividualDecays[counter],
                                                                            beta_rec,beta_mean_tof,beta_average,decay_z);
        h_doppler_total_individual_decays->Fill(dopplerEnergyIndividualDecays[counter]);
        counter++;
      }
      //____________________________________________________
      //Making the gamma-gamma coincidencs
      if(crystalMult==2)  {
        if (dopplerEnergy[0]>=dopplerEnergy[1]) 
          h_doppler_gamma_gamma->Fill(dopplerEnergy[0],dopplerEnergy[1]);
        else h_doppler_gamma_gamma->Fill(dopplerEnergy[1],dopplerEnergy[0]);
      }
      //____________________________________________________
      // The variables have to be reset.
      crystalMult   = 0;
      energySum     = 0.;
      energyMax     = 0.;
      for(int i =0;i<maxDecays;i++)  {
        dopplerEnergyIndividualDecays[i] = 0.0;
        energySumIndividualDecays[i] = 0.0;
        energyMaxIndividualDecays[i] = 0.0;
      }
      for(int i =0;i<NUMBEROFSHOGUNCRYSTALS;i++){
        energy[i]         = 0.;
        time[i]           = 0.;
        dopplerEnergy[i]  = 0.;
      }
    }
  }
  
  //---------------------------------------------------------------
  // Writing into the root file:
  // 
  rootfile->cd();
  h_crystal_mult->Write();
  h_energy->Write();
  h_doppler->Write();
  
  h_doppler_total->Write();
  h_doppler_total_individual_decays->Write();

  return 0;
} 
