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
  double vz2 = pos_det_z - target_z - decay_z;                 //cout<<"vy2: "<<vz2<<endl;
 
  double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
          
  // The theta angle:
  double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);//cout<<"Theta lab: "<<theta_lab<<endl;
  // the beta value:
  double beta = beta_average + (beta_rec-beta_mean_tof);         //cout<<"beta: "<<beta<<endl;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}
void GetAverageInteractionPoint();
const int numberOfCrystals = 182;
double fi_average[3][numberOfCrystals] = {{0.0}};
int    fi_interactions[numberOfCrystals] = {0};
float dali2_pos[3][numberOfCrystals];
float dali2_pos_ave[3][numberOfCrystals];

void Dali2Reconstructor(){ 
  char root_in[200];
  char root_out[200];
  char temp[200];
  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
  float beta_average     = 0.0;        // The average value at the moment of decay for the correction
  float beta_mean_tof    = 0.0;        // The mean beta from Tof
  float decay_z          = 0.0;        // The lifetime moves the average decay point along the beam axis
  int fi_option          = 0;
 
  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/Dali2Reconstructor.in","r");
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
      if(reduction_factor<1) return;
    }
    else if(strcmp(temp,"FIFIND")==0)  {
      fscanf(fin,"%i",&fi_option); 
      printf("%s %i \n",temp,fi_option);
    }

    else if(strcmp(temp,"END")==0) break;
    else {
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
  float p0[3],pvb[3],pva[3];
  float energy_p;
  float beta_b,beta_r,beta_a,halflife,decay_time_after_interaction;
  float g0[3],gv[3];
  float e_rest,e_doppler;
  float theta_gamma_rest, theta_gamma_lab;
  float energy_vertex_stored;
  int kind_target;
  //***********************************************
  float total_event_num;//energy_notcor;
  int gamma_det_type;  // 1 for dali, 2 for grape, 3 for sgt 
  float pos_det[3];     //reconstructed HI position from detector after target
  float ver_rec[3];
  float beta_rec;
  //-----------------------------------------------
  // For the Dali2 array
  int dali2_flag[numberOfCrystals];
  float dali2_fi[3][numberOfCrystals];
  float dali2_energy_not_cor[numberOfCrystals];
  float dali2_time[numberOfCrystals] = {0.};
  // Input of the resolutions:
  int dali2_en_res_opt = 0;
  float dali2_en_res[2];
  float dali2_time_res[2];
  float pos_det_at_target_res;
  float pos_det_after_target_res; // Pos Resolution in mm FWHM!
  float beta_res; //Resolution of beta in FWHM!

  //End, same variables
  //----------------------------------------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------------------------------
  char tempname[100];

  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();

  //Testing the position
  TH2F *h_gamma_vector_xy = new TH2F("gamma_vector_xy","",200,-1.1,1.1,200,-1.1,1.1);

  //Creating Beta and gamma lab spectra!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TH2F *h_beta_before_beta_after = new TH2F("beta_before_beta_after","beta_before_beta_after",200,0.50,0.65,200,0.50,0.65);
  TH2F *h_beta_real_beta_before = new TH2F("beta_real_beta_before","beta_real_beta_before",200,0.50,0.65,200,0.50,0.65);
  TH2F *h_beta_real_beta_after = new TH2F("beta_real_beta_after","beta_real_beta_after",200,0.50,0.65,200,0.50,0.65);
  
  TH1F *h_crystal_mult;
  TH1F *h_doppler;
  TH1F *h_doppler_ave_pos;
  TH1F *h_doppler_mult[20];
  TH1F *h_doppler_crystal[numberOfCrystals];
  TH1F *h_doppler_ring[17];
  TH2F *h_doppler_max_ring;
  TH1F *h_energy;
  TH1F *h_energy_crystal[numberOfCrystals];
  TH2F *h_energy_beta_real;

  TH2F *h_beta_after_doppler;
  TH2F *h_crystal_fired_doppler;

  
  //Creating dali spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Crystal multiplicity,ok
  h_crystal_mult = new TH1F("crystal_Mult","crystal_mult",numberOfCrystals,0,numberOfCrystals);  
  // All, OK
  h_doppler = new TH1F("doppler","doppler",numBin,firstBin,lastBin);
  
  // According to the multiplicity
  for(int i=0;i<20;i++)  {
    sprintf(tempname,"doppler_mult[%i]",i);
    h_doppler_mult[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  } 

  // The spectra crystalwise, OK
  for(int i=0;i<numberOfCrystals;i++)  {
    sprintf(tempname,"doppler_crystal[%i]",i);
    h_doppler_crystal[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }  
  //Ringwise:
  for(int i=0;i<17;i++)  {
    sprintf(tempname,"doppler_ring[%i]",i);
    h_doppler_ring[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }
  //The ring with the most energy:
  //for(int i=0;i<15;i++)
  //{
  sprintf(tempname,"doppler_max_ring");
  h_doppler_max_ring = new TH2F(tempname,tempname,numBin,firstBin,lastBin,15,0,15);
  //}

  //Not doppler corrected, OK
  h_energy = new TH1F("energy","energy",numBin,firstBin,lastBin);
  // The spectra crystalwise,OK
  for(int i=0;i<numberOfCrystals;i++)  {
    sprintf(tempname,"energy_crystal[%i]",i);
    h_energy_crystal[i] = new TH1F(tempname,tempname,numBin,firstBin,lastBin);
  }
  //comparing beta REAL with dali energy:
  h_energy_beta_real = new TH2F("energy_beta_real","energy_beta_real",numBin,firstBin,lastBin,200,0.40,0.45);

  //comparing beta after the target with doppler corrected energy
  h_beta_after_doppler = new TH2F("beta_after_doppler","",200,0.5,0.64,numBin,firstBin,lastBin);

  //comparing beta after the target with doppler corrected energy
  h_crystal_fired_doppler = new TH2F("crystal_fired_doppler","",numberOfCrystals,0,numberOfCrystals,numBin,firstBin,lastBin);

  //______________________________________________
  if(fi_option ==1)  {
    FILE *fInAve  = fopen("./output/Dali2AverageInteractionPoint.out","r");
    int counter = 0;
    float x,y,z;
    int counts;
    while(!feof(fInAve) || counter<numberOfCrystals-1){
      fscanf(fInAve,"%i %i %f %f %f",&counter,&counts,&x,&y,&z); 
      //counter++;                                            
      dali2_pos_ave[0][counter] = x;
      dali2_pos_ave[1][counter] = y;
      dali2_pos_ave[2][counter] = z;
      cout<<counter<<" "<<dali2_pos_ave[0][counter]<<" "<<dali2_pos_ave[1][counter]<<" "<<dali2_pos_ave[2][counter]<<endl;
    }
    h_doppler_ave_pos = new TH1F("doppler_ave_pos","doppler_ave_pos",numBin,firstBin,lastBin);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile *infile = new TFile(root_in,"READ");
  TTree *t = (TTree*)infile->Get("ObservedEvents");

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
  //So far identical to the ROOT-Tree of step one
  //The new stuff:
  //Setting the type of the gamma detector
  t->SetBranchAddress("GammaDetType",&gamma_det_type);
  //-------------------------------
  //DALI:
  t->SetBranchAddress("DALI2Flag",dali2_flag);                                  // ID of which det recorded a gamma in the event
  t->SetBranchAddress("DALI2EnergyNotCor",dali2_energy_not_cor);
  t->SetBranchAddress("DALI2Time",dali2_time);
  if(fi_option == 1)  {
    t->SetBranchAddress("DALI2FI",dali2_fi);
  }
  //-------------------------------
  //Position detector after targetpos_det_after_target_res
  t->SetBranchAddress("PosDet",pos_det);                                        // Reconstructed position from a det after the secondary target
  t->SetBranchAddress("VertexReconstructed",ver_rec);                           // Reconstructed vertex postion, 
 
  //reconstructed beta 
  t->SetBranchAddress("BetaReconstructed",&beta_rec);                           // Reconstructed beta, including resolution

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
  //DALI2:
  tHeader->SetBranchAddress("DALI2EnResOpt",&dali2_en_res_opt);
  tHeader->SetBranchAddress("DALI2EnRes",dali2_en_res);
  tHeader->SetBranchAddress("DALI2TimeRes",dali2_time_res);
  tHeader->SetBranchAddress("DALI2Pos",dali2_pos);
  tHeader->SetBranchAddress("BetaResolution",&beta_res);  
  tHeader->SetBranchAddress("PosDetAtTargetRes",&pos_det_at_target_res); 
  tHeader->SetBranchAddress("PosDetAfterTargetRes",&pos_det_after_target_res);
  tHeader->GetEntry(0);
  //-----------------------------------------------------------------------------------------------------------
  //Finished reading the Root file

  Int_t nentries = (Int_t)t->GetEntries()/reduction_factor;
  cout <<"Entries: "<< nentries << endl;

  int crystalMult = 0;
  //Dividing old dali2 into 16 rings of equal z-position:
  //int dali2_ring_low[16]={0,6,12,20,32,45,56,70,80,90,104,116,128,140,148,154};
  //int dali2_ring_high[16]={5,11,19,31,44,55,69,79,89,103,115,127,139,147,153,159};
  //Dividing new dali2 into 17 rings of equal theta angle:
  int ringLow[17] = {0, 6,13,20,30,40,52,66,80,94,108,122,136,148,156,164,170};
  int ringHigh[17]= {5,12,19,29,39,51,65,79,93,107,121,135,147,155,163,169,181};
  int ringEnergyMax = -1;

  float energy[numberOfCrystals]                 = {0.0};
  float time[numberOfCrystals]                   = {0.0};
  float dopplerEnergy[numberOfCrystals]          = {0.0};
 
  float energySum = 0.0;
  float energyMax = -999.9;
   
  float posThatsit[3] = {0.0};
 
  int eventNumberNext = 0;

  cout<<"Starting the analysis"<<endl;
  for(int iii=0;iii<nentries;iii++)  {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    if(iii<nentries-1)  {
      t->GetEntry(iii+1);
      eventNumberNext = evnum;
    }

    t->GetEntry(iii);

    //Fill the gamma xy spectrum
    h_gamma_vector_xy->Fill(gv[0],gv[1]);

    // Fill the beta spectra:
    h_beta_before_beta_after->Fill(beta_b,beta_a);
    h_beta_real_beta_before->Fill(beta_r,beta_b);
    h_beta_real_beta_after->Fill(beta_r,beta_a);

    //-------------------------------------------------------------
    // Starting with the gamma ray detector spectra:
    for(int nnn=0;nnn<numberOfCrystals;nnn++)  {
      if(dali2_flag[nnn]==1)   {
        energy[nnn] = energy[nnn] + dali2_energy_not_cor[nnn];
      }
    }
    if(eventNumberNext!=evnum)  {                                               // The last observed decay, we can start to fill the spectra
      for(int nnn=0;nnn<numberOfCrystals;nnn++)  {
        if(energy[nnn]>0)  {
       	  // Filling the energy spectra without Doppler-correction        
	  h_energy->Fill(energy[nnn]);
	  h_energy_crystal[nnn]->Fill(energy[nnn]);
	  h_energy_beta_real->Fill(energy[nnn],beta_r);
   	  //Performing the Doppler correction,
	  //This one is detectorwise!
	  float correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                                   energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
	  h_doppler->Fill(correctedEnergyCrystal);
          h_doppler_crystal[nnn]->Fill(correctedEnergyCrystal);
          h_crystal_fired_doppler->Fill(nnn,correctedEnergyCrystal);
          
	  dopplerEnergy[crystalMult] = energy[nnn];
          crystalMult++;
	  //Filling Dali ringwise:
	  for(int ppp=0;ppp<17;ppp++)  {
	    if(nnn>=ringLow[ppp] && nnn<=ringHigh[ppp]) h_doppler_ring[ppp]->Fill(correctedEnergyCrystal);
	  }
	  //See how addback procedures work
	  if(energy[nnn]>energyMax)  {
	    energyMax = energy[nnn];
	    posThatsit[0]  = dali2_pos[0][nnn];
            posThatsit[1]  = dali2_pos[1][nnn];
            posThatsit[2]  = dali2_pos[2][nnn];
	    for(int ppp=0;ppp<17;ppp++) {if(nnn>=ringLow[ppp] && nnn<=ringHigh[ppp]) ringEnergyMax = ppp;}
	  }
          //________________________________________________
          // Only if almost all the energy was released within one crystal the average position is determined:
          if(fi_option == 1 && correctedEnergyCrystal>= e_rest/1.1)  {
            for(int jjj=0;jjj<3;jjj++)  {
              fi_average[jjj][nnn] = fi_average[jjj][nnn] + dali2_fi[jjj][nnn]; 
              if(jjj==0) fi_interactions[nnn]++;
            }
          }
          //_________________________________________________
          // Checking the effect of getting the average position:
          if(fi_option ==1)  {
            float correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos_ave[0][nnn],dali2_pos_ave[1][nnn],dali2_pos_ave[2][nnn],
                                                                     ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                                     energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
            h_doppler_ave_pos->Fill(correctedEnergyCrystal);
          }
	}
      }
      h_crystal_mult->Fill(crystalMult);
      //if(crystalMult>0){
      for(int i = 1;i<=20;i++){
        if(crystalMult==i){	 
          for(int j = 0;j<i;j++){
            h_doppler_mult[i]->Fill(dopplerEnergy[j]);
          }
        }
      }
      //}
      // The variables have to be reset.
      ringEnergyMax = -1;
      crystalMult   = 0;
      for(int i =0;i<numberOfCrystals;i++){
        energy[i]         = 0.;
        time[i]           = 0.;
        dopplerEnergy[i]  = 0.;
      }
    }
  }
  //_______________________________________________________________
  // Getting the average point:
  if(fi_option == 1) GetAverageInteractionPoint();

  //---------------------------------------------------------------
  // Writing into the root file:
  // BETA
  rootfile->cd();
  h_gamma_vector_xy->Write();
  h_beta_before_beta_after->Write();
  h_beta_real_beta_before->Write();
  h_beta_real_beta_after->Write();
  //---------------------------------------------------------------
  // Dali----------------------------------------------------------
  //---------------------------------------------------------------
  h_crystal_mult->Write();
  for(int i=0;i<20;i++){
    h_doppler_mult[i]->Write();
  }
  h_doppler->Write();
  if(fi_option ==1) h_doppler_ave_pos->Write();
  for(int i=0;i<numberOfCrystals;i++){
    h_doppler_crystal[i]->Write();
  }   
  h_doppler_max_ring->Write();
  h_energy->Write();
  for(int i=0;i<numberOfCrystals;i++){
    h_energy_crystal[i]->Write();
  }
  h_energy_beta_real->Write();
  h_beta_after_doppler->Write();
  h_crystal_fired_doppler->Write();
}

void GetAverageInteractionPoint(){
  cout<<" Getting the average FI point for all the DALI2 crystals..."<<endl;

  FILE *fOutAve  = fopen("./output/Dali2AverageInteractionPoint.out","w");
  FILE *fOutDiff = fopen("./output/Dali2InteractionDifferenceToCenterOfGravity.out","w");
  float difference[3];
  float distance1, distance2;
  float thetaDifference;
  for( int i = 0;i<numberOfCrystals;i++)  {
    for(int jjj=0;jjj<3;jjj++)  {
      fi_average[jjj][i] = fi_average[jjj][i]/fi_interactions[i];
      difference[jjj] = fi_average[jjj][i] - dali2_pos[jjj][i];
    }
    // The average interaction point:
    fprintf(fOutAve,"%4i %4i %10.2f %10.2f %10.2f \n", i, fi_interactions[i], fi_average[0][i], fi_average[1][i], fi_average[2][i]); 
    cout<<i<<" "<<fi_interactions[i]<<" "<<fi_average[0][i]<<" "<<fi_average[1][i]<<" "<<fi_average[2][i];
    
    distance1 = sqrt(dali2_pos[0][i]*dali2_pos[0][i] + dali2_pos[1][i]*dali2_pos[1][i] + dali2_pos[2][i]*dali2_pos[2][i]);
    distance2 = sqrt(fi_average[0][i]*fi_average[0][i] + fi_average[1][i]*fi_average[1][i] + fi_average[2][i]*fi_average[2][i]);

    thetaDifference = acos(dali2_pos[0][i]*fi_average[0][i] + dali2_pos[1][i]*fi_average[1][i] + dali2_pos[2][i]*fi_average[2][i]
                           /distance1/distance2);

    // Comparison to the center of gravity points:
    fprintf(fOutDiff,"%4i %10.2f %10.2f %10.2f %10.2f\n", i, difference[0],difference[1],difference[2],thetaDifference); 
    cout<<i<<" "<<difference[0]<<" "<<difference[1]<<" "<<difference[2]<<thetaDifference<<endl;
  }
  fclose(fOutAve); 
  fclose(fOutDiff);
}
