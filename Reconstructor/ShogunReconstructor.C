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
#include "TVector3.h" 

using namespace std;
using namespace MyNamespace;

double thetaLast;

//__________________________________________________________________________________________________
bool IncludeAddbackTable(float det[3][2])  {

  float distance = TMath::Sqrt(TMath::Power(det[0][0]-det[0][1],2) +
                               TMath::Power(det[1][0]-det[1][1],2) + 
                               TMath::Power(det[2][0]-det[2][1],2));

  //cout<<"Distance: "<<distance<<endl;
  if( distance > maxAddbackDistance ) return false;
  else return true;
}

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
  thetaLast = theta_lab * 180./3.14159;                        //the last calculated theta value in degrees;
  // the beta value:
  double beta = beta_average + (beta_rec-beta_mean_tof);         //cout<<"beta: "<<beta<<endl;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

//_________________________________________________________________________________________________
float GetDopplerCorrectedEnergySimple(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                               float gamma_energy,float beta,float decay_z){ 

  TVector3 v1(gamma_det_x,gamma_det_y,(gamma_det_z-decay_z));

  double theta_lab = v1.Theta();
  thetaLast = theta_lab * 180./3.14159;                        //the last calculated theta value in degrees;
  
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

void GetAverageInteractionPoint();

const  int maxDecays = 5;
double fi_average[3][NUMBEROFSHOGUNCRYSTALS] = {{0.0}};
int    fi_interactions[NUMBEROFSHOGUNCRYSTALS] = {0};
float  Shogun_pos_ave[3][NUMBEROFSHOGUNCRYSTALS];

//int main(int argc,char** argv){ 
  int ShogunReconstructor () {
  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
  float beta_average     = 0.0;        // The average value at the moment of decay for the correction
  float beta_mean_tof    = 0.0;        // The mean beta from Tof
  float decay_z          = 0.0;        // The lifetime moves the average decay point along the beam axis
  int fi_option          = 0;
 
  float energyThreshold = 0;                // Sets the threshold for observed lab energies

  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/ShogunReconstructor.in","r");
  while(!feof(fin)){
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
    else if(strcmp(temp,"ADDBACK")==0)  {
      fscanf(fin,"%i %f",&shogun_addback_opt,&maxAddbackDistance); 
      printf("%s %i %f\n",temp,shogun_addback_opt,maxAddbackDistance);
    }
    else if(strcmp(temp,"ENERGYTHRESHOLD")==0)  {
      fscanf(fin,"%f",&energyThreshold); 
      printf("%s %f\n",temp,energyThreshold);
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

  TH1F *h_crystal_mult;

  TH1F *h_doppler;
  TH1F *h_doppler_total;

  TH1F *h_energy;
  TH1F *h_energy_total;

  TH1F *h_doppler_total_individual_decays; 

  TH2F *h_crystal_fired_doppler;
  TH2F *h_crystal_fired_doppler_mult1;
  TH2F *h_crystal_fired_doppler_addback;
  
  TH2F *h_crystal_fired_energy;
  TH2F *h_crystal_fired_energy_addback;

  TH2F *h_doppler_gamma_gamma;

  //Creating SHOGUN spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  h_crystal_mult = new TH1F("crystal_Mult","crystal_mult",NUMBEROFSHOGUNCRYSTALS,
                            0,NUMBEROFSHOGUNCRYSTALS);  
 
  h_doppler = new TH1F("doppler","doppler",numBin,firstBin,lastBin);

  // Taking the sum energy and making the Doppler correction with the crystal that has the highest energy:
  h_doppler_total = new TH1F("doppler_total","doppler_total",numBin,firstBin,lastBin);
  h_doppler_total_individual_decays = new TH1F("doppler_total_individual_decays","doppler_total_individual_decays",numBin,firstBin,lastBin);

  h_energy = new TH1F("energy","energy",numBin,firstBin,lastBin);  
  h_energy_total = new TH1F("energy_total","energy_total",numBin,firstBin,lastBin); 
  
  // Comparing crystal fired with doppler corrected energy
  h_crystal_fired_doppler = new TH2F("crystal_fired_doppler","",NUMBEROFSHOGUNCRYSTALS,0,NUMBEROFSHOGUNCRYSTALS,numBin,firstBin,lastBin);
  h_crystal_fired_doppler_mult1 = new TH2F("crystal_fired_doppler_mult1","",NUMBEROFSHOGUNCRYSTALS,0,NUMBEROFSHOGUNCRYSTALS,numBin,firstBin,lastBin);
  h_crystal_fired_energy = new TH2F("crystal_fired_energy","",NUMBEROFSHOGUNCRYSTALS,0,NUMBEROFSHOGUNCRYSTALS,numBin,firstBin,lastBin);
   
  //______________________________________________
  if(fi_option==2)  {
    FILE *fInAve  = fopen("./output/AverageInteractionPoint.out","r");
    int counter = 0;
    float x,y,z;
    int counts;
    while(!feof(fInAve) && counter<NUMBEROFSHOGUNCRYSTALS){
      fscanf(fInAve,"%i %i %f %f %f",&counter,&counts,&x,&y,&z); 
      //counter++;                                            
      Shogun_pos_ave[0][counter] = x;
      Shogun_pos_ave[1][counter] = y;
      Shogun_pos_ave[2][counter] = z;
      cout<<counter<<" "<<Shogun_pos_ave[0][counter]<<" "<<Shogun_pos_ave[1][counter]<<" "<<Shogun_pos_ave[2][counter]<<endl;
    }
  }

  //For gamma-gamma analysis:
  h_doppler_gamma_gamma = new TH2F("doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);
  
  h_crystal_fired_doppler_addback = new TH2F("crystal_fired_doppler_addback","",NUMBEROFSHOGUNCRYSTALS,0,NUMBEROFSHOGUNCRYSTALS,numBin,firstBin,lastBin);
  h_crystal_fired_energy_addback = new TH2F("crystal_fired_energy_addback","",NUMBEROFSHOGUNCRYSTALS,0,NUMBEROFSHOGUNCRYSTALS,numBin,firstBin,lastBin);

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
  //make an addback-table
  //________________________________________________________
  if(shogun_addback_opt==1) {
    FILE *fAddbackTableIn  = fopen("./output/ShogunAddbackTable.out","w");
    float dummy[3][2];
    bool inTable;
    int counter = 0;
    for(int i=0;i<NUMBEROFSHOGUNCRYSTALS;i++) {
      fprintf(fAddbackTableIn," %i",i);
      for(int j=0;j<NUMBEROFSHOGUNCRYSTALS;j++) {
        if(j==i) continue;

        for(int k=0;k<3;k++) {
          if(fi_option==2) {
            dummy[k][0] = Shogun_pos_ave[k][i];
            dummy[k][1] = Shogun_pos_ave[k][j];
          }
          else {
            dummy[k][0] = Shogun_pos[k][i];
            dummy[k][1] = Shogun_pos[k][j];
          }
        }
        inTable = IncludeAddbackTable(dummy);  
        if(inTable && counter< NUMBEROFSHOGUNADDBACKCRYSTALS) {
          fprintf(fAddbackTableIn," %i",j);
          shogunAddbackTable[i][counter] = j;
          counter++;
        }
        if(counter == NUMBEROFSHOGUNADDBACKCRYSTALS)  { //Too many detectors 
          cout<<"You have to increase the variable NUMBEROFSHOGUNADDBACKCRYSTALS!!!"<<endl;
          abort();
        }
      }
      counter = 0;
      fprintf(fAddbackTableIn," \n");
    }
    fclose(fAddbackTableIn);
  }
  

  nentries = (Int_t)t->GetEntries()/reduction_factor;
  cout <<"Entries: "<< nentries << endl;
  
  int crystalMult = 0;

  float energy[NUMBEROFSHOGUNCRYSTALS]                 = {0.0};  // Energy of single crystal. NOT summed over individual decays from a single event
  float time[NUMBEROFSHOGUNCRYSTALS]                   = {0.0};
  float dopplerEnergy[NUMBEROFSHOGUNCRYSTALS]          = {0.0};  // Doppler corrected energy. summed over individual decays.
  float labEnergy[NUMBEROFSHOGUNCRYSTALS]              = {0.0};  // Same as energy[], but summed.
  
  int crystalFired[NUMBEROFSHOGUNCRYSTALS]             = {0};    // The crystal that fired the unsorted energy
  bool crystalUsedForAddback[NUMBEROFSHOGUNCRYSTALS]   = {false};// This flag ist different and gives information if the individual crystals have
  float dopplerEnergyIndividualDecays[maxDecays]       = {0.0};  // Every individual gamma decay branch
  
  float energySum = 0.0;
  float energySumIndividualDecays[maxDecays] = {0.0};
  float energyMax = -999.9;
  float energyMaxIndividualDecays[maxDecays] = {-999.};

  float posThatsIt[3] = {0.0};
  float posThatsItIndividualDecays[maxDecays][3] = {{0.0}};
 
  float AdBkEnergy = 0.0;

  float dummyEnergy = 0;
  int highestCrystal = -1;

  int idThatsIt = -1;
  int eventNumberNext = 0;

  cout<<"Starting the analysis"<<endl;
  for(int iii=0;iii<nentries;iii++)  {
    if(iii%1000==0){
      std::cout << "Event: " << iii <<", "<< (100.*iii/nentries) <<"% of events done\r" <<std::flush;
    }
    if(iii<nentries-1)  {
      t->GetEntry(iii+1);
     
      eventNumberNext = evnum;
    }
    t->GetEntry(iii);
    
    //-------------------------------------------------------------
    // Starting with the gamma ray detector spectra:
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
        if(energy[nnn]>energyThreshold)  {
       	  
          h_crystal_fired_energy->Fill(nnn,energy[nnn]);
   	  h_energy->Fill(energy[nnn]);
          //_____________________________________________________________
          //Performing the Doppler correction,
          float correctedEnergyCrystal;

          if(fi_option !=2){
            // correctedEnergyCrystal = GetDopplerCorrectedEnergy(Shogun_pos[0][nnn],Shogun_pos[1][nnn],Shogun_pos[2][nnn],
            //                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
            //                                                  energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
            correctedEnergyCrystal = GetDopplerCorrectedEnergySimple(Shogun_pos[0][nnn],Shogun_pos[1][nnn],Shogun_pos[2][nnn],
                                                                     energy[nnn],beta_average,decay_z);
          }
     
          else{
            // correctedEnergyCrystal = GetDopplerCorrectedEnergy(Shogun_pos_ave[0][nnn],Shogun_pos_ave[1][nnn],Shogun_pos_ave[2][nnn],
            //                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                               //energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z); 
            correctedEnergyCrystal = GetDopplerCorrectedEnergySimple(Shogun_pos_ave[0][nnn],Shogun_pos_ave[1][nnn],Shogun_pos_ave[2][nnn], 
                                                                     energy[nnn],beta_average,decay_z);                                                   
          }

	  dopplerEnergy[crystalMult] = correctedEnergyCrystal;
          labEnergy[crystalMult] = energy[nnn];

          h_doppler->Fill(correctedEnergyCrystal);

          h_crystal_fired_doppler->Fill(nnn,correctedEnergyCrystal);
          
          crystalFired[crystalMult] = nnn;

          crystalMult++;
	  
	  //See how addback procedures work
	  if(energy[nnn]>energyMax)  {
	    energyMax = energy[nnn];
            if(fi_option==2)  {
              posThatsIt[0]  = Shogun_pos_ave[0][nnn];
              posThatsIt[1]  = Shogun_pos_ave[1][nnn];
              posThatsIt[2]  = Shogun_pos_ave[2][nnn];
            }
            else{
              posThatsIt[0]  = Shogun_pos[0][nnn];
              posThatsIt[1]  = Shogun_pos[1][nnn];
              posThatsIt[2]  = Shogun_pos[2][nnn];
            }
            idThatsIt = nnn;
	  }
          //________________________________________________
          // Only if almost all the energy was released within one crystal the average position is determined:
          if(fi_option == 1 && correctedEnergyCrystal>= e_rest/1.1)  {
            for(int jjj=0;jjj<3;jjj++)  {
              fi_average[jjj][nnn] = fi_average[jjj][nnn] + Shogun_fi[jjj][nnn]; 
              if(jjj==0) fi_interactions[nnn]++;
            }
          }
	}
      }
      h_crystal_mult->Fill(crystalMult);
      if(crystalMult==1) h_crystal_fired_doppler_mult1->Fill(crystalFired[0],dopplerEnergy[0]);
     
      // ___________________________________________________
      // Add-back procedure:
      // It includes also events for which only one detector saw a gammaray!!!
      // Thus -> Addback + singles

      if(crystalMult>=1 && shogun_addback_opt==1)  {
        for(int l =0;l<crystalMult;l++) {
          
          if(crystalUsedForAddback[crystalFired[l]]==true) continue;
          
          float highestEnergy = labEnergy[l];
          dummyEnergy = highestEnergy;
          
          highestCrystal = crystalFired[l];
                   
          crystalUsedForAddback[crystalFired[l]]=true; 

          //Of the remaining crystals, check which one has the highest energy
           for(int i = 0;i<NUMBEROFSHOGUNADDBACKCRYSTALS;i++)  {
             if(crystalUsedForAddback[shogunAddbackTable[crystalFired[l]][i]]==true) continue;
             if(shogunAddbackTable[crystalFired[l]][i]>0)  {
               if(highestEnergy < energy[shogunAddbackTable[crystalFired[l]][i]]) {
                 highestEnergy = energy[shogunAddbackTable[crystalFired[l]][i]];
                 highestCrystal = shogunAddbackTable[crystalFired[l]][i];
               }
               dummyEnergy = dummyEnergy + energy[shogunAddbackTable[crystalFired[l]][i]];
               crystalUsedForAddback[shogunAddbackTable[crystalFired[l]][i]]=true;
             }
           }
           float correctedEnergyCrystal = 0.0;

           if(fi_option !=2){
             //correctedEnergyCrystal = GetDopplerCorrectedEnergy(Shogun_pos[0][highestCrystal],Shogun_pos[1][highestCrystal],Shogun_pos[2][highestCrystal],
             //                                                          ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
             //dummyEnergy,beta_rec,beta_mean_tof,beta_average,decay_z);
             correctedEnergyCrystal = GetDopplerCorrectedEnergySimple(Shogun_pos[0][highestCrystal],Shogun_pos[1][highestCrystal],
                                                                           Shogun_pos[2][highestCrystal], 
                                                                           dummyEnergy,beta_average,decay_z);    
           }
           else {
             //correctedEnergyCrystal = GetDopplerCorrectedEnergy(Shogun_pos_ave[0][highestCrystal],Shogun_pos_ave[1][highestCrystal],Shogun_pos_ave[2][highestCrystal],
             //                                                          ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
             //dummyEnergy,beta_rec,beta_mean_tof,beta_average,decay_z);
             correctedEnergyCrystal = GetDopplerCorrectedEnergySimple(Shogun_pos_ave[0][highestCrystal],Shogun_pos_ave[1][highestCrystal],
                                                                            Shogun_pos_ave[2][highestCrystal], 
                                                                            dummyEnergy,beta_average,decay_z); 
           }
           h_crystal_fired_doppler_addback->Fill(idThatsIt,correctedEnergyCrystal);
           h_crystal_fired_energy_addback->Fill(idThatsIt,dummyEnergy);      
           
           dummyEnergy=0.;
         
          //for(int ppp=0;ppp<11;ppp++)  {
          //  h_doppler_addback_with_threshold[0][ppp]->Fill(dummyEnergyWithThreshold[ppp]);
          //}
        }
      }
      //}
      //____________________________________________________
      // Making a Doppler correction with the Sum of the Energy detected in the array
      // The Position is taken from the crystal that recorded the highest energy
      // float correctedEnergyTotal = GetDopplerCorrectedEnergy(posThatsIt[0],posThatsIt[1],posThatsIt[2],
                                                             //ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                             //energySum,beta_rec,beta_mean_tof,beta_average,decay_z);
                                                             
      float correctedEnergyTotal = GetDopplerCorrectedEnergySimple(posThatsIt[0],posThatsIt[1],posThatsIt[2],
                                                                   energySum,beta_average,decay_z);
                                     
      h_doppler_total->Fill(correctedEnergyTotal);
      h_energy_total->Fill(energySum);

      int counter = 0;

      // Treating the decays individually:
      while(energySumIndividualDecays[counter]>0 && counter<maxDecays)  {
        //        dopplerEnergyIndividualDecays[counter] =  GetDopplerCorrectedEnergy(posThatsItIndividualDecays[counter][0],
        //                                                                    posThatsItIndividualDecays[counter][1],
        //                                                                    posThatsItIndividualDecays[counter][2],
        //                                                                    ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
        //                                                                    energySumIndividualDecays[counter],
        //                                                                   beta_rec,beta_mean_tof,beta_average,decay_z);
        dopplerEnergyIndividualDecays[counter] =  GetDopplerCorrectedEnergySimple(posThatsItIndividualDecays[counter][0],
                                                                            posThatsItIndividualDecays[counter][1],
                                                                            posThatsItIndividualDecays[counter][2],
                                                                            energySumIndividualDecays[counter],beta_average,decay_z);

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
        labEnergy[i]      = 0.;
        crystalFired[i]   = -1;
        crystalUsedForAddback[i] = false;
      }
    }
  }
  //_______________________________________________________________
  // Getting the average point:
  if(fi_option == 1) GetAverageInteractionPoint();

  //---------------------------------------------------------------
  // Writing into the root file:
  cout<<endl;
  cout<<"Writing to the root file..."<<endl;
  rootfile->cd();
    //---------------------------------------------------------------
  // Shogun--------------------------------------------------------
  //---------------------------------------------------------------
  h_crystal_mult->Write();
  
  h_doppler->Write();
  h_doppler_total->Write();
  h_doppler_total_individual_decays->Write();

  h_energy->Write();
  h_energy_total->Write();

  h_crystal_fired_doppler->Write();
  h_crystal_fired_doppler_mult1->Write();
  h_crystal_fired_energy->Write();

  h_doppler_gamma_gamma->Write();
 
  h_crystal_fired_doppler_addback->Write();
  h_crystal_fired_energy_addback->Write();
   
  rootfile->Close();
  infile->Close();
  return 0;
}

void GetAverageInteractionPoint()  {
  cout<<" Getting the average FI point for all the SHOGUN crystals..."<<endl;

  FILE *fOutAve  = fopen("./output/AverageInteractionPoint.out","w");
  FILE *fOutDiff = fopen("./output/InteractionDifferenceToCenterOfGravity.out","w");
  float difference[3];
  float distance1, distance2;
  
  for( int i = 0;i<NUMBEROFSHOGUNCRYSTALS;i++)  {
    for(int jjj=0;jjj<3;jjj++)  {
      fi_average[jjj][i] = fi_average[jjj][i]/fi_interactions[i];
      difference[jjj] = fi_average[jjj][i] - Shogun_pos[jjj][i];
    }
    // The average interaction point:
    fprintf(fOutAve,"%4i %4i %10.2f %10.2f %10.2f \n", i, fi_interactions[i], fi_average[0][i], fi_average[1][i], fi_average[2][i]); 
    cout<<i<<" "<<fi_interactions[i]<<" "<<fi_average[0][i]<<" "<<fi_average[1][i]<<" "<<fi_average[2][i];
    
    distance1 = sqrt(Shogun_pos[0][i]*Shogun_pos[0][i] + Shogun_pos[1][i]*Shogun_pos[1][i] + Shogun_pos[2][i]*Shogun_pos[2][i]);
    distance2 = sqrt(fi_average[0][i]*fi_average[0][i] + fi_average[1][i]*fi_average[1][i] + fi_average[2][i]*fi_average[2][i]);

    // Comparison to the center of gravity points:
    fprintf(fOutDiff,"%4i %10.2f %10.2f %10.2f\n", i, difference[0],difference[1],difference[2]); 
    cout<<i<<" "<<difference[0]<<" "<<difference[1]<<" "<<difference[2]<<endl;

    fi_interactions[i]=0;
  }
  fclose(fOutAve); 
  fclose(fOutDiff);
}
