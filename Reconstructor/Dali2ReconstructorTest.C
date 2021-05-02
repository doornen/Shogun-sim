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
  // the beta value:
  double beta = beta_average + (beta_rec-beta_mean_tof);         //cout<<"beta: "<<beta<<endl;
  //double beta = beta_average;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

void GetAverageInteractionPoint();


const  int maxDecays = 20;
double fi_average[3][NUMBEROFDALI2CRYSTALS] = {{0.0}};
int    fi_interactions[NUMBEROFDALI2CRYSTALS] = {0};
float  dali2_pos_ave[3][NUMBEROFDALI2CRYSTALS];

//int main(int argc,char** argv){ 
int Dali2Reconstructor()  {
  
  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
   
  float beta_average     = 0.0;             // The average value at the moment of decay for the correction
  float beta_mean_tof    = 0.0;             // The mean beta from Tof
  float decay_z          = 0.0;             // The lifetime moves the average decay point along the beam axis
 
  int trigger_opt        = 0;
  int trigger_counter[11][3]= {{0}};        // Counts how many events were seen above a certain thresshold (i*50 keV)
  bool trigger_flag[11][3]  = {{false}};    // To prevent double triggering 0=all,1=backward,2=forward angles.


  float energyThreshold = 0;                // Sets the threshold for DALI2 energies
  
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
      if(reduction_factor<1) return 0;
    }
    else if(strcmp(temp,"FIFIND")==0)  {
      fscanf(fin,"%i",&dali2_fi_opt); 
      printf("%s %i \n",temp,dali2_fi_opt);
    }
    else if(strcmp(temp,"ADDBACK")==0)  {
      fscanf(fin,"%i %f",&dali2_addback_opt,&maxAddbackDistance); 
      printf("%s %i %f\n",temp,dali2_addback_opt,maxAddbackDistance);
    }
    else if(strcmp(temp,"TRIGGER")==0)  {
      fscanf(fin,"%i",&trigger_opt); 
      printf("%s %i\n",temp,trigger_opt);
    }
    else if(strcmp(temp,"ENERGYTHRESHOLD")==0)  {
      fscanf(fin,"%f",&energyThreshold); 
      printf("%s %f\n",temp,energyThreshold);
    }
    //else if(strcmp(temp,"DETECTORLIMIT")==0)  {
    //  fscanf(fin,"%i",&detectorLimit); 
    //  printf("%s %i\n",temp,detectorLimit);
    //}
    else if(strcmp(temp,"END")==0) break;
    else {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      return 0;
    }
  }
  //--------------------------------------------------------------------------------------------------
  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();

  TH1F *h_trigger_prob[3];
  
  TH1F *h_crystal_mult;

  TH1F *h_doppler;

  TH1F *h_doppler_total;
  TH1F *h_doppler_total_individual_decays; 
    
  TH2F *h_crystal_fired_doppler;
  TH2F *h_crystal_fired_doppler_addback;//Include also add-back crystals
  TH2F *h_crystal_fired_energy;

  TH2F *h_doppler_gamma_gamma;
  
  // Spectra of the trigger probablility, also separated into backward and forward angles:
  if(trigger_opt == 1)  {
    for(int i=0;i<3;i++)  {
      sprintf(temp,"h_trigger_prob[%i]",i);
      h_trigger_prob[i] = new TH1F(temp,temp,501,-0.5,500.5);
    }
  }

  //Creating dali spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Crystal multiplicity,ok
  h_crystal_mult = new TH1F("crystal_Mult","crystal_mult",NUMBEROFDALI2CRYSTALS,0,NUMBEROFDALI2CRYSTALS);  
  // All, OK
  h_doppler = new TH1F("doppler","doppler",numBin,firstBin,lastBin);
    
  // Taking the sum energy and making the Doppler correction with the crystal that has the highest energy:
  h_doppler_total = new TH1F("doppler_total","doppler_total",numBin,firstBin,lastBin);
  h_doppler_total_individual_decays = new TH1F("doppler_total_individual_decays","doppler_total_individual_decays",numBin,firstBin,lastBin);
  
  h_crystal_fired_doppler = new TH2F("crystal_fired_doppler","",NUMBEROFDALI2CRYSTALS,0,
                                     NUMBEROFDALI2CRYSTALS,numBin,firstBin,lastBin);
  h_crystal_fired_energy = new TH2F("crystal_fired_energy","",NUMBEROFDALI2CRYSTALS,0,
                                    NUMBEROFDALI2CRYSTALS,numBin,firstBin,lastBin);
  //______________________________________________
  if(dali2_fi_opt==2)  {
    FILE *fInAve  = fopen("./output/AverageInteractionPoint.out","r");
    int counter = 0;
    float x,y,z;
    int counts;
    int dummy;
    float dummy2;
    while(!feof(fInAve) && counter<NUMBEROFDALI2CRYSTALS){
      fscanf(fInAve,"%f %f %f %i %i %f %i %i",&x,&y,&z,&counter,&counts,&dummy2,&dummy,&dummy); 
      dali2_pos_ave[0][counter] = x;
      dali2_pos_ave[1][counter] = y;
      dali2_pos_ave[2][counter] = z;
      cout<<counter<<" "<<dali2_pos_ave[0][counter]<<" "<<dali2_pos_ave[1][counter]<<" "<<dali2_pos_ave[2][counter]<<endl;
      //counter++;                                            
    }
  }

  // For gamma-gamma analysis:
  h_doppler_gamma_gamma = new TH2F("doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);

  // Testing the addback procedure with cluster grouping
  if(dali2_addback_opt==1)  {
    h_crystal_fired_doppler_addback = new TH2F("crystal_fired_doppler_addback","",NUMBEROFDALI2CRYSTALS,0,
                                               NUMBEROFDALI2CRYSTALS,numBin,firstBin,lastBin);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  infile = new TFile(root_in,"READ");
  t = (TTree*)infile->Get("ObservedEvents");
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
  t->SetBranchAddress("EnergyVertex",&energy_vertex_stored);                    // Energy of fragment at fragmentation
  t->SetBranchAddress("BGGamma",&bgGamma);                                      // To know if the gamma is from bg or a good event
  //So far identical to the ROOT-Tree of step one
  //The new stuff for step 2:
  //Setting the type of the gamma detector
  t->SetBranchAddress("GammaDetType",&gamma_det_type);
  //-------------------------------
  //DALI:
  t->SetBranchAddress("DALI2Flag",dali2_flag);                                  // ID of which det recorded a gamma in the event
  t->SetBranchAddress("DALI2EnergyNotCor",dali2_energy_not_cor);
  t->SetBranchAddress("DALI2Time",dali2_time);
  if(dali2_fi_opt >= 1)  {
    t->SetBranchAddress("DALI2FI",dali2_fi);
  }
  //-------------------------------
  //Position detector after targetpos_det_after_target_res
  t->SetBranchAddress("PosDet",pos_det);                                        // Reconstructed position after the secondary target
  t->SetBranchAddress("VertexReconstructed",ver_rec);                           // Reconstructed vertex postion, 
 
  //reconstructed beta 
  t->SetBranchAddress("BetaReconstructed",&beta_rec);                           // Reconstructed beta, including resolution

  //------------------------------------------------------------------------------------------------------
  //Going to the Header file, which carries the information of constant values.
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

  // Before I start the analysis, I make an addback-table:
  //________________________________________________________
  if(dali2_addback_opt==1) {
    FILE *fAddbackTableIn  = fopen("./output/AddbackTable.out","w");
    float dummy[3][2];
    bool inTable;
    int counter = 0;
    for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++) {
      fprintf(fAddbackTableIn," %i",i);
      for(int j=0;j<NUMBEROFDALI2CRYSTALS;j++) {
        if(j==i) continue;
        for(int k=0;k<3;k++) {
          if(dali2_fi_opt==2) {
            dummy[k][0] = dali2_pos_ave[k][i];
            dummy[k][1] = dali2_pos_ave[k][j];
          }
          else {
            dummy[k][0] = dali2_pos[k][i];
            dummy[k][1] = dali2_pos[k][j];
          }
        }
        inTable = IncludeAddbackTable(dummy);  
        if(inTable && counter< NUMBEROFDALI2ADDBACKCRYSTALS) {
          fprintf(fAddbackTableIn," %i",j);
          addbackTable[i][counter] = j;
          counter++;
        }
        if(counter == NUMBEROFDALI2ADDBACKCRYSTALS)  { //Too many detectors 
          cout<<"You have to increase the variable NUMBEROFDALI2ADDBACKCRYSTALS!!!"<<endl;
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
  
  float energy[NUMBEROFDALI2CRYSTALS]                 = {0.0};  // Unsorted energy, filled according to the loop from 0 to the last crystal
  float time[NUMBEROFDALI2CRYSTALS]                   = {0.0};  // same for the time
  float dopplerEnergy[NUMBEROFDALI2CRYSTALS]          = {0.0};  // Unsorted dopplerEnergy
  float labEnergy[NUMBEROFDALI2CRYSTALS]              = {0.0};  // Unsorted labEnergy

  int crystalFired[NUMBEROFDALI2CRYSTALS]             = {0};    // The crystal that fired the unsorted energy
  bool crystalUsedForAddback[NUMBEROFDALI2CRYSTALS]   = {false};// This flag ist different and gives information if the individual crystals have
  // been used for addback already
  float dopplerEnergyIndividualDecays[maxDecays]      = {0.0};  // Every individual gamma decay branch
 
  float energySum = 0.0;
  float energySumIndividualDecays[maxDecays] = {0.0};
  float energyMax = -999.9;
  float energyMaxIndividualDecays[maxDecays] = {-999.};

  float posThatsIt[3] = {0.0};
  float posThatsItIndividualDecays[maxDecays][3] = {{0.0}};
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
    for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++)  {
      if(dali2_flag[nnn]==1)   {
        energy[nnn] += dali2_energy_not_cor[nnn];
        energySum   +=  dali2_energy_not_cor[nnn];
        if(maxDecays>decayIDevnum) {
          energySumIndividualDecays[(int)decayIDevnum] = energySumIndividualDecays[(int)decayIDevnum] + dali2_energy_not_cor[nnn];
          // Getting the right position for the Doppler correction of individual decays:
          if(dali2_energy_not_cor[nnn]>energyMaxIndividualDecays[(int)decayIDevnum])  {
	    energyMaxIndividualDecays[(int)decayIDevnum] = dali2_energy_not_cor[nnn];
	    posThatsItIndividualDecays[(int)decayIDevnum][0]  = dali2_pos[0][nnn];
            posThatsItIndividualDecays[(int)decayIDevnum][1]  = dali2_pos[1][nnn];
            posThatsItIndividualDecays[(int)decayIDevnum][2]  = dali2_pos[2][nnn];
	  }
        }
        //Getting the trigger counter:
        if(trigger_opt==1)
          for(int j=0;j<11;j++)  {
            if(energy[nnn]>(j*50))  {
              for(int k=0;k<3;k++)
                if(trigger_flag[j][k]==false){
                  if(k==0) {
                    trigger_counter[j][k]++;
                    trigger_flag[j][k]=true;
                  }
                  if(k==1&&dali2_pos[2][nnn]<0)  {
                    trigger_counter[j][k]++;
                    trigger_flag[j][k]=true;
                  }
                  if(k==2&&dali2_pos[2][nnn]>=0)  {
                    trigger_counter[j][k]++;
                    trigger_flag[j][k]=true;
                  }
                }
            }
          }   
      }
    }
    if(eventNumberNext!=evnum)  {                 // The last observed decay we can start to fill the spectra
      for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++)  {
        if(energy[nnn]>energyThreshold)  {
       	  // Filling the energy spectra without Doppler-correction        
          h_crystal_fired_energy->Fill(nnn,energy[nnn]);
          
          //_____________________________________________________________
          //Performing the Doppler correction,
	  //This one is detectorwise!
          float correctedEnergyCrystal;
          
          //dali2_fi_opt==2 I use the average FI point from the simulation
          if(dali2_fi_opt !=2)
            correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                               ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                               energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
          else
            correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos_ave[0][nnn],dali2_pos_ave[1][nnn],dali2_pos_ave[2][nnn],
                                                               ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                               energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z); 
  
          
	  dopplerEnergy[crystalMult] = correctedEnergyCrystal;
          labEnergy[crystalMult] = energy[nnn];

          h_doppler->Fill(correctedEnergyCrystal);
          
          h_crystal_fired_doppler->Fill(nnn,correctedEnergyCrystal);
          
          crystalFired[crystalMult] = nnn;
          
          //if(energy[nnn]>150)
          crystalMult++;
	  
	  //See how addback procedures work
	  if(energy[nnn]>energyMax)  {
	    energyMax = energy[nnn];
            if(dali2_fi_opt==2)  {
              posThatsIt[0]  = dali2_pos_ave[0][nnn];
              posThatsIt[1]  = dali2_pos_ave[1][nnn];
              posThatsIt[2]  = dali2_pos_ave[2][nnn];
            }
            else {
              posThatsIt[0]  = dali2_pos[0][nnn];
              posThatsIt[1]  = dali2_pos[1][nnn];
              posThatsIt[2]  = dali2_pos[2][nnn];
            }

            idThatsIt      = nnn;
	  }
          //________________________________________________
          // Only if almost all the energy was released within one crystal the average position is determined:
          if(dali2_fi_opt == 1 && correctedEnergyCrystal>= e_rest/1.1)  {
            for(int jjj=0;jjj<3;jjj++)  {
              fi_average[jjj][nnn] = fi_average[jjj][nnn] + dali2_fi[jjj][nnn]; 
              if(jjj==0) fi_interactions[nnn]++;
            }
          }
	}
      }
      h_crystal_mult->Fill(crystalMult);
      //SortEnergies();           //Sort the observed energies from highest to lowest!!!!

      // ___________________________________________________
      // Add-back procedure:
      // It includes also events for which only one detector saw a gammaray!!!
      // Thus -> Addback + singles
      if(crystalMult>=1 && dali2_addback_opt==1)  {
        for(int l =0;l<crystalMult;l++) {
          float highestEnergy = 0.;
          int highestCrystal = -1;
          for(int i = 0;i<crystalMult;i++)  {  
            if(crystalUsedForAddback[crystalFired[i]]==true) continue;
            if(highestEnergy < labEnergy[i]) {
              highestEnergy = labEnergy[i];
              highestCrystal = crystalFired[i];
            }
          }

          if(highestCrystal<0) break;

          float dummyEnergy = highestEnergy;
          crystalUsedForAddback[highestCrystal]=true; 
           
          for(int j = 0;j<crystalMult;j++)  {
            if(crystalUsedForAddback[crystalFired[j]]==false && dopplerEnergy[j]>0)  {
              for(int k = 0;k<NUMBEROFDALI2ADDBACKCRYSTALS;k++) {
                if(addbackTable[highestCrystal][k]<0) continue;
                
                if(crystalFired[j] == addbackTable[highestCrystal][k])  {
                  crystalUsedForAddback[crystalFired[j]]=true;
                  dummyEnergy = dummyEnergy + labEnergy[j];
                }
              }
            }
          }

          float correctedEnergyCrystal = GetDopplerCorrectedEnergy(posThatsIt[0],posThatsIt[1],posThatsIt[2],
                                                                   ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                                   dummyEnergy,beta_rec,beta_mean_tof,beta_average,decay_z);
          
          h_crystal_fired_doppler_addback->Fill(idThatsIt,correctedEnergyCrystal);
          //for(int ppp=0;ppp<11;ppp++)  {
          //  h_doppler_addback_with_threshold[0][ppp]->Fill(dummyEnergyWithThreshold[ppp]);
          //}
        }
      }
      //____________________________________________________
      // Making a Doppler correction with the Sum of the Energy detected in the array
      // The Position is taken from the crystal that recorded the highest energy
      float correctedEnergyTotal = GetDopplerCorrectedEnergy(posThatsIt[0],posThatsIt[1],posThatsIt[2],
                                                             ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                             energySum,beta_rec,beta_mean_tof,beta_average,decay_z);
      h_doppler_total->Fill(correctedEnergyTotal);
      int counter = 0;

      //____________________________________________________
      // Treating decays individually:
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
      if(crystalMult>=2)  {
        for(int i = 0;i<crystalMult;i++)
          for(int j = i+1;j<crystalMult;j++){
            if (dopplerEnergy[i]>=dopplerEnergy[j]) 
              h_doppler_gamma_gamma->Fill(dopplerEnergy[i],dopplerEnergy[j]);
            else h_doppler_gamma_gamma->Fill(dopplerEnergy[j],dopplerEnergy[i]);
          }
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
      for(int i =0;i<NUMBEROFDALI2CRYSTALS;i++){
        energy[i]         = 0.;
        time[i]           = 0.;
        dopplerEnergy[i]  = 0.;
        labEnergy[i]      = 0.;
        crystalFired[i]   = -1;
        crystalUsedForAddback[i] = false;
      }
      //Getting the trigger counter:
      if(trigger_opt==1)
        for(int j=0;j<11;j++)  
          for(int k=0;k<3;k++)
            trigger_flag[j][k]=false;
    }
  } // End of loop through events
  //_______________________________________________________________
  // Getting the average first interaction point:
  if(dali2_fi_opt == 1) GetAverageInteractionPoint();

  //---------------------------------------------------------------
  cout<<endl;
  cout<<"Writing to the root file..."<<endl;
  rootfile->cd();
  //_______________________________________________________________
  // Filling the trigger probability spectra:
  if(trigger_opt == 1)  {
    for(int k=0;k<3;k++)  {
      for(int j=0;j<11;j++)  {
        h_trigger_prob[k]->SetBinContent((j*50+k+1),(trigger_counter[j][k]/evnum));
      }
      h_trigger_prob[k]->Write();
    }
  }
  //---------------------------------------------------------------
  // DALI----------------------------------------------------------
  //---------------------------------------------------------------
  h_crystal_mult->Write();
  h_doppler->Write();
  
  h_doppler_total->Write();
  h_doppler_total_individual_decays->Write();
  
  h_crystal_fired_doppler->Write();
  h_crystal_fired_energy->Write();
  
  h_doppler_gamma_gamma->Write();
  if(dali2_addback_opt==1){
    h_crystal_fired_doppler_addback->Write();
  }
  rootfile->Close();
  infile->Close();
  return 0;
}

void GetAverageInteractionPoint()  {
  cout<<" Getting the average FI point for all the DALI2 crystals..."<<endl;

  FILE *fOutAve  = fopen("./output/AverageInteractionPoint.out","w");
  FILE *fOutDiff = fopen("./output/InteractionDifferenceToCenterOfGravity.out","w");
  float difference[3];
  float distance1, distance2;
  
  for( int i = 0;i<NUMBEROFDALI2CRYSTALS;i++)  {
    for(int jjj=0;jjj<3;jjj++)  {
      fi_average[jjj][i] = fi_average[jjj][i]/fi_interactions[i];
      difference[jjj] = fi_average[jjj][i] - dali2_pos[jjj][i];
    }

    distance1 = sqrt(dali2_pos[0][i]*dali2_pos[0][i] + dali2_pos[1][i]*dali2_pos[1][i] + dali2_pos[2][i]*dali2_pos[2][i]);
    distance2 = sqrt(fi_average[0][i]*fi_average[0][i] + fi_average[1][i]*fi_average[1][i] + fi_average[2][i]*fi_average[2][i]);
  

    TVector3 thetaangle(fi_average[0][i],fi_average[1][i],fi_average[2][i]);
    TVector3 thetaanglecg(dali2_pos[0][i],dali2_pos[1][i],dali2_pos[2][i]);
    // The average interaction point:
    fprintf(fOutAve,"%10.2f %10.2f %10.2f %4i %4i %10.2f 0 0\n",fi_average[0][i], fi_average[1][i], fi_average[2][i],i, fi_interactions[i], 180*thetaangle.Theta()/3.14159); 
    cout<<i<<" "<<fi_interactions[i]<<" "<<fi_average[0][i]<<" "<<fi_average[1][i]<<" "<<fi_average[2][i];
    
    // Comparison to the center of gravity points:
    fprintf(fOutDiff,"%4i %10.2f %10.2f %10.2f %10.3f\n", i, difference[0],difference[1],difference[2],180*thetaangle.Angle(thetaanglecg)/3.14159); 
    cout<<" "<<i<<" "<<difference[0]<<" "<<difference[1]<<" "<<difference[2]<<endl;
  
    fi_interactions[i]=0;
  }
  fclose(fOutAve); 
  fclose(fOutDiff);
}
