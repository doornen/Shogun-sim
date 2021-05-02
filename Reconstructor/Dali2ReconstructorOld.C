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
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

void GetAverageInteractionPoint();

 
const  int maxDecays = 20;
double fi_average[3][NUMBEROFDALI2CRYSTALS] = {{0.0}};
int    fi_interactions[NUMBEROFDALI2CRYSTALS] = {0};
float  dali2_pos_ave[3][NUMBEROFDALI2CRYSTALS];

//int main(int argc,char** argv){ 
int Dali2ReconstructorOld()  {
  
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
  int detectorLimit = 40;                   // From which crystal ID detectors are defined as forward  
 
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
    else if(strcmp(temp,"DETECTORLIMIT")==0)  {
      fscanf(fin,"%i",&detectorLimit); 
      printf("%s %i\n",temp,detectorLimit);
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

  //Testing the position
  //TH2F *h_gamma_vector_xy = new TH2F("gamma_vector_xy","",200,-1.1,1.1,200,-1.1,1.1);

  //Creating Beta and gamma lab spectra!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //TH2F *h_beta_before_beta_after = new TH2F("beta_before_beta_after","beta_before_beta_after",200,0.50,0.65,200,0.50,0.65);
  //TH2F *h_beta_real_beta_before = new TH2F("beta_real_beta_before","beta_real_beta_before",200,0.50,0.65,200,0.50,0.65);
  //TH2F *h_beta_real_beta_after = new TH2F("beta_real_beta_after","beta_real_beta_after",200,0.50,0.65,200,0.50,0.65);

  TH1F *h_trigger_prob[3];
  
  TH1F *h_crystal_mult;

  TH1F *h_doppler;

  TH1F *h_doppler_forward_angles;
  TH1F *h_doppler_backward_angles;

  TH1F *h_doppler_simple;

  TH1F *h_doppler_mult[20];

  TH1F *h_doppler_mult1and2;
  //TH1F *h_doppler_crystal[NUMBEROFDALI2CRYSTALS];
  TH1F *h_doppler_layer[17];

  TH1F *h_doppler_total;
  TH1F *h_doppler_total_individual_decays; 
  //TH2F *h_doppler_max_ring;

  TH1F *h_energy;
  TH1F *h_energy_forward_angles;
  TH1F *h_energy_backward_angles;
  //TH1F *h_energy_crystal[NUMBEROFDALI2CRYSTALS];
  //TH2F *h_energy_beta_real;

  //TH2F *h_beta_after_doppler;
  TH2F *h_crystal_fired_doppler;
  TH2F *h_crystal_fired_doppler_addback;//Include also add-back crystals
  TH2F *h_crystal_fired_energy;

  TH1F *h_doppler_ave_pos;
  
  TH2F *h_doppler_gamma_gamma;

  TH1F *h_doppler_addback[3];  // 0: all, 1: backward,2: forward angles
  //TH1F *h_doppler_addback_with_thresold[3][11];  // Threshold from 0 to 500 keV in steps of 50 keV

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
  h_doppler_forward_angles = new TH1F("doppler_forward_angles","doppler_forward_angles",numBin,firstBin,lastBin);
  h_doppler_backward_angles = new TH1F("doppler_backward_angles","doppler_backward_angles",numBin,firstBin,lastBin);
  h_doppler_simple = new TH1F("doppler_simple","doppler_simple",numBin,firstBin,lastBin);

  h_doppler_mult1and2 = new TH1F("doppler_mult1and2","doppler_mult1and2",numBin,firstBin,lastBin);
  
  // According to the multiplicity
  for(int i=0;i<20;i++)  {
    sprintf(temp,"doppler_mult[%i]",i);
    h_doppler_mult[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  } 

  // The spectra crystalwise, OK
  //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
  //  sprintf(temp,"doppler_crystal[%i]",i);
  //  h_doppler_crystal[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  //}  
  // By layer:
  for(int i=0;i<17;i++)  {
    sprintf(temp,"doppler_layer[%i]",i);
    h_doppler_layer[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  }
  // Taking the sum energy and making the Doppler correction with the crystal that has the highest energy:
  h_doppler_total = new TH1F("doppler_total","doppler_total",numBin,firstBin,lastBin);
  h_doppler_total_individual_decays = new TH1F("doppler_total_individual_decays","doppler_total_individual_decays",numBin,firstBin,lastBin);
  // The ring with the most energy:
  //for(int i=0;i<15;i++)
  //{
  //sprintf(temp,"doppler_max_ring");
  //h_doppler_max_ring = new TH2F(temp,temp,numBin,firstBin,lastBin,15,0,15);
  //}

  // Not doppler corrected, OK
  h_energy = new TH1F("energy","energy",numBin,firstBin,lastBin);
  h_energy_forward_angles = new TH1F("energy_forward_angles","energy_forward_angles",numBin,firstBin,lastBin);
  h_energy_backward_angles = new TH1F("energy_backward_angles","energy_backward_angles",numBin,firstBin,lastBin);
  // The spectra crystalwise,OK
  //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
  //  sprintf(temp,"energy_crystal[%i]",i);
  //  h_energy_crystal[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  //}
  // Comparing beta REAL with dali energy:
  //h_energy_beta_real = new TH2F("energy_beta_real","energy_beta_real",numBin,firstBin,lastBin,200,0.40,0.45);

  // Comparing beta after the target with doppler corrected energy
  //h_beta_after_doppler = new TH2F("beta_after_doppler","",200,0.5,0.64,numBin,firstBin,lastBin);

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
    h_doppler_ave_pos = new TH1F("doppler_ave_pos","doppler_ave_pos",numBin,firstBin,lastBin);
  }

  // For gamma-gamma analysis:
  h_doppler_gamma_gamma = new TH2F("doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);

  // Testing the addback procedure with cluster grouping
  if(dali2_addback_opt==1)  {
    for(int i=0;i<3;i++)  {
      sprintf(temp,"h_doppler_addback[%i]",i);
      h_doppler_addback[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
      //for(int j=0;j<11;j++)  {
      //  sprintf(temp,"h_doppler_addback_with_threshold[%i][%i]",i,j);
      //  h_doppler_addback_with_threshold[i][j] = new TH1F(temp,temp,numBin,firstBin,lastBin);
      //}
    }
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
  if(dali2_addback_opt==1)  {
    FILE *fAddbackTableIn  = fopen("./output/AddbackTable.out","w");
    float dummy[3][2];
    bool inTable;
    int counter = 0;
    for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
      fprintf(fAddbackTableIn," %i",i);
      for(int j=i+1;j<NUMBEROFDALI2CRYSTALS;j++)  {
        for(int k=0;k<3;k++) {
          dummy[k][0] = dali2_pos[k][i];
          dummy[k][1] = dali2_pos[k][j];
          //cout<<"dali2_pos: "<<dali2_pos[k][j]<<endl;
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
  //Dividing old dali2 into 16 rings of equal z-position:
  //int dali2_layer_low[16]={0,6,12,20,32,45,56,70,80,90,104,116,128,140,148,154};
  //int dali2_layer_high[16]={5,11,19,31,44,55,69,79,89,103,115,127,139,147,153,159};
  //Dividing new dali2. into 17 rings of equal theta angle:
  //int ringLow[17] = {0, 6,13,20,30,40,52,66,80,94,108,122,136,148,156,164,170};
  //int ringHigh[17]= {5,12,19,29,39,51,65,79,93,107,121,135,147,155,163,169,181};
  // Dali2.2 "pineapple" which includes the forward wall
  int layerLow[17] = {0, 6,12,20,30,40,52,66,80, 94,108,122,199,199,199,199,199};
  int layerHigh[17]= {5,11,19,29,39,51,65,79,93,107,121,185,199,199,199,199,199};
  
  int ringEnergyMax = -1;

  float energy[NUMBEROFDALI2CRYSTALS]                 = {0.0};  // Unsorted energy, filled according to the loop from 0 to the last crystal
  float time[NUMBEROFDALI2CRYSTALS]                   = {0.0};  // same for the time
  float dopplerEnergy[NUMBEROFDALI2CRYSTALS]          = {0.0};  // Unsorted dopplerEnergy
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
    //if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    if(iii%1000==0){
      std::cout << "Event: " << iii <<", "<< (100.*iii/nentries) <<"% of events done\r" <<std::flush;
    }
    if(iii<nentries-1)  {
      t->GetEntry(iii+1);
      
      eventNumberNext = evnum;
    }
    t->GetEntry(iii);
    //Fill the gamma xy spectrum
    //h_gamma_vector_xy->Fill(gv[0],gv[1]);
    
    // Fill the beta spectra:
    //h_beta_before_beta_after->Fill(beta_b,beta_a);
    //h_beta_real_beta_before->Fill(beta_r,beta_b);
    //h_beta_real_beta_after->Fill(beta_r,beta_a);
    
    //-------------------------------------------------------------
    // Starting with the gamma ray detector spectra:
    for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++)  {
      //removing crystals that have a bad energy resolution.
      //NP0811-RIBF70
      //       if(nnn==0  || nnn==1  || nnn==2  || nnn==3 || nnn==5  || nnn==6  || nnn==39 || nnn==63 || 
      //     nnn==123 || nnn==124|| nnn==125||nnn==126|| nnn==127|| nnn==130|| nnn==132) continue;
      //NP0702-RIBF32exp2010
      //if(nnn==0  || nnn==9  || nnn==132 || nnn==134) continue;
      //Efficiency for Upgrade:
      //if(nnn==0  || nnn==1  || nnn==2  || nnn==3 || nnn==5  || nnn==6  || nnn==39 || nnn==63 || nnn==65 || 
      //nnn==123 || nnn==124|| nnn==125||nnn==126|| nnn==127||  nnn==130|| nnn==132) continue;
      if(dali2_flag[nnn]==1)   {
        energy[nnn] += dali2_energy_not_cor[nnn];
        //if(nnn==180) cout<<"Energy : "<<energy[nnn]<<endl;
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
        //if(nnn==23 || nnn==69) continue;  //removing crystals that have a bad energy resolutoin.
        if(energy[nnn]>energyThreshold)  {
       	  // Filling the energy spectra without Doppler-correction        
	  h_energy->Fill(energy[nnn]);
            
          if(nnn>=detectorLimit)  
            h_energy_forward_angles->Fill(energy[nnn]);
          else h_energy_backward_angles->Fill(energy[nnn]);
	  //h_energy_crystal[nnn]->Fill(energy[nnn]);
	  //h_energy_beta_real->Fill(energy[nnn],beta_r);
          h_crystal_fired_energy->Fill(nnn,energy[nnn]);
          
          //_____________________________________________________________
          //Performing the Doppler correction,
	  //This one is detectorwise!
          float correctedEnergyCrystal;
          float correctedEnergySimple;
          //dali2_fi_opt==2 I use the average FI point from the simulation
          if(dali2_fi_opt !=2)
            correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                               ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                               energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
          else
            correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos_ave[0][nnn],dali2_pos_ave[1][nnn],dali2_pos_ave[2][nnn],
                                                               ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                               energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z); 
  
          if(dali2_fi_opt !=2)
            correctedEnergySimple = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                              0,0,0,0,0,10000, 
                                                              energy[nnn],0,0,beta_average,decay_z);
          else
            correctedEnergySimple = GetDopplerCorrectedEnergy(dali2_pos_ave[0][nnn],dali2_pos_ave[1][nnn],dali2_pos_ave[2][nnn],
                                                              0,0,0,0,0,10000, 
                                                              energy[nnn],0,0,beta_average,0);
          
	  dopplerEnergy[crystalMult] = correctedEnergyCrystal;
          h_doppler->Fill(correctedEnergyCrystal);
          // Checking the difference between forward and backward angles.
          if(nnn>=detectorLimit) 
            h_doppler_forward_angles->Fill(correctedEnergyCrystal);
          else h_doppler_backward_angles->Fill(correctedEnergyCrystal);
          
          h_doppler_simple->Fill(correctedEnergySimple);

          //h_doppler_crystal[nnn]->Fill(correctedEnergyCrystal);
          h_crystal_fired_doppler->Fill(nnn,correctedEnergyCrystal);
          
	  dopplerEnergy[crystalMult] = correctedEnergyCrystal;
          crystalFired[crystalMult] = nnn;
          //Set the threshold
          //if(energy[nnn]>150)
          crystalMult++;
	  //Filling Dali layer wise:
	  for(int ppp=0;ppp<17;ppp++)  {
	    if(nnn>=layerLow[ppp] && nnn<=layerHigh[ppp]) h_doppler_layer[ppp]->Fill(correctedEnergyCrystal);
	  }
	  //See how addback procedures work
	  if(energy[nnn]>energyMax)  {
	    energyMax = energy[nnn];
	    posThatsIt[0]  = dali2_pos[0][nnn];
            posThatsIt[1]  = dali2_pos[1][nnn];
            posThatsIt[2]  = dali2_pos[2][nnn];
            idThatsIt      = nnn;
	    for(int ppp=0;ppp<17;ppp++) {if(nnn>=layerLow[ppp] && nnn<=layerHigh[ppp]) ringEnergyMax = ppp;}
	  }
          //________________________________________________
          // Only if almost all the energy was released within one crystal the average position is determined:
          if(dali2_fi_opt == 1 && correctedEnergyCrystal>= e_rest/1.1)  {
            for(int jjj=0;jjj<3;jjj++)  {
              fi_average[jjj][nnn] = fi_average[jjj][nnn] + dali2_fi[jjj][nnn]; 
              if(jjj==0) fi_interactions[nnn]++;
            }
          }
          //_________________________________________________
          // Checking the effect of getting the average position:
          if(dali2_fi_opt == 2)  {
            h_doppler_ave_pos->Fill(correctedEnergyCrystal);
          }
	}
      }
      h_crystal_mult->Fill(crystalMult);
      //SortEnergies();           //Sort the observed energies from highest to lowest!!!!

      if(crystalMult>0)  {
        for(int i = 1;i<20;i++)  {
          if(crystalMult==i)  {	 
            for(int j = 0;j<i;j++)  {
              h_doppler_mult[i]->Fill(dopplerEnergy[j]);
              if(crystalMult==1 || crystalMult==2)h_doppler_mult1and2->Fill(dopplerEnergy[j]);
            }
          }
        }
      }
      // ___________________________________________________
      // Checking realistic add-back procedure:
      // It includes also events for which only one detector saw a gammaray!!!
      // Thus -> Addback + singles
      if(crystalMult>=1 && dali2_addback_opt==1)  {
        for(int i = 0;i<crystalMult;i++)  {  
          if(crystalUsedForAddback[crystalFired[i]]==true) continue;
          float dummyEnergy = dopplerEnergy[i];
          //float dummyEnergyWithThreshold[11];
          //for(int ppp=0;ppp<11;ppp++)  {
          //  dummyEnergyWithThreshold[ppp]= dopplerEnergy[i];
          //}

          crystalUsedForAddback[crystalFired[i]]=true; 
           
          for(int j = i;j<crystalMult;j++)  {
            if(crystalUsedForAddback[crystalFired[j]]==false&& dopplerEnergy[j]>0)  {
              for(int k = 0;k<NUMBEROFDALI2ADDBACKCRYSTALS;k++) {
                if(crystalFired[j] == addbackTable[crystalFired[i]][k])  {
          
                  crystalUsedForAddback[crystalFired[j]]=true;
                  dummyEnergy = dummyEnergy + dopplerEnergy[j];
                  //for(int ppp=0;ppp<11;ppp++)  {
                  //  if(dopplerEnergy[i]>=ppp*50&& dopplerEnergy[j]>=ppp*50) 
                  //    dummyEnergyWithThreshold[ppp] = dummyEnergyWithThreshold[ppp] +  dopplerEnergy[j];
                  //}
                }
              }
            }
          }
          h_doppler_addback[0]->Fill(dummyEnergy);
          h_crystal_fired_doppler_addback->Fill(idThatsIt,dummyEnergy);
          //for(int ppp=0;ppp<11;ppp++)  {
          //  h_doppler_addback_with_threshold[0][ppp]->Fill(dummyEnergyWithThreshold[ppp]);
          //}
          if(posThatsIt[2]<0) {
            h_doppler_addback[1]->Fill(dummyEnergy);
          }
          else {
            h_doppler_addback[2]->Fill(dummyEnergy);
          }   
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
      ringEnergyMax = -1;
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
  // BETA, GAMMA
  //h_gamma_vector_xy->Write();
  //h_beta_before_beta_after->Write();
  //h_beta_real_beta_before->Write();
  //h_beta_real_beta_after->Write();
  //---------------------------------------------------------------
  // DALI----------------------------------------------------------
  //---------------------------------------------------------------
  h_crystal_mult->Write();
  h_doppler->Write();
  h_doppler_forward_angles->Write();
  h_doppler_backward_angles->Write();
  h_doppler_simple->Write();
  h_doppler_mult1and2->Write();

  for(int i=0;i<15;i++){
    h_doppler_mult[i]->Write();
  }
  //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++){
  //  h_doppler_crystal[i]->Write();
  //}
  for(int i=0;i<17;i++)  {
    h_doppler_layer[i]->Write();
  }
  h_doppler_total->Write();
  h_doppler_total_individual_decays->Write();
  //h_doppler_max_ring->Write();
  h_energy->Write();
  h_energy_forward_angles->Write();
  h_energy_backward_angles->Write();
  //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++){
  //  h_energy_crystal[i]->Write();
  //}
  //h_energy_beta_real->Write();
  //h_beta_after_doppler->Write();
  h_crystal_fired_doppler->Write();
  h_crystal_fired_energy->Write();
  if(dali2_fi_opt >1) 
    h_doppler_ave_pos->Write();
  h_doppler_gamma_gamma->Write();
  if(dali2_addback_opt==1){
    h_crystal_fired_doppler_addback->Write();
    for(int i=0;i<3;i++){
      h_doppler_addback[i]->Write();
      //      for(int ppp=0;ppp<11;ppp++)  {
      //  h_doppler_addback_with_threshold[i][ppp]->Write();
      //}
    }
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
