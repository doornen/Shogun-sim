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
#include "TH3F.h"
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




void GetAverageInteractionPoint()  {
  cout<<" Getting the average FI point for all the DALI2 crystals..."<<endl;

  FILE *fOutAve  = fopen("../../output/AverageInteractionPoint.out","w");
  FILE *fOutDiff = fopen("../../output/InteractionDifferenceToCenterOfGravity.out","w");
  float difference[3];
  float distance1, distance2;
  float thetaDifference;
  for( int i = 0;i<NUMBEROFDALI2CRYSTALS;i++)  {
    for(int jjj=0;jjj<3;jjj++)  {
      fi_average[jjj][i] = fi_average[jjj][i]/fi_interactions[i];
      difference[jjj] = fi_average[jjj][i] - dali2_pos[jjj][i];
    }


    distance1 = sqrt(dali2_pos[0][i]*dali2_pos[0][i] + dali2_pos[1][i]*dali2_pos[1][i] + dali2_pos[2][i]*dali2_pos[2][i]);
    distance2 = sqrt(fi_average[0][i]*fi_average[0][i] + fi_average[1][i]*fi_average[1][i] + fi_average[2][i]*fi_average[2][i]);
    thetaDifference = acos(dali2_pos[0][i]*fi_average[0][i] + dali2_pos[1][i]*fi_average[1][i] + dali2_pos[2][i]*fi_average[2][i]
                           /distance1/distance2);

    TVector3 thetaangle(fi_average[0][i],fi_average[1][i],fi_average[2][i]);   
    // The average interaction point:
    fprintf(fOutAve,"%10.2f %10.2f %10.2f %4i %4i %10.2f 0 0\n",fi_average[0][i], fi_average[1][i], fi_average[2][i],i, fi_interactions[i], 180*thetaangle.Theta()/3.14159); 
//     cout<<i<<" "<<fi_interactions[i]<<" "<<fi_average[0][i]<<" "<<fi_average[1][i]<<" "<<fi_average[2][i];
    
    // Comparison to the center of gravity points:
    fprintf(fOutDiff,"%4i %10.2f %10.2f %10.2f %10.3f\n", i, difference[0],difference[1],difference[2],thetaDifference); 
//     cout<<" "<<i<<" "<<difference[0]<<" "<<difference[1]<<" "<<difference[2]<<endl;
  
    fi_interactions[i]=0;
  }
  fclose(fOutAve); 
  fclose(fOutDiff);
} //end of GetAverageInteractionPoint


bool IncludeAddbackTable(float det[3][2])  {

  float distance = TMath::Sqrt(TMath::Power(det[0][0]-det[0][1],2) +
                               TMath::Power(det[1][0]-det[1][1],2) + 
                               TMath::Power(det[2][0]-det[2][1],2));

  //cout<<"Distance: "<<distance<<endl;
  if( distance > maxAddbackDistance ) return false;
  else return true;
} //end of IncludeAddbackTable





void ReadInputFile(){
  
  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/RikenReconstructorDLB.in","r");
  while(!feof(fin))  {
    fscanf(fin,"%s ",temp);  
    if(strcmp(temp,"INPUTFILE")==0)  {
      fscanf(fin,"%300s",&root_in[0]);
    }
    else if(strcmp(temp,"OUTPUTFILE")==0)  {
      fscanf(fin,"%300s ",&root_out[0]);
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
      if(reduction_factor<1) {
        cout << "Error: Reduction factor must not be smaller than 1! Aborting." << endl;
        abort();
//        return 0;
      }
    }
//     else if(strcmp(temp,"DALI2FIFIND")==0)  {
//       fscanf(fin,"%i",&dali2_fi_option); 
//       printf("%s %i \n",temp,dali2_fi_option);
//     }
    else if(strcmp(temp,"DALI2INCLUDE")==0) {
      fscanf(fin,"%i ",&dali2_opt); 
      printf("%s %i \n",temp,dali2_opt);
    }
    else if(strcmp(temp,"LABR3INCLUDE")==0)  {
      fscanf(fin,"%i ",&LaBr3Opt); 
      printf("%s %i \n",temp,LaBr3Opt);
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
      cout<<"Could not read your input keyword '" << temp << "'. Aborting program."<<endl; 
      abort();
    }
  }
  
  
//  return 0;
} // end of ReadInputFile







void ResetValuesAll(){
  
  
  ringEnergyMax = -1;
      crystalMultDali2   = 0;
      dali2_energySum     = 0.;
      dali2_energyMax     = 0.;
      dali2_doppler_total = 0.0;
//       dali2_doppler_addback= 0.0;
      daliAddbackMul=0;
      for(int i=0; i<20; i++){dali2_doppler_addback[i]=0.0;}
      for(int i =0;i<maxDecays;i++)  {
//         dali2_crystal_id[i]=-1;
        dali2_crystalFired[i]   = -1;
        IndividualDecaysMul  = 0;
        dopplerEnergyIndividualDecays[i] = 0.0;
        energySumIndividualDecays[i] = 0.0;
        energyMaxIndividualDecays[i] = 0.0;
        posThatsItDetId[i]=0; //1=dali2, 6=labr3
      }
      
      
      for(int i =0;i<NUMBEROFDALI2CRYSTALS;i++){
        dali2_energy[i]         = 0.;
        dali2_time[i]           = 0.;
        dali2_dopplerEnergy[i]  = 0.;
//         dali2_crystalFired[i]   = -1;
        crystalUsedForAddback[i] = false;
      }
      //Getting the trigger counter:
        if(trigger_opt==1)
          for(int j=0;j<11;j++)  
            for(int k=0;k<3;k++)
              trigger_flag[j][k]=false;
  
  crystalMultLaBr3  =0;
  labr3_energySum  =0.0;
//   labr3_energyMax  = 0.0;
  for(int i=0;i<NUMBEROFLABR3ARRAYCRYSTALS;i++)  {
  labr3_energy[i]    = 0.0;  // Unsorted energy, filled according to the loop from 0 to the last crystal
  labr3_dopplerEnergy[i]  = 0.0;  // Unsorted dopplerEnergy
  labr3_crystalFired[i]  = 0;
  }
  
  for(int i=0; i<maxDecays; i++){
  labr3_crystal_id[i]  =-1;
  labr3_crystal_energy[i]  =0.0;
  labr3_crystal_dopplerEnergy[i]=0.0;
  }
  
  
  
} // end of ResetValuesAll






















//int RikenReconstructorDLB()  {
Int_t main()  {
  // int main(int argc,char** argv)  { 
  
  //pschrock: just for testing
  //use this option only for setups with 3 LaBr crystals
  bool threeGammaPlot=false;
  
  
  //Reading the input file, which can change the values
  ReadInputFile();
  
  
  //________________________________________________________
  // Check dali2 options before starting the analysis
  if(dali2_opt==0 && dali2_addback_opt==1){
    cout << "Dali-2 addback requested, but Dali-2 is not included. " << endl;
    cout << "Dali-2 addback will be deactivated!!!" << endl;
    dali2_addback_opt=0;
  }
  if(dali2_opt==0 && dali2_fi_opt==1){
    cout << "Dali-2 first interaction option requested, but Dali-2 is not included." << endl;
    cout << "Dali-2 first interaction will be deactivated!!!" << endl;
    dali2_fi_opt=0;
  }
  if(dali2_opt==0 && trigger_opt==1){
    cout << "Dali-2 trigger option requested, but Dali-2 is not included." << endl;
    cout << "Dali-2 trigger will be deactivated!!!" << endl;
    trigger_opt=0;
  }
  
  
  //--------------------------------------------------------------------------------------------------
  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();

  
  TH1F *h_trigger_prob[3];
  
  TH1F *h_dali2_crystal_mult;

  TH1F *h_dali2_doppler;

  TH1F *h_dali2_doppler_forward_angles;
  TH1F *h_dali2_doppler_backward_angles;

  TH1F *h_dali2_doppler_simple;

  TH1F *h_dali2_doppler_mult[20];

  TH1F *h_dali2_doppler_mult1and2;
  //TH1F *h_doppler_crystal[NUMBEROFDALI2CRYSTALS];
  TH1F *h_dali2_doppler_layer[17];

  TH1F *h_dali2_doppler_total;
  TH1F *h_doppler_total_individual_decays; 
  //TH2F *h_doppler_max_ring;

  TH1F *h_dali2_energy;
  TH1F *h_dali2_energy_forward_angles;
  TH1F *h_dali2_energy_backward_angles;
  //TH1F *h_dali2_energy_crystal[NUMBEROFDALI2CRYSTALS];
  //TH2F *h_dali2_energy_beta_real;

  //TH2F *h_beta_after_doppler;
  TH2F *h_dali2_crystal_fired_doppler;
  TH2F *h_dali2_crystal_fired_energy;

  TH1F *h_dali2_doppler_ave_pos;
  
  TH2F *h_dali2_doppler_gamma_gamma;
  TH2F *h_dali2_addback_gamma_gamma;

  TH1F *h_dali2_doppler_addback[3];  // 0: all, 1: backward,2: forward angles
  //TH1F *h_doppler_addback_with_thresold[3][11];  // Threshold from 0 to 500 keV in steps of 50 keV

  // Spectra of the trigger probablility, also separated into backward and forward angles:
  if(trigger_opt == 1)  {
    for(int i=0;i<3;i++)  {
      sprintf(temp,"h_trigger_prob[%i]",i);
      h_trigger_prob[i] = new TH1F(temp,temp,501,-0.5,500.5);
    }
  }

  //Creating dali spectra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Crystal multiplicity,ok
  h_dali2_crystal_mult = new TH1F("dali2_crystal_Mult","dali2_crystal_mult",NUMBEROFDALI2CRYSTALS,0,NUMBEROFDALI2CRYSTALS);  
  // All, OK
  h_dali2_doppler = new TH1F("dali2_doppler","dali2_doppler",numBin,firstBin,lastBin);
  h_dali2_doppler_forward_angles = new TH1F("dali2_doppler_forward_angles","dali2_doppler_forward_angles",numBin,firstBin,lastBin);
  h_dali2_doppler_backward_angles = new TH1F("dali2_doppler_backward_angles","dali2_doppler_backward_angles",numBin,firstBin,lastBin);
  h_dali2_doppler_simple = new TH1F("dali2_doppler_simple","dali2_doppler_simple",numBin,firstBin,lastBin);

  h_dali2_doppler_mult1and2 = new TH1F("dali2_doppler_mult1and2","dali2_doppler_mult1and2",numBin,firstBin,lastBin);
  
  // According to the multiplicity
  for(int i=0;i<20;i++)  {
    sprintf(temp,"dali2_doppler_mult[%i]",i);
    h_dali2_doppler_mult[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  } 

  // The spectra crystalwise, OK
  //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
  //  sprintf(temp,"doppler_crystal[%i]",i);
  //  h_doppler_crystal[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  //}  
  // By layer:
  for(int i=0;i<17;i++)  {
    sprintf(temp,"dali2_doppler_layer[%i]",i);
    h_dali2_doppler_layer[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  }
  // Taking the sum energy and making the Doppler correction with the crystal that has the highest energy:
  h_dali2_doppler_total = new TH1F("dali2_doppler_total","dali2_doppler_total",numBin,firstBin,lastBin);
  h_doppler_total_individual_decays = new TH1F("doppler_total_individual_decays","doppler_total_individual_decays",numBin,firstBin,lastBin);
  // The ring with the most energy:
  //for(int i=0;i<15;i++)
  //{
  //sprintf(temp,"doppler_max_ring");
  //h_doppler_max_ring = new TH2F(temp,temp,numBin,firstBin,lastBin,15,0,15);
  //}

  // Not doppler corrected, OK
  h_dali2_energy = new TH1F("dali2_energy","dali2_energy",numBin,firstBin,lastBin);
  h_dali2_energy_forward_angles = new TH1F("dali2_energy_forward_angles","dali2_energy_forward_angles",numBin,firstBin,lastBin);
  h_dali2_energy_backward_angles = new TH1F("dali2_energy_backward_angles","dali2_energy_backward_angles",numBin,firstBin,lastBin);
  // The spectra crystalwise,OK
  //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
  //  sprintf(temp,"energy_crystal[%i]",i);
  //  h_energy_crystal[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  //}
  // Comparing beta REAL with dali energy:
  //h_energy_beta_real = new TH2F("energy_beta_real","energy_beta_real",numBin,firstBin,lastBin,200,0.40,0.45);

  // Comparing beta after the target with doppler corrected energy
  //h_beta_after_doppler = new TH2F("beta_after_doppler","",200,0.5,0.64,numBin,firstBin,lastBin);

  h_dali2_crystal_fired_doppler = new TH2F("dali2_crystal_fired_doppler","",NUMBEROFDALI2CRYSTALS,0,
                                     NUMBEROFDALI2CRYSTALS,numBin,firstBin,lastBin);
  h_dali2_crystal_fired_energy = new TH2F("dali2_crystal_fired_energy","",NUMBEROFDALI2CRYSTALS,0,
                                    NUMBEROFDALI2CRYSTALS,numBin,firstBin,lastBin);
  
  
  
  //LaBr3 histograms
  TH1F* h_labr3_crystal_mult = new TH1F("labr3_crystal_Mult","labr3_crystal_mult",NUMBEROFLABR3ARRAYCRYSTALS,0,NUMBEROFLABR3ARRAYCRYSTALS);  
  TH1F* h_labr3_energy = new TH1F("labr3_energy","labr3_energy",2.0*numBin,firstBin,2.0*lastBin);
  TH2F* h_labr3_crystal_fired_energy = new TH2F("labr3_crystal_fired_energy","",NUMBEROFLABR3ARRAYCRYSTALS,0,NUMBEROFLABR3ARRAYCRYSTALS,2.0*numBin,firstBin,2.0*lastBin);
  TH2F* h_labr3_crystal_fired_doppler= new TH2F("labr3_crystal_fired_doppler","",NUMBEROFLABR3ARRAYCRYSTALS,0,NUMBEROFLABR3ARRAYCRYSTALS,numBin,firstBin,lastBin);
  
  TH1F* h_labr3_doppler_total = new TH1F("labr3_doppler_total","labr3_doppler_total",numBin,firstBin,lastBin);
  
  
  

  // For gamma-gamma coincidences
  
  h_dali2_doppler_gamma_gamma = new TH2F("dali2_doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);
//   TH2F* h_labr3_doppler_gamma_gamma = new TH2F("labr3_doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);
  h_dali2_addback_gamma_gamma = new TH2F("h_dali2_addback_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);
  
  TH1F* h_labr3_doppler_gamma[4];
  TH2F* h_labr3_doppler_gamma_gamma[4];
  TH2F* h_labr3_doppler_gamma_gamma_crystalMul3 = new TH2F("h_labr3_doppler_gamma_gamma_crystalMul3","h_labr3_doppler_gamma_gamma_crystalMul3",numBin,firstBin,lastBin,numBin/2,firstBin,lastBin/2.0);
  
  for(int i=0; i<4; i++){
    sprintf(temp,"labr3_doppler_gamma[%i]",i);
    h_labr3_doppler_gamma[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
    
    sprintf(temp,"labr3_doppler_gamma_gamma[%i]",i);
     h_labr3_doppler_gamma_gamma[i] = new TH2F(temp,temp,numBin,firstBin,lastBin,numBin/2,firstBin,lastBin/2.0);
     
  }
   
  TH3F* h3_labr;
  TH3F* h3_labr_cuts;
  if(threeGammaPlot){
    //h3_labr=new TH3F("h3_labr","h3_labr",numBin,firstBin,lastBin,numBin,firstBin,lastBin,numBin,firstBin,lastBin);
    h3_labr=new TH3F("h3_labr","h3_labr",200,firstBin,lastBin,200,firstBin,lastBin,200,firstBin,lastBin);
    h3_labr_cuts=new TH3F("h3_labr_cuts","h3_labr_cuts",200,firstBin,lastBin,200,firstBin,lastBin,200,firstBin,lastBin);
  }
  
  TH2F* h_labr3_dali2_doppler_gamma_gamma = new TH2F("labr3_dali2_doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);
  TH2F* h_labr3_dali2Total_doppler_gamma_gamma = new TH2F("h_labr3_dali2Total_doppler_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);
  TH2F* h_labr3_dali2addback_gamma_gamma = new TH2F("h_labr3_dali2addback_gamma_gamma","",numBin,firstBin,lastBin,numBin,firstBin,lastBin);
  
  
  
  
  //______________________________________________
  if(dali2_fi_opt==2)  {
    FILE *fInAve  = fopen("output/AverageInteractionPoint.out","r");
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
    h_dali2_doppler_ave_pos = new TH1F("dali2_doppler_ave_pos","dali2_doppler_ave_pos",numBin,firstBin,lastBin);
  }

  // Testing the addback procedure with cluster grouping
  if(dali2_addback_opt==1)  {
    for(int i=0;i<3;i++)  {
      sprintf(temp,"h_dali2_doppler_addback[%i]",i);
      h_dali2_doppler_addback[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
      //for(int j=0;j<11;j++)  {
      //  sprintf(temp,"h_doppler_addback_with_threshold[%i][%i]",i,j);
      //  h_doppler_addback_with_threshold[i][j] = new TH1F(temp,temp,numBin,firstBin,lastBin);
      //}
    }
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
  t->SetBranchAddress("VGamma",gv);                                    // Gamma vector
  t->SetBranchAddress("EGammaRest",&e_rest);                            // Energy at rest
  t->SetBranchAddress("EGammaDoppler",&e_doppler);                              // Theta of doppler boosted gamma
  t->SetBranchAddress("ThetaGammaRest",&theta_gamma_rest);
  t->SetBranchAddress("ThetaGammaLab",&theta_gamma_lab);
  t->SetBranchAddress("EnergyVertex",&energy_vertex_stored);                    // Energy of fragment at fragmentation
  t->SetBranchAddress("BGGamma",&bgGamma);                                      // To know if the gamma is from bg or a good event
  t->SetBranchAddress("ObsGammaMul",&gammaMul);                                 // Multiplicity information from simulation: how many gammas were emmittet
  //So far identical to the ROOT-Tree of step one
  //The new stuff for step 2:
  //Setting the type of the gamma detector
  t->SetBranchAddress("GammaDetType",&gamma_det_type);
  //-------------------------------
  //DALI:
  if(dali2_opt==1){
  t->SetBranchAddress("DALI2Flag",dali2_flag);                                  // ID of which det recorded a gamma in the event
  t->SetBranchAddress("DALI2EnergyNotCor",dali2_energy_not_cor);
  t->SetBranchAddress("DALI2Time",dali2_time);
  }
  if(dali2_fi_opt >= 1)  {
    t->SetBranchAddress("DALI2FI",dali2_fi);
  }
  //LaBr3
  if(LaBr3Opt==1)  {
    t->SetBranchAddress("LaBr3Flag",labr3_flag);
    t->SetBranchAddress("LaBr3EnergyNotCor",labr3_energy_not_cor);
    t->SetBranchAddress("LaBr3Time",labr3_time);
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
  tHeader->SetBranchAddress("ZFragment",&z_f);                                  // Fragment z
  tHeader->SetBranchAddress("ChargeFragment",&charge_f);                        // Fragment charge
  tHeader->SetBranchAddress("TargetKind",&kind_target);
  tHeader->SetBranchAddress("TargetThicknessCM",&thickness);
  tHeader->SetBranchAddress("TotalEventNumber",&total_event_num);
  tHeader->SetBranchAddress("ThetaLow",&theta_low);
  tHeader->SetBranchAddress("ThetaHigh",&theta_high);
  //DALI2:
  if(dali2_opt==1){
  tHeader->SetBranchAddress("DALI2EnResOpt",&dali2_en_res_opt);
  tHeader->SetBranchAddress("DALI2EnRes",dali2_en_res);
  tHeader->SetBranchAddress("DALI2TimeRes",dali2_time_res);
  tHeader->SetBranchAddress("DALI2Pos",dali2_pos);
  }
  //LaBr
  if(LaBr3Opt==1)  {
    tHeader->SetBranchAddress("LaBr3EnResOpt",&LaBr3EnResOpt); //not used in this Reconstructor
    tHeader->SetBranchAddress("LaBr3EnRes",LaBr3EnRes);
    tHeader->SetBranchAddress("LaBr3TimeRes",LaBr3TimeRes);
    tHeader->SetBranchAddress("LaBr3Pos",labr3_pos);
  }
  tHeader->SetBranchAddress("BetaResolution",&beta_res);  
  tHeader->SetBranchAddress("PosDetAtTargetRes",&pos_det_at_target_res); 
  tHeader->SetBranchAddress("PosDetAfterTargetRes",&pos_det_after_target_res);
  tHeader->GetEntry(0);
  //-----------------------------------------------------------------------------------------------------------
  //Finished reading the Root file

  
  //create output root tree
  //add new branches with Dali2 and LaBr3 data
  rootfile->cd();
  
  //make a copy of the root tree from EventBuilder
  t2Header = tHeader->CloneTree();
  t2 = t->CloneTree(0);
  t2->SetName("ReconstructedEvents");
  
  
  if(dali2_opt==1){
    t2->Branch("RecDali2CrystalMul", &crystalMultDali2, "RecDali2CrystalMul/I");
    t2->Branch("RecDali2CrystalId", dali2_crystalFired, "RecDali2CrystalId[20]/I");
    t2->Branch("RecDali2CrystalEnergyDopplerCorrected", dali2_dopplerEnergy, "RecDali2EnergyDopplerCorrected[20]/F");
    t2->Branch("RecDali2EnergyDopplerCorrectedTotal", &dali2_doppler_total, "RecDali2EnergyDopplerCorrectedTotal/F");
    
  }
  if(dali2_addback_opt==1){
//     t2->Branch("RecDali2EnergyDopplerCorrectedAddback", &dali2_doppler_addback, "RecDali2EnergyDopplerCorrectedAddback/F");
    t2->Branch("RecDali2AddbackMul", &daliAddbackMul, "RecDali2AddbackMul/I");
    t2->Branch("RecDali2EnergyDopplerCorrectedAddback", dali2_doppler_addback, "RecDali2EnergyDopplerCorrectedAddback[20]/F");
  
  }
  
  if(LaBr3Opt==1)  {
    
    //zero suppressed LaBr3 energies:
    t2->Branch("RecLaBr3CrystalMul", &crystalMultLaBr3,"RecLaBr3CrystalMul/I");
    t2->Branch("RecLaBr3CrystalId", labr3_crystal_id,"RecLaBr3CrystalId[20]/I");
    t2->Branch("RecLaBr3CrystalEnergyNotCorrected", labr3_crystal_energy,"RecLaBr3CrystalEnergyNotCorrected[20]/F");
    t2->Branch("RecLaBr3CrystalEnergyDopplerCorrected", labr3_crystal_dopplerEnergy,"RecLaBr3CrystalEnergyDopplerCorrected[20]/F");
    
    t2->Branch("RecLaBr3EnergyDopplerCorrected", labr3_dopplerEnergy,"RecLaBr3EnergyDopplerCorrected[8]/F");
    
    
  }
  
  //individual decays, independent of the chosen detectors
  t2->Branch("RecIndividualDecaysMul", &IndividualDecaysMul,"RecIndividualDecaysMul/I");
  t2->Branch("RecIndividualDecaysDetId", posThatsItDetId,"RecIndividualDecaysDetId[20]/I");
  t2->Branch("RecIndividualDecaysDopplerEnergy", dopplerEnergyIndividualDecays, "RecIndividualDecaysDopplerEnergy[20]/F");
  
  
  
  
  // Before I start the analysis, I make an addback-table:
  //________________________________________________________
  if(dali2_addback_opt==1)  {
    FILE *fAddbackTableIn  = fopen("output/AddbackTable.out","w");
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
          cout<<"You have to increase the variable NUMBEROFDALI2ADDBACKCRYSTALS in Globals.hh!!!"<<endl;
          abort();
        }
      }
      counter = 0;
      fprintf(fAddbackTableIn," \n");
    }
    fclose(fAddbackTableIn);
  } //dali2 addback
  
  
  
  nentries = (Int_t)t->GetEntries()/reduction_factor;
  cout <<"Entries to process: "<< nentries << " (reduction factor " << reduction_factor << ")" << endl;

  //variables defined in MyNamespace.hh

  cout<<"Starting the analysis"<<endl;
  for(int iii=0;iii<nentries;iii++)  {
    if((iii+1)%100000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    
    if(iii<nentries-1)  {
      t->GetEntry(iii+1);
      
      eventNumberNext = evnum;
    }
    t->GetEntry(iii);
    
    
    
    //-------------------------------------------------------------
    // Starting with the gamma ray detector spectra:
    if(dali2_opt==1){
    for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++)  {
      //removing crystals that have a bad energy resolution.
      //NP0801-RIBF70
      //       if(nnn==0  || nnn==1  || nnn==2  || nnn==3 || nnn==5  || nnn==6  || nnn==39 || nnn==63 || 
      //     nnn==123 || nnn==124|| nnn==125||nnn==126|| nnn==127|| nnn==130|| nnn==132) continue;
      
      if(dali2_flag[nnn]==1)   {
        dali2_energy[nnn] += dali2_energy_not_cor[nnn];
        //if(nnn==180) cout<<"Energy : "<<energy[nnn]<<endl;
        dali2_energySum   +=  dali2_energy_not_cor[nnn];
        if(maxDecays>decayIDevnum) {
          energySumIndividualDecays[(int)decayIDevnum] = energySumIndividualDecays[(int)decayIDevnum] + dali2_energy_not_cor[nnn];
          // Getting the right position for the Doppler correction of individual decays:
          // The hit with the highest energy entry is expected to be the correct one
          if(dali2_energy_not_cor[nnn]>energyMaxIndividualDecays[(int)decayIDevnum])  {
            energyMaxIndividualDecays[(int)decayIDevnum]    = dali2_energy_not_cor[nnn];
            posThatsItIndividualDecays[(int)decayIDevnum][0]= dali2_pos[0][nnn];
            posThatsItIndividualDecays[(int)decayIDevnum][1]= dali2_pos[1][nnn];
            posThatsItIndividualDecays[(int)decayIDevnum][2]= dali2_pos[2][nnn];
            posThatsItDetId[(int)decayIDevnum]=1;
          }
        }else{
          cout << "There are " << decayIDevnum+1 << " decays registered, but analysis supports only " << maxDecays << "! Please increase 'maxDecays' in MyNamespace.hh!" << endl;
          abort();
        }
        //Getting the trigger counter:
        if(trigger_opt==1)
          for(int j=0;j<11;j++)  {
            if(dali2_energy[nnn]>(j*500))  {
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
      }//dali2 flag 1
    }//loop over dali2 crystals
    }//dali2opt
    
    
    //the same for LaBr3
    if(LaBr3Opt==1){
      for(int nnn=0;nnn<NUMBEROFLABR3ARRAYCRYSTALS;nnn++){
        if(labr3_flag[nnn]==1){
            // bool newCrystal=false;
            // if(labr3_energy[nnn]==0.0){newCrystal=true;}
          
            labr3_energy[nnn] += labr3_energy_not_cor[nnn];
            labr3_energySum   += labr3_energy_not_cor[nnn];
            
            if(maxDecays>decayIDevnum) {
              energySumIndividualDecays[(int)decayIDevnum] = energySumIndividualDecays[(int)decayIDevnum] + labr3_energy_not_cor[nnn];
              // Getting the right position for the Doppler correction of individual decays:
              // The hit with the highest energy entry is expected to be the correct one
              if(labr3_energy_not_cor[nnn]>energyMaxIndividualDecays[(int)decayIDevnum])  {
              energyMaxIndividualDecays[(int)decayIDevnum] = labr3_energy_not_cor[nnn];
              posThatsItIndividualDecays[(int)decayIDevnum][0]  = labr3_pos[0][nnn];
              posThatsItIndividualDecays[(int)decayIDevnum][1]  = labr3_pos[1][nnn];
              posThatsItIndividualDecays[(int)decayIDevnum][2]  = labr3_pos[2][nnn];
              posThatsItDetId[(int)decayIDevnum]=6;
            }
            }else{
              cout << "There are " << decayIDevnum+1 << " decays registered, but analysis supports only " << maxDecays << "! Please increase 'maxDecays' in MyNamespace.hh!" << endl;
              abort();
            }
      
      
          } //LaBr3 flag
      
      } //loop over LaBr3 crystals
    } //LaBr3 opt 1
    
    
    
    
    
    
    if(eventNumberNext!=evnum) {
    // The last observed decay - filling the spectra can be started
      
      
      if(dali2_opt==1){
      for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++)  {
        //if(nnn==23 || nnn==69) continue;  //removing crystals that have a bad energy resolution.
        if(dali2_energy[nnn]>energyThreshold)  {
           // Filling the energy spectra without Doppler-correction
          h_dali2_energy->Fill(dali2_energy[nnn]);
            
          if(nnn>=detectorLimit) {
            h_dali2_energy_forward_angles->Fill(dali2_energy[nnn]);
          }
          else {h_dali2_energy_backward_angles->Fill(dali2_energy[nnn]);}
          //h_dali2_energy_crystal[nnn]->Fill(dali2_energy[nnn]);
          //h_dali2_energy_beta_real->Fill(dali2_energy[nnn],beta_r);
          h_dali2_crystal_fired_energy->Fill(nnn,dali2_energy[nnn]);
          
          //_____________________________________________________________
          //Performing the Doppler correction,
          //This one is detectorwise!
          float correctedEnergyCrystal;
          float correctedEnergySimple;
          //dali2_fi_opt==2 I use the average FI point from the simulation
          if(dali2_fi_opt !=2){
            correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                               ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                               dali2_energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
          }
          else{
            correctedEnergyCrystal = GetDopplerCorrectedEnergy(dali2_pos_ave[0][nnn],dali2_pos_ave[1][nnn],dali2_pos_ave[2][nnn],
                                                               ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                               dali2_energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z); 
          }
          if(dali2_fi_opt !=2) {
            correctedEnergySimple = GetDopplerCorrectedEnergy(dali2_pos[0][nnn],dali2_pos[1][nnn],dali2_pos[2][nnn],
                                                              0,0,0,0,0,10000, 
                                                              dali2_energy[nnn],0,0,beta_average,decay_z);
          }
          else{
            correctedEnergySimple = GetDopplerCorrectedEnergy(dali2_pos_ave[0][nnn],dali2_pos_ave[1][nnn],dali2_pos_ave[2][nnn],
                                                              0,0,0,0,0,10000, 
                                                              dali2_energy[nnn],0,0,beta_average,0);
          }
          
          dali2_dopplerEnergy[crystalMultDali2] = correctedEnergyCrystal;
          h_dali2_doppler->Fill(correctedEnergyCrystal);
          // Checking the difference between forward and backward angles.
          if(nnn>=detectorLimit) {
            h_dali2_doppler_forward_angles->Fill(correctedEnergyCrystal);
          }
          else {h_dali2_doppler_backward_angles->Fill(correctedEnergyCrystal);}
          
          h_dali2_doppler_simple->Fill(correctedEnergySimple);

          //h_doppler_crystal[nnn]->Fill(correctedEnergyCrystal);
          h_dali2_crystal_fired_doppler->Fill(nnn,correctedEnergyCrystal);
          
          dali2_dopplerEnergy[crystalMultDali2] = correctedEnergyCrystal;
          dali2_crystalFired[crystalMultDali2] = nnn;
          
          // dali2_crystal_id[crystalMultDali2] = nnn+1; //start counting from 1 for the root tree
          
          //Set the threshold
          //if(dali2_energy[nnn]>150)
          crystalMultDali2++;
          //Filling Dali layer wise:
          for(int ppp=0;ppp<17;ppp++)  {
            if(nnn>=layerLow[ppp] && nnn<=layerHigh[ppp]) {h_dali2_doppler_layer[ppp]->Fill(correctedEnergyCrystal);}
          }
          //See how addback procedures work
          if(dali2_energy[nnn]>dali2_energyMax)  {
            dali2_energyMax = dali2_energy[nnn];
            posThatsItDali2[0]  = dali2_pos[0][nnn];
                  posThatsItDali2[1]  = dali2_pos[1][nnn];
                  posThatsItDali2[2]  = dali2_pos[2][nnn];
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
                  h_dali2_doppler_ave_pos->Fill(correctedEnergyCrystal);
                }
        } //dali energy > energy threshold
      } //loop over dali2 crystals
      
      h_dali2_crystal_mult->Fill(crystalMultDali2);
      //SortEnergies();           //Sort the observed energies from highest to lowest!!!!

      if(crystalMultDali2>0)  {
        for(int i = 1;i<20;i++)  {
          if(crystalMultDali2==i)  {   
            for(int j = 0;j<i;j++)  {
              h_dali2_doppler_mult[i]->Fill(dali2_dopplerEnergy[j]);
              if(crystalMultDali2==1 || crystalMultDali2==2) {h_dali2_doppler_mult1and2->Fill(dali2_dopplerEnergy[j]);}
            }
          }
        }
      }
      
      } //end of dali2opt
      
      
      
      
      // ___________________________________________________
      // Checking realistic add-back procedure:
      // It includes also events for which only one detector saw a gammaray!!!
      // Thus -> Addback + singles
      if(crystalMultDali2>=1 && dali2_addback_opt==1)  {
        for(int i = 0;i<crystalMultDali2;i++)  {  
          if(crystalUsedForAddback[dali2_crystalFired[i]]==true) continue;
          float dummyEnergy = dali2_dopplerEnergy[i];
          //float dummyEnergyWithThreshold[11];
          //for(int ppp=0;ppp<11;ppp++)  {
          //  dummyEnergyWithThreshold[ppp]= dali2_dopplerEnergy[i];
          //}

          crystalUsedForAddback[dali2_crystalFired[i]]=true; 
           
          for(int j = i;j<crystalMultDali2;j++){
            if(crystalUsedForAddback[dali2_crystalFired[j]]==false&& dali2_dopplerEnergy[j]>0){
              for(int k = 0;k<NUMBEROFDALI2ADDBACKCRYSTALS;k++) {
                if(dali2_crystalFired[j] == addbackTable[dali2_crystalFired[i]][k]){
          
                  crystalUsedForAddback[dali2_crystalFired[j]]=true;
                  dummyEnergy = dummyEnergy + dali2_dopplerEnergy[j];
                  //for(int ppp=0;ppp<11;ppp++)  {
                  //  if(dali2_dopplerEnergy[i]>=ppp*50&& dali2_dopplerEnergy[j]>=ppp*50) 
                  //    dummyEnergyWithThreshold[ppp] = dummyEnergyWithThreshold[ppp] +  dali2_dopplerEnergy[j];
                  //}
                }
              }
            }
          }
          h_dali2_doppler_addback[0]->Fill(dummyEnergy);
          
          
          dali2_doppler_addback[daliAddbackMul]=dummyEnergy;
          daliAddbackMul++;
          
          //for(int ppp=0;ppp<11;ppp++)  {
          //  h_doppler_addback_with_threshold[0][ppp]->Fill(dummyEnergyWithThreshold[ppp]);
          //}
          if(posThatsItDali2[2]<0) {
            h_dali2_doppler_addback[1]->Fill(dummyEnergy);
          }
          else {
            h_dali2_doppler_addback[2]->Fill(dummyEnergy);
          }   
        }
      } // end of addback
      
      
      //____________________________________________________
      // Making a Doppler correction with the Sum of the Energy detected in the array
      // The Position is taken from the crystal that recorded the highest energy
      if(dali2_opt==1){
        float correctedEnergyTotal = GetDopplerCorrectedEnergy(posThatsItDali2[0],posThatsItDali2[1],posThatsItDali2[2],
                      ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                      dali2_energySum,beta_rec,beta_mean_tof,beta_average,decay_z);
        h_dali2_doppler_total->Fill(correctedEnergyTotal);
        dali2_doppler_total=correctedEnergyTotal;
      }
      
      
      
  if(LaBr3Opt==1){
    
    //first: create a zero suppressed array
    crystalMultLaBr3=0;
    for(int ii=0; ii<NUMBEROFLABR3ARRAYCRYSTALS; ii++){
      //zero suppressed array:
      if(labr3_energy[ii]>0.0){
        labr3_crystal_id[crystalMultLaBr3]=ii+1;
        labr3_crystal_energy[crystalMultLaBr3] += labr3_energy[ii];
        
        crystalMultLaBr3++;
      }
    }
    
    
    h_labr3_crystal_mult->Fill(crystalMultLaBr3);
    
    float correctedEnergyCrystal=0.0;
    
    //sum over all hits in LaBr3
    for(int nnn=0;nnn<crystalMultLaBr3;nnn++){
      correctedEnergyCrystal=0.0;
      
      if(labr3_crystal_energy[nnn]>energyThreshold)  {
        
        // Filling the energy spectra without Doppler-correction
        h_labr3_energy->Fill(labr3_crystal_energy[nnn]);
        h_labr3_crystal_fired_energy->Fill(nnn,labr3_crystal_energy[nnn]);
        
        //make doppler correction crystalwise
        if(labr3_crystal_energy[nnn]>energyThreshold){
          correctedEnergyCrystal = GetDopplerCorrectedEnergy(  labr3_pos[0][labr3_crystal_id[nnn]-1],
                        labr3_pos[1][labr3_crystal_id[nnn]-1], 
                        labr3_pos[2][labr3_crystal_id[nnn]-1],
                  ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                  labr3_crystal_energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
        }
        
        h_labr3_doppler_gamma[0]->Fill(correctedEnergyCrystal);
        if(gammaMul>0 && gammaMul<4){
          h_labr3_doppler_gamma[gammaMul]->Fill(correctedEnergyCrystal);
        }
        h_labr3_crystal_fired_doppler->Fill(nnn,correctedEnergyCrystal);
        labr3_crystal_dopplerEnergy[nnn]=correctedEnergyCrystal;
        labr3_dopplerEnergy[labr3_crystal_id[nnn]-1]=correctedEnergyCrystal; //not zero suppressed
        
        
      
      }
    } //end of loop over LaBr3 crystals
    
    
//     for(int nnn=0;nnn<NUMBEROFLABR3ARRAYCRYSTALS;nnn++){
//     correctedEnergyCrystal=0.0;
//     
//     if(labr3_energy[nnn]>energyThreshold)  {
//       
//       // Filling the energy spectra without Doppler-correction
//       h_labr3_energy->Fill(labr3_energy[nnn]);
//       h_labr3_crystal_fired_energy->Fill(nnn,labr3_energy[nnn]);
//       
//       //make doppler correction crystalwise
//       if(labr3_energy[nnn]>energyThreshold){
//         correctedEnergyCrystal = GetDopplerCorrectedEnergy(labr3_pos[0][nnn],labr3_pos[1][nnn],labr3_pos[2][nnn],
//                 ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
//                 labr3_energy[nnn],beta_rec,beta_mean_tof,beta_average,decay_z);
//       }
//       
//       h_labr3_crystal_fired_doppler->Fill(nnn,correctedEnergyCrystal);
//       labr3_dopplerEnergy[nnn]=correctedEnergyCrystal;
//     
//     }
//     } //end of loop over LaBr3 crystals
//     
    
    
    
    
    
    //correct the sum energy in LaBr
    if(labr3_energySum>0.0){
      float correctedEnergyTotal = GetDopplerCorrectedEnergy(posThatsItLaBr3[0],posThatsItLaBr3[1],posThatsItLaBr3[2],
                    ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                    labr3_energySum,beta_rec,beta_mean_tof,beta_average,decay_z);
      h_labr3_doppler_total->Fill(correctedEnergyTotal);
      
    }
    
    
    
  } //end of LaBr3
      
      
      
      
      
      
      
      
      
      int counter = 0;

      //____________________________________________________
      // Treating decays individually:
//       if(dali2_opt==1){
      while(energySumIndividualDecays[counter]>0 && counter<maxDecays)  {
        dopplerEnergyIndividualDecays[counter] =  GetDopplerCorrectedEnergy(posThatsItIndividualDecays[counter][0],
                                                                            posThatsItIndividualDecays[counter][1],
                                                                            posThatsItIndividualDecays[counter][2],
                                                                            ver_rec[0],ver_rec[1],ver_rec[2],pos_det[0],pos_det[1],pos_det[2], 
                                                                            energySumIndividualDecays[counter],
                                                                            beta_rec,beta_mean_tof,beta_average,decay_z);
        h_doppler_total_individual_decays->Fill(dopplerEnergyIndividualDecays[counter]);
        
        
        
        counter++;
        IndividualDecaysMul=counter;
      }
//       }
      //____________________________________________________
      //Making the gamma-gamma coincidencs
      //dali2-dali2
      if(crystalMultDali2>=2)  {
        for(int i = 0;i<crystalMultDali2;i++){
          for(int j = i+1;j<crystalMultDali2;j++){
            if (dali2_dopplerEnergy[i]>=dali2_dopplerEnergy[j]) {
              h_dali2_doppler_gamma_gamma->Fill(dali2_dopplerEnergy[i],dali2_dopplerEnergy[j]);
            }else {h_dali2_doppler_gamma_gamma->Fill(dali2_dopplerEnergy[j],dali2_dopplerEnergy[i]);}
          }
        }
      }
      //if addback was done
      if(dali2_addback_opt==1){
        for(int i = 0;i<daliAddbackMul;i++){
          for(int j = i+1;j<daliAddbackMul;j++){
            if(dali2_doppler_addback[i]>=dali2_doppler_addback[j]){
              h_dali2_addback_gamma_gamma->Fill(dali2_doppler_addback[i],dali2_doppler_addback[j]);
            }
            else{
              h_dali2_addback_gamma_gamma->Fill(dali2_doppler_addback[j],dali2_doppler_addback[i]);
            }
            
            
            
          }
        }
      }
      
      //labr3-labr3
      double en1=0.0, en2=0.0;
      if(crystalMultLaBr3>=2)  {
        for(int i = 0;i<crystalMultLaBr3;i++){
          for(int j = i+1;j<crystalMultLaBr3;j++){
            if (labr3_crystal_dopplerEnergy[i] >= labr3_crystal_dopplerEnergy[j]) {
              en1=labr3_crystal_dopplerEnergy[i];
              en2=labr3_crystal_dopplerEnergy[j];
            }else{
              en1=labr3_crystal_dopplerEnergy[j];
              en2=labr3_crystal_dopplerEnergy[i];
            }
            h_labr3_doppler_gamma_gamma[0]->Fill(en1,en2);
            
            if(crystalMultLaBr3==3 && labr3_crystal_dopplerEnergy[0]>500.0 && labr3_crystal_dopplerEnergy[1]>500.0 && labr3_crystal_dopplerEnergy[2]>500.0){
              h_labr3_doppler_gamma_gamma_crystalMul3->Fill(en1,en2);
            }
            
            if(gammaMul>0 && gammaMul<4){
              h_labr3_doppler_gamma_gamma[gammaMul]->Fill(en1,en2);
            }
          }
        }
      }
      //labr3-dali2crystal
      if(crystalMultLaBr3>=1 && crystalMultDali2>=1)  {
        for(int i = 0;i<crystalMultLaBr3;i++){
          for(int j = 0;j<crystalMultDali2;j++){
//             if (labr3_crystal_dopplerEnergy[i] >= dali2_dopplerEnergy[j]) {
              h_labr3_dali2_doppler_gamma_gamma->Fill(labr3_crystal_dopplerEnergy[i],dali2_dopplerEnergy[j]);
//             }else {h_labr3_dali2_doppler_gamma_gamma->Fill(dali2_dopplerEnergy[j],labr3_crystal_dopplerEnergy[i]);}
          }
          
        }
      }
      //labr3-dali2crystal
      if(crystalMultLaBr3>=1 && crystalMultDali2>=1)  {
        for(int i = 0;i<crystalMultLaBr3;i++){
          
          h_labr3_dali2Total_doppler_gamma_gamma->Fill(labr3_crystal_dopplerEnergy[i],dali2_doppler_total);
        }
        
        
      
        
      }
      
      //labr 3 gamma coincidence cube:
      //use this only for testsetup with maximum 3 LaBr detectors
      if(threeGammaPlot && crystalMultLaBr3==3){
        Float_t energyLabr[3]={0.0};
        for(int i=0; i<3; i++){
          if(labr3_crystal_id[i]<4){
            energyLabr[labr3_crystal_id[i]-1]=labr3_crystal_dopplerEnergy[i];
            //printf("  LaBr crystal %i energy %f\n",labr3_crystal_id[i]-1,labr3_crystal_dopplerEnergy[i]);
          }else{
            //continue;
            printf("Warning! LaBr crystal id larger than 3! 3 Gamma Coincidence plot not correct!\n");
          }
        }
        
        h3_labr->Fill(energyLabr[0],energyLabr[1],energyLabr[2]);
        if(energyLabr[0]>485 && energyLabr[0]<535 &&
          energyLabr[1]>485 && energyLabr[1]<535 &&
          energyLabr[2]>1235){ // && energyLabr[2]<
          h3_labr_cuts->Fill(energyLabr[0],energyLabr[1],energyLabr[2]);
        }
      }
      
      
      //filling the root tree
      t2->Fill();
      
      
      //____________________________________________________
      // The variables have to be reset.
      ResetValuesAll();
      
      
      
    } //end of evnum!=nextevnum
  } // End of loop through events
  //_______________________________________________________________
  // Getting the average first interaction point:
  if(dali2_fi_opt == 1) {GetAverageInteractionPoint();}

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
  
  //---------------------------------------------------------------
  // DALI----------------------------------------------------------
  //---------------------------------------------------------------
  
  if(dali2_opt==1){
    //printf("Writing dali plots\n");
    h_dali2_crystal_mult->Write();
    h_dali2_doppler->Write();
    h_dali2_doppler_forward_angles->Write();
    h_dali2_doppler_backward_angles->Write();
    h_dali2_doppler_simple->Write();
    h_dali2_doppler_mult1and2->Write();
    
    for(int i=0;i<15;i++){
      h_dali2_doppler_mult[i]->Write();
    }
    //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++){
    //  h_doppler_crystal[i]->Write();
    //}
    for(int i=0;i<17;i++)  {
      h_dali2_doppler_layer[i]->Write();
    }
    h_dali2_doppler_total->Write();
    
    //h_doppler_max_ring->Write();
    h_dali2_energy->Write();
    h_dali2_energy_forward_angles->Write();
    h_dali2_energy_backward_angles->Write();
    //for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++){
    //  h_dali2_energy_crystal[i]->Write();
    //}
    //h_dali2_energy_beta_real->Write();
    //h_beta_after_doppler->Write();
    h_dali2_crystal_fired_doppler->Write();
    h_dali2_crystal_fired_energy->Write();
    if(dali2_fi_opt >1) {h_dali2_doppler_ave_pos->Write();}
    h_dali2_doppler_gamma_gamma->Write();
  }
  if(dali2_addback_opt==1){
    for(int i=0;i<3;i++){
      h_dali2_doppler_addback[i]->Write();
      //      for(int ppp=0;ppp<11;ppp++)  {
      //  h_doppler_addback_with_threshold[i][ppp]->Write();
      //}
    }
    h_dali2_addback_gamma_gamma->Write();
  }
  
  
  //LaBr3 histograms
  
  if(LaBr3Opt==1){
    //printf("Writing LaBr plots\n");
    
    h_labr3_energy->Write();
    h_labr3_crystal_mult->Write();
    h_labr3_crystal_fired_energy->Write();
    h_labr3_doppler_total->Write();
    h_labr3_crystal_fired_doppler->Write();
    
    //printf("  h_labr3_doppler_gamma_gamma_crystalMul3\n");
    h_labr3_doppler_gamma_gamma_crystalMul3->Write();
    
    for(int i=0;i<4;i++){
      
      if(i==0){sprintf(temp,"LaBr3 Doppler corrected Energy, GammaMul=all");}
      else{sprintf(temp,"LaBr3 Doppler corrected Energy, GammaMul=%i",i);}
      //set titles
      h_labr3_doppler_gamma[i]->SetTitle(temp);
      h_labr3_doppler_gamma_gamma[i]->SetTitle(temp);
      
      //set axis titles
      h_labr3_doppler_gamma[i]->GetXaxis()->SetTitle("Energy / keV");
      
      h_labr3_doppler_gamma_gamma[i]->GetXaxis()->SetTitle("Energy / keV");
      h_labr3_doppler_gamma_gamma[i]->GetYaxis()->SetTitle("Energy / keV");
      
      //printf("  h_labr3_doppler_gamma\n");
      h_labr3_doppler_gamma[i]->Write();
      //printf("h_labr3_doppler_gamma_gamma\n");  
      h_labr3_doppler_gamma_gamma[i]->Write();
    
    }
    
    if(threeGammaPlot){
      printf("Saving LaBr 3 gamma coincidence plot\n");
      h3_labr->GetXaxis()->SetTitle("LaBr 1");
      h3_labr->GetXaxis()->CenterTitle(true);
      h3_labr->GetXaxis()->SetTitleOffset(2);
      h3_labr->GetYaxis()->SetTitle("LaBr 2");
      h3_labr->GetYaxis()->CenterTitle(true);
      h3_labr->GetYaxis()->SetTitleOffset(2);
      h3_labr->GetZaxis()->SetTitle("LaBr 3");
      h3_labr->GetZaxis()->CenterTitle(true);
      h3_labr->GetZaxis()->SetTitleOffset(1.5);
      h3_labr->Write();
      
      h3_labr_cuts->GetXaxis()->SetTitle("LaBr 1");
      h3_labr_cuts->GetXaxis()->CenterTitle(true);
      h3_labr_cuts->GetXaxis()->SetTitleOffset(2);
      h3_labr_cuts->GetYaxis()->SetTitle("LaBr 2");
      h3_labr_cuts->GetYaxis()->CenterTitle(true);
      h3_labr_cuts->GetYaxis()->SetTitleOffset(2);
      h3_labr_cuts->GetZaxis()->SetTitle("LaBr 3");
      h3_labr_cuts->GetZaxis()->CenterTitle(true);
      h3_labr_cuts->GetZaxis()->SetTitleOffset(1.5);
      h3_labr_cuts->Write();
      printf("Saved!\n");
      
    }
   
    
    //write labr-dali coincidence histogram
    
    if(dali2_opt==1){
      //printf("Writing coincidence plots\n"); 
      h_labr3_dali2_doppler_gamma_gamma->GetXaxis()->SetTitle("LaBr3 Doppler Cor. Energy / keV");
      h_labr3_dali2_doppler_gamma_gamma->GetYaxis()->SetTitle("Dali2 Doppler Cor. Energy / keV");
      h_labr3_dali2_doppler_gamma_gamma->Write();
      
      h_labr3_dali2Total_doppler_gamma_gamma->GetXaxis()->SetTitle("LaBr3 Doppler Cor. Energy / keV");
      h_labr3_dali2Total_doppler_gamma_gamma->GetYaxis()->SetTitle("Dali2 Total Doppler Cor. Energy / keV");
      h_labr3_dali2Total_doppler_gamma_gamma->Write();
    }
  
    if(dali2_addback_opt==1){
      //printf("  with dali addback\n"); 
      h_labr3_dali2addback_gamma_gamma->GetXaxis()->SetTitle("LaBr3 Doppler Cor. Energy / keV");
      h_labr3_dali2addback_gamma_gamma->GetYaxis()->SetTitle("Dali2 Doppler Cor. Energy with addback / keV");
      h_labr3_dali2addback_gamma_gamma->Write();
    }
    
  }
  
  
  
  //overall histograms
  //printf("h_doppler_total_individual_decays\n");
  h_doppler_total_individual_decays->Write();
  
  printf("Writing tree\n");
  t2->Write();
  
  printf("Close\n");
  rootfile->Close();
  infile->Close();
  return 0;
  
  
  
} // end of RikenReconstructorDLB / main


