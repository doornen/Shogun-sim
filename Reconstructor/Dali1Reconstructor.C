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

//_________________________________________________________________________________________________
float GetDopplerCorrectedEnergy(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                                float target_x,float target_y,float target_z,
                                float pos_det_x, float pos_det_y, float pos_det_z, 
                                float gamma_energy, float beta_rec,float beta_mean_tof,float beta_average,float decay_z)
{
  float vx1 = gamma_det_x - target_x;                                      //cout<<"vx1: "<<vx1<<endl;
  float vy1 = gamma_det_y - target_y;                                      //cout<<"vy1: "<<vy1<<endl;
  float vz1 = gamma_det_z - target_z - decay_z;                            //cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  float vx2 = pos_det_x - target_x;                                        //cout<<"vx2: "<<vx2<<endl;
  float vy2 = pos_det_y - target_y;                                        //cout<<"vy2: "<<vy2<<endl;
  float vz2 = pos_det_z - target_z - decay_z;                              //cout<<"vy2: "<<vz2<<endl;
 
  double l2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2);
          
  // The theta angle:
  double theta_lab = acos((vx1*vx2 + vy1*vy2 + vz1*vz2)/l1/l2);             //cout<<"Theta lab: "<<theta_lab<<endl;
  // the beta value:
  float beta                     = beta_average + (beta_rec-beta_mean_tof); //cout<<"beta: "<<beta<<endl;
  float doppler_corrected_energy = gamma_energy*(1-beta*cos(theta_lab))/sqrt(1.0-beta*beta);
  return doppler_corrected_energy;
}

//_________________________________________________________________________________________________
float GetThetaDetector(float gamma_det_x,float gamma_det_y,float gamma_det_z)
{
  double l1 = sqrt(gamma_det_x*gamma_det_x + gamma_det_y*gamma_det_y + gamma_det_z*gamma_det_z);

  // The theta angle:
  float theta_lab = 180. * acos(gamma_det_z/l1) / 3.14159;//cout<<"Theta lab: "<<theta_lab<<endl;
  return theta_lab;
}

//___________________________________________________________________________________________________
float GetTOFGammaNs(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                    float target_x,float target_y,float target_z)
{
  float vx1 = gamma_det_x - target_x;                         //cout<<"vx1: "<<vx1<<endl;
  float vy1 = gamma_det_y - target_y;                         //cout<<"vy1: "<<vy1<<endl;
  float vz1 = gamma_det_z - target_z;                         //cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  return  l1/29.97925;
} 

//_________________________________________________________________________________________________
void SortGammaEnergies(float gammaEnergy[])
{
  float dummyGammaEnergy[10];
  int rankOfGamma[10]={0};
  for(int i=0;i<10;i++) dummyGammaEnergy[i] = gammaEnergy[i];
  for(int i=0;i<10;i++)
  {
    if(gammaEnergy[i]>0)
    for(int j=0;j<10;j++){if(gammaEnergy[i]>dummyGammaEnergy[j]) rankOfGamma[i]++;}
  }
  for(int i=0;i<10;i++) gammaEnergy[i] = dummyGammaEnergy[rankOfGamma[i]];
}

//_________________________________________________________________________________________________
void Dali1Reconstructor()
{ 
  const int numberOfdali3Crystals = 28;
  float eventNumberPrevious = 0.;
  float dopplerCorrectedGammaEnergy[numberOfdali3Crystals];
  float tofCorrectedMeasuredTime[numberOfdali3Crystals];  
  float thetaDetectorLab[numberOfdali3Crystals];
  for(int i=0;i<numberOfdali3Crystals;i++)
  {
    tofCorrectedMeasuredTime[numberOfdali3Crystals] = -999.;
  }
  char root_in[200];
  char root_out[200];
  char temp[200];
  int numBin             = 1000;
  float firstBin         = 0;
  float lastBin          = 4000.0;
  float reduction_factor = 1.;
  float beta_average     = 0.4295;        // The average value at the moment of decay for the correction
  float beta_mean_tof    = 0.4295;        // The mean beta from Tof
  float decay_z          = 0.0;           // The lifetime of the excited state moves the average decay point along the beam axis
  float bin_reduction    = 1.;            // By which factor the bining of the 2D spectra shall be reduced
  //-----------------------------------------------------------------------------------------------
  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/Dali1Reconstructor.in","r");
  while(!feof(fin))
  {
    fscanf(fin,"%s ",temp); 
    if(strcmp(temp,"INPUTFILE")==0)
    {
      fscanf(fin,"%s",&root_in);
      printf("%s %s \n",temp,root_in);
    }
    else if(strcmp(temp,"OUTPUTFILE")==0)  
    {
      fscanf(fin,"%s ",&root_out); 
      printf("%s %s \n",temp,root_out);
    }
    else if(strcmp(temp,"SPECTRABINANDRANGE")==0)
    {
      fscanf(fin,"%i %f %f",&numBin,&firstBin,&lastBin);
      printf("%s %i %f %f\n",temp, numBin,firstBin,lastBin);
    }
    else if(strcmp(temp,"BETADOPPLERAVERAGE")==0)
    {
      fscanf(fin,"%f",&beta_average); 
      printf("%s %f \n",temp,beta_average);
    }
    else if(strcmp(temp,"BETATOFAVERAGE")==0)
    {
      fscanf(fin,"%f",&beta_mean_tof); 
      printf("%s %f \n",temp,beta_mean_tof);
    }
    else if(strcmp(temp,"DECAYPOSITIONZ")==0)
    {
      fscanf(fin,"%f",&decay_z); 
      printf("%s %f \n",temp,decay_z);
    }
    else if(strcmp(temp,"STATISTICSREDUCTIONFACTOR")==0)
    {
      fscanf(fin,"%f",&reduction_factor); 
      printf("%s %f \n",temp,reduction_factor);
      if(reduction_factor<1) return;
    }
    else if(strcmp(temp,"END")==0) break;
    else 
    {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      abort();
    }
  }

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// These variables are the same as for the eventbuilder 
float theta_low, theta_high; 
int mass,z,charge, mass_f, z_f, charge_f;
float thickness; 
float evnum;
float decayIDevnum;
float x_p0,y_p0,z_p0,x_pvb,y_pvb,z_pvb,x_pva,y_pva,z_pva;
float energy_p;
float beta_b,beta_r,beta_a,halflife,decay_time_after_interaction;
float x_g0,y_g0,z_g0,x_gv,y_gv,z_gv;
float e_rest,e_doppler;
float theta_gamma_rest, theta_gamma_lab;
float energy_vertex_stored;
int kind_target;
//***********************************************
float total_event_num;
int gamma_det_type;   
float pos_det_x,pos_det_y,pos_det_z;  
float ver_rec_x, ver_rec_y, ver_rec_z;
float beta_rec;
//-----------------------------------------------
//For the dali3 array
int dali3_flag[numberOfdali3Crystals];
float dali3_x[numberOfdali3Crystals];
float dali3_y[numberOfdali3Crystals];
float dali3_z[numberOfdali3Crystals];
float dali3_energy_not_cor[numberOfdali3Crystals];
float dali3_time[numberOfdali3Crystals];
//Input of the resolutions:
int dali3_en_res_opt;
float dali3_en_res[2];
float dali3_time_res[2];
float pos_det_at_target_res;
float pos_det_after_target_res; // Pos Resolution in mm FWHM!
float beta_res; //Resolution of beta in FWHM!
//End, same variables 
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------    

  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();
  
  TH1F *h_dali3_gamma_fold;
  TH1F *h_dali3_doppler_cor_single;
  TH1F *h_dali3_doppler_cor;
  TH1F *h_dali3_energy_sum;

  TH1F *h_dali3_doppler_det[numberOfdali3Crystals];            // doppler correction of the single crystals

  // Crystal foldiplicity, actually it is the fold
  h_dali3_gamma_fold = new TH1F("dali3_Crystal_Fold","dali3_Crystal_Fold",100,-0.5,99.5);

  //Single hits
  h_dali3_doppler_cor_single = new TH1F("dali3_doppler_cor_single","dali3_doppler_cor_single",numBin,firstBin,lastBin);
  
  // All, This includes addback events. So far no routine for addback developed, so I take the energy released in the entire array
  h_dali3_doppler_cor = new TH1F("h_dali3_doppler_cor","h_dali3_doppler_cor",numBin,firstBin,lastBin);
  h_dali3_energy_sum  = new TH1F("h_dali3_energy_sum","h_dali3_energy_sum",numBin,firstBin,lastBin);

  for(int i=0;i<numberOfdali3Crystals;i++)
  {
    sprintf(temp,"h_dali3_doppler_det[%i]",i);
    h_dali3_doppler_det[i] = new TH1F(temp,temp,numBin,firstBin,lastBin);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile *infile = new TFile(root_in,"READ");
  TTree *t = (TTree*)infile->Get("ObservedEvents");

  t->SetBranchAddress("EventNumber",&evnum);
  t->SetBranchAddress("DecayIDOfEventNumber",&decayIDevnum);
  t->SetBranchAddress("X_vertex",&x_p0);                      // X Position at the fragmentation point
  t->SetBranchAddress("Y_vertex",&y_p0);                      // Y Position at the fragmentation point
  t->SetBranchAddress("Z_vertex",&z_p0);                      // Z Position at the fragmentation point
  t->SetBranchAddress("XV_projectile_before_target",&x_pvb);  // Normalized Vector of beam before the target
  t->SetBranchAddress("YV_projectile_before_target",&y_pvb);
  t->SetBranchAddress("ZV_projectile_before_target",&z_pvb);
  t->SetBranchAddress("XV_projectile_after_target",&x_pva);   // Normalized Vector of beam after the target
  t->SetBranchAddress("YV_projectile_after_target",&y_pva);
  t->SetBranchAddress("ZV_projectile_after_target",&z_pva);
  t->SetBranchAddress("Energy_projectile",&energy_p);         // Energy of beam before the target in MeV/u
  t->SetBranchAddress("Beta_before_target",&beta_b);          // Beta before the target
  t->SetBranchAddress("Beta_real",&beta_r);                   // Beta during deexcitation	
  t->SetBranchAddress("Beta_after_target",&beta_a);           // Beta After Target
  t->SetBranchAddress("Halflife",&halflife);                  // Halflife
  t->SetBranchAddress("Decay_Time_after_interaction",&decay_time_after_interaction);
  t->SetBranchAddress("X_vertex_gamma",&x_g0);                // X Position at the gamma emmittance point
  t->SetBranchAddress("Y_vertex_gamma",&y_g0);                // Y Position at the gamma emmittance point
  t->SetBranchAddress("Z_vertex_gamma",&z_g0);                // Z Position at the gamma emmittance point
  t->SetBranchAddress("XV_gamma",&x_gv);                      // Gamma vector
  t->SetBranchAddress("YV_gamma",&y_gv);
  t->SetBranchAddress("ZV_gamma",&z_gv);
  t->SetBranchAddress("E_gamma_rest",&e_rest);		      // Energy at rest
  t->SetBranchAddress("E_gamma_Doppler",&e_doppler);          // Theta of doppler boosted gamma
  t->SetBranchAddress("Theta_gamma_rest",&theta_gamma_rest);
  t->SetBranchAddress("Theta_gamma_lab",&theta_gamma_lab);
  t->SetBranchAddress("Energy_Vertex",&energy_vertex_stored); // Energy of fragment at fragmentation 
  t->SetBranchAddress("Gamma_det_type",&gamma_det_type);
  t->SetBranchAddress("DALI3_flag",dali3_flag);               // ID of which det recorded a gamma in the event
  t->SetBranchAddress("DALI3_energy_not_cor",dali3_energy_not_cor);
  t->SetBranchAddress("DALI3_time",dali3_time);
  t->SetBranchAddress("pos_det_rec_X",&pos_det_x);            // Reconstructed position from a det after the secondary target
  t->SetBranchAddress("pos_det_rec_Y",&pos_det_y);            // Including the resolution of the detectors
  t->SetBranchAddress("pos_det_rec_Z",&pos_det_z);
  t->SetBranchAddress("Vertex_Reconstructed_X",&ver_rec_x);   // Reconstructed vertex postion, 
  t->SetBranchAddress("Vertex_Reconstructed_Y",&ver_rec_y);   // Including resolution
  t->SetBranchAddress("Vertex_Reconstructed_Z",&ver_rec_z);    
  t->SetBranchAddress("Beta_Reconstructed",&beta_rec);        // Reconstructed beta, including resolution

  //-----------------------------------------------------------------------------------------------
  //Going to the Header file, which has the information of constant values.
  infile->cd();
  TTree *tHeader = (TTree*)infile->Get("Header");
  tHeader->SetBranchAddress("Mass",&mass);                    // Beam mass
  tHeader->SetBranchAddress("Z",&z);                          // Beam z	
  tHeader->SetBranchAddress("Charge",&charge);                // Beam charge	
  tHeader->SetBranchAddress("Mass_Fragment",&mass_f);         // Fragment mass
  tHeader->SetBranchAddress("Z_Fragment",&z_f);	              // Fragment z
  tHeader->SetBranchAddress("Charge_Fragment",&charge_f);     // Fragment charge
  tHeader->SetBranchAddress("TargetKind",&kind_target);
  tHeader->SetBranchAddress("TargetThicknessCM",&thickness);
  tHeader->SetBranchAddress("Total_Event_Number",&total_event_num);
  tHeader->SetBranchAddress("ThetaLow",&theta_low);
  tHeader->SetBranchAddress("ThetaHigh",&theta_high);
  tHeader->SetBranchAddress("DALI3_en_res_opt",&dali3_en_res_opt);
  tHeader->SetBranchAddress("DALI3_en_res",dali3_en_res);
  tHeader->SetBranchAddress("DALI3_time_res",dali3_time_res);
  tHeader->SetBranchAddress("DALI3_x",dali3_x);
  tHeader->SetBranchAddress("DALI3_y",dali3_y);
  tHeader->SetBranchAddress("DALI3_z",dali3_z);
  tHeader->SetBranchAddress("Beta_Resolution",&beta_res); 
  tHeader->SetBranchAddress("Pos_Det_at_Target_res",&pos_det_at_target_res); 
  tHeader->SetBranchAddress("Pos_Det_After_Target_Res",&pos_det_after_target_res);
  tHeader->GetEntry(0);
  //-----------------------------------------------------------------------------------------------

  //Int_t nentries = (Int_t)(t->GetEntries()/reduction_factor);
  Int_t nentries = (Int_t)(t->GetEntries());
  cout <<"Entries: "<< nentries << endl;

  //Variables that are determined event by event:
  int dali3_gamma_fold    = 0;
  float dali3_energy_sum  = 0.;
  float corrected_energy  = 0.;
  float dali3_energy_max  = -999.9;
  
  float dali3_x_thatsit   = 0.;
  float dali3_y_thatsit   = 0.;
  float dali3_z_thatsit   = 0.;
  float gamma_tof         = 0.;

  int observedDecayNumber = 0;
//_________________________________________________________________________________________________
  cout<<"Starting the analysis"<<endl;
  for(int iii=0;iii<nentries;iii++)
  {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    t->GetEntry(iii);

    //____________________________________________________________________
    //---------------------------------------------
    // Starting with the dali3 analysis:-----------
    //---------------------------------------------
    //First looping through all the crystals and see, what we have.
    for(int nnn=0;nnn<numberOfdali3Crystals;nnn++)
    {
      //cout<<" dali flag:"<<dali3_flag[nnn]<<endl;
      if(dali3_flag[nnn]==1.0) 
      {
        dali3_gamma_fold++;
        //cout<<"gamma fold true"<<iii<<endl;
        thetaDetectorLab[dali3_gamma_fold-1] = GetThetaDetector(dali3_x[nnn],dali3_y[nnn],dali3_z[nnn]);
         
        //summed energy
        dali3_energy_sum = dali3_energy_sum+dali3_energy_not_cor[nnn];
        //cout<<"dali3_energy_sum: "<< dali3_energy_sum<<endl;
        // getting the MI
        if(dali3_energy_not_cor[nnn]>dali3_energy_max)
        {
          dali3_energy_max = dali3_energy_not_cor[nnn];
          dali3_x_thatsit  = dali3_x[nnn];
          dali3_y_thatsit  = dali3_y[nnn];
          dali3_z_thatsit  = dali3_z[nnn];
        }
      }
    }
    //____________________________________________________________________
    //Performing the analysis if 
    if(dali3_gamma_fold>0)
    {
      //Performing the Doppler-correction from the angles
      //of the detector that registered the highest gamma-ray energy.
      corrected_energy = GetDopplerCorrectedEnergy(dali3_x_thatsit,dali3_y_thatsit,dali3_z_thatsit,
                               ver_rec_x,ver_rec_y,ver_rec_z, 
                               pos_det_x,pos_det_y,pos_det_z, 
                               dali3_energy_sum,beta_rec,beta_mean_tof,beta_average,decay_z);

      if(evnum==eventNumberPrevious) dopplerCorrectedGammaEnergy[observedDecayNumber+1] = corrected_energy;

      h_dali3_doppler_cor->Fill(corrected_energy);
      //_________________________________________
      // Checking the spectra of the first 10 detectors:
      for(int i=0;i<numberOfdali3Crystals;i++)
      {
        if(dali3_flag[i]==1.0)
        {
          h_dali3_doppler_det[i]->Fill(corrected_energy);
        } 
      }
      //_________________________________________
      // Filling the spectra accorging to the theta angle of the MI
      float dummyangle = GetThetaDetector(dali3_x_thatsit,dali3_y_thatsit,dali3_z_thatsit);
      h_dali3_gamma_fold->Fill(dali3_gamma_fold);

      h_dali3_energy_sum          ->Fill(dali3_energy_sum);
      if(dali3_gamma_fold==1) 
        h_dali3_doppler_cor_single->Fill(corrected_energy);
    }
    //___________________________________________
    // The variables have to be reset.
    dali3_gamma_fold = 0;
    dali3_energy_sum = 0.0;
    dali3_energy_max = -999.9;
    observedDecayNumber++;
    if(evnum != eventNumberPrevious)
    { 
      observedDecayNumber = 0;
      for(int i =0;i<20;i++)
      {
        dopplerCorrectedGammaEnergy[i] = 0.;
        tofCorrectedMeasuredTime[i]    = -999.;
      }
      dopplerCorrectedGammaEnergy[observedDecayNumber] = corrected_energy;
    }
    eventNumberPrevious = evnum;
  }
  //______________________________________________________________________
  //Writing into the root file:
  rootfile->cd();
  h_dali3_gamma_fold                                                                   ->Write();
  h_dali3_doppler_cor_single                                                           ->Write();
  h_dali3_doppler_cor                                                                  ->Write();
  h_dali3_energy_sum                                                                   ->Write();
  for(int i=0;i<numberOfdali3Crystals;i++)h_dali3_doppler_det[i]                       ->Write(); 
}
