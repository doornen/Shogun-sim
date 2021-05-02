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
                                float gamma_energy, float beta_rec,float beta_mean_tof,float beta_average,float decay_z)  {
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
float GetThetaDetector(float gamma_det_x,float gamma_det_y,float gamma_det_z)  {
  double l1 = sqrt(gamma_det_x*gamma_det_x + gamma_det_y*gamma_det_y + gamma_det_z*gamma_det_z);

  // The theta angle:
  float theta_lab = 180. * acos(gamma_det_z/l1) / 3.14159;//cout<<"Theta lab: "<<theta_lab<<endl;
  return theta_lab;
}

//___________________________________________________________________________________________________
float GetTOFGammaNs(float gamma_det_x,float gamma_det_y,float gamma_det_z,
                    float target_x,float target_y,float target_z)  {
  float vx1 = gamma_det_x - target_x;                         //cout<<"vx1: "<<vx1<<endl;
  float vy1 = gamma_det_y - target_y;                         //cout<<"vy1: "<<vy1<<endl;
  float vz1 = gamma_det_z - target_z;                         //cout<<"vz1: "<<vz1<<endl;

  double l1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1);

  return  l1/29.97925;
} 

//_________________________________________________________________________________________________
void SortGammaEnergies(float gammaEnergy[])  {
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
//void Dali3Reconstructor()  { 
int main(int argc,char** argv) {
  const int numberOfDali3Crystals = 1001;
  float eventNumberPrevious = 0.;
  float dopplerCorrectedGammaEnergy[20];
  float tofCorrectedMeasuredTime[20];  
  float thetaDetectorLab[20];
  for(int i=0;i<20;i++)  {
    tofCorrectedMeasuredTime[20] = -999.;
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
  int opt_gamma_gamma    = 0;             // Option to include gamma-gamma matrices in the analysis. May blow up the memory use!
  float bin_reduction    = 1.;            // By which factor the bining of the 2D spectra shall be reduced
  int opt_energy_gate    = 0;             // Option to set a gate on specific energies.
  //-----------------------------------------------------------------------------------------------
  //Gates for the analysis:
  float e_rest_gate[20]       = {0.};
  float energy_gate_lower[20] = {0.};
  float energy_gate_upper[20] = {0.};
  //-----------------------------------------------------------------------------------------------
  //Reading the input file, which can change the values
  FILE *fin = fopen("./input/Dali3Reconstructor.in","r");
  while(!feof(fin))  {
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
      if(reduction_factor<1) return;
    }
    else if(strcmp(temp,"INCLUDEGAMMAGAMMA")==0)  {
      fscanf(fin,"%i %f",&opt_gamma_gamma,&bin_reduction); 
      if(bin_reduction<1.) bin_reduction = 1.;
      printf("%s %i %f\n",temp,opt_gamma_gamma,bin_reduction);
    }
    else if(strcmp(temp,"ENERGYGATE")==0)  {
      char temp2[200];
      fscanf(fin,"%i %s",&opt_energy_gate,&temp2); 
      printf("%s %i %s\n",temp,opt_energy_gate,temp2);
      if(opt_energy_gate==1)
      {
        FILE *fGateIn = fopen(temp2,"r");
        float dummy1 = 0.,dummy2 = 0.,dummy3 = 0.;
        for(int i=0;!feof(fGateIn)&& i<20;i++)  {
          fscanf(fGateIn,"%f %f %f",&dummy1,&dummy2,&dummy3);
          if(dummy1==-1.){energy_gate_lower[i]=0.;energy_gate_upper[i]=0.; break;}
          e_rest_gate[i]       = dummy1;
          energy_gate_lower[i] = dummy2;
          energy_gate_upper[i] = dummy3;
        }
      }
    }
    else if(strcmp(temp,"END")==0) break;
    else {
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
float p0[3],pvb[3],pva[3];
float energy_p;
float beta_b,beta_r,beta_a,halflife,decay_time_after_interaction;
float g0[3],gv[3];
float e_rest,e_doppler;
float theta_gamma_rest, theta_gamma_lab;
float energy_vertex_stored;
int kind_target;
//***********************************************
float total_event_num;
int gamma_det_type;   
float pos_det[3];  
float ver_rec[3];
float beta_rec;
//-----------------------------------------------
//For the Dali3 array
int dali3_flag[numberOfDali3Crystals];
float dali3_pos[3][numberOfDali3Crystals];
float dali3_energy_not_cor[numberOfDali3Crystals];
//float dali3_SizeX,dali3_SizeY,dali3_SizeZ;
float dali3_time[numberOfDali3Crystals];
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

  //Creating Beta and gamma lab spectra!!!!!!!!!!!!!!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TH2F *h_beta_before_beta_after    = new TH2F("beta_before_beta_after","beta_before_beta_after",200,0.50,0.65,200,0.50,0.65);
  TH2F *h_beta_real_beta_before     = new TH2F("beta_real_beta_before","beta_real_beta_before",200,0.50,0.65,200,0.50,0.65);
  TH2F *h_beta_real_beta_after      = new TH2F("beta_real_beta_after","beta_real_beta_after",200,0.50,0.65,200,0.50,0.65);
  
  TH1F *h_number_of_observed_decays = new TH1F("number_of_observed_decays","number_of_observed_decays",100,-0.5,99.5);
  TH2F *h_decayId_vs_energy_sum     = new TH2F("decayId_vs_energy_sum","decayId_vs_energy_sum",10,0,10,1000,0,5000);
  
  TH1F *h_gamma_fold;
  TH1F *h_gamma_fold_gamma_gate[10];  //gamma fold with gates on the energy
  TH1F *h_doppler_mult[10];
  TH1F *h_doppler;
  TH1F *h_doppler_vs_detector[]; 
  TH1F *h_energy;
  TH1F *h_dali3_time;
  TH1F *h_dali3_time_minus_tof;
  TH1F *h_dali3_time_plus_decay_time;
  TH1F *h_dali3_time_plus_decay_time_minus_tof;
  TH1F *h_dali3_time_plus_decay_time_minus_tof_energy_gated[20];
  TH1F *h_dali3_time_plus_decay_time_minus_tof_energy_and_angle_gated[9];
  TH2F *h_dali3_gamma_gamma[10];

  TH1F *h_dali3_gamma_fold_1_in_coupled_detector;
  TH2F *h_dali3_gamma_vs_angle;
  // Crystal foldiplicity, actually it is the fold
  h_dali3_gamma_fold = new TH1F("crystal_Fold","crystal_Fold",100,-0.5,99.5);
  for(int i=0;i<10;i++) {
    sprintf(temp,"gamma_fold_gamma_gate[%i]",i);
    h_dali3_gamma_fold_gamma_gate[i] = new TH1F(temp,temp,100,-0.5,99.5);
  }

   h_dali3_doppler_cor_single = new TH1F("doppler_cor_mult","doppler_cor_single",numBin,firstBin,lastBin);
  
  // All, This includes addback events. So far no routine for addback developed, so I take the energy released in the entire array
  h_dali3_doppler_cor = new TH1F("h_dali3_doppler_cor","h_dali3_doppler_cor",numBin,firstBin,lastBin);
  for(int i=0;i<8;i++)
  {
    sprintf(temp,"h_dali3_doppler_cor_theta_angle[%i]",i);
    h_dali3_doppler_cor_theta_angle[i]= new TH1F(temp,temp,numBin,firstBin,lastBin);
  }

  h_dali3_energy_sum  = new TH1F("h_dali3_energy_sum","h_dali3_energy_sum",numBin,firstBin,lastBin);

  //---------------------------------------
  //Time Investigations
  h_dali3_time                           = new TH1F("h_dali3_time","",10000,-0.01,99.99);
  h_dali3_time_minus_tof                 = new TH1F("h_dali3_time_minus_tof","",11000,-10.01,99.99);
  h_dali3_time_plus_decay_time           = new TH1F("h_dali3_time_plus_decay_time","",11000,-10.01,99.99);
  h_dali3_time_plus_decay_time_minus_tof = new TH1F("h_dali3_time_plus_decay_time_minus_tof","",11000,-10.01,99.99);
  
  if(opt_gamma_gamma==1)
  for(int i=0;i<10;i++)
  {
    sprintf(temp,"h_dali3_gamma_gamma[%i]",i);
    h_dali3_gamma_gamma[i]= new TH2F(temp,"",(int)(numBin/bin_reduction),firstBin,lastBin,(int)(numBin/bin_reduction),firstBin,lastBin);
  }
  for(int i=0;i<20;i++)
  {
    sprintf(temp,"h_dali3_time_plus_decay_time_minus_tof_energy_gated[%i]",i);
    h_dali3_time_plus_decay_time_minus_tof_energy_gated[i] = new TH1F(temp,"",11000,-10,100);
  }  
  for(int i=0;i<9;i++)
  {
    sprintf(temp,"h_dali3_time_plus_decay_time_minus_tof_energy_and_angle_gated[%i]",i);
    h_dali3_time_plus_decay_time_minus_tof_energy_and_angle_gated[i] = new TH1F(temp,"",11000,-10,100);
  } 

  h_dali3_gamma_fold_1_in_coupled_detector = new TH1F("h_dali3_gamma_fold_1_in_coupled_detector","h_dali3_gamma_fold_1_in_coupled_detector",numBin,firstBin,lastBin);

  h_dali3_gamma_vs_angle = new TH2F("h_dali3_gamma_vs_angle","h_dali3_gamma_vs_angle",numBin,firstBin,lastBin,180,0,180);

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
  //tHeader->SetBranchAddress("DALI3_SizeX",&dali3_SizeX);
  //tHeader->SetBranchAddress("DALI3_SizeY",&dali3_SizeY);
  //tHeader->SetBranchAddress("DALI3_SizeZ",&dali3_SizeZ);
  tHeader->SetBranchAddress("Beta_Resolution",&beta_res); 
  tHeader->SetBranchAddress("Pos_Det_at_Target_res",&pos_det_at_target_res); 
  tHeader->SetBranchAddress("Pos_Det_After_Target_Res",&pos_det_after_target_res);
  tHeader->GetEntry(0);
  //-----------------------------------------------------------------------------------------------

  Int_t nentries = (Int_t)(t->GetEntries()/reduction_factor);
  //Int_t nentries = t->GetEntries();
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
  //bool det_fired[1000]    = {false};

  int observedDecayNumber = 0;
//_________________________________________________________________________________________________
  cout<<"Starting the analysis"<<endl;
  for(int iii=0;iii<nentries;iii++)
  {
    if((iii+1)%1000 == 0) cout << iii+1 << "/" << nentries << " events DONE!" << endl;
    t->GetEntry(iii);

    // Fill the beta spectra:
    h_beta_before_beta_after->Fill(beta_b,beta_a);
    h_beta_real_beta_before ->Fill(beta_r,beta_b);
    h_beta_real_beta_after  ->Fill(beta_r,beta_a);

    //____________________________________________________________________
    //---------------------------------------------
    // Starting with the DALI3 analysis:-----------
    //---------------------------------------------
    //First looping through all the crystals and see, what we have.
    for(int nnn=0;nnn<numberOfDali3Crystals;nnn++)
    {
      if(dali3_flag[nnn]==1.0) 
      {
        dali3_gamma_fold++;
        
        thetaDetectorLab[dali3_gamma_fold-1] = GetThetaDetector(dali3_x[nnn],dali3_y[nnn],dali3_z[nnn]);

        //time investigations:
        gamma_tof                                    = GetTOFGammaNs(dali3_x[nnn],dali3_y[nnn],dali3_z[nnn],ver_rec_x,ver_rec_y,ver_rec_z);
        tofCorrectedMeasuredTime[dali3_gamma_fold-1] = dali3_time[nnn]+decay_time_after_interaction-gamma_tof;
        //Filling the array of measured times:
        h_dali3_time                          ->Fill(dali3_time[nnn]);
        h_dali3_time_minus_tof                ->Fill((dali3_time[nnn]-gamma_tof));
        h_dali3_time_plus_decay_time          ->Fill((dali3_time[nnn]+decay_time_after_interaction));
        h_dali3_time_plus_decay_time_minus_tof->Fill((dali3_time[nnn]+decay_time_after_interaction-gamma_tof));
         
        //summed energy
        dali3_energy_sum = dali3_energy_sum+dali3_energy_not_cor[nnn];
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
      // Filling the spectra accorging to the theta angle of the MI
      float dummyangle = GetThetaDetector(dali3_x_thatsit,dali3_y_thatsit,dali3_z_thatsit);
      for(int i=0;i<8;i++)
      {
        if(dummyangle>(10+20*i) && dummyangle<(10+20*(i+1)))
        {
          h_dali3_doppler_cor_theta_angle[i]->Fill(corrected_energy);
        }
      }
      //_________________________________________
      //Some time checks
      for(int i=0;i<20;i++)
      {
        if(energy_gate_lower[i]<corrected_energy && energy_gate_upper[i]>corrected_energy)
        {
          for(int j=0;j<20;j++)
          {
            if(tofCorrectedMeasuredTime[j]>-999.)
            {
              h_dali3_time_plus_decay_time_minus_tof_energy_gated[i]->Fill(tofCorrectedMeasuredTime[j]);
              for(int k=0;k<9;k++)
              {
                if(thetaDetectorLab[j]>20*k && thetaDetectorLab[j]<=20*k+20)
                {
                  h_dali3_time_plus_decay_time_minus_tof_energy_and_angle_gated[k]->Fill(tofCorrectedMeasuredTime[j]);
                }
              }
            }
          }
        }
      } 
      //_________________________________________
      //The gamma_fold gated on different energies:
      h_dali3_gamma_fold->Fill(dali3_gamma_fold);
      for(int i =0;i<4;i++)
      {
        if(e_rest>=(e_rest_gate[i]-1) && e_rest<=(e_rest_gate[i]+1))
        {    
          h_dali3_gamma_fold_gamma_gate[i]->Fill(dali3_gamma_fold);
          if( corrected_energy>energy_gate_lower[i] && corrected_energy<energy_gate_upper[i])
            h_dali3_gamma_fold_gamma_gate[i+4]->Fill(dali3_gamma_fold);
        }
      }
      h_dali3_energy_sum          ->Fill(dali3_energy_sum);
      h_decayId_vs_energy_sum     ->Fill(decayIDevnum,dali3_energy_sum);
      if(dali3_gamma_fold==1) 
        h_dali3_doppler_cor_single->Fill(corrected_energy);
      //_________________________________________
      //Energy versus angle:
      h_dali3_gamma_vs_angle->Fill(corrected_energy,dummyangle);
    }
    //___________________________________________
    // The variables have to be reset.
    dali3_gamma_fold = 0;
    dali3_energy_sum = 0.0;
    dali3_energy_max = -999.9;
    observedDecayNumber++;
    if(evnum != eventNumberPrevious)
    { 
      h_number_of_observed_decays->Fill(observedDecayNumber);
      //Filling the gamma-gamma spectra:
      //Sorting in the right order
      if(opt_gamma_gamma==1)
      for(int i =1;i<observedDecayNumber;i++)
      {
        if(dopplerCorrectedGammaEnergy[0]>=dopplerCorrectedGammaEnergy[i])
        {
	  h_dali3_gamma_gamma[i]->Fill(dopplerCorrectedGammaEnergy[0],dopplerCorrectedGammaEnergy[i]);
          h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[0],dopplerCorrectedGammaEnergy[i]);
        }
        else
        {
          h_dali3_gamma_gamma[i]->Fill(dopplerCorrectedGammaEnergy[i],dopplerCorrectedGammaEnergy[0]);
          h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[i],dopplerCorrectedGammaEnergy[0]);
        }
        if(i>1)
        {
          if(dopplerCorrectedGammaEnergy[1]>=dopplerCorrectedGammaEnergy[i])
               h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[1],dopplerCorrectedGammaEnergy[i]);
          else h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[i],dopplerCorrectedGammaEnergy[1]);
        }
        if(i>2)
        { 
          if(dopplerCorrectedGammaEnergy[2]>=dopplerCorrectedGammaEnergy[i])
               h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[2],dopplerCorrectedGammaEnergy[i]);
          else h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[i],dopplerCorrectedGammaEnergy[2]);
        }
        if(i>3)
          if(dopplerCorrectedGammaEnergy[3]>=dopplerCorrectedGammaEnergy[i])
               h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[3],dopplerCorrectedGammaEnergy[i]);
          else h_dali3_gamma_gamma[0]->Fill(dopplerCorrectedGammaEnergy[i],dopplerCorrectedGammaEnergy[3]);
      } 
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
  //BETA
  rootfile->cd();
  h_beta_before_beta_after                                                             ->Write();
  h_beta_real_beta_before                                                              ->Write();
  h_beta_real_beta_after                                                               ->Write();    
  //---------------------------------------------------------------
  //Dali3-----------------------------------------------------------
  //---------------------------------------------------------------
  h_number_of_observed_decays                                                          ->Write();
  h_dali3_gamma_fold                                                                   ->Write();
  for(int i=0;i<10;i++)h_dali3_gamma_fold_gamma_gate[i]                                ->Write();
  h_dali3_doppler_cor_single                                                           ->Write();
  h_dali3_doppler_cor                                                                  ->Write();
  for(int i=0;i<5;i++)h_dali3_doppler_cor_theta_angle[i]                               ->Write();
  h_dali3_energy_sum                                                                   ->Write();
  h_dali3_time                                                                         ->Write();
  h_dali3_time_minus_tof                                                               ->Write();
  h_dali3_time_plus_decay_time                                                         ->Write();
  h_dali3_time_plus_decay_time_minus_tof                                               ->Write();
  for(int i=0;i<20;i++)h_dali3_time_plus_decay_time_minus_tof_energy_gated[i]          ->Write(); 
  for(int i=0;i<9;i++) h_dali3_time_plus_decay_time_minus_tof_energy_and_angle_gated[i]->Write(); 
  if(opt_gamma_gamma==1){for(int i=0;i<10;i++)h_dali3_gamma_gamma[i]                   ->Write();
  h_dali3_gamma_vs_angle                                                               ->Write();} 
}
