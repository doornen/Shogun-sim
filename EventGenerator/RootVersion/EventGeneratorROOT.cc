#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom2.h"
#include "TVector3.h"
#include "TMath.h"

#include <iostream>
#define PI 3.14159265


using namespace std;
using namespace TMath;

                                                                               // Defining all the variables
                                                                               // Initialzing the variables that need an initial value
                                                                               // order: First what kind of particle, second its specific parameter
int    gBeamMass                   = 1;
int    gBeamZet                    = 1;
int    gBeamCharge                 = 1;
float  gBeamEnergyMean             = 100.;
float  gBeamEnergySigma            = 1.;                                       // Energy of beam A MeV
float  gBeamX                      = 0.; 
float  gBeamY                      = 0.; 
float  gBeamXSigma                 = 0.1;
float  gBeamYSigma                 = 0.1;
float  gBeamTheta                  = 0.;
float  gBeamThetaSigma             = 0.1;
float  gBeamPhiMin                 = 0.;
float  gBeamPhiMax                 = 360.;
float  gBeamInteractionPointX      = 0.;
float  gBeamInteractionPointY      = 0.;
float  gBeamInteractionPointZ      = 0.0;
float  gBeamVecBeforeX             = 0.0; 
float  gBeamVecBeforeY             = 0.0;
float  gBeamVecBeforeZ             = 0.;
float  gBeamEnergyVertexStored     = 0.;
float  gBeamEnergy                 = 0.;
//_______________________________________________
int    gFragmentMass               = 1;
int    gFragmentZet                = 1;
int    gFragmentCharge             = 1; 
float  gFragmentVecAfterX          = 0.;
float  gFragmentVecAfterY          = 0.;
float  gFragmentVecAfterZ          = 0.;
float  gFragmentEnergyAfterTarget  = 0.;
//_______________________________________________
int    gTargetKind                 = 2;                                        // target kind 1 Au 2 Be
float  gTargetSizeX                = 1.;                                       // Dimension of the target in cm for X,Y and mg/cm^2 for Z
float  gTargetSizeY                = 1.;
float  gTargetThickness            = 1.;   
double gTargetDensity[7]           = {0.};
int    gTargetBroadeningOption     = 0;                                        // Option to take also the angular broadening of the beam due to fragmentation
float  gTargetBroadeningThetaSigma = 0.0;
//_______________________________________________
char   gTemp[200]                  = "dummy";
char   gGammaIn[200]               = "dummy.in";
char   gRootOut[200]               = "dummy.root";
char   gTempParticle[200]          = "dummy";
//_______________________________________________
float  gGammaThetaLow              = 0.;
float  gGammaThetaHigh             = 0.;                                       // Angle covered
float  gGammaDecayPointX           = 0.;
float  gGammaDecayPointY           = 0.;
float  gGammaDecayPointZ           = 0.;
float  gGammaVecX                  = 0.;
float  gGammaVecY                  = 0.;
float  gGammaVecZ                  = 0.;
float  gGammaEnergyRest            = 0.;
float  gGammaEnergyDoppler         = 0.;
float  gGammaThetaRest             = 0.;
float  gGammaThetaLab              = 0.;
//_______________________________________________
int    gBorrelOption               = 0;
float  gBindingEnergyNucleon       = 8.;
int    gGoldhaberOption            = 0;
float  gGoldhaberSigma0            = 90.;
float  gDefaultCutValue            = 0.1;
//_______________________________________________
double gSumHomo                    = 0.;
double gSumHomoUpToNow[18001]      = {0.};
//_______________________________________________
float  gEvNum                      = 0.;
float  gDecayIDEvNum               = 0.;
int    gNumEventTotal              = 10000; 
int    gDeltaMass                  = 0;
int    gDeltaZet                   = 0;            
float  gBetaBefore                 = 0.;
float  gBetaReal                   = 0;
float  gBetaAfter                  = 0.;
float  gHalflife                   = 0.;
float  gDecayTimeAfterInteraction  = 0.;
//_______________________________________________
int    gdEdXTableOption             = 0;
char   gdEdXTableInputBeam[200]     = "./dEdXTables/dummy.in";
char   gdEdXTableInputFragment[200] = "./dEdXTables/dummy.in";
float  gBeamdEdX[2][50]             = {{0.0}};                                 // The first value gives the energy, the second the actual value;
float  gFragmentdEdX[2][50]         = {{0.0}};
//_______________________________________________
struct Level
{
  int ID;
  float excitationProbability;
  float beginOfTotalExcitationProbability;
  float endOfTotalExcitationProbability;
  float energy;
  double halflife;
  int numberOfDecayBranches;
  int decayIntoLevelId[5];
  float decayIntoLevelProbability[5];
  float totalDecayProbability;
  float beginOfDecayProbability[5];
  float endOfDecayProbability[5];
} SLevel[100];
//_______________________________________________
float  gTotalLevelExcitationProbability = 0.;
int    gNumberOfLevels                  = 0; 
int    gNumberOfExcitedLevels           = 0; 
//_______________________________________________                              // The ROOT-Tree
TTree *t;
TTree *tHeader;                                                                // The header for information that doesn't change.
TRandom2 *fRandom;                                                             // The generator for random numbers
//_______________________________________________
void ReadInputFile();
void ReadGammaFile();
void RunSimulation();
void ReaddEdXTable();
double EnergyDetermination();
TVector3 GetVector(TVector3 vectorIn, float theta);
//_________________________________________________________________________________________________
int main(int argc,char** argv) 
{
  ReadInputFile();                                                             // Reading the input file 
  ReadGammaFile();                                                             // Reading the gamma-ray decay file

  gFragmentZet    = gBeamZet    - gDeltaZet;                                   // Calculating the fragment mass and charge
  gFragmentMass   = gBeamMass   - gDeltaMass;
  gFragmentCharge = gBeamCharge - gDeltaZet;

  gTargetDensity[0] = 0.0;                                                     // Vacuum
  gTargetDensity[1] = 19.3;                                                    // Au
  gTargetDensity[2] = 1.848;                                                   // Be
  gTargetDensity[3] = 1.8;                                                     // C
  gTargetDensity[4] = 7.874;                                                   // Fe
  gTargetDensity[5] = 11.35;                                                   // Pb
  gTargetDensity[6] = 0.07099;                                                 // LH2
  gTargetThickness = gTargetThickness/gTargetDensity[gTargetKind]/1000.0;      // Calculating the target thickness in cm!

  gSumHomo = 0.0;                                                              // Making a uniform theta distribution
  for(int m2=(int)(gGammaThetaLow*100);m2<(int)(gGammaThetaHigh*100+1);m2++)
  {					
    double angleHomo    = m2/100.0;
    gSumHomo            = gSumHomo + 1.0 * Sin((angleHomo*PI/180.));
    gSumHomoUpToNow[m2] = gSumHomo;
  }

  TFile *rootFile = new TFile(gRootOut,"RECREATE");
  rootFile->cd();

  t = new TTree("Events","Events");
  t->Branch("EventNumber",&gEvNum,"evnum/F");
  t->Branch("DecayIDOfEventNumber",&gDecayIDEvNum,"decayIDevnum/F");
  t->Branch("X_vertex",&gBeamInteractionPointX,"x_p0/F");                      // X Position at the fragmentation point
  t->Branch("Y_vertex",&gBeamInteractionPointY,"y_p0/F");                                      
  t->Branch("Z_vertex",&gBeamInteractionPointZ,"z_p0/F");                                      
  t->Branch("XV_projectile_before_target",&gBeamVecBeforeX,"x_pvb/F");         // Normalized Vector of beam before the target
  t->Branch("YV_projectile_before_target",&gBeamVecBeforeY,"y_pvb/F");
  t->Branch("ZV_projectile_before_target",&gBeamVecBeforeZ,"z_pvb/F");
  t->Branch("XV_projectile_after_target",&gFragmentVecAfterX,"x_pva/F");       // Normalized Vector of beam after the target
  t->Branch("YV_projectile_after_target",&gFragmentVecAfterY,"y_pva/F");
  t->Branch("ZV_projectile_after_target",&gFragmentVecAfterZ,"z_pva/F");
  t->Branch("Energy_projectile",&gBeamEnergy,"energy_p/F");                    // Energy of beam before the target in MeV/u
  t->Branch("Beta_before_target",&gBetaBefore,"beta_b/F");                     // Beta before the target
  t->Branch("Beta_real",&gBetaReal,"beta_r/F");                                // Beta at deexcitation	
  t->Branch("Beta_after_target",&gBetaAfter,"beta_a/F");                       // Beta After Target
  t->Branch("Halflife",&gHalflife,"halflife/F");                               // Halflife
  t->Branch("Decay_Time_after_interaction",&gDecayTimeAfterInteraction,"decay_time_after_interaction/F");
  t->Branch("X_vertex_gamma",&gGammaDecayPointX,"x_g0/F");                     // X Position at the gamma emmittance point
  t->Branch("Y_vertex_gamma",&gGammaDecayPointY,"y_g0/F");                                
  t->Branch("Z_vertex_gamma",&gGammaDecayPointZ,"z_g0/F");                               
  t->Branch("XV_gamma",&gGammaVecX,"x_gv/F");		                       // Gamma vector
  t->Branch("YV_gamma",&gGammaVecY,"y_gv/F");
  t->Branch("ZV_gamma",&gGammaVecZ,"z_gv/F");
  t->Branch("E_gamma_rest",&gGammaEnergyRest,"e_rest/F");		       // Energy at rest
  t->Branch("E_gamma_Doppler",&gGammaEnergyDoppler,"e_doppler/F");             // Theta of doppler boosted gamma
  t->Branch("Theta_gamma_rest",&gGammaThetaRest,"theta_gamma_rest/F");
  t->Branch("Theta_gamma_lab",&gGammaThetaLab,"theta_gamma_lab/F");
  t->Branch("Energy_Vertex",&gBeamEnergyVertexStored,"energy_vertex_stored/F");// Energy of fragment at fragmentation
                                                                               
  tHeader  = new TTree("Header","Header");                                     // The header tree gets the information that is not changed.
  tHeader->Branch("Mass",&gBeamMass,"mass/I");                                 // Beam mass
  tHeader->Branch("Z",&gBeamZet,"z/I");                                        // Beam z	
  tHeader->Branch("Charge",&gBeamCharge,"charge/I");                           // Beam charge	
  tHeader->Branch("Mass_Fragment",&gFragmentMass,"mass_f/I");                  // Fragment mass
  tHeader->Branch("Z_Fragment",&gFragmentZet,"z_f/I");	                       // Fragment z
  tHeader->Branch("Charge_Fragment",&gFragmentCharge,"charge_f/I");            // Fragment charge
  tHeader->Branch("TargetKind",&gTargetKind,"kind_target/I");
  tHeader->Branch("TargetThicknessCM",&gTargetThickness,"thickness/F");
  tHeader->Branch("ThetaLow",&gGammaThetaLow,"theta_low/F");                   // Theta angle covered by simulation
  tHeader->Branch("ThetaHigh",&gGammaThetaHigh,"theta_high/F");
  tHeader->Fill();
 
  fRandom = new TRandom2();                                                    // The generator for random numbers
  RunSimulation();                                                             // Starting the simulation

  tHeader->Write();                                                            // Writing the the header information to the root tree
  t->Write();

  return 0;
}

//_________________________________________________________________________________________________
void ReadInputFile()                                                           // The Input files are of free type.
{
  FILE *fin = fopen("../input/EventGenerator.in","r");
  while(!feof(fin))
  {
    fscanf(fin,"%s ",gTemp); 
    if(strcmp(gTemp,"BEAMISOTOPE")==0)          
    {
      fscanf(fin,"%i %i %i ",&gBeamMass,&gBeamZet,&gBeamCharge); 
      printf("%s %i %i %i\n",gTemp,gBeamMass,gBeamZet,gBeamCharge);
    }
    else if(strcmp(gTemp,"BEAMENERGY")==0)
    {
      fscanf(fin,"%f %f",&gBeamEnergyMean,&gBeamEnergySigma); 
      gBeamEnergySigma = gBeamEnergySigma/2.35;
      printf("%s %f %f \n",gTemp,gBeamEnergyMean,gBeamEnergySigma);
    }
    else if(strcmp(gTemp,"BEAMPOSITION")==0)    
    {
      fscanf(fin,"%f %f %f %f",&gBeamX,&gBeamXSigma,&gBeamY,&gBeamYSigma);
      gBeamXSigma = gBeamXSigma/2.35;gBeamYSigma = gBeamYSigma/2.35;
      printf("%s %f %f %f %f \n",gTemp,gBeamX,gBeamXSigma,gBeamY,gBeamYSigma);
    }
    else if(strcmp(gTemp,"BEAMANGLE")==0)       
    {
      fscanf(fin,"%f %f %f %f",&gBeamTheta,&gBeamThetaSigma,&gBeamPhiMin,&gBeamPhiMax);
      gBeamThetaSigma = gBeamThetaSigma/2.35;
      printf("%s %f %f %f %f \n",gTemp,gBeamTheta,gBeamThetaSigma,gBeamPhiMin,gBeamPhiMax);
    }
    else if(strcmp(gTemp,"TARGET")==0)   
    {
      fscanf(fin,"%i %f %f %f",&gTargetKind,&gTargetSizeX,&gTargetSizeY,&gTargetThickness);
      if(gTargetKind<0 && gTargetKind>6)
      {
        cout<<"Could not read your input keyword. Aborting program."<<endl; 
        abort();
      }
      printf("%s %i %f %f %f \n",gTemp,gTargetKind,gTargetSizeX,gTargetSizeY,gTargetThickness);
    }   
    else if(strcmp(gTemp,"TARGETANGULARBROADENING")==0)
    {
      fscanf(fin,"%i %f",&gTargetBroadeningOption,&gTargetBroadeningThetaSigma);
      gTargetBroadeningThetaSigma = gTargetBroadeningThetaSigma/2.35;
      printf("%s %i %f \n",gTemp,gTargetBroadeningOption,gTargetBroadeningThetaSigma);
    }
    else if(strcmp(gTemp,"MASSCHANGE")==0)      
    {
      fscanf(fin,"%i %i",&gDeltaMass,&gDeltaZet);
      printf("%s %i %i \n",gTemp,gDeltaMass,gDeltaZet);
    }
    else if(strcmp(gTemp,"BORREL")==0)    
    {
      fscanf(fin,"%i %f",&gBorrelOption,&gBindingEnergyNucleon);
      printf("%s %i %f \n",gTemp,gBorrelOption,gBindingEnergyNucleon);
    }
    else if (strcmp(gTemp,"GOLDHABER")==0)
    {
      fscanf(fin,"%i %f",&gGoldhaberOption,&gGoldhaberSigma0);
      printf("%s %i %f \n",gTemp,gGoldhaberOption,gGoldhaberSigma0);
    }
    else if(strcmp(gTemp,"GAMMAINPUT")==0)  
    {
      fscanf(fin,"%s",&gGammaIn);
      printf("%s %s \n",gTemp,gGammaIn);
    }
    else if(strcmp(gTemp,"THETARANGE")==0)      
    {
      fscanf(fin,"%f %f",&gGammaThetaLow,&gGammaThetaHigh);
      printf("%s %f %f \n",gTemp,gGammaThetaLow,gGammaThetaHigh);
    }
    else if(strcmp(gTemp,"NUMBEROFEVENTS")==0)  
    {
      fscanf(fin,"%i ",&gNumEventTotal);
      printf("%s %i\n",gTemp,gNumEventTotal);
    }
    else if(strcmp(gTemp,"DEFAULTCUTVALUE")==0)  
    {
      fscanf(fin,"%f ",&gDefaultCutValue);
      printf("%s %f\n",gTemp,gDefaultCutValue);
    }
    else if(strcmp(gTemp,"OUTPUTFILE")==0)
    {
      fscanf(fin,"%s ",&gRootOut); 
      printf("%s %s \n",gTemp,gRootOut);
    }
    else if(strcmp(gTemp,"DEDXTABLE")==0)
    {
      fscanf(fin,"%i %s %s",&gdEdXTableOption,&gdEdXTableInputBeam,&gdEdXTableInputFragment); 
      printf("%s %i %s %s\n",gTemp,gdEdXTableOption,gdEdXTableInputBeam,gdEdXTableInputFragment);
      if(gdEdXTableOption==1) ReaddEdXTable();
    }
    else if(strcmp(gTemp,"END")==0) break;
    else 
    {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      abort();
    }
  }
  fclose(fin);
}

//_________________________________________________________________________________________________
void ReadGammaFile()
{
  for(int i=0;i<100;i++)                                                        // At first, initializing the struct of SLevel[100]
  {
    SLevel[i].ID = -1;
    SLevel[i].excitationProbability             = 0.;
    SLevel[i].beginOfTotalExcitationProbability = 0.;
    SLevel[i].endOfTotalExcitationProbability   = 0.;
    SLevel[i].energy                            = 0.;
    SLevel[i].halflife                          = 0.;
    SLevel[i].numberOfDecayBranches             = 0;
    SLevel[i].totalDecayProbability             = 0.;
    for(int j=0;j<5;j++)
    {
      SLevel[i].decayIntoLevelId[j]             = 0;
      SLevel[i].decayIntoLevelProbability[j]    = 0.;
      SLevel[i].beginOfDecayProbability[j]      = 0.;
      SLevel[i].endOfDecayProbability[j]        = 0.;
    }
  }
  int   levelID,decayLevelID; 
  float levelExcitationProbability, levelEnergy, levelHalflife;
  float branchRatio;
  
  FILE *gammaIn = fopen(gGammaIn,"r");
  while(!feof(gammaIn))
  {
    fscanf(gammaIn,"%s ",&gTemp);
    if(strcmp(gTemp,"LEVEL")==0) 
    {
      fscanf(gammaIn,"%i %f %f %f",&levelID,&levelExcitationProbability,&levelEnergy,&levelHalflife);
      if(levelID < 0 || levelExcitationProbability < 0. || levelEnergy < 0.)   // Check if the input values are greater than zero
      {
        cout<<"At least one of your LEVEL input values is smaller than zero. Aborting program."<<endl; 
        abort();
      }
      if(SLevel[levelID].ID != -1)                                             // Check if the level has been assigned already
      {
        cout<<"This LEVEL has been assigned already. Aborting program."<<endl; 
        abort();
      }
      SLevel[levelID].ID = levelID;
      SLevel[levelID].excitationProbability = levelExcitationProbability;
                                                                               // Determine the range of this level within the total excitation probabilty
                                                                               // To be used later for the determination of the initial state of excitation
      SLevel[levelID].beginOfTotalExcitationProbability = gTotalLevelExcitationProbability;
      gTotalLevelExcitationProbability = gTotalLevelExcitationProbability + levelExcitationProbability;
      
      cout<<"gTotalLevelExcitationProbability: "<<gTotalLevelExcitationProbability<<endl;
      
      SLevel[levelID].endOfTotalExcitationProbability = gTotalLevelExcitationProbability;
      SLevel[levelID].energy=levelEnergy;
      SLevel[levelID].halflife=levelHalflife;
      gNumberOfLevels++;
      if(levelEnergy>0.) gNumberOfExcitedLevels++;
    } 
    else if(strcmp(gTemp,"DECAY")==0) 
    {
      fscanf(gammaIn,"%i %i %f",&levelID,&decayLevelID,&branchRatio);
      int branchID = SLevel[levelID].numberOfDecayBranches;                    // Setting the maximum of decay branches to five:
      if(branchID>4)
      {
        cout<<"This LEVEL has already five decay branches. Aborting program."<<endl; 
        abort();
      }
      SLevel[levelID].decayIntoLevelId[branchID] = decayLevelID;
      SLevel[levelID].decayIntoLevelProbability[branchID] = branchRatio;
                                                                               // Determine the range of this decay within all decay branches
                                                                               // To be used later for the determination of the decay branch
      SLevel[levelID].beginOfDecayProbability[branchID] = SLevel[levelID].totalDecayProbability;
      SLevel[levelID].totalDecayProbability = SLevel[levelID].totalDecayProbability + branchRatio;
      cout<<" Total Decay Probability of Level "<<SLevel[levelID].ID<<": "<< SLevel[levelID].totalDecayProbability<<endl;
      SLevel[levelID].endOfDecayProbability[branchID] = SLevel[levelID].totalDecayProbability;
      SLevel[levelID].numberOfDecayBranches++;
    }
    else if(strcmp(gTemp,"END")==0) break;
    else 
    {
      cout<<"Could not read your input keyword of the gamma-ray decay file. Aborting program."<<endl; 
      abort();
    }
  }
  fclose(gammaIn);
} 

//_________________________________________________________________________________________________
void ReaddEdXTable()
{
  cout<<"Building the energy loss table..."<<endl;
  float dummy1,dummy2;
  int i=0,j=0;
  FILE *tableInBeam = fopen(gdEdXTableInputBeam,"r");
  while(!feof(tableInBeam))
  {
    fscanf(tableInBeam,"%s ",gTemp); 
    if(strcmp(gTemp,"dEdX")==0)   
    {
      fscanf(tableInBeam,"%f %f",&dummy1,&dummy2);
      gBeamdEdX[0][i] = dummy1;
      gBeamdEdX[1][i] = dummy2;
      cout<< dummy1<<" "<<dummy2<<endl;
      i++;
    }
    if(strcmp(gTemp,"END")==0)  break;
  }
  fclose(tableInBeam);
  
  FILE *tableInFragment = fopen(gdEdXTableInputFragment,"r");
  while(!feof(tableInFragment))
  {
    fscanf(tableInFragment,"%s ",gTemp); 
    if(strcmp(gTemp,"dEdX")==0)   
    {
      fscanf(tableInFragment,"%f %f",&dummy1,&dummy2);
      gFragmentdEdX[0][j] = dummy1;
      gFragmentdEdX[1][j] = dummy2;
      cout<< dummy1<<" "<<dummy2<<endl;
      j++;
    }
    if(strcmp(gTemp,"END")==0)  break;  
  }
  fclose(tableInFragment);
}
                            
//_________________________________________________________________________________________________
void RunSimulation()
{
  cout<<"Running the simulation"<<endl;

  double gammaTemp;
  int i=0;
  while(i<gNumEventTotal)
  {
    if(i%10000==0) cout << i << "/" << gNumEventTotal << "   DONE!" << endl;
    int restartLoop = 0;

    gDecayTimeAfterInteraction = 0.;                                           // Reseting the decay ID and the time of the event to zero;
    gDecayIDEvNum              = 0.;
    double beamOriginX         = 999.9;
    double beamOriginY         = 999.9;
    double beamOriginZ         = 999.9;
                                                                               // Defining the beam's starting position
    while(beamOriginX>gTargetSizeX || beamOriginX<(-1.0*gTargetSizeX) || beamOriginY>gTargetSizeY || beamOriginY<(-1.*gTargetSizeY))
    {
      beamOriginX = fRandom->Gaus(gBeamX,gBeamXSigma);
      beamOriginY = fRandom->Gaus(gBeamY,gBeamYSigma);
      beamOriginZ = -1.0*gTargetThickness/2.0;                                 // cout<<"Beam origin Z: "<<beam_origin_z<<endl;
    }
    double theta = -10.0;                                                      // Defining the beam's starting vector.
    double phi   = -10.0;	
    while(theta< 0.0 || theta > 180.0 || phi < 0.0 || phi>360.0)
    {  
      theta = fRandom->Gaus(gBeamTheta, gBeamThetaSigma);                      // Gaussian distribution from input values
      phi   = gBeamPhiMin + fRandom->Rndm() * (gBeamPhiMax-gBeamPhiMin);  
    }
    phi   = phi * PI/180.;                                                     // cout<<"phi = "<<phi<<endl;
    theta = theta * PI/180.;                                                   // cout<<"theta = "<<theta<<endl;
                                                                               // Normalized vector of the projectile before the target					
    gBeamVecBeforeX = Sin(theta)*Cos(phi);                                     // cout<<"x_pvb: "<<gBeamVecBeforeX<<endl;
    gBeamVecBeforeY = Sin(theta)*Sin(phi);                                     // cout<<"y_pvb: "<<gBeamVecBeforeY<<endl;
    gBeamVecBeforeZ = Cos(theta);                                              // cout<<"z_pvb: "<<gBeamVecBeforeZ<<endl;
    
    gBeamEnergy     = fRandom->Gaus(gBeamEnergyMean,gBeamEnergySigma);         // Energy of the projectile
    gammaTemp       = (931.494+gBeamEnergy)/931.494;                           // Gamma before the target
    gBetaBefore     = sqrt(1.0 - 1.0/gammaTemp/gammaTemp);                     // Beta before the target
    
    double energyTotal = gBeamEnergy * (double)gBeamMass;                      // Total beam energy before the target
 
    //___________________________________________
    //*****************Starting Beam*************
    //___________________________________________ 
    // First, getting the fragmentation point to know how far the energy loss has
    // to be calculated for the secondary beam
    gBeamInteractionPointX = beamOriginX;	
    gBeamInteractionPointY = beamOriginY;
    gBeamInteractionPointZ = -1.0*gTargetThickness/2.0 + fRandom->Rndm() * gTargetThickness;

    double energyVertex;
    double energyHalfStep = (gBeamdEdX[0][0] - gBeamdEdX[0][1])/2;
       
    for(int j=0;j<50;j++)                                                      // Have to calculate the energy loss:
    {
      if(gBeamdEdX[0][j]+energyHalfStep<gBeamEnergy && j>0 && gBeamdEdX[0][j-1]+energyHalfStep>gBeamEnergy)
      {
        float dummyEnergy1  = energyTotal;                                     // cout<<"condition fullfilled: "<<j<<endl;
        float dummyEnergy2  = gBeamdEdX[0][j] + energyHalfStep;
        float dummyPosition = beamOriginZ;
        int dummyCounter    = j;
        while(dummyPosition<gBeamInteractionPointZ)
        {
          dummyPosition = dummyPosition + 0.1*gDefaultCutValue;                // Convert into cm because beamdEdX is given in cm
          dummyEnergy1  = dummyEnergy1 - gBeamdEdX[1][dummyCounter]*0.1*gDefaultCutValue * gTargetDensity[gTargetKind]*1000;
                                                                               // cout<<"dummyEnergy1: "<<dummyEnergy1<<endl;
          if(dummyEnergy1/(double)gBeamMass < dummyEnergy2)
          {
            dummyCounter++;
            dummyEnergy2 = gBeamdEdX[0][dummyCounter] + energyHalfStep;
          }
        }
        energyVertex = dummyEnergy1;
      }
    }
    //_____________________________________________________________________________________________
    // Energy of fragment at fragmentation point:
    // Calculating the energy spread due to goldhaber
    if(gGoldhaberOption==1)
    {	  
      double goldhaberSigma   = Sqrt(gGoldhaberSigma0 * gGoldhaberSigma0 * gFragmentMass * (gBeamMass-gFragmentMass)/(gBeamMass-1)); 
      double momentumSigma    = fRandom->Gaus(0,goldhaberSigma); 	
   
      gammaTemp               = (931.494+energyVertex/gBeamMass) / 931.494;
      float betaTemp          = Sqrt(1.0-1.0/gammaTemp/gammaTemp);
      double energySigma      = momentumSigma * gammaTemp * betaTemp;

      gBeamEnergyVertexStored = energyVertex * gFragmentMass/gBeamMass + energySigma;
      energyVertex            = gBeamEnergyVertexStored;                       // Energy at fragmentation point
    }
    else
    {
      gBeamEnergyVertexStored = energyVertex * gFragmentMass/gBeamMass;
      energyVertex            = gBeamEnergyVertexStored;
    }
	   
    double betaVertex;
    if(gBorrelOption==1)
    {                                                                          // Fragment velocity change
      double velocityChange   = Sqrt(1-((gBindingEnergyNucleon*(gBeamMass-gFragmentMass))/(energyVertex)));  
                                                                               // Relative change of velocity (Borrel et al.)
      gammaTemp               = (931.494 + energyVertex/gFragmentMass)/931.494;
                                                                               // Initial gamma value
      betaVertex              = velocityChange * Sqrt(1.0 - 1.0/gammaTemp/gammaTemp);              
                                                                               // Changed absolute beta
      gammaTemp               = Sqrt(1/(1 - betaVertex * betaVertex));  // Changed gamma value
      energyVertex            = (931.494 * gammaTemp - 931.494)*gFragmentMass; // Changed energy from gamma value
      gBeamEnergyVertexStored = energyVertex;                                  // cout<<"Energy vertex including Borrel: "<<energy_vertex<<endl;  
    }
    else
    {
      gammaTemp               = (931.494 + energyVertex/gFragmentMass)/931.494;
      betaVertex              = Sqrt(1.0 - 1.0/gammaTemp/gammaTemp);
    }
    //End energy at fragmentation point
    //---------------------------------
    //---------------------------------------------------------
    // Broadening in angles due to fragmentation in the target:
    double theta2 = -10.0;
    double phi2   = -10.0;  
    if(gTargetBroadeningOption==1)
    {
      while(theta2< 0.0 || theta2>180.0 || phi2<0.0 || phi2>360.0)
      {  
        theta2 = fRandom->Gaus(0.0,gTargetBroadeningThetaSigma);               // cout<<"ftheta2 = "<<ftheta2<<endl;
        phi2   = fRandom->Rndm() * 360.;                                       // cout<<"fphi2 = "<<fphi2<<endl;
      }

      phi2   = phi2 * PI/180.;                                                 // cout<<"fphi2 = "<<fphi2<<endl;
      theta2 = theta2 * PI/180.;                                               // cout<<"ftheta2 = "<<ftheta2<<endl;
                                                                               // vector after the fragmentation
      gFragmentVecAfterX = Cos(phi2) * Tan(theta2) + gBeamVecBeforeX;   					
      gFragmentVecAfterY = Sin(phi2) * Tan(theta2) + gBeamVecBeforeY;
      gFragmentVecAfterZ = Cos(theta2) * Tan(theta2) + gBeamVecBeforeZ; 
                                                                               // cout<<"z_pva = "<<gFragmentVecAfterZ<<endl;

      double length1 = sqrt(gFragmentVecAfterX*gFragmentVecAfterX + gFragmentVecAfterY*gFragmentVecAfterY + gFragmentVecAfterZ*gFragmentVecAfterZ); 
                                                                               // cout<<"length1 = "<<length1<<endl;

      gFragmentVecAfterX = gFragmentVecAfterX/length1;                         // cout<<"x_pva = "<<x_pva<<endl;
      gFragmentVecAfterY = gFragmentVecAfterY/length1;                         // cout<<"y_pva = "<<y_pva<<endl;
      gFragmentVecAfterZ = gFragmentVecAfterZ/length1;                         // cout<<"z_pva = "<<z_pva<<endl;
    }
    else
    {
      gFragmentVecAfterX = gBeamVecBeforeX;
      gFragmentVecAfterY = gBeamVecBeforeY;
      gFragmentVecAfterZ = gBeamVecBeforeZ;
    }
    
    // Moving the decay point already to the point of fragmentation
    // the later "real" decay point will be added by the time they need to decay times velocity and vector
    gGammaDecayPointX = gBeamInteractionPointX;
    gGammaDecayPointY = gBeamInteractionPointY;
    gGammaDecayPointZ = gBeamInteractionPointZ;

    // Selecting the populated excitation level
    float randNumber       = fRandom->Rndm() * gTotalLevelExcitationProbability;
    int populatedLevelID   = 0;
    int decayIntoLevelID   = 0;
    float excitationEnergy = 1.;
    float decayEnergy      = 0.;

    double halflife, meanlife, lambda, decayTime, decayTimeAfterInteraction;
    for(int j=0;j<gNumberOfLevels;j++)
    {
      if(randNumber>=SLevel[j].beginOfTotalExcitationProbability && randNumber<SLevel[j].endOfTotalExcitationProbability)
      populatedLevelID = j; 
    }
    // cout<<"Populated Level:"<<populatedLevelID<<endl;
    // Ok, we have the initialy excited level. Now we have to determine the decay pattern
    while(excitationEnergy != 0. && SLevel[populatedLevelID].numberOfDecayBranches >0 && SLevel[populatedLevelID].totalDecayProbability>0.)
    {
      randNumber = fRandom->Rndm() * SLevel[populatedLevelID].totalDecayProbability;
      for(int j=0;j<SLevel[populatedLevelID].numberOfDecayBranches;j++)
      {  
        if(randNumber>=SLevel[populatedLevelID].beginOfDecayProbability[j] && randNumber<SLevel[populatedLevelID].endOfDecayProbability[j])
        { 
          decayIntoLevelID = SLevel[populatedLevelID].decayIntoLevelId[j];
          decayEnergy      = SLevel[populatedLevelID].energy - SLevel[decayIntoLevelID].energy;
          if(decayEnergy<0.){cout<<"Decay energy smaller than zero. Aborting program."<<endl; abort();}
        }
      }
  
      // The beam has to be shot up to the decay point. For simplicity, I assume that gamma_temp is constant during that time,
      // which causes errors of about 1% or so, depending on the beta spread in the target.
      halflife                  = SLevel[populatedLevelID].halflife;
      meanlife                  = halflife*0.001/log(2.0)*gammaTemp;
      lambda                    = 1.0/meanlife;
      decayTime                 = fRandom->Exp(1/lambda);                      // Only the decay time for the level!!!
      decayTimeAfterInteraction = decayTimeAfterInteraction + decayTime;       // Decay time from reaction point
        
      float energyDeexcite;                                                    // cout<<"Starting Beam"<<endl;
      //_________________________________________________________________________________________
      float energyHalfStep = (gFragmentdEdX[0][0] - gFragmentdEdX[0][1])/2;    // Have to calculate the energy loss myself:
      for(int j=0;j<50;j++)
      {
        if(gFragmentdEdX[0][j]+energyHalfStep<energyVertex/gFragmentMass && j>0 && gFragmentdEdX[0][j-1]+energyHalfStep>energyVertex/gFragmentMass)
        {
          float dummyEnergy1  = energyVertex;
          float dummyEnergy2  = gFragmentdEdX[0][j] - energyHalfStep;
          float dummyTime     = 0.;
          float dummyPosition = gGammaDecayPointZ;
          int dummyCounter    = j;
          double dummyGamma   = 0.;
          double dummyBeta    = 0.;
          while(dummyTime<decayTimeAfterInteraction)
          {
            if(dummyPosition>gTargetThickness/2.0) break;
            dummyGamma        = (931.494 + dummyEnergy1/gFragmentMass)/931.494;
            dummyBeta         = Sqrt(1.0 - 1.0/dummyGamma/dummyGamma);
            dummyTime         = dummyTime + 0.1*gDefaultCutValue/(dummyBeta*29.97925);  //have to convert into cm because beamdEdX is given in MeV/cm**2
            dummyEnergy1      = dummyEnergy1 - gFragmentdEdX[1][dummyCounter]*0.1*gDefaultCutValue * gTargetDensity[gTargetKind]*1000;
            dummyPosition     = dummyPosition + 0.1*gDefaultCutValue;            
            gGammaDecayPointZ = dummyPosition;
            if(dummyEnergy1/gFragmentMass < dummyEnergy2)
            {
              dummyCounter++;
              dummyEnergy2 = gFragmentdEdX[0][dummyCounter] + energyHalfStep;
            }
          }
          if(dummyTime<decayTimeAfterInteraction)
          {
            gGammaDecayPointZ = gGammaDecayPointZ + (decayTimeAfterInteraction - dummyTime) * (dummyBeta*29.97925);
          }
          energyDeexcite = dummyEnergy1;
        }
      }
      gammaTemp = (931.494+energyDeexcite/gFragmentMass)/931.494;    // Getting beta at deexcitation time.
      gBetaReal = Sqrt(1.0 - 1.0/gammaTemp/gammaTemp);
    }
    //_________________________________________________________________________________________
    //------------------------------------------ 
    // Now comes the deexcitation part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //------------------------------------------
      
    // Determine the theta distribution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double ranHomo = fRandom->Rndm() * gSumHomo;
    for(int mmm=1;mmm<18001;mmm++)
    {
      if(ranHomo>=gSumHomoUpToNow[mmm-1] && ranHomo<=gSumHomoUpToNow[mmm])
      {
        gGammaThetaRest = mmm/100.0/180.0*PI;
      }
    }
    
    gGammaThetaLab  = ACos((Cos(gGammaThetaRest)+gBetaReal)/(1.0+gBetaReal*Cos(gGammaThetaRest)));
    TVector3 vecIn(gFragmentVecAfterX,gFragmentVecAfterY,gFragmentVecAfterZ);
    TVector3 vecOut = GetVector(vecIn, gGammaThetaLab);

    if(vecOut.X()==-999.9 && vecOut.Y()==-999.9 && vecOut.Z()==-999.9)
    {
      restartLoop = 1; 
      continue;
    }
    gGammaVecX = vecOut.X();
    gGammaVecY = vecOut.Y();
    gGammaVecZ = vecOut.Z();
		
    gFragmentEnergyAfterTarget = 0;                                           // The energy after the target
    energyHalfStep = (gFragmentdEdX[0][0] - gFragmentdEdX[0][1])/2;
    for(int j=0;j<50;j++)
    {
      if(gFragmentdEdX[0][j]+energyHalfStep<energyVertex/gFragmentMass && j>0 && gFragmentdEdX[0][j-1]+energyHalfStep>energyVertex/gFragmentMass)
      {
        // cout<<"condition fullfilled: "<<j<<endl;
        float dummyEnergy1  = energyVertex;
        float dummyEnergy2  = gFragmentdEdX[0][j] - energyHalfStep;
        float dummyPosition = gGammaDecayPointZ;
        int dummyCounter    = j;
        while(dummyPosition<gTargetThickness/2.0)
        {
          dummyEnergy1 = dummyEnergy1 - gFragmentdEdX[1][dummyCounter]*0.1*gDefaultCutValue * gTargetDensity[gTargetKind]*1000;
          dummyPosition = dummyPosition + 0.1*gDefaultCutValue;            
          if(dummyEnergy1/gFragmentMass < dummyEnergy2)
          {
            dummyCounter++;
            dummyEnergy2 = gFragmentdEdX[0][dummyCounter] + energyHalfStep;
          }
        }
        gFragmentEnergyAfterTarget = dummyEnergy1;
      }
    }
    //_________________________________________________________________________________________
    gammaTemp           = (931.494+gFragmentEnergyAfterTarget/gFragmentMass)/931.494;
    gBetaAfter          = Sqrt(1.0 - 1.0/gammaTemp/gammaTemp);	    
    gGammaEnergyRest    = decayEnergy;
    gGammaEnergyDoppler = gGammaEnergyRest *(Sqrt(1.0-gBetaReal*gBetaReal)/(1.0-gBetaReal*Cos(gGammaThetaLab)));
	        
    gEvNum = i;
	
    t->Fill();                                                                 // cout<<"Filling root tree"<<endl;

    //Setting the new excitation energy and populated level:
    excitationEnergy = SLevel[decayIntoLevelID].energy;
    populatedLevelID = SLevel[decayIntoLevelID].ID;                   
    gDecayIDEvNum++;   
    if(restartLoop==1) continue; 
    i = i + 1;
  } // end while loop through events  
}

TVector3 GetVector(TVector3 vectorIn, float theta)
{
  double x1     = vectorIn.X();
  double y1     = vectorIn.Y();
  double z1     = vectorIn.Z();
  double l1     = sqrt(x1*x1 + y1*y1 + z1*z1);
  double theta1 = ACos(z1/l1);
  double phi1   = 0.0;

  if(Sqrt((x1/l1/Sin(theta1))*(x1/l1/Sin(theta1)))>1.0) 
  {
    x1     = 0.0;
    y1     = 0.0;
    z1     = 1.0;
    theta1 = 0.0;
    phi1   = 0.0;
  }
  if(Sqrt((x1/l1/Sin(theta1))*(x1/l1/Sin(theta1)))<=1.0) 
  {
    phi1 = ACos(x1/l1/Sin(theta1));
  }

  TVector3 rotatedPos;

  double phi    = fRandom->Rndm()*360.0 * PI/180.;
  rotatedPos.RotateY(theta);  
  rotatedPos.RotateZ(phi);
  rotatedPos.RotateY(theta1);  
  rotatedPos.RotateZ(phi1);
	      
  return rotatedPos;
}

double EnergyDetermination()
{
return 0.;
}


