#ifndef MYNAMESPACE_H
#define MYNAMESPACE_H
//#include "G4RunManager.hh"
//#include "RunAction.hh"
//#include "DetectorConstruction.hh"
//#include "PrimaryGeneratorAction.hh"
//#include "SteppingAction.hh"
//#include "G4UImanager.hh"
//#include "PhysicsList.hh"
//
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"

// Defining namespace for all the global variables
namespace MyNamespace  {

  //Defining all the variables
  //Initialzing the variables that need an initial value
  int mass                = 1;
  int z                   = 1;
  int charge              = 1;
  int numevent_total      = 10000;
  int mass_f              = 1;
  int z_f                 = 1;
  int charge_f            = 1; 
  int kind_target         = 2;            // kind_target 1 Au 2 Be.... 
  int deltamass           = 0;
  int deltazet            = 0;            

  float energy_beam_mean  = 0.;
  float energy_beam_sigma = 0.;           // Energy_of_beam A MeV
  float target_X          = 3.8;           // dimension of the target in cm for X,Y and mg/cm^2 for Z
  float target_Y          = 3.8;
  float thickness         = 0.001;   

  char temp[300];
  char gamma_in[300]      = "dummy.in";
  char root_in[300]       = "dummyIn.root";
  char root_out[300]      = "dummy.root";
  char tempparticle[300];

  float theta_low         = 0.;
  float theta_high        = 180.;           // Angle covered
  float x_beam            = 0.; 
  float y_beam            = 0.; 
  float x_beam_sigma      = 0.0;
  float y_beam_sigma      = 0.0;
  float theta_beam        = 0.;
  float theta_beam_sigma  = 0.0;
  float phi_beam_min      = 0.;
  float phi_beam_max      = 360.;
  int target_broadening_option = 0;        // Option to take also the angular broadening of the beam due to reaction
  float theta_target_broadening_sigma = 0.0;

  int borrel_option       = 0;
  int goldhaber_option    = 0;
  float sigma0_goldhaber  = 90.;
  float defaultCutValue   = 0.001;
  int degraderOption      = 0;             // Option to use a degrader after the taget for RDDS measurements
  int atomicBGOption      = 0;             // Option to insert atomic bg
  //_______________________________________________
  double sum_dist;
  double sum_dist_uptonow[18001];
  //_______________________________________________
  float evnum;
  float decayIDevnum;
  float p0[3],pvb[3],pva[3];
  float energy_p;
  float beta_b,beta_r,beta_a,halflife,decay_time_after_interaction;
  float g0[3],gv[3];
  float e_rest,e_doppler;
  float theta_gamma_rest, theta_gamma_lab;
  float energy_vertex_stored;
  float binding_energy_nucleon;
  double density[9];
  //_______________________________________________
  int dEdXTableOption          = 0;
  char dEdXTableInputBeam[300] = "./dEdXTables/dummy.in";
  char dEdXTableInputFragment[300] = "./dEdXTables/dummy.in";
  float beamdEdX[2][50]        = {{0.0}}; // The first value gives the projectile energy in A MeV, the second dEdX in MeV/mgcm^2;
  float fragmentdEdX[2][50]    = {{0.0}};
  //_______________________________________________
  struct Level  {
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
  float gTotalLevelExcitationProbability = 0.;
  int gNumberOfLevels                    = 0; 
  int gNumberOfExcitedLevels             = 0; 
  //_______________________________________________
  int fDistributionType = 0;
  //_______________________________________________
  // The ROOT-Tree
  TFile *infile;                                                                  // The ROOT-Tree Reading
  TTree *t;
  TTree *tHeader;  // The header for information that doesn't change.

  TFile *outfile;                                                                 // The ROOT-Tree Writing
  TTree *t2;
  TTree *t2Header;
  //_______________________________________________
  float atomicCrossSection       = 0.0;               // The atomic cross section given in mbarn
  bool  bgGamma                  = false;             // To know if the simulated gamma stems from the background or not
  float totalCrossSection        = 0.0;               // The total cross-section of atomic bg gamma or from the excitation in mbarn
  char  bgFileName[300]          = "bgFileName.root"; // The location of the bg file name, which must be a root file!!!!
  char  bgSpectrumName[300]      = "bgSpectrumName";  // The name of the spectrum in the bg root file
  TFile *bgFile;                                      // The root file of the bg
  TH2F  *h_atomicBackground;                          // The histogram of the bg
  float gAverageNumberOfGammas   = 1;                 // Average number of gammas per incident particle. Can be changed in the input file
 
  //_______________________________________________
  // Variables needed only from builder on:


  float total_event_num;  // Put in the header of the root tree to know, how many incident events were simulated
  int gamma_det_type = 0; // RIKEN:1 for dali2, 2 for grape, 3 for sgt, 4 for Shogun   
  float pos_det[3]={0.,0.,200.};
  float ver_rec[3];
  float beta_rec;

  int stq_opt     = 0;
  int target_holder_opt = 0;
  int beam_pipe_opt     = 0;

  int nentries;                                                                   // Number of events in the input file
  float dummy=0.0;                                                                // Used very often for quick calculations;

  // RIKEN:
  //-----------------------------------------------
  //For the Dali2 array
  //const int numberOfDali2Crystals = NUMBEROFDALI2CRYSTALS; 

  int dali2_flag[NUMBEROFDALI2CRYSTALS]= {0};
  float dali2_pos[3][NUMBEROFDALI2CRYSTALS];
  float dali2_fi[3][NUMBEROFDALI2CRYSTALS];
  float dali2_energy_not_cor[NUMBEROFDALI2CRYSTALS] = {0.};
  float dali2_time[NUMBEROFDALI2CRYSTALS] = {0.};
  //-----------------------------------------------
  //For the Shogun array using LaBr3
  int Shogun_flag[NUMBEROFSHOGUNCRYSTALS] = {-1};
  float Shogun_pos[3][NUMBEROFSHOGUNCRYSTALS] = {{0.}};
  float Shogun_fi[3][NUMBEROFSHOGUNCRYSTALS];
  float Shogun_energy_not_cor[NUMBEROFSHOGUNCRYSTALS] = {0.};
  float Shogun_time[NUMBEROFSHOGUNCRYSTALS] = {0.};
  int   Shogun_number_of_crystals = 0;
  float ShogunHousingThicknessX=0.,ShogunHousingThicknessY=0.,ShogunHousingThicknessZ=0.;
  float ShogunMgOThicknessX=0.,ShogunMgOThicknessY=0.,ShogunMgOThicknessZ=0.;
  //-----------------------------------------------
  //For the Grape array
  //pos,energy [xx][xx][0] is for the core
  int grape_flag[NUMBEROFGRAPEDETECTORS], grape_crystal_flag[NUMBEROFGRAPEDETECTORS][2][10];
  float grape_pos[3][NUMBEROFGRAPEDETECTORS][2][10];
  float grape_energy_not_cor[NUMBEROFGRAPEDETECTORS][2][10] = {{{0.}}};
  float grape_time[NUMBEROFGRAPEDETECTORS][2][10] = {{{0.}}};
  //----------------------------------------------
  //For the strip ge telescope
  //----------------------------------------------
  int sgt_flag[10], sgt_coax_crystal_flag[10], sgt_planar_strip_flag[10][26];     // Fictive central contact plus the 25 strips
  float sgt_pos[3][10][27];                                                       // First the coax, then the fictive central planar contact, then the strips
  float sgt_energy_not_cor[10][27] = {{0.}}; 
  float sgt_time[10][27] = {{0.}};
  //----------------------------------------------
  //For the sphere:
  float sphere_inner_radius = 30.;
  float sphere_outer_radius = 36.;
  float sphere_fi[3]={0.};                                                        // fi = first interaction
  float sphere_energy_not_cor = 0.;

  //----------------------------------------------
  // For the LaBr3 array of huge detector size:
  int LaBr3Flag[NUMBEROFLABR3ARRAYCRYSTALS] = {-1};
  float LaBr3Pos[3][NUMBEROFLABR3ARRAYCRYSTALS] = {{0.}};
  float LaBr3EnergyNotCor[NUMBEROFLABR3ARRAYCRYSTALS] = {0.};
  float LaBr3Time[NUMBEROFLABR3ARRAYCRYSTALS] = {0.};
  int   LaBr3NumberOfCrystals = 0;
  float LaBr3HousingThicknessFront=0.,LaBr3HousingThicknessSide=0;
  float LaBr3InsulationThicknessFront=0.,LaBr3InsulationThicknessSide=0;
  //________________________________________________
  int dali2_opt          = 0;                                                     // Option, which detector system is simulated
  float z_pos_shift= 0;                                                           // Shift everything along the Z axis in order to account for the target thickness.  
  int dali2_fi_opt       = 0; 
  int Shogun_opt         = 0;
  int Shogun_fi_opt      = 0;
  int Shogun_housing_opt = 0;
  int Shogun_MgO_opt     = 0;
  int grape_opt          = 0;
  int sgt_opt            = 0;
  int sphere_opt         = 0;
  int LaBr3Opt           = 0;
  int LaBr3HousingOpt    = 0;
  int LaBr3InsulationOpt = 0;
  int collimator_opt     = 0;                                                     // Option to insert a collimator
  int en_det_after_target_opt = 0;                                                // Option for energy detector after the target
  // Input of the resolutions:
  int dali2_en_res_opt           = 2;                                             // Option for the type of resolution 
  int Shogun_en_res_opt          = 2;                                             // In principle this parameter is not necessary anymore but still present 
  int grape_en_res_opt           = 2;                                             // for compability reasons... :(
  int sgt_en_res_opt             = 2;
  int LaBr3EnResOpt              = 2;
  int dali2_en_res_ind           = 0;                                             // Option to give indiviual resolution values
  
  float dali2_en_res[2]          = {0.0,0.5};
  float Shogun_en_res[2]         = {0.};
  float grape_en_res[2]          = {0.0};
  float sgt_en_res[2]            = {0.0};
  float LaBr3EnRes[2]            = {0.0};
  float dali2_time_res[2]        = {0.0};
  float Shogun_time_res[2]       = {0.0};
  float grape_time_res[2]        = {0.0};
  float sgt_time_res[2]          = {0.0};
  float LaBr3TimeRes[2]          = {0.0};
  float pos_det_at_target_res    = 0.0;                                           // Pos Resolution at the target position in mm FWHM!
  float pos_det_after_target_res = 0.0;                                           // Pos Resolution  after the target in mm FWHM
  float beta_res                 = 0.00;                                         // Resolution of beta in FWHM!
  
  //For the Absorber:
  float ShieldRadius      = 0.0;
  float PbShieldThickness = 0.0;
  float SnShieldThickness = 0.0;
 
  //For the addback
  int shogun_addback_opt          = 0;
  int dali2_addback_opt          = 0;                                             // Use Cluster addback in Reconstructor
  float maxAddbackDistance       = 10;
  int addbackTable[NUMBEROFDALI2CRYSTALS][NUMBEROFDALI2ADDBACKCRYSTALS] = {{0}};
  int shogunAddbackTable[NUMBEROFSHOGUNCRYSTALS][NUMBEROFSHOGUNADDBACKCRYSTALS] = {{0}};

  // RISING:
  //------------------------------------

  //For the Cluster array
  int cluster_flag[NUMBEROFCLUSTERDETECTORS][7]             = {{0}};
  float cluster_pos[3][NUMBEROFCLUSTERDETECTORS][7]         = {{{0.}}};
  float cluster_fi[3][NUMBEROFCLUSTERDETECTORS][7]         = {{{0.}}};
  float cluster_energy_not_cor[NUMBEROFCLUSTERDETECTORS][7] = {{0.}};
  float cluster_time[NUMBEROFCLUSTERDETECTORS][7]           = {{0.}};
  //-----------------------------------------------
  //For the Miniball array
  int miniball_flag[NUMBEROFMINIBALLDETECTORS][3][7]             = {{{0}}};
  float miniball_pos[3][NUMBEROFMINIBALLDETECTORS][3][7]         = {{{{0.}}}};
  float miniball_fi[3][NUMBEROFMINIBALLDETECTORS][3][7]         = {{{{0.}}}};
  float miniball_energy_not_cor[NUMBEROFMINIBALLDETECTORS][3][7] = {{{0.}}};
  float miniball_time[NUMBEROFMINIBALLDETECTORS][3][7]           = {{{0.}}};
  //-----------------------------------------------
  //For the Hector array
  int hector_flag[NUMBEROFHECTORDETECTORS]             = {0};
  float hector_pos[3][NUMBEROFHECTORDETECTORS]         = {{0.}};
  float hector_energy_not_cor[NUMBEROFHECTORDETECTORS] = {0.};
  float hector_time[NUMBEROFHECTORDETECTORS]           = {0.};
  //-----------------------------------------------
  //Option, which detector system is simulated
  int cluster_opt             = 0;
  int miniball_opt            = 0;
  int hector_opt              = 0;

  // Option to determine the first interaction points.
  int cluster_fi_opt          = 0;
  int miniball_fi_opt         = 0;
  //Input of the resolutions:
  //Option for the type of resolution 1==linear increase, 2== more realistic increase with sqrt function
  int cluster_en_res_opt  = 0;
  int miniball_en_res_opt = 0;
  int hector_en_res_opt   = 0;
  float cluster_en_res[2];
  float miniball_en_res[2];
  float hector_en_res[2];
  float cluster_time_res[2];
  float miniball_time_res[2];
  float hector_time_res[2];
 
  float cluster_thickness_Pb,cluster_thickness_Sn,cluster_thickness_Al;           // Absorber thicknesses
  float miniball_thickness_Pb,miniball_thickness_Sn,miniball_thickness_Al;
  float hector_thickness_Pb,hector_thickness_Sn,hector_thickness_Al;

  double bp_inside  = 7.0;
  double bp_outside = 7.5;
}
#endif
