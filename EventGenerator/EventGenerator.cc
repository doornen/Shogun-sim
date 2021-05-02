// Some of the variables may not be necessary to be written into the tree.
// These are:
// beta_b or energy_p: one value can be calculated from the other
// theta_gamma_rest: just a check of the gamma-ray angular distribution
// theta_gamma_lab: just a check of the gamma-ray angular distribution
// energy_vertex_stored: not a measurable quantity and no influence on the Doppler correction

#include "Globals.hh"
#include "MyNamespace.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
  #include "G4VisExecutive.hh"
#endif

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include <iostream>
#include "G4ios.hh"

using namespace std;
using namespace MyNamespace;

//_______________________________________________
G4RunManager           *runManager;
RunAction              *run;
DetectorConstruction   *det;
PrimaryGeneratorAction *kin;
SteppingAction         *steppingaction;
G4UImanager            *UI;
PhysicsList            *phys;

//_______________________________________________
void ReadInputFile();
void ReadGammaFile();
void RunSimulation();
void ReaddEdXTable();
void ReadAtomicBG();
G4ThreeVector GetNextVector(G4ThreeVector vector_in, float theta);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc,char** argv) {

  // Reading the input file 
  ReadInputFile();

  // Read the atomic bg, if the option is turned on
  if(atomicBGOption==1) ReadAtomicBG();

  // Reading the gamma-ray decay file
  ReadGammaFile();

  // Calculating the fragment mass and charge:
  z_f      = z      - deltazet;
  mass_f   = mass   - deltamass;
  charge_f = charge - deltazet;

  // Some density of some material:
  density[0] = 0.0;     // Vacuum
  density[1] = 19.3;    // Au
  density[2] = 1.848;   // Be
  density[3] = 1.82;/// 2009.3.11: density changed 2.2;     // C
  density[4] = 7.874;   // Fe
  density[5] = 11.35;   // Pb
  density[6] = 0.07099; // LH2
  density[7] = 6.52;    // Zr
  density[8] = 0.95;    // CH2
  
  // Calculating the target thickness in cm!
  thickness = thickness/density[kind_target]/1000.0;

  // Making a theta distribution
  ///* 
  sum_dist = 0.0;
  //distribution 0: uniform; 2200: J,M 2,2-> 0,0; 2100: J,M 2,1->0,0 2000: J,M 2,0 -> 0,0
  // The zet axis is the beam axis
  for(int m2=(int)(theta_low*100);m2<(int)(theta_high*100+1);m2++){					
    G4double angle  = m2/100.0*degree;
    if(fDistributionType == 2000)  
      sum_dist = sum_dist + 1.0*sin(angle/rad)
        *( sin(angle/rad)*sin(angle/rad)*cos(angle/rad)*cos(angle/rad));
    else if(fDistributionType == 2100)
      sum_dist = sum_dist + 1.0*sin(angle/rad)
        *(1 - 3*cos(angle/rad)*cos(angle/rad) + 4*cos(angle/rad)*cos(angle/rad)*cos(angle/rad)*cos(angle/rad));
    else if (fDistributionType == 2200)
      sum_dist = sum_dist + 1.0*sin(angle/rad)
        *(1 - cos(angle/rad)*cos(angle/rad)*cos(angle/rad)*cos(angle/rad));
    // Making it uniform:
    else sum_dist = sum_dist + 1.0*sin(angle/rad);
    sum_dist_uptonow[m2] = sum_dist;
  }
  
  // Choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed;
  seed = time(0);  CLHEP::HepRandom::setTheSeed(seed);

  // My Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    
  // Construct the default run manager
  runManager = new G4RunManager;

#ifdef G4VIS_USE
   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();
#endif   

  runManager->SetUserInitialization(det  = new DetectorConstruction);
  det       ->UpdateGeometry();
  runManager->SetUserInitialization(phys = new PhysicsList(defaultCutValue));
  runManager->SetUserAction(kin = new PrimaryGeneratorAction(det));
  
  // Set user action classes
  runManager->SetUserAction(run = new RunAction(det,phys,kin)); 
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new TrackingAction(run));
  runManager->SetUserAction(steppingaction = new SteppingAction(det,run));

  // Get the pointer to the User Interface manager 
  // G4UImanager *UI = G4UImanager::GetUIpointer();
  UI = G4UImanager::GetUIpointer();
 
  // Setting target material and thickness in cm
  det->SetTargetSize(target_X,target_Y,thickness);
  if(kind_target==1) det->SetTargetMaterial("Au");
  if(kind_target==2) det->SetTargetMaterial("Be");
  if(kind_target==3) det->SetTargetMaterial("C");
  if(kind_target==4) det->SetTargetMaterial("Fe");
  if(kind_target==5) det->SetTargetMaterial("Pb");
  if(kind_target==6) det->SetTargetMaterial("LH2");
  if(kind_target==7) det->SetTargetMaterial("Zr");
  if(kind_target==8) det->SetTargetMaterial("CH2");  

  det->UpdateGeometry();

  if(argc==1){
  // Define (G)UI terminal for interactive mode  
    cout<<"Interactive mode"<<endl;
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = 0; 
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    

    UI->ApplyCommand("/control/execute vis.mac");    
    session->SessionStart();
    delete session;
  }
  else  {
  // Batch mode
    cout<<"Batch mode"<<endl;
    G4String command  = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  UI->ApplyCommand("/vis/scene/notifyHandlers");
  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 0");

  TFile *rootfile = new TFile(root_out,"RECREATE");
  rootfile->cd();

  t = new TTree("Events","Events");
  t->Branch("EventNumber",&evnum,"evnum/F");
  t->Branch("DecayIDOfEventNumber",&decayIDevnum,"decayIDevnum/F");
  t->Branch("ProjectileVertex",p0,"p0[3]/F");                                   // Position at the reaction point
  t->Branch("VProjectileBeforeTarget",pvb,"pvb[3]/F");                          // Normalized Vector of beam before the target
  t->Branch("VProjectileAfterTarget",pva,"pva[3]/F");                           // Normalized Vector of beam after the target
  t->Branch("EnergyProjectile",&energy_p,"energy_p/F");                         // Energy of beam before the target in MeV/u
  t->Branch("BetaBeforeTarget",&beta_b,"beta_b/F");                             // Beta before the target
  t->Branch("BetaReal",&beta_r,"beta_r/F");                                     // Beta at deexcitation	
  t->Branch("BetaAfterTarget",&beta_a,"beta_a/F");                              // Beta After Target
  t->Branch("Halflife",&halflife,"halflife/F");                                 // Halflife
  t->Branch("DecayTimeAfterInteraction",&decay_time_after_interaction,"decay_time_after_interaction/F");
  t->Branch("VertexGamma",g0,"g0[3]/F");                                        // Position at the gamma emmittance point
  t->Branch("VGamma",gv,"gv[3]/F");		                                // Gamma vector
  t->Branch("EGammaRest",&e_rest,"e_rest/F");		                        // Energy at rest
  t->Branch("EGammaDoppler",&e_doppler,"e_doppler/F");                          // Theta of doppler boosted gamma
  t->Branch("ThetaGammaRest",&theta_gamma_rest,"theta_gamma_rest/F");
  t->Branch("ThetaGammaLab",&theta_gamma_lab,"theta_gamma_lab/F");
  t->Branch("EnergyVertex",&energy_vertex_stored,"energy_vertex_stored/F");     // Energy of projectile at reaction
  t->Branch("BGGamma",&bgGamma,"bgGamma/O");                                    // To know if the gamma is from bg or a good event
  // The header tree gets the information that is not changed.
  tHeader  = new TTree("Header","Header");
  tHeader->Branch("Mass",&mass,"mass/I");                                       // Beam mass
  tHeader->Branch("Z",&z,"z/I");                                                // Beam z	
  tHeader->Branch("Charge",&charge,"charge/I");                                 // Beam charge	
  tHeader->Branch("MassFragment",&mass_f,"mass_f/I");                           // Fragment mass
  tHeader->Branch("ZFragment",&z_f,"z_f/I");	                                // Fragment z
  tHeader->Branch("ChargeFragment",&charge_f,"charge_f/I");                     // Fragment charge
  tHeader->Branch("TargetKind",&kind_target,"kind_target/I");
  tHeader->Branch("TargetThicknessCM",&thickness,"thickness/F");
  tHeader->Branch("ThetaLow",&theta_low,"theta_low/F");                         // Theta angle covered by simulation
  tHeader->Branch("ThetaHigh",&theta_high,"theta_high/F");
  tHeader->Fill();

  UI->ApplyCommand("/testem/phys/addPhysics standard");
  UI->ApplyCommand("/run/initialize");

  sprintf(tempparticle,"/gun/ion %i %i %i",z,mass,charge);
  UI->ApplyCommand("/gun/particle ion");
  UI->ApplyCommand(tempparticle);
 
  // Starting the simulation:
  RunSimulation();
  // Writing stuff to file and spectra:
  tHeader->Write();
  t->Write();

  UI->ApplyCommand("/vis/viewer/update");
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager; 
  
  rootfile->Close();

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void ReadInputFile()  {
  // The Input files are of free type.
  FILE *fin = fopen("./input/EventGenerator.in","r");
  while(!feof(fin))  {
    fscanf(fin,"%s ",temp); 
    if(strcmp(temp,"BEAMISOTOPE")==0)  {
      fscanf(fin,"%i %i %i ",&mass,&z,&charge); 
      printf("%s %i %i %i\n",temp,mass,z,charge);
    }
    else if(strcmp(temp,"BEAMENERGY")==0)  {
      fscanf(fin,"%f %f",&energy_beam_mean,&energy_beam_sigma); 
      energy_beam_sigma = energy_beam_sigma/2.35;
      printf("%s %f %f \n",temp,energy_beam_mean,energy_beam_sigma);
    }
    else if(strcmp(temp,"BEAMPOSITION")==0)  {
      fscanf(fin,"%f %f %f %f",&x_beam,&x_beam_sigma,&y_beam,&y_beam_sigma);
      x_beam_sigma = x_beam_sigma/2.35;y_beam_sigma = y_beam_sigma/2.35;
      printf("%s %f %f %f %f \n",temp,x_beam,x_beam_sigma,y_beam,y_beam_sigma);
    }
    else if(strcmp(temp,"BEAMANGLE")==0)  {
      fscanf(fin,"%f %f %f %f",&theta_beam,&theta_beam_sigma,&phi_beam_min,&phi_beam_max);
      theta_beam_sigma=theta_beam_sigma/2.35;
      printf("%s %f %f %f %f \n",temp,theta_beam,theta_beam_sigma,phi_beam_min,phi_beam_max);
    }
    else if(strcmp(temp,"TARGET")==0)  {
      fscanf(fin,"%i %f %f %f",&kind_target,&target_X,&target_Y,&thickness);
      if(kind_target<0 || kind_target>8)  {
        cout<<"Your Target doesn't exist. Aborting program."<<endl; 
        abort();
      }
      printf("%s %i %f %f %f \n",temp,kind_target,target_X,target_Y,thickness);
    }   
    else if(strcmp(temp,"TARGETANGULARBROADENING")==0)  {
      fscanf(fin,"%i %f",&target_broadening_option,&theta_target_broadening_sigma);
      theta_target_broadening_sigma = theta_target_broadening_sigma/2.35;
      printf("%s %i %f \n",temp,target_broadening_option,theta_target_broadening_sigma);
    }
    else if(strcmp(temp,"MASSCHANGE")==0)  {
      fscanf(fin,"%i %i",&deltamass,&deltazet);
      printf("%s %i %i \n",temp,deltamass,deltazet);
    }
    else if(strcmp(temp,"BORREL")==0)  {
      fscanf(fin,"%i %f",&borrel_option,&binding_energy_nucleon);
      printf("%s %i %f \n",temp,borrel_option,binding_energy_nucleon);
    }
    else if (strcmp(temp,"GOLDHABER")==0)  {
      fscanf(fin,"%i %f",&goldhaber_option,&sigma0_goldhaber);
      printf("%s %i %f \n",temp,goldhaber_option,sigma0_goldhaber);
    }
    else if(strcmp(temp,"GAMMAINPUT")==0)  {
      fscanf(fin,"%s",&gamma_in);
      printf("%s %s \n",temp,gamma_in);
    }
    else if(strcmp(temp,"THETARANGE")==0)  {
      fscanf(fin,"%f %f",&theta_low,&theta_high);
      printf("%s %f %f \n",temp,theta_low,theta_high);
    }
    else if(strcmp(temp,"NUMBEROFEVENTS")==0)  {
      fscanf(fin,"%i ",&numevent_total);
      printf("%s %i\n",temp,numevent_total);
    }
    else if(strcmp(temp,"DEFAULTCUTVALUE")==0)  {
      fscanf(fin,"%f ",&defaultCutValue);
      printf("%s %f\n",temp,defaultCutValue);
    }
    else if(strcmp(temp,"OUTPUTFILE")==0)  {
      fscanf(fin,"%s ",&root_out); 
      printf("%s %s \n",temp,root_out);
    }
    else if(strcmp(temp,"DEDXTABLE")==0)  {
      fscanf(fin,"%i %s %s",&dEdXTableOption,&dEdXTableInputBeam,&dEdXTableInputFragment); 
      printf("%s %i %s %s\n",temp,dEdXTableOption,dEdXTableInputBeam,dEdXTableInputFragment);
      if(dEdXTableOption==1) ReaddEdXTable();
    }
    else if(strcmp(temp,"DISTRIBUTIONTYPE")==0)  {
      fscanf(fin,"%i ",&fDistributionType); 
      printf("%s %i \n",temp,fDistributionType);
    }
    else if(strcmp(temp,"DEGRADER")==0)  {//Option to use a degrader after the target for RDDS measurements
      fscanf(fin,"%i ",&degraderOption);  //Not yet implemented
      printf("%s %i \n",temp,degraderOption);
    }
    else if(strcmp(temp,"ATOMICBG")==0)  {// Option to simulate the atomic background from a file
      fscanf(fin,"%i %s %s %f %f",&atomicBGOption,&bgFileName,&bgSpectrumName,&atomicCrossSection,&gAverageNumberOfGammas); 
      printf("%s %i %s %s %f %f\n",temp,atomicBGOption,bgFileName,bgSpectrumName,atomicCrossSection,gAverageNumberOfGammas);
      if(atomicBGOption==1)
         totalCrossSection = atomicCrossSection; // The gamma-ray cross-section will be added later.
    }
    else if(strcmp(temp,"END")==0) break;
    else {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      abort();
    }
  }
  fclose(fin);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ReadGammaFile()  {
  // At first, initializing the struct of SLevel[100]
  cout<<"Reading Gamma input file"<<endl;
  for(int i=0;i<100;i++)  {
    SLevel[i].ID = -1;
    SLevel[i].excitationProbability             = 0.;
    SLevel[i].beginOfTotalExcitationProbability = 0.;
    SLevel[i].endOfTotalExcitationProbability   = 0.;
    SLevel[i].energy                            = 0.;
    SLevel[i].halflife                          = 0.;
    SLevel[i].numberOfDecayBranches             = 0;
    SLevel[i].totalDecayProbability             = 0.;
    for(int j=0;j<5;j++)  {
      SLevel[i].decayIntoLevelId[j]             = 0;
      SLevel[i].decayIntoLevelProbability[j]    = 0.;
      SLevel[i].beginOfDecayProbability[j]      = 0.;
      SLevel[i].endOfDecayProbability[j]        = 0.;
    }
  }
  int levelID,decayLevelID; 
  float levelExcitationProbability, levelEnergy, levelHalflife;
  float branchRatio;
  
  FILE *gammaIn = fopen(gamma_in,"r");
  while(!feof(gammaIn))  {
    fscanf(gammaIn,"%s ",&temp);
    if(strcmp(temp,"LEVEL")==0)  {
      fscanf(gammaIn,"%i %f %f %f",&levelID,&levelExcitationProbability,&levelEnergy,&levelHalflife);
      // Check if the input values are greater than zero
      if(levelID < 0 || levelExcitationProbability < 0. || levelEnergy < 0.){
        cout<<"At least one of your LEVEL input values is smaller than zero. Aborting program."<<endl; 
        abort();
      }
      // Check if the level has been assigned already
      if(SLevel[levelID].ID != -1)  {
        cout<<"This LEVEL has been assigned already. Aborting program."<<endl; 
        abort();
      }
      SLevel[levelID].ID=levelID;
      SLevel[levelID].excitationProbability=levelExcitationProbability;
      totalCrossSection = totalCrossSection + levelExcitationProbability;  // adding the excitation of the nucleus and the atomic bg.
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
    else if(strcmp(temp,"DECAY")==0)  {
      fscanf(gammaIn,"%i %i %f",&levelID,&decayLevelID,&branchRatio);
      // Setting the maximum of decay branches to five:
      int branchID = SLevel[levelID].numberOfDecayBranches;
      if(branchID>4)  {
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
    else if(strcmp(temp,"END")==0) break;
    else  {
      cout<<"Could not read your input keyword of the gamma-ray decay file. Aborting program."<<endl; 
      abort();
    }
  }
  fclose(gammaIn);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ReaddEdXTable()  {
  G4cout<<"Building the energy loss table..."<<G4endl;
  float dummy1,dummy2;
  int i=0,j=0;
  FILE *tableInBeam = fopen(dEdXTableInputBeam,"r");
  while(!feof(tableInBeam)&&i<50)  {
    fscanf(tableInBeam,"%s ",temp); 
    if(strcmp(temp,"dEdX")==0)  {
      fscanf(tableInBeam,"%f %f",&dummy1,&dummy2);
      beamdEdX[0][i] = dummy1;
      beamdEdX[1][i] = dummy2;
      cout<< dummy1<<" "<<dummy2<<endl;
      i++;
    }
    if(strcmp(temp,"END")==0)  break;
  }
  fclose(tableInBeam);
  
  FILE *tableInFragment = fopen(dEdXTableInputFragment,"r");
  while(!feof(tableInFragment)&&j<50)  {
    fscanf(tableInFragment,"%s ",temp); 
    if(strcmp(temp,"dEdX")==0){
      fscanf(tableInFragment,"%f %f",&dummy1,&dummy2);
      fragmentdEdX[0][j] = dummy1;
      fragmentdEdX[1][j] = dummy2;
      cout<< dummy1<<" "<<dummy2<<endl;
      j++;
    }
    if(strcmp(temp,"END")==0)  break;  
  }
  fclose(tableInFragment);
}
         
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ReadAtomicBG()  {
  G4cout<<"Reading the atomic background"<<G4endl;

  // Get the background spectrum:
  bgFile = new TFile(bgFileName,"READ");
  bgFile->GetObject(bgSpectrumName,h_atomicBackground);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunSimulation()  {
  G4cout<<"Running the simulation"<<G4endl;

  double gamma_temp;
  int i=0;
  TRandom *randomGammas =  new TRandom;

  // Starting to loop through the incident particle events:
  while(i<numevent_total)  {
    //if(i%1000==0) cout << i << "/" << numevent_total << "   DONE!" << endl;
    if(i%1000==0){
      std::cout << "Event: " << i <<", "<< (100.*i/numevent_total) <<"% of events done\r" <<std::flush;
    }

    int restart_loop = 0;
    
    // Reseting the decay ID and the time of the event to zero;
    decay_time_after_interaction = 0.;
    decayIDevnum                 = 0.;
    G4double beam_origin_x       = 999.9;
    G4double beam_origin_y       = 999.9;
    G4double beam_origin_z       = 999.9;

    // Defining the beam's starting position
    while(beam_origin_x>(0.5*target_X) || beam_origin_x<(-0.5*target_X) || beam_origin_y>(0.5*target_Y) || beam_origin_y<(-0.5*target_Y))  {
      beam_origin_x = G4RandGauss::shoot(x_beam,x_beam_sigma);
      beam_origin_y = G4RandGauss::shoot(y_beam,y_beam_sigma);
      beam_origin_z = -1.0*thickness/2.0;                      // cout<<"Beam origin Z: "<<beam_origin_z<<endl;
    }
    // Defining the beam's starting vector.
    G4double ftheta = -10.0;
    G4double fphi   = -10.0;	
    while(ftheta< 0.0 || ftheta > 180.0 || fphi < 0.0 || fphi>360.0)  {  
      // Gaussian distribution from input values
      ftheta = G4RandGauss::shoot(theta_beam, theta_beam_sigma);
      fphi   = CLHEP::RandFlat::shoot(phi_beam_min, phi_beam_max);  
    }
    fphi   = fphi*degree;           // cout<<"fphi = "<<fphi<<endl;
    ftheta = ftheta*degree;         // cout<<"ftheta = "<<ftheta<<endl;
    
    // Normalized vector before the target					
    pvb[0] = sin(ftheta)*cos(fphi);  // cout<<"x_pvb: "<<x_pvb<<endl;
    pvb[1] = sin(ftheta)*sin(fphi);  // cout<<"y_pvb: "<<y_pvb<<endl;
    pvb[2] = cos(ftheta);            // cout<<"z_pvb: "<<z_pvb<<endl;
    
    // Energy of the projectile
    energy_p   = G4RandGauss::shoot(energy_beam_mean,energy_beam_sigma); 
    gamma_temp = (931.494+energy_p)/931.494;            // Gamma before the target
    beta_b     = sqrt(1.0 - 1.0/gamma_temp/gamma_temp); // Beta before the target
    
    // Total beam energy before the target
    double energy_total = energy_p * (double)mass;        

    // Assigning the initial beam values
    kin->GetParticleGun()->SetParticleEnergy(energy_total*MeV);
    kin->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(pvb[0]*m,pvb[1]*m,pvb[2]*m));
    kin->GetParticleGun()->SetParticlePosition(G4ThreeVector(beam_origin_x*cm,beam_origin_y*cm,beam_origin_z*cm));
   
    //------------------------------
    //****Starting Beam*************
    //------------------------------ 
    // First, getting the fragmentation point to know how far the energy loss has
    // to be calculated for the secondary beam
    p0[2] = CLHEP::RandFlat::shoot(-1.0*thickness/2.0,1.0*thickness/2.0);  // G4cout<<"Z_Point of fragmentation: "<<z_p0<<endl;
    p0[0] = beam_origin_x;	
    p0[1] = beam_origin_y;
    steppingaction->SetInteractionPoint(p0[0]*cm,p0[1]*cm,p0[2]*cm);
   
    // Reseting the target and gamma excitation flag 
    steppingaction->StartBeam();
    // Shooting the beam;
    if(dEdXTableOption==0 && energy_total>0.) runManager->BeamOn(1);  // Let GEANT4 only simulate if no dEdX tables are provided

    if(steppingaction->GetFlagTarget()==1 || dEdXTableOption==1 || energy_beam_mean==0)  {
      double energy_vertex=0.;
      if(dEdXTableOption==0 && energy_total>0.) energy_vertex = steppingaction->GetTotalEnergy(); 
      //___________________________________________________________________________________________
      else if(dEdXTableOption==1)  {
        float energyHalfStep = (beamdEdX[0][0] - beamdEdX[0][1])/2;
        // Have to calculate the energy loss myself:
        for(int j=0;j<50;j++){
          if(beamdEdX[0][j]+energyHalfStep>energy_p && 
             beamdEdX[0][j]-energyHalfStep<energy_p)  {
            //cout<<"condition fullfilled: "<<j<<endl;
            float dummyEnergy1  = energy_total;
            float dummyEnergy2  = beamdEdX[0][j+1] + energyHalfStep;
            float dummyPosition = beam_origin_z;
            int dummyCounter    = j;
            while(dummyPosition<p0[2])  {
              // cout<<"Dummy position: "<<dummyPosition<<endl;
              dummyPosition = dummyPosition + 0.1*defaultCutValue;  // Convert into cm because beamdEdX is given in cm
              dummyEnergy1  = dummyEnergy1 - beamdEdX[1][dummyCounter]*0.1*defaultCutValue * density[kind_target]*1000;
              // cout<<"dummyEnergy1: "<<dummyEnergy1<<endl;
              if(dummyEnergy1/(double)mass < dummyEnergy2){
                dummyCounter++;
                dummyEnergy2 = beamdEdX[0][dummyCounter+1] + energyHalfStep;
              }
            }
            energy_vertex = dummyEnergy1;
          }
        }
      }
      else energy_vertex = 0.0;
      //___________________________________________________________________________________________
      // Energy of fragment at fragmentation point:
      // Calculating the energy spread due to goldhaber
      if(goldhaber_option==1 && energy_vertex >0.)  {	  
        double sigma_goldhaber = sqrt(sigma0_goldhaber * sigma0_goldhaber * mass_f * (mass-mass_f)/(mass-1)); 
        double sigma_momentum  = G4RandGauss::shoot(0,sigma_goldhaber); 	
        
        gamma_temp             = (931.494+energy_vertex/mass) / 931.494;
        float beta_temp        = sqrt(1.0-1.0/gamma_temp/gamma_temp);
        double sigma_energy    = sigma_momentum * gamma_temp * beta_temp;
        
        energy_vertex_stored   = energy_vertex * mass_f/mass + sigma_energy; // Energy of fragment at fragmentation
        energy_vertex          = energy_vertex_stored;                       // cout<<"Energy vertex of fragment: "<<energy_vertex<<endl;
      }
      else  {
        energy_vertex_stored   = energy_vertex * mass_f/mass;
        energy_vertex          = energy_vertex_stored;
      }
      
      double beta_vertex;
      if(borrel_option==1 && energy_vertex>0.)  {
        // Fragment velocity change
	double VelocityChange = sqrt(1 - ((binding_energy_nucleon*(mass-mass_f))/(energy_vertex)));  // Relative change of velocity (Borrel et al.)
	gamma_temp            = (931.494 + energy_vertex/mass_f)/931.494;	                     // Initial gamma value
	beta_vertex           = VelocityChange * sqrt(1.0 - 1.0/gamma_temp/gamma_temp);              // Changed absolute beta
        gamma_temp            = sqrt(1/(1 - beta_vertex * beta_vertex));                             // Changed gamma value
	energy_vertex         = (931.494 * gamma_temp - 931.494)*mass_f;	                     // Changed energy from gamma value
	energy_vertex_stored  = energy_vertex;                                                       // cout<<"Energy vertex including Borrel: "<<energy_vertex<<endl;  
      }
      else{
        gamma_temp  = (931.494 + energy_vertex/mass_f)/931.494;
	beta_vertex = sqrt(1.0 - 1.0/gamma_temp/gamma_temp);
      }
      //End energy at fragmentation point
      //---------------------------------
      //---------------------------------------------------------
      // Broadening in angles due to fragmentation in the target:
      G4double ftheta2 = -10.0;
      G4double fphi2   = -10.0;  
      if(target_broadening_option==1 && energy_vertex>0.) {
        while(ftheta2< 0.0 || ftheta2>180.0 || fphi2<0.0 || fphi2>360.0)  {  
          ftheta2 = G4RandGauss::shoot(0.0,theta_target_broadening_sigma); // cout<<"ftheta2 = "<<ftheta2<<endl;
          fphi2   = CLHEP::RandFlat::shoot(0.0,360.0);                     // cout<<"fphi2 = "<<fphi2<<endl;
        }
        
	fphi2   = fphi2*degree;   // cout<<"fphi2 = "<<fphi2<<endl;
	ftheta2 = ftheta2*degree; // cout<<"ftheta2 = "<<ftheta2<<endl;
        
        pva[0] = cos(fphi2) * tan(ftheta2) + pvb[0];   // vector after the fragmentation					
        pva[1] = sin(fphi2) * tan(ftheta2) + pvb[1];
        pva[2] = cos(ftheta2) * tan(ftheta2) + pvb[2]; // cout<<"z_pva = "<<z_pva<<endl;
        
        double length1 = sqrt(pva[0]*pva[0] + pva[1]*pva[1] + pva[2]*pva[2]); // cout<<"length1 = "<<length1<<endl;
        
        for(int j=0;j<3;j++) pva[j] = pva[j]/length1; // cout<<"x_pva = "<<x_pva<<endl;
        // End Broadening in angles due to fragmentation in the target
        //-----------------------------------------------------------
      }
      else {
        for(int j=0;j<3;j++){ pva[j] = pvb[j];}       
      }
      // Changing tempparticle to the fragment!!!!!!!!!!!!!!!!!!!
      if(deltazet!=0  || deltamass!=0) {
        sprintf(tempparticle,"/gun/ion %i %i %i",z_f,mass_f,charge_f);  
        UI->ApplyCommand("/gun/particle ion");
        UI->ApplyCommand(tempparticle);
      } 	 
      kin->GetParticleGun()->SetParticleEnergy(energy_vertex*MeV);
      kin->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(pva[0]*m,pva[1]*m,pva[2]*m));
      kin->GetParticleGun()->SetParticlePosition(G4ThreeVector(p0[0]*cm,p0[1]*cm,p0[2]*cm));		
      
      // Not really necessary
      steppingaction->SetInteractionPoint(9999.0,9999.0,9999.0);
      
      // Moving the decay point already to the point of fragmentation
      // the later "real" decay point will be added by the time they need to decay times velocity and vector
      for(int j=0;j<3;j++){ g0[j] = p0[j];}       
      
      int numberOfGammas;           // The number of gammas produced in one incident event according to a Poisson distribution
      bool goodEvent;               // To know if the simulated gamma comes from bg or from "good" excitation
      bool goodGammaShot = false;   // To make sure that an excitation occures only ones
      if(atomicBGOption!=0)  {
        //First determine the number of gammas excited
        numberOfGammas = randomGammas->Poisson((double)gAverageNumberOfGammas);
      }
      else numberOfGammas = 1;
      
      // Looping through the number of produced excitations per incident event:
      for(int gammasShot=0;gammasShot<numberOfGammas;gammasShot++)  {
        if(atomicBGOption==0)  {
          goodEvent = 1;
          bgGamma = false;
        }
        else {
          float dummy = CLHEP::RandFlat::shoot(0.0,totalCrossSection);
          if(dummy<= atomicCrossSection)  {
            goodEvent = 0;
            bgGamma = true;
          }
          else  { 
            goodEvent = 1;
            bgGamma = false;
          }
        }
        
        if(goodEvent==1 && goodGammaShot==false)   {
          // Selecting the populated excitation level
          float randNumber       = CLHEP::RandFlat::shoot(0.0,gTotalLevelExcitationProbability);
          int populatedLevelID   = 0;
          int decayIntoLevelID   = 0;
          float excitationEnergy = 1.;
          float decayEnergy      = 0.;
          
          G4double meanlife, lambda, decayTime;
          for(int j=0;j<gNumberOfLevels;j++){
            if(randNumber>=SLevel[j].beginOfTotalExcitationProbability && 
               randNumber<SLevel[j].endOfTotalExcitationProbability)
              populatedLevelID = j; 
          }
          // cout<<"Populated Level:"<<populatedLevelID<<endl;
          //Ok, we have the initialy excited level. Now we have to determine the decay pattern
          while(excitationEnergy != 0. && SLevel[populatedLevelID].numberOfDecayBranches >0 && 
                SLevel[populatedLevelID].totalDecayProbability>0.)  {
            randNumber = CLHEP::RandFlat::shoot(0.0,SLevel[populatedLevelID].totalDecayProbability);
            for(int j=0;j<SLevel[populatedLevelID].numberOfDecayBranches;j++)  {  
              if(randNumber>=SLevel[populatedLevelID].beginOfDecayProbability[j] && 
                 randNumber<SLevel[populatedLevelID].endOfDecayProbability[j])  { 
                decayIntoLevelID = SLevel[populatedLevelID].decayIntoLevelId[j];
                decayEnergy      = SLevel[populatedLevelID].energy - SLevel[decayIntoLevelID].energy;
              
                if(decayEnergy<0.)  {
                  cout<<"Decay energy smaller than zero. Aborting program."<<endl; 
                  abort();
                }
              }
            } 
            
            // The beam has to be shot up to the decay point. For simplicity, I assume that gamma_temp is constant during that time,
            // which causes errors of about 1% or so, depending on the beta spread in the target.
            halflife                     = SLevel[populatedLevelID].halflife;
            meanlife                     = halflife*0.001/log(2.0)*gamma_temp*ns;
            lambda                       = 1.0/(meanlife/ns);
            decayTime                    = CLHEP::RandExponential::shoot(1/lambda)*ns; // Only the decay time for the level!!!
            decay_time_after_interaction = decay_time_after_interaction + decayTime;   // Decay time from reaction point
            
            steppingaction->SetTime(0.);
            steppingaction->SetDecayTime(decayTime);
            steppingaction->SetAfterGammaExcitation(TRUE);
            // Setting the position before starting the beam to get the moved position within the target
            steppingaction->InitialParticlePosition(G4ThreeVector(g0[0]*cm,g0[1]*cm,g0[2]*cm));
            
            float energy_deexcite=0.;
            // cout<<"Starting Beam"<<endl;
            if(dEdXTableOption==0 && energy_vertex>0.)  {
              runManager->BeamOn(1);
              // Getting energy at deexcitation time.
              
              if(steppingaction->GetFlagDecayWithinTarget()==TRUE) 
                energy_deexcite = steppingaction->GetTotalEnergyAtPointOfGammaDecay()/MeV;
              else energy_deexcite = steppingaction->GetTotalEnergyAfterTarget()/MeV;
              
              // Reseting the hi energy for the next gamma decay:
              kin->GetParticleGun()->SetParticleEnergy(energy_deexcite*MeV);
              
              // Getting beta at deexcitation time.
              gamma_temp = (931.494+energy_deexcite/mass_f)/931.494;
              beta_r     = sqrt(1.0 - 1.0/gamma_temp/gamma_temp);
              
              // Getting position at deexcitation time.
              g0[0] = g0[0] + steppingaction->GetMovedTrackInTarget().getX()/cm; // cout<<"x_g0 :"<<x_g0<<endl;
              g0[1] = g0[1] + steppingaction->GetMovedTrackInTarget().getY()/cm; // cout<<"y_g0 :"<<y_g0<<endl;
              g0[2] = g0[2] + steppingaction->GetMovedTrackInTarget().getZ()/cm; // cout<<"z_g0 :"<<z_g0<<endl;
              
              // Reseting the hi position for the next gamma decay:
              kin->GetParticleGun()->SetParticlePosition(G4ThreeVector(g0[0]*cm,g0[1]*cm,g0[2]*cm));
            }
            //_________________________________________________________________________________________
            else if(dEdXTableOption==1 && energy_vertex>0.)  {
              float energyHalfStep = (fragmentdEdX[0][0] - fragmentdEdX[0][1])/2;
              // Have to calculate the energy loss myself:
              for(int j=0;j<50;j++)  {
                if(fragmentdEdX[0][j]+ energyHalfStep >energy_vertex/mass_f &&  
                   fragmentdEdX[0][j]- energyHalfStep <energy_vertex/mass_f)  {
                  float dummyEnergy1  = energy_vertex;
                  float dummyEnergy2  = fragmentdEdX[0][j+1] + energyHalfStep;
                  float dummyTime     = 0.;
                  float dummyPosition = p0[2];
                  int dummyCounter    = j;
                  double dummyGamma   = 0.;
                  double dummyBeta    = 0.;
                  while(dummyTime<decay_time_after_interaction)  {
                    if(dummyPosition>thickness/2.0) break;
                    dummyGamma    = (931.494 + dummyEnergy1/mass_f)/931.494;
                    dummyBeta     = sqrt(1.0 - 1.0/dummyGamma/dummyGamma);
                    dummyTime     = dummyTime + 0.1*defaultCutValue/(dummyBeta*29.97925)/ns;  //have to convert into cm because beamdEdX is given in MeV/cm**2
                    dummyEnergy1  = dummyEnergy1 - fragmentdEdX[1][dummyCounter]*0.1*defaultCutValue * density[kind_target]*1000;
                    dummyPosition = dummyPosition + 0.1*defaultCutValue;            
                    g0[2]          = dummyPosition;
                    if(dummyEnergy1/mass_f < dummyEnergy2){
                      dummyCounter++;
                      dummyEnergy2 = fragmentdEdX[0][dummyCounter+1] + energyHalfStep;
                    }
                  }
                  if(dummyTime<decay_time_after_interaction){
                    g0[2] = g0[2] + (decay_time_after_interaction - dummyTime)/ns * (dummyBeta*29.97925);
                    // cout<<" Never true"<<endl;
                  }
                  energy_deexcite = dummyEnergy1;
                }
              }
              // Getting beta at deexcitation time.
              gamma_temp = (931.494+energy_deexcite/mass_f)/931.494;
              beta_r     = sqrt(1.0 - 1.0/gamma_temp/gamma_temp);
            }
            else beta_r = 0.;
            //_________________________________________________________________________________________
            
            
            //------------------------------------------ 
            //Now comes the deexcitation part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //------------------------------------------
            
            // Determine the theta distribution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            double ran_dist = CLHEP::RandFlat::shoot(0.0,sum_dist);
            for(int mmm=1;mmm<18001;mmm++)  {
              if(ran_dist>=sum_dist_uptonow[mmm-1] && ran_dist<=sum_dist_uptonow[mmm])  {
                theta_gamma_rest = mmm/100.0/180.0*3.14159;
                break;
              }
            }
            
            theta_gamma_lab       = acos((cos(theta_gamma_rest)+beta_r)/(1.0+beta_r*cos(theta_gamma_rest)));
            G4ThreeVector vec_in  = G4ThreeVector(pva[0], pva[1], pva[2]);
            G4ThreeVector vec_out = GetNextVector(vec_in, theta_gamma_lab);
            
            if(vec_out.getX()==-999.9 && vec_out.getY()==-999.9 && vec_out.getZ()==-999.9){
              restart_loop = 1; 
              continue;
             }
            gv[0] = vec_out.getX();
            gv[1] = vec_out.getY();
            gv[2] = vec_out.getZ();
            
            // The energy after the target		
            double energy_after_target = 0;
            if(dEdXTableOption==0) energy_after_target = steppingaction->GetTotalEnergyAfterTarget()/MeV;
            //_________________________________________________________________________________________
            else if(dEdXTableOption==1)  {
              float energyHalfStep = (fragmentdEdX[0][0] - fragmentdEdX[0][1])/2;
              // Have to calculate the energy loss myself:
              for(int j=0;j<50;j++)  {
                if(fragmentdEdX[0][j] + energyHalfStep >energy_vertex/mass_f && 
                   fragmentdEdX[0][j] - energyHalfStep <energy_vertex/mass_f)  {
                  // cout<<"condition fullfilled: "<<j<<endl;
                  float dummyEnergy1  = energy_vertex;
                  float dummyEnergy2  = fragmentdEdX[0][j+1] + energyHalfStep;
                  float dummyPosition = p0[2];
                  int dummyCounter    = j;
                  while(dummyPosition<thickness/2.0)  {
                    dummyEnergy1 = dummyEnergy1 - fragmentdEdX[1][dummyCounter]*0.1*defaultCutValue * density[kind_target]*1000;
                    dummyPosition = dummyPosition + 0.1*defaultCutValue;            
                    if(dummyEnergy1/mass_f < dummyEnergy2)  {
                      dummyCounter++;
                      dummyEnergy2 = fragmentdEdX[0][dummyCounter+1] + energyHalfStep;
                    }
                  }
                  energy_after_target = dummyEnergy1;
                }
              }
            }
            //_________________________________________________________________________________________
            
            gamma_temp = (931.494+energy_after_target/mass_f)/931.494;
            beta_a     = sqrt(1.0 - 1.0/gamma_temp/gamma_temp);

            if(dEdXTableOption==0 && energy_vertex>0.)
              if(steppingaction->GetFlagDecayWithinTarget()==FALSE)  {
                G4double timeThrougTarget = steppingaction->GetTime();
                g0[0] = g0[0] + (decayTime-timeThrougTarget)/ns*29.97925*beta_r*pva[0];
                g0[1] = g0[1] + (decayTime-timeThrougTarget)/ns*29.97925*beta_r*pva[1];
                g0[2] = g0[2] + (decayTime-timeThrougTarget)/ns*29.97925*beta_r*pva[2];
                // Reseting the hi position for the next gamma decay:
                kin->GetParticleGun()->SetParticlePosition(G4ThreeVector(g0[0]*cm,g0[1]*cm,g0[2]*cm));
              }

            e_rest    = decayEnergy;
            e_doppler = e_rest*(sqrt(1.0-beta_r*beta_r)/(1.0-beta_r*cos(theta_gamma_lab)));
            
            evnum = i;
            
            t->Fill(); // cout<<"Filling root tree"<<endl;
            
            //Setting the new excitation energy and populated level:
            excitationEnergy = SLevel[decayIntoLevelID].energy;
            populatedLevelID = SLevel[decayIntoLevelID].ID;                   
            decayIDevnum++;   
          }// end while loop through decay scheme
          goodGammaShot=true;
        }
        else if(goodEvent==0)  {
          //What variables do I have to fill?
          double atomicBGEnery;
          double atomicBGTheta;
          h_atomicBackground->GetRandom2(atomicBGTheta,atomicBGEnery);
          float atomicBGPhi = CLHEP::RandFlat::shoot(0., 360.);
          //have to transfer it to the correct coordinate system:

          atomicBGTheta = atomicBGTheta*degree;
          atomicBGPhi = atomicBGPhi*degree;
          theta_gamma_rest = atomicBGTheta;
          theta_gamma_lab = atomicBGTheta;

          gv[0] = sin(atomicBGTheta)*cos(atomicBGPhi);  
          gv[1] = sin(atomicBGTheta)*sin(atomicBGPhi);  
          gv[2] = cos(atomicBGTheta);           

          e_rest = atomicBGEnery;
          e_doppler = e_rest;
          halflife = 0;
          beta_r = 0; 
          beta_a = 0;
          decayIDevnum = -1;
          evnum = i;
          //Filling the tree:
          t->Fill(); 
          //Resetting the Id to zero for possible good decay
          decayIDevnum  = 0.;
        }
        else continue;         
      }// end the loop throug the number of gammas
    }// end if
    if(restart_loop==1) continue; 
    i = i + 1;
  } // end while loop through events  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......        
G4ThreeVector GetNextVector(G4ThreeVector vector_in, float theta)  {
  double x1     = vector_in.getX();
  double y1     = vector_in.getY();
  double z1     = vector_in.getZ();
  double l1     = sqrt(x1*x1 + y1*y1 + z1*z1);
  double theta1 = acos(z1/l1);
  double phi1   = 0.0;

  if(sqrt((x1/l1/sin(theta1))*(x1/l1/sin(theta1)))>1.0)  {
    x1     = 0.0;
    y1     = 0.0;
    z1     = 1.0;
    theta1 = 0.0;
    phi1   = 0.0;
  }
  if(sqrt((x1/l1/sin(theta1))*(x1/l1/sin(theta1)))<=1.0)  {
    phi1 = acos(x1/l1/sin(theta1));
  }

  G4RotationMatrix rot3D;
  G4ThreeVector rotatedPos;

  G4double g4phi    = CLHEP::RandFlat::shoot(0.0,360.0)*degree;
  G4double g4theta  = theta*rad;
  G4double g4theta1 = theta1*rad;
  G4double g4phi1   = phi1*rad;

  rot3D.set(0, 0, 0);
  rot3D.rotateY(g4theta);  
  rot3D.rotateZ(g4phi);
  rot3D.rotateY(g4theta1);  
  rot3D.rotateZ(g4phi1);
	      
  rotatedPos = rot3D(G4ThreeVector(0.0,0.0,1.0));
  return rotatedPos;
}
