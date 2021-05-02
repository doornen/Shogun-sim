#include "Globals.hh"
#include "MyNamespace.hh"
 
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"

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
#include "TMath.h"
#include "G4ios.hh"

#include "CLHEP/Random/Randomize.h"

#include <iostream>

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
void FillRootTreeHeader();
void RunSimulation();
void SetDetDefaultValues();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc,char** argv)  {

  cout<<"Entering main function ......................."<<endl;

  ReadInputFile();                                                              // Reading the input file
 
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);                      // Choose the Random engine
  G4long seed;
  seed = time(0);  CLHEP::HepRandom::setTheSeed(seed);


  G4VSteppingVerbose::SetInstance(new SteppingVerbose);                         // My Verbose output class
    
  runManager = new G4RunManager;                                                // Construct the default run manager

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif   

  runManager->SetUserInitialization(det = new DetectorConstruction);            // Set mandatory initialization classes
  runManager->SetUserInitialization(phys = new PhysicsList);

  runManager->SetUserAction(kin = new PrimaryGeneratorAction(det));             // Set user action classes
  runManager->SetUserAction(run = new RunAction(det,phys,kin));                 
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new TrackingAction(run));
  runManager->SetUserAction(steppingaction = new SteppingAction(det,run));

  UI = G4UImanager::GetUIpointer();
 
  det->SetTargetSize((double)target_X,(double)target_Y,(double)thickness); 
  det->SetDetectorInclude(dali2_opt,Shogun_opt,grape_opt,sgt_opt,sphere_opt,LaBr3Opt);
  det->SetBeamPipeInclude(beam_pipe_opt,(20.*bp_inside),(20.*bp_outside)); //the radius is given in cm, but its made with diam. in mm
  det->SetSTQInclude(stq_opt);
  det->SetTargetHolderInclude(target_holder_opt);
  det->SetCollimatorInclude(collimator_opt);
 
  UI->ApplyCommand("/testem/phys/addPhysics standard");
  UI->ApplyCommand("/run/initialize");
 
  if(dali2_opt==1)  {
    // The absorbers are now part of the beam pipe!!!!
    // Was never fully implemented anyway ;)
    //det->fDali2Array->SetPbAbsorberThickness(dali2_thickness_Pb);
    //det->fDali2Array->SetSnAbsorberThickness(dali2_thickness_Sn);
    //det->fDali2Array->SetAlAbsorberThickness(dali2_thickness_Al);
    det->fDali2Array->SetEnergyResolution(dali2_en_res_opt,dali2_en_res[0],dali2_en_res[1]);  // The default energy resolution for dali
    det->fDali2Array->SetEnergyResolutionInd(dali2_en_res_ind);             // Option to set Energy resolution individually
    det->fDali2Array->SetTimeResolution(dali2_time_res[0],dali2_time_res[1]);
    det->fDali2Array->SetZPosShift(z_pos_shift);
    det->fDali2Array->CreateArrayXYZPsiThetaPhi();
    det->fDali2Array->SetFIPointOption(dali2_fi_opt);
  }
  if(Shogun_opt==1)  {
    if(Shogun_housing_opt==1)det->fShogunArray->SetHousingThickness(ShogunHousingThicknessX,ShogunHousingThicknessY,ShogunHousingThicknessZ);
    if(Shogun_MgO_opt==1)det->fShogunArray->SetMgOThickness(ShogunMgOThicknessX,ShogunMgOThicknessY,ShogunMgOThicknessZ);
    det->fShogunArray->SetEnergyResolution(Shogun_en_res_opt,Shogun_en_res[0],Shogun_en_res[1]);
    det->fShogunArray->SetTimeResolution(Shogun_time_res[0],Shogun_time_res[1]);
    det->fShogunArray->SetZPosShift(z_pos_shift);
    det->fShogunArray->CreateArrayXYZPsiThetaPhi();
    det->fShogunArray->SetFIPointOption(Shogun_fi_opt);
  } 
  if(grape_opt==1)  {
    //det->fGrapeArray->SetPbAbsorberThickness(grape_thickness_Pb);
    //det->fGrapeArray->SetSnAbsorberThickness(grape_thickness_Sn);
    //det->fGrapeArray->SetAlAbsorberThickness(grape_thickness_Al);
    det->fGrapeArray->SetEnergyResolution(grape_en_res_opt,grape_en_res[0],grape_en_res[1]);
    det->fGrapeArray->SetTimeResolution(grape_time_res[0],grape_time_res[1]);
    det->fGrapeArray->CreateArrayXYZPsiThetaPhi();
  }
  if(sgt_opt==1)  {
    //det->fSGTArray->SetPbAbsorberThickness(sgt_thickness_Pb);
    //det->fSGTArray->SetSnAbsorberThickness(sgt_thickness_Sn);
    //det->fSGTArray->SetAlAbsorberThickness(sgt_thickness_Al);
    det->fSGTArray->SetEnergyResolution(sgt_en_res_opt,sgt_en_res[0],sgt_en_res[1]);
    det->fSGTArray->SetTimeResolution(sgt_time_res[0],sgt_time_res[1]);
    det->fSGTArray->CreateArrayXYZPsiThetaPhi();
  }
  if(sphere_opt==1)  {
    det->fSphere->SetDimensions(sphere_outer_radius,sphere_inner_radius);
    det->fSphere->CreateSphere();
  }
  if(LaBr3Opt==1)  {
    if(LaBr3HousingOpt==1)det->fLaBr3Array->SetHousingThickness(LaBr3HousingThicknessFront,LaBr3HousingThicknessSide);
    if(LaBr3InsulationOpt==1)det->fLaBr3Array->SetInsulationThickness(LaBr3InsulationThicknessFront,LaBr3InsulationThicknessSide);
    det->fLaBr3Array->SetEnergyResolution(LaBr3EnResOpt,LaBr3EnRes[0],LaBr3EnRes[1]);
    det->fLaBr3Array->SetTimeResolution(LaBr3TimeRes[0],LaBr3TimeRes[1]);
    det->fLaBr3Array->CreateArrayXYZPsiThetaPhi();
  }

  if(collimator_opt==1)
    det->fCollimator->CreateCollimator();

  if(beam_pipe_opt==1)  {
    det->fBeamPipe->SetShieldRadius(ShieldRadius);
    det->fBeamPipe->SetPbShieldThickness(PbShieldThickness);
    det->fBeamPipe->SetSnShieldThickness(SnShieldThickness);
    det->fBeamPipe->ConstructShield();
  }

  if(kind_target==0) det->SetTargetMaterial("Air");
  if(kind_target==1) det->SetTargetMaterial("Au");
  if(kind_target==2) det->SetTargetMaterial("Be");
  if(kind_target==3) det->SetTargetMaterial("C");
  if(kind_target==4) det->SetTargetMaterial("Fe");
  if(kind_target==5) det->SetTargetMaterial("Pb");
  if(kind_target==6) det->SetTargetMaterial("LH2");
  if(kind_target==7) det->SetTargetMaterial("Zr");
  if(kind_target==8) det->SetTargetMaterial("CH2");

  if(target_holder_opt==1)
    det->fTargetHolder->CreateTargetHolder();

  cout<<"creating output ROOT-tree"<<endl;
  outfile = new TFile(root_out,"RECREATE");
  outfile->cd();

  t2 = new TTree("ObservedEvents","ObservedEvents");
  t2->Branch("EventNumber",&evnum,"evnum/F");
  t2->Branch("DecayIDOfEventNumber",&decayIDevnum,"decayIDevnum/F");
  t2->Branch("ProjectileVertex",p0,"p0[3]/F");                                  // Position at the fragmentation point
  t2->Branch("VProjectileBeforeTarget",pvb,"pvb[3]/F");                         // Normalized Vector of beam before the target
  t2->Branch("VProjectileAfterTarget",pva,"pva[3]/F");                          // Normalized Vector of beam after the target
  t2->Branch("EnergyProjectile",&energy_p,"energy_p/F");                        // Energy of beam before the target in MeV/u
  t2->Branch("BetaBeforeTarget",&beta_b,"beta_b/F");                            // Beta before the target
  t2->Branch("BetaReal",&beta_r,"beta_r/F");                                    // Beta during deexcitation	
  t2->Branch("BetaAfterTarget",&beta_a,"beta_a/F");                             // Beta After Target
  t2->Branch("Halflife",&halflife,"halflife/F");                                // Halflife
  t2->Branch("DecayTimeAfterInteraction",&decay_time_after_interaction,"decay_time_after_interaction/F");
  t2->Branch("VertexGamma",g0,"g0[3]/F");                                       // Position at the gamma emmittance point
  t2->Branch("VGamma",gv,"gv[3]/F");		                                // Gamma vector
  t2->Branch("EGammaRest",&e_rest,"e_rest/F");		                        // Energy at rest
  t2->Branch("EGammaDoppler",&e_doppler,"e_doppler/F");                         // Theta of doppler boosted gamma
  t2->Branch("ThetaGammaRest",&theta_gamma_rest,"theta_gamma_rest/F");
  t2->Branch("ThetaGammaLab",&theta_gamma_lab,"theta_gamma_lab/F");
  t2->Branch("EnergyVertex",&energy_vertex_stored,"energy_vertex_stored/F");    // Energy of fragment at fragmentation
  t2->Branch("BGGamma",&bgGamma,"bgGamma/O");                                    // To know if the gamma is from bg or a good event
  // So far identical to the ROOT-Tree of step one
  //The new stuff:
  //Setting the type of the gamma detector, in order to distinguish the det. types
  t2->Branch("GammaDetType",&gamma_det_type,"gamma_det_type/I");
  //-------------------------------
  //DALI:
  if(dali2_opt==1)  {
    sprintf(temp,"dali2_flag[%i]/I",NUMBEROFDALI2CRYSTALS);
    t2->Branch("DALI2Flag",dali2_flag,temp);
    sprintf(temp,"dali2_energy_not_cor[%i]/F",NUMBEROFDALI2CRYSTALS);
    t2->Branch("DALI2EnergyNotCor",dali2_energy_not_cor,temp);
    sprintf(temp,"dali2_time[%i]/F",NUMBEROFDALI2CRYSTALS);
    t2->Branch("DALI2Time",dali2_time,temp);
    if(dali2_fi_opt==1)  {
      sprintf(temp,"dali2_fi[3][%i]/F",NUMBEROFDALI2CRYSTALS);
      t2->Branch("DALI2FI",dali2_fi,temp);
    }
  }
  //-------------------------------
  //Shogun:
  if(Shogun_opt==1)  {
    sprintf(temp,"Shogun_flag[%i]/I",NUMBEROFSHOGUNCRYSTALS);
    t2->Branch("ShogunFlag",Shogun_flag,temp);
    sprintf(temp,"Shogun_energy_not_cor[%i]/F",NUMBEROFSHOGUNCRYSTALS);
    t2->Branch("ShogunEnergyNotCor",Shogun_energy_not_cor,temp);
    sprintf(temp,"Shogun_time[%i]/F",NUMBEROFSHOGUNCRYSTALS);
    t2->Branch("ShogunTime",Shogun_time,temp);
    if(Shogun_fi_opt==1)  {
      sprintf(temp,"Shogun_fi[3][%i]/F",NUMBEROFSHOGUNCRYSTALS);
      t2->Branch("ShogunFI",Shogun_fi,temp);
    }
  }
  //-----------------------------------------------------------
  //GRAPE:
  if(grape_opt==1)  {
    t2->Branch("GRAPEFlag",grape_flag,"grape_flag[18]/I");
    t2->Branch("GRAPECrystalFlag",grape_crystal_flag,"grape_crystal_flag[18][2][10]/I");
    t2->Branch("GRAPEEnergyNotCor",grape_energy_not_cor,"grape_energy_not_cor[18][2][10]/F");
    t2->Branch("GRAPETime",grape_time,"grape_time[18][2][10]/F");
  }
  //-----------------------------------------------------------
  //SGT
  if(sgt_opt==1)  {
    t2->Branch("SGTFlag",sgt_flag,"sgt_flag[10]/I");
    t2->Branch("SGTCoaxCrystalFlag",sgt_coax_crystal_flag,"sgt_coax_crystal_flag[10]/I");
    t2->Branch("SGTPlanarStripFlag",sgt_planar_strip_flag,"sgt_planar_strip_flag[10][26]/I");// Central plus 25 strips
    t2->Branch("SGTEnergyNotCor",sgt_energy_not_cor,"sgt_energy_not_cor[10][27]/F");
    t2->Branch("SGTTime",sgt_time,"sgt_time[10][27]/F");
  }
  //Sphere
  if(sphere_opt==1)  {
    //fi: First interaction
    t2->Branch("SphereFI",sphere_fi,"sphere_fi[3]/F");                          
    t2->Branch("sphereEnergyNotCor",&sphere_energy_not_cor,"sphere_energy_not_cor/F");
  }
  //LaBr3Array:
  if(LaBr3Opt==1)  {
    sprintf(temp,"LaBr3Flag[%i]/I",NUMBEROFLABR3ARRAYCRYSTALS);
    t2->Branch("LaBr3Flag",LaBr3Flag,temp);
    sprintf(temp,"LaBr3EnergyNotCor[%i]/F",NUMBEROFLABR3ARRAYCRYSTALS);
    t2->Branch("LaBr3EnergyNotCor",LaBr3EnergyNotCor,temp);
    sprintf(temp,"LaBr3Time[%i]/F",NUMBEROFLABR3ARRAYCRYSTALS);
    t2->Branch("LaBr3Time",LaBr3Time,temp);
  }
  //-------------------------------
  t2->Branch("PosDet",pos_det,"pos_det[3]/F");                                  // Reconstructed after the secondary target
  t2->Branch("VertexReconstructed",ver_rec,"ver_rec[3]/F");                     // Reconstructed vertex position, 
  t2->Branch("BetaReconstructed",&beta_rec,"beta_rec/F");                       // Reconstructed beta, including resolution

  //Making the header file with information that is needed only once
  outfile->cd();
  t2Header = new TTree("Header","Header");
  
  t2Header->Branch("Mass",&mass,"mass/I");                                      // Beam mass
  t2Header->Branch("Z",&z,"z/I");                                               // Beam z	
  t2Header->Branch("Charge",&charge,"charge/I");                                // Beam charge	
  t2Header->Branch("MassFragment",&mass_f,"mass_f/I");                          // Fragment mass
  t2Header->Branch("ZFragment",&z_f,"z_f/I");	                                // Fragment z
  t2Header->Branch("ChargeFragment",&charge_f,"charge_f/I");                    // Fragment charge
  t2Header->Branch("TargetKind",&kind_target,"kind_target/I");
  t2Header->Branch("TargetThicknessCM",&thickness,"thickness/F");
  t2Header->Branch("TotalEventNumber",&total_event_num,"total_event_num/F");
  t2Header->Branch("ThetaLow",&theta_low,"theta_low/F");
  t2Header->Branch("ThetaHigh",&theta_high,"theta_high/F");

  //DALI:
  if(dali2_opt==1)  {
    t2Header->Branch("DALI2EnResOpt",&dali2_en_res_opt,"dali2_en_res_opt/I");
    t2Header->Branch("DALI2EnRes",dali2_en_res,"dali2_en_res[2]/F");
    t2Header->Branch("DALI2TimeRes",dali2_time_res,"dali2_time_res[2]/F");
    sprintf(temp,"dali2_Pos[3][%i]/F",NUMBEROFDALI2CRYSTALS);
    t2Header->Branch("DALI2Pos",dali2_pos,temp);
  }
  //Shogun
  if(Shogun_opt==1)  {
    t2Header->Branch("ShogunEnResOpt",&Shogun_en_res_opt,"Shogun_en_res_opt/I");
    t2Header->Branch("ShogunEnRes",Shogun_en_res,"Shogun_en_res[2]/F");
    t2Header->Branch("ShogunTimeRes",Shogun_time_res,"Shogun_time_res[2]/F");
    sprintf(temp,"ShogunPos[3][%i]/F",NUMBEROFSHOGUNCRYSTALS);
    t2Header->Branch("ShogunPos",Shogun_pos,temp);
  }
  //GRAPE:
  if(grape_opt==1)  {
    t2Header->Branch("GRAPEEnResOpt",&grape_en_res_opt,"grape_en_res_opt/I");
    t2Header->Branch("GRAPEEnRes",grape_en_res,"grape_en_res[2]/F");
    t2Header->Branch("GRAPETimeRes",grape_time_res,"grape_time_res[2]/F");
    sprintf(temp,"grape_pos[3][%i][2][10]/F",NUMBEROFGRAPEDETECTORS);
    t2Header->Branch("GRAPEPos",grape_pos,temp);
  }
  //SGT
  if(sgt_opt==1)  {
    t2Header->Branch("SGTEnResOpt",&sgt_en_res_opt,"sgt_en_res_opt/I");
    t2Header->Branch("SGTEnRes",sgt_en_res,"sgt_en_res[2]/F");
    t2Header->Branch("SGTTimeRes",sgt_time_res,"sgt_time_res[2]/F");
    sprintf(temp,"sgt_pos[3][%i][27]/F",NUMBEROFSGTDETECTORS);
    t2Header->Branch("SGTPos",sgt_pos,temp);                   // First the coax, then the fictive planar central, then the strips
  }
  //LaBr3
  if(LaBr3Opt==1)  {
    t2Header->Branch("LaBr3EnResOpt",&LaBr3EnResOpt,"LaBr3EnResOpt/I");
    t2Header->Branch("LaBr3EnRes",LaBr3EnRes,"LaBr3EnRes[2]/F");
    t2Header->Branch("LaBr3TimeRes",LaBr3TimeRes,"LaBr3TimeRes[2]/F");
    sprintf(temp,"LaBr3Pos[3][%i]/F",NUMBEROFLABR3ARRAYCRYSTALS);
    t2Header->Branch("LaBr3Pos",LaBr3Pos,temp);
  }
  t2Header->Branch("BetaResolution",&beta_res,"beta_res/F");                    // Beta and position resolution
  t2Header->Branch("PosDetAtTargetRes",&pos_det_at_target_res,"pos_det_at_target_res/F"); 
  t2Header->Branch("PosDetAfterTargetRes",&pos_det_after_target_res,"pos_det_after_target_res/F"); 

  if(argc==1)  {                                                                // Define (G)UI terminal for interactive mode   
    G4UIsession * session = 0;                                                  // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif    

    UI->ApplyCommand("/control/execute vis.mac");    
    session->SessionStart();
    delete session;
  }
  else  {                                                                       // Batch mode 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }
  UI->ApplyCommand("/vis/scene/notifyHandlers");

  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 0");

  total_event_num = nentries;
  RunSimulation();
  cout<<"Finished the simulation"<<endl;                                        // Writing stuff to file and spectra:
  FillRootTreeHeader();
  outfile->cd();
  t2Header->Write();
  t2->Write();
  cout<<"Everything written to the Root-trees"<<endl;

  UI->ApplyCommand("/vis/viewer/update");
#ifdef G4VIS_USE
  delete visManager;
#endif  

  outfile->Close();
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void ReadInputFile()  {
  float dummy1, dummy2;
  FILE *fin = fopen("./input/EventBuilder.in","r");
  // All input resolutions must be given in FWHM!!!!!!!!!!!!!!!!!!!
  while(!feof(fin)){
    fscanf(fin,"%s ",temp); 
    if(strcmp(temp,"INPUTFILE")==0)  {
      fscanf(fin,"%s ",&root_in); printf("%s %s \n",temp,root_in);
    }
    else if(strcmp(temp,"OUTPUTFILE")==0)  {
      fscanf(fin,"%s ",&root_out); 
      printf("%s %s \n",temp,root_out);
    }
    else if(strcmp(temp,"SHIELD")==0)  { // The inner radius of the shield!!!!
      fscanf(fin,"%f %f %f",&ShieldRadius,&PbShieldThickness,&SnShieldThickness); 
      printf("%s %f %f %f\n",temp,ShieldRadius,PbShieldThickness,SnShieldThickness);
    }
    else if(strcmp(temp,"DALI2INCLUDE")==0)  {
      fscanf(fin,"%i ",&dali2_opt); 
      printf("%s %i \n",temp,dali2_opt);
    }
    else if(strcmp(temp,"ZPOSSHIFT")==0)  {
      fscanf(fin,"%f ",&z_pos_shift); 
      printf("%s %f \n",temp,z_pos_shift);
    }
    else if(strcmp(temp,"DALI2FIINCLUDE")==0)  {
      fscanf(fin,"%i ",&dali2_fi_opt); 
      printf("%s %i \n",temp,dali2_fi_opt);
    }
    else if(strcmp(temp,"SHOGUNINCLUDE")==0)  {
      fscanf(fin,"%i ",&Shogun_opt); 
      printf("%s %i \n",temp,Shogun_opt);
    }
    else if(strcmp(temp,"SHOGUNFIINCLUDE")==0)  {
      fscanf(fin,"%i ",&Shogun_fi_opt); 
      printf("%s %i \n",temp,Shogun_fi_opt);
    }
    else if(strcmp(temp,"SHOGUNHOUSINGTHICKNESSXYZ")==0)  {
      fscanf(fin,"%i %f %f %f",&Shogun_housing_opt,&ShogunHousingThicknessX,&ShogunHousingThicknessY,&ShogunHousingThicknessZ); 
      printf("%s %i %f %f %f\n",temp,Shogun_housing_opt,ShogunHousingThicknessX,ShogunHousingThicknessY,ShogunHousingThicknessZ);
    }
    else if(strcmp(temp,"SHOGUNMGOTHICKNESSXYZ")==0)  {
      fscanf(fin,"%i %f %f %f",&Shogun_MgO_opt,&ShogunMgOThicknessX,&ShogunMgOThicknessY,&ShogunMgOThicknessZ); 
      printf("%s %i %f %f %f\n",temp,Shogun_MgO_opt,ShogunMgOThicknessX,ShogunMgOThicknessY,ShogunMgOThicknessZ);
    }
    else if(strcmp(temp,"LABR3INCLUDE")==0)  {
      fscanf(fin,"%i ",&LaBr3Opt); 
      printf("%s %i \n",temp,LaBr3Opt);
    }
    else if(strcmp(temp,"LABR3HOUSINGTHICKNESS")==0)  {
      fscanf(fin,"%i %f %f",&LaBr3HousingOpt,&LaBr3HousingThicknessFront,&LaBr3HousingThicknessSide); 
      printf("%s %i %f %f\n",temp,LaBr3HousingOpt,LaBr3HousingThicknessFront,LaBr3HousingThicknessSide);
    }
    else if(strcmp(temp,"LABR3INSULATIONTHICKNESS")==0)  {
      fscanf(fin,"%i %f %f",&LaBr3InsulationOpt,&LaBr3InsulationThicknessFront,&LaBr3InsulationThicknessSide); 
      printf("%s %i %f %f\n",temp,LaBr3InsulationOpt,LaBr3InsulationThicknessFront,LaBr3InsulationThicknessSide);
    }
    else if(strcmp(temp,"GRAPEINCLUDE")==0)  {
      fscanf(fin,"%i ",&grape_opt); 
      printf("%s %i \n",temp,grape_opt);
    }
    else if(strcmp(temp,"SGTINCLUDE")==0)  {
      fscanf(fin,"%i ",&sgt_opt); 
      printf("%s %i \n",temp,sgt_opt);
    }
    else if(strcmp(temp,"SPHEREINCLUDE")==0)  {
      fscanf(fin,"%i %f %f",&sphere_opt,&sphere_outer_radius,&sphere_inner_radius); 
      printf("%s %i %f %f\n",temp,sphere_opt,sphere_outer_radius,sphere_inner_radius);
    }
    else if(strcmp(temp,"DALI2ENERGYRESOLUTION")==0)  {  
      fscanf(fin,"%i %f %f",&dali2_en_res_opt,&dummy1,&dummy2);
      dali2_en_res[0] = dummy1/2.35; dali2_en_res[1] = dummy2;
      printf("%s %i %f %f\n",temp,dali2_en_res_opt,dali2_en_res[0],dali2_en_res[1]);
    }
    else if(strcmp(temp,"DALI2ENERGYRESOLUTIONINDIVIDUAL")==0)  { //2011-11-07:Pieter:was sigma -> changed to FWHM 
      fscanf(fin,"%i",&dali2_en_res_ind);
      printf("%s %i\n",temp,dali2_en_res_ind);
    }
    else if(strcmp(temp,"SHOGUNENERGYRESOLUTION")==0)  {
      fscanf(fin,"%i %f %f",&Shogun_en_res_opt,&dummy1,&dummy2); 
      Shogun_en_res[0] = dummy1/2.35; Shogun_en_res[1] = dummy2;
      printf("%s %i %f %f\n",temp,Shogun_en_res_opt,Shogun_en_res[0],Shogun_en_res[1]);
    }
    else if(strcmp(temp,"GRAPEENERGYRESOLUTION")==0)  {
      fscanf(fin,"%i %f %f",&grape_en_res_opt,&dummy1,&dummy2); 
      grape_en_res[0] = dummy1/2.35; grape_en_res[1] = dummy2;
      printf("%s %i %f %f\n",temp,grape_en_res_opt,grape_en_res[0],grape_en_res[1]);
    }
    else if(strcmp(temp,"SGTENERGYRESOLUTION")==0)  {
      fscanf(fin,"%i %f %f",&sgt_en_res_opt,&dummy1,&dummy2);
      sgt_en_res[0] = dummy1/2.35; sgt_en_res[1] = dummy2;
      printf("%s %i %f %f\n",temp,sgt_en_res_opt,sgt_en_res[0],sgt_en_res[1]);
    }
    else if(strcmp(temp,"LABR3ENERGYRESOLUTION")==0)  {
      fscanf(fin,"%i %f %f",&LaBr3EnResOpt,&dummy1,&dummy2); 
      LaBr3EnRes[0] = dummy1/2.35; LaBr3EnRes[1] = dummy2;
      printf("%s %i %f %f\n",temp,Shogun_en_res_opt,LaBr3EnRes[0],LaBr3EnRes[1]);
    }
    else if(strcmp(temp,"DALI2TIMERESOLUTION")==0)  {
      fscanf(fin,"%f %f",&dummy1,&dummy2);
      dali2_time_res[0] = dummy1/2.35;dali2_time_res[1] = dummy2;
      printf("%s %f %f\n",temp,dali2_time_res[0],dali2_time_res[1]);
    }
    else if(strcmp(temp,"SHOGUNTIMERESOLUTION")==0)  {
      fscanf(fin,"%f %f",&dummy1,&dummy2); 
      Shogun_time_res[0] = dummy1/2.35; Shogun_time_res[1] = dummy2/2.35;
      printf("%s %f %f\n",temp,Shogun_time_res[0],Shogun_time_res[1]);
    }
    else if(strcmp(temp,"GRAPETIMERESOLUTION")==0)  {
      fscanf(fin,"%f %f",&dummy1,&dummy2);
      grape_time_res[0] = dummy1/2.35;grape_time_res[1] = dummy2/2.35;
      printf("%s %f %f\n",temp,grape_time_res[0],grape_time_res[1]);
    }
    else if(strcmp(temp,"SGTTIMERESOLUTION")==0)  {
      fscanf(fin,"%f %f",&dummy1,&dummy2);
      sgt_time_res[0] = dummy1/2.35;sgt_time_res[1] = dummy2/2.35;
      printf("%s %f %f\n",temp,sgt_time_res[0],sgt_time_res[1]);
    }
    else if(strcmp(temp,"LABR3TIMERESOLUTION")==0)  {
      fscanf(fin,"%f %f",&dummy1,&dummy2); 
      LaBr3TimeRes[0] = dummy1/2.35; LaBr3TimeRes[1] = dummy2/2.35;
      printf("%s %f %f\n",temp,LaBr3TimeRes[0],LaBr3TimeRes[1]);
    }
    else if(strcmp(temp,"POSDETECTORONTARGETRESOLUTION")==0)  {
      fscanf(fin,"%f",&pos_det_at_target_res);
      pos_det_at_target_res = pos_det_at_target_res/2.35;
      printf("%s %f\n",temp,pos_det_at_target_res);
    }
    else if(strcmp(temp,"ENERGYDETECTORAFTERTARGETINCLUDE")==0)  {
      fscanf(fin,"%i",&en_det_after_target_opt);
      printf("%s %i\n",temp,en_det_after_target_opt);
    }
    else if(strcmp(temp,"POSDETECTORAFTERTARGETDISTANCE")==0)  {
      fscanf(fin,"%f",&pos_det[2]); 
      printf("%s %f\n",temp,pos_det[2]);
    }
    else if(strcmp(temp,"POSDETECTORAFTERTARGETRESOLUTION")==0)  {
      fscanf(fin,"%f",&pos_det_after_target_res);
      pos_det_after_target_res = pos_det_after_target_res/2.35;
      printf("%s %f\n",temp,pos_det_after_target_res);
    }
    else if(strcmp(temp,"BETARESOLUTION")==0)  {
      fscanf(fin,"%f",&beta_res); 
      beta_res = beta_res/2.35;
      printf("%s %f\n",temp,beta_res);
    }
    else if(strcmp(temp,"BEAMPIPEINCLUDE")==0)  {
      fscanf(fin,"%i %lf %lf",&beam_pipe_opt,&bp_inside,&bp_outside);
      printf("%s %i %lf %lf \n",temp,beam_pipe_opt,bp_inside, bp_outside);
    }
    else if(strcmp(temp,"TARGETHOLDERINCLUDE")==0)  {
      fscanf(fin,"%i",&target_holder_opt); 
      printf("%s %i\n",temp,target_holder_opt);
    }
    else if(strcmp(temp,"STQINCLUDE")==0)  {
      fscanf(fin,"%i",&stq_opt); 
      printf("%s %i\n",temp,stq_opt);
    }
    else if(strcmp(temp,"COLLIMATORINCLUDE")==0)  {
      fscanf(fin,"%i",&collimator_opt); 
      printf("%s %i\n",temp,collimator_opt);
    }
    else if(strcmp(temp,"END")==0) break;
    //else if(strcmp(temp[0],"#")==0) break;
    else  {
      cout<<"Could not read your input keyword. Aborting program."<<endl; 
      abort();
    }
  }
  fclose(fin);

  infile = new TFile(root_in,"READ");
  infile->cd();
  tHeader = (TTree*)infile->Get("Header");                                      // Reading the Header data:
  tHeader->SetBranchAddress("Mass",&mass);                                      // Beam mass
  tHeader->SetBranchAddress("Z",&z);                                            // Beam z	
  tHeader->SetBranchAddress("Charge",&charge);                                  // Beam charge	
  tHeader->SetBranchAddress("MassFragment",&mass_f);                            // Fragment mass
  tHeader->SetBranchAddress("ZFragment",&z_f);	                                // Fragment z
  tHeader->SetBranchAddress("ChargeFragment",&charge_f);                        // Fragment charge
  tHeader->SetBranchAddress("TargetKind",&kind_target);
  tHeader->SetBranchAddress("TargetThicknessCM",&thickness);
  tHeader->SetBranchAddress("ThetaLow",&theta_low);                             // Theta angle covered by simulation
  tHeader->SetBranchAddress("ThetaHigh",&theta_high);
  tHeader->GetEntry(0);
  cout<<"the thickness is: "<<thickness<<endl; 

  t = (TTree*)infile->Get("Events");                                            // Setting the address of event-by-event data.
  t->SetBranchAddress("EventNumber",&evnum);  
  t->SetBranchAddress("DecayIDOfEventNumber",&decayIDevnum);
  t->SetBranchAddress("ProjectileVertex",p0);                                   // X Position at the fragmentation point
  t->SetBranchAddress("VProjectileBeforeTarget",pvb);                           // Normalized Vector of beam before the target
  t->SetBranchAddress("VProjectileAfterTarget",pva);                            // Normalized Vector of beam after the target
  t->SetBranchAddress("EnergyProjectile",&energy_p);                            // Energy of beam before the target in MeV/u
  t->SetBranchAddress("BetaBeforeTarget",&beta_b);                              // Beta before the target
  t->SetBranchAddress("BetaReal",&beta_r);                                      // Beta during deexcitation	
  t->SetBranchAddress("BetaAfterTarget",&beta_a);                               // Beta After Target
  t->SetBranchAddress("Halflife",&halflife);                                    // Halflife
  t->SetBranchAddress("DecayTimeAfterInteraction",&decay_time_after_interaction);
  t->SetBranchAddress("VertexGamma",g0);                                        // Position at the gamma emmittance point
  t->SetBranchAddress("VGamma",gv);                                             // Gamma vector
  t->SetBranchAddress("EGammaRest",&e_rest);	                                // Energy of gamma at rest
  t->SetBranchAddress("EGammaDoppler",&e_doppler);                              // Energy of doppler boosted gamma
  t->SetBranchAddress("ThetaGammaRest",&theta_gamma_rest);                      // Theta of gamma at rest
  t->SetBranchAddress("ThetaGammaLab",&theta_gamma_lab);                        // Theta of doppler boosted gamma
  t->SetBranchAddress("EnergyVertex",&energy_vertex_stored);                    // Energy of fragment at fragmentation
  t->SetBranchAddress("BGGamma",&bgGamma);                                      // To know if the gamma is from bg or a good event
  nentries = (Int_t)t->GetEntries();
  cout << nentries << endl;
 
  t->GetEntry(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunSimulation()  {
  G4cout<<"Running the simulation"<<G4endl;
  SetDetDefaultValues();
  
  for(int i=0;i<nentries;i++)  {
    //if(i%1000==0) cout << i << "/" << nentries << "   DONE!" << endl;
    if(i%1000==0){
      std::cout << "Event: " << i <<", "<< (100.*i/nentries) <<"% of events done\r" <<std::flush;
    }

    //Only needed to get the position (and maybe tof-deltaE-E which is not included at the moment)
    if(en_det_after_target_opt==1)  {
      sprintf(tempparticle,"/gun/ion %i %i %i",(int)z_f,(int)mass_f,(int)charge_f);
      UI->ApplyCommand("/gun/particle ion");
      UI->ApplyCommand(tempparticle);

      kin->GetParticleGun()->SetParticleEnergy(energy_vertex_stored*MeV);
      kin->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(pva[0]*m, pva[1]*m, pva[2]*m));
      kin->GetParticleGun()->SetParticlePosition(G4ThreeVector(p0[0]*cm, p0[1]*cm, p0[2]*cm));
      
      runManager->BeamOn(1);
    }  
    
    if(steppingaction->GetHIEnDetFlag()==1 || en_det_after_target_opt==0)  {
      // Making the analysis only if position information (or tof-deltaE-E) 
      //cout<<"PosDetFlag true"<<endl;
      //dummy     = pos_det_after_target_res;
      // reconstsruction of the pos det point
      //pos_det_x = G4RandGauss::shoot(steppingaction->GetPosDetX()/cm,dummy);
      //pos_det_y = G4RandGauss::shoot(steppingaction->GetPosDetY()/cm,dummy);
      //pos_det_z = 140.0;   //Has to be set according to the distance in cm downstream of the target
      // Also need to get some points in X,Y if the pos det option is turned off

      pos_det[0] = G4RandGauss::shoot((p0[0]+pos_det[2]*pva[0]),pos_det_after_target_res);
      pos_det[1] = G4RandGauss::shoot((p0[1]+pos_det[2]*pva[1]),pos_det_after_target_res);
      
      if(en_det_after_target_opt==1)  {
        //energy, time resolution of tof-deltaE-E can be placed here  
        //dummy       = steppingaction->GetDetXYZDeltaE()/MeV*xyz_res*0.01/2.35;
        //xyz_delta_e = G4RandGauss::shoot(steppingaction->GetDetXYZDeltaE()/MeV,dummy);
        //dummy       = steppingaction->GetDetXYZE()/MeV*xyz_res*0.01/2.35;
        //xyz_e      = G4RandGauss::shoot(steppingaction->GetDETXYZE()/MeV,dummy);
      }
      else  {
        //xyz_e=-999.0;
        //xyz_delta_e =-999.0;
      }

      ver_rec[0] = G4RandGauss::shoot(p0[0],pos_det_at_target_res);             // Reconstruction of target point
      ver_rec[1] = G4RandGauss::shoot(p0[1],pos_det_at_target_res);
      ver_rec[2] = 0.0;

      dummy     = beta_res*beta_b;                                              // Reconstruction of beta
      beta_rec  = G4RandGauss::shoot(beta_b,dummy);
			      
      kin->GetParticleGun()->SetParticleDefinition(G4Gamma::GammaDefinition()); // Changing to the gamma particle:
      kin->GetParticleGun()->SetParticleEnergy(e_doppler*keV);
      kin->GetParticleGun()->SetParticleMomentumDirection(G4ThreeVector(gv[0]*m, gv[1]*m, gv[2]*m));
      kin->GetParticleGun()->SetParticlePosition(G4ThreeVector(g0[0]*cm,g0[1]*cm,g0[2]*cm));
      steppingaction       ->StartBeam();
      steppingaction       ->SetGammaShot(1);
      runManager           ->BeamOn(1);
      steppingaction       ->SetGammaShot(0);
   	      
      //---------------------------------------------------------------------------------------
      if(dali2_opt==1)  {                                                       // Checking first the DALI2 array
        det->fDali2Array->DetermineCrystalMult();
        if(det->fDali2Array->GetCrystalMult() > 0)  {
          for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++)  {
	    if(det->fDali2Array->GetCrystalFlag(nnn)==true && det->fDali2Array->GetCrystalEnergy(nnn) > 0.0)  {
              gamma_det_type            = 1;
              dali2_energy_not_cor[nnn] = det->fDali2Array->GetCrystalMeasuredEnergy(nnn);
              dali2_time[nnn]           = det->fDali2Array->GetCrystalMeasuredTime(nnn);
              dali2_flag[nnn]           = 1;
              if(dali2_fi_opt==1)  {
                dali2_fi[0][nnn] = det->fDali2Array->GetFIX(nnn);                //cout<<"Dali2_fi[0] = "<<dali2_fi[0][nnn]<<endl;
                dali2_fi[1][nnn] = det->fDali2Array->GetFIY(nnn);                //cout<<"Dali2_fi[1] = "<<dali2_fi[1][nnn]<<endl;
                dali2_fi[2][nnn] = det->fDali2Array->GetFIZ(nnn);                //cout<<"Dali2_fi[2] = "<<dali2_fi[2][nnn]<<endl;
              }
	    }
	  }
          t2->Fill();
        }
      }
      //---------------------------------------------------------------------------------------
      if(Shogun_opt==1)  {                                                       // Checking the Shogun array
        det->fShogunArray->DetermineCrystalMult();
        if(det->fShogunArray->GetCrystalMult() > 0) {
          for(int nnn=0;nnn<NUMBEROFSHOGUNCRYSTALS;nnn++) {
	    if(det->fShogunArray->GetCrystalFlag(nnn)==true && det->fShogunArray->GetCrystalEnergy(nnn) > 0.0){
              gamma_det_type            = 4;
              Shogun_energy_not_cor[nnn] = det->fShogunArray->GetCrystalMeasuredEnergy(nnn);
              Shogun_time[nnn]           = det->fShogunArray->GetCrystalMeasuredTime(nnn);
              Shogun_flag[nnn]           = 1;
              if(Shogun_fi_opt==1)  {
                Shogun_fi[0][nnn] = det->fShogunArray->GetFIX(nnn); 
                Shogun_fi[1][nnn] = det->fShogunArray->GetFIY(nnn); 
                Shogun_fi[2][nnn] = det->fShogunArray->GetFIZ(nnn); 
              }
	    }
	  }
          t2->Fill();
        }
      }
      //---------------------------------------------------------------------------------------
      if(LaBr3Opt==1)  {                                                       // Checking the array of large size LaBr3 detectors
        det->fLaBr3Array->DetermineCrystalMult();
        if(det->fLaBr3Array->GetCrystalMult() > 0){
          for(int nnn=0;nnn<NUMBEROFLABR3ARRAYCRYSTALS;nnn++){
	    if(det->fLaBr3Array->GetCrystalFlag(nnn)==true && det->fLaBr3Array->GetCrystalEnergy(nnn) > 0.0){
              gamma_det_type             = 6;
              LaBr3EnergyNotCor[nnn]     = det->fLaBr3Array->GetCrystalMeasuredEnergy(nnn);
              LaBr3Time[nnn]             = det->fLaBr3Array->GetCrystalMeasuredTime(nnn);
              LaBr3Flag[nnn]             = 1;
	    }
	  }
          t2->Fill();
        }
      }

      //------------------------------------------------------------------------------------------
      //reminder of the variables [xx][xx][0] is for the entire crystal:
      //float grape_flag[18], grape_crystal_flag[18][2][10]; 
      //float grape_energy_not_cor[18][2][10];
      if(grape_opt==1)  {                                                       // Grape start
        det->fGrapeArray->DetermineCrystalMult();
        if(det->fGrapeArray->GetCrystalMult()>0)  {
          for(int nnn=0;nnn<NUMBEROFGRAPEDETECTORS;nnn++)  {                                        // Detectors 
            for(int ppp=0;ppp<2;ppp++)  {                                       // Crystal by crystal (two per detector) 
              if(det->fGrapeArray->GetCrystalFlag(nnn,ppp,0)==true && det->fGrapeArray->GetCrystalEnergy(nnn,ppp,0) > 0.0)  {
                gamma_det_type                    = 2;
                grape_flag[nnn]                   = 1;
                grape_crystal_flag[nnn][ppp][0]   = 1;
                grape_energy_not_cor[nnn][ppp][0] = det->fGrapeArray->GetCrystalMeasuredEnergy(nnn,ppp,0);
                grape_time[nnn][ppp][0]           = det->fGrapeArray->GetCrystalMeasuredTime(nnn,ppp,0);
                for(int qqq=1;qqq<10;qqq++)  {                                  // Going to the segments 
                  if(det->fGrapeArray->GetCrystalFlag(nnn,ppp,qqq)==1 && det->fGrapeArray->GetCrystalEnergy(nnn,ppp,qqq) > 0.0)  {
                    grape_crystal_flag[nnn][ppp][qqq]   = 1;
                    grape_energy_not_cor[nnn][ppp][qqq] = det->fGrapeArray->GetCrystalMeasuredEnergy(nnn,ppp,qqq);
                    grape_time[nnn][ppp][qqq]           = det->fGrapeArray->GetCrystalMeasuredTime(nnn,ppp,qqq);
                  }
                }
              }
            } 
          }  
          t2->Fill();
        }
      }
      //--------------------------------------------------------------------------------------
      //Reminder of the variables:
      //float sgt_flag[1], sgt_coax_crystal_flag[1], sgt_planar_strip_flag[1][26]; //fictive central contactplus the 25 strips
      //float sgt_energy_not_cor[1][27]; //first the coax, then the planar, then the strips
      if(sgt_opt==1) {                                                          // Starting now with SGT
        det->fSGTArray->DetermineDetectorMult();
        if(det->fSGTArray->GetDetectorMult()>0)  {
          for(int nnn=0;nnn<NUMBEROFSGTDETECTORS;nnn++) {                                          // Only four detectors 
            if(det->fSGTArray->GetDetectorFlag(nnn)==true && 
               det->fSGTArray->GetDetectorEnergy(nnn)>0.0)  {                   // Coax
              gamma_det_type = 3;
              sgt_flag[nnn]  = 1;
              if(det->fSGTArray->GetCoaxFlag(nnn)==true && det->fSGTArray->GetCoaxEnergy(nnn)>0.0)  {
                sgt_coax_crystal_flag[nnn] = 1;
	        sgt_energy_not_cor[nnn][0] = det->fSGTArray->GetCoaxMeasuredEnergy(nnn);
                sgt_time[nnn][0]           = det->fSGTArray->GetCoaxMeasuredTime(nnn);
              }
              for(int ppp=0;ppp<26;ppp++)  {
                if(det->fSGTArray->GetPlanarFlag(nnn,ppp)==true && 
                   det->fSGTArray->GetPlanarEnergy(nnn,ppp)>0.0)  {             // Planar
                
                  sgt_planar_strip_flag[nnn][ppp] = 1;
	          sgt_energy_not_cor[nnn][ppp+1]  = det->fSGTArray->GetPlanarMeasuredEnergy(nnn,ppp);
                  sgt_time[nnn][ppp+1]            = det->fSGTArray->GetPlanarMeasuredTime(nnn,ppp);
                }
              }
            }
          }
          t2->Fill();
        }
      }
      //--------------------------------------------------------------------------------------------
      if(sphere_opt==1) {                                                       // Sphere
        if(det->fSphere->GetFlagFirstInteraction()==true)  {
          gamma_det_type        = 5;

          //This doesn't work. I don't know why. Better use the gamma vector from the generator.
          sphere_fi[0]         = det->fSphere->GetFirstInteractionX();
          sphere_fi[1]         = det->fSphere->GetFirstInteractionY();
          sphere_fi[2]         = det->fSphere->GetFirstInteractionZ();
          sphere_energy_not_cor = det->fSphere->GetEnergy();
          t2->Fill();
        }
      }
    }
    if(i<nentries-1) {
      SetDetDefaultValues();
      t->GetEntry(i+1);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SetDetDefaultValues() {
  int i,j,k;
  if(dali2_opt==1)
    for(i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
        dali2_flag[i]           = 0;
        dali2_energy_not_cor[i] = 0.;
        dali2_time[i]           = 0.;
        if(dali2_fi_opt==1)  {
          for(j=0;j<3;j++)dali2_fi[j][i] = -999.0;
        }
      }
  if(Shogun_opt==1)
    for(i=0;i<NUMBEROFSHOGUNCRYSTALS;i++)  {
      Shogun_flag[i]           = 0;
      Shogun_energy_not_cor[i] = 0.;
      Shogun_time[i]           = 0.;
    }
  if(grape_opt==1)
    for(i=0;i<NUMBEROFGRAPEDETECTORS;i++)  {
      grape_flag[i] = 0;
      for(j=0;j<2;j++)  {
        for(k=0;k<10;k++)  {
          grape_crystal_flag[i][j][k]   = 0;
          grape_energy_not_cor[i][j][k] = 0.;
          grape_time[i][j][k]           = 0.;
        }
      }
    } 
  if(sgt_opt==1)
    for(i=0;i<NUMBEROFSGTDETECTORS;i++)  {
      sgt_flag[i]=sgt_coax_crystal_flag[i]   = 0;
      for(j=0;j<27;j++)  {
        if(j<26) sgt_planar_strip_flag[i][j] = 0;                                 // Fictive central contactplus the 25 strips
        sgt_energy_not_cor[i][j]             = 0.;                                // First the coax, then the planar, then the strips
        sgt_time[i][j]                       = 0.;
      }
    }
  if(sphere_opt==1)  {                                                          // For the Sphere
    sphere_fi[0]          = 0.;
    sphere_fi[1]          = 0.;
    sphere_fi[2]          = 0.;
    sphere_energy_not_cor = 0.;
  }
  if(LaBr3Opt==1)
    for(i=0;i<NUMBEROFLABR3ARRAYCRYSTALS;i++)  {
      LaBr3Flag[i]         = 0;
      LaBr3EnergyNotCor[i] = 0.;
      LaBr3Time[i]         = 0.;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FillRootTreeHeader()  {
  // Filling the header Tree 
  if(dali2_opt==1)  {                                                           // DALI2
    for(int nnn=0;nnn<NUMBEROFDALI2CRYSTALS;nnn++)  {
      dali2_pos[0][nnn] = det->fDali2Array->GetPosXCrystal(nnn)/cm;
      dali2_pos[1][nnn] = det->fDali2Array->GetPosYCrystal(nnn)/cm;
      dali2_pos[2][nnn] = det->fDali2Array->GetPosZCrystal(nnn)/cm;
    }
  }
  if(Shogun_opt==1) {                                                            // Shogun
    for(int nnn=0;nnn<NUMBEROFSHOGUNCRYSTALS;nnn++)  {
      Shogun_pos[0][nnn] = det->fShogunArray->GetPosXCrystal(nnn)/cm;
      Shogun_pos[1][nnn] = det->fShogunArray->GetPosYCrystal(nnn)/cm;
      Shogun_pos[2][nnn] = det->fShogunArray->GetPosZCrystal(nnn)/cm;
    }
  }
  if(grape_opt==1)  {                                                           // Grape
    for(int nnn=0;nnn<NUMBEROFGRAPEDETECTORS;nnn++)  {
      for(int ppp=0;ppp<2;ppp++)  {
        for(int qqq=0;qqq<10;qqq++)  {
          grape_pos[0][nnn][ppp][qqq] = det->fGrapeArray->GetPosXCrystal(nnn,ppp,qqq)/cm; 
          grape_pos[1][nnn][ppp][qqq] = det->fGrapeArray->GetPosYCrystal(nnn,ppp,qqq)/cm;
          grape_pos[2][nnn][ppp][qqq] = det->fGrapeArray->GetPosZCrystal(nnn,ppp,qqq)/cm;
        }
      }
    }
  }
  if(sgt_opt==1) {                                                               // SGT
    for(int nnn=0;nnn<NUMBEROFSGTDETECTORS;nnn++)  {
      for(int ppp=0;ppp<30;ppp++)  {
        sgt_pos[0][nnn][ppp] = det->fSGTArray->GetPosXCrystal(nnn,ppp)/cm;
        sgt_pos[1][nnn][ppp] = det->fSGTArray->GetPosYCrystal(nnn,ppp)/cm;
        sgt_pos[2][nnn][ppp] = det->fSGTArray->GetPosZCrystal(nnn,ppp)/cm;
      }
    }
  }
  if(LaBr3Opt==1) {   
    for(int nnn=0;nnn<NUMBEROFLABR3ARRAYCRYSTALS;nnn++)  {
      LaBr3Pos[0][nnn] = det->fLaBr3Array->GetPosXCrystal(nnn)/cm;
      LaBr3Pos[1][nnn] = det->fLaBr3Array->GetPosYCrystal(nnn)/cm;
      LaBr3Pos[2][nnn] = det->fLaBr3Array->GetPosZCrystal(nnn)/cm;
    }
  }
  t2Header->Fill();
} 
