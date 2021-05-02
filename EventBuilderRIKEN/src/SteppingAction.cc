#include "Globals.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(DetectorConstruction* det,RunAction* RuAct)
:fDetector(det),fRunAction(RuAct){   
  fTargetHolderHit = 0;
  fBeamPipeHit     = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction()  { 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* aStep){  
  G4double edep = aStep->GetTotalEnergyDeposit()/keV;
  if (edep <= 1.e-10) return;
 
  G4Track * theTrack = aStep->GetTrack(); 
  //float totalEnergy = theTrack->GetKineticEnergy();
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
  // Checking if the particle is a nucleus:
  if(particle->GetParticleType()=="nucleus"){
    // Hit on Pos Det----------------------------------------------
    //if(theTrack->GetVolume()->GetCopyNo() >= ???? && theTrack->GetVolume()->GetCopyNo() <= ?????   
    //&& theTrack->GetVolume()->GetName() == "????" && theTrack->GetTrackID()==1 && posDetFlag==0)
    //G4cout<<"Hit on Pos Det"<<G4endl;
    //posDetFlag++;
    //PosDetFlag[(int)theTrack->GetVolume()->GetCopyNo() - ????] = 1;
    //PosDetX = -1.0*theTrack->GetPosition().x();
    //PosDetY = theTrack->GetPosition().y();
    //PosDetZ = theTrack->GetPosition().z();
  }
  // Checking if the particle shot is a gamma:
  if(fGammaShot==1){
    int trackVolume = theTrack->GetVolume()->GetCopyNo();
    //___________________________________________
    // Dali2
    if(trackVolume >=100100 && trackVolume < (100000+NUMBEROFDALI2CRYSTALS*100+2)){
      if(trackVolume%100==0){
        int dummy = ((trackVolume-100000)/100)-1;
        if(fDetector->fDali2Array->GetFIPointOption()==1){
          //if( fDetector->fDali2Array->GetCrystalFlag(dummy)==false)
          //  fDetector->fDali2Array->SetFITime(dummy,10000);
          //G4cout<<"Got hit"<<G4endl;
          double newTime = theTrack->GetGlobalTime();
          double oldTime = fDetector->fDali2Array->GetFITime(dummy);
          if(newTime<oldTime ||  fDetector->fDali2Array->GetCrystalFlag(dummy)==false) {
            fDetector->fDali2Array->SetFIX(dummy,(float)theTrack->GetPosition().x()/cm); 
            fDetector->fDali2Array->SetFIY(dummy,(float)theTrack->GetPosition().y()/cm);
            fDetector->fDali2Array->SetFIZ(dummy,(float)theTrack->GetPosition().z()/cm);
            fDetector->fDali2Array->SetFITime(dummy,newTime);
          }
        } 
        fDetector->fDali2Array->SetCrystalFlagTrue(dummy);
        fDetector->fDali2Array->AddCrystalEnergy(dummy,edep);
      }
    }
    //___________________________________________
    // Shogun
    if(trackVolume >=400000 && trackVolume < 405000){
      int dummy = trackVolume-400000; 
      //Setting the time only if the flag was false!
      if(fDetector->fShogunArray->GetCrystalFlag(dummy)==false){
         fDetector->fShogunArray->SetCrystalTime(dummy,theTrack->GetGlobalTime());
      //G4cout <<"Time: "<<theTrack->GetGlobalTime()<<G4endl;
      }
      if(fDetector->fShogunArray->GetFIPointOption()==1)  {
        //if( fDetector->fDali2Array->GetCrystalFlag(dummy)==false)
        //  fDetector->fDali2Array->SetFITime(dummy,10000);
        //G4cout<<"Got hit"<<G4endl;
        double newTime = theTrack->GetGlobalTime();
        double oldTime = fDetector->fShogunArray->GetFITime(dummy);
        if(newTime<oldTime ||  fDetector->fShogunArray->GetCrystalFlag(dummy)==false) {
          fDetector->fShogunArray->SetFIX(dummy,(float)theTrack->GetPosition().x()/cm); 
          fDetector->fShogunArray->SetFIY(dummy,(float)theTrack->GetPosition().y()/cm);
          fDetector->fShogunArray->SetFIZ(dummy,(float)theTrack->GetPosition().z()/cm);
          fDetector->fShogunArray->SetFITime(dummy,newTime);
        }
      } 
      fDetector->fShogunArray->SetCrystalFlagTrue(dummy);
      fDetector->fShogunArray->AddCrystalEnergy(dummy,edep);
    }
    //___________________________________________
    // LaBr3Array
    if(trackVolume >=600000 && trackVolume < (600000+NUMBEROFLABR3ARRAYCRYSTALS))  {
      int dummy = trackVolume-600000; 
      //Setting the time only if the flag was false!
      if(fDetector->fLaBr3Array->GetCrystalFlag(dummy)==false){
         fDetector->fLaBr3Array->SetCrystalTime(dummy,theTrack->GetGlobalTime());
      //G4cout <<"Time: "<<theTrack->GetGlobalTime()<<G4endl;
      }
      fDetector->fLaBr3Array->SetCrystalFlagTrue(dummy);
      fDetector->fLaBr3Array->AddCrystalEnergy(dummy,edep);
    }
    //___________________________________________
    if(trackVolume >=200100 && trackVolume <= 202000)
    {
      //the copy number is 2xxxyz, with xxx=detector number(0-17),y=crystal(0-1),z=segment(0-8)
      // for all other material belonging to the detector, y>=20 
      // reminder of the variables:
      // grapeCrystalFlag[20][2][20];
      //  grapeCrystalEnergyGamma[20][2][20];
      if(trackVolume%100<20)
      {
        //G4cout <<"TrackVolume: "<<trackVolume<<G4endl;
        int seg = trackVolume%10;
        int cry = ((trackVolume-seg)/10)%10;
        int det = (trackVolume-10*cry-seg-200000)/100 -1;
        //G4cout <<"det: "<<det<<G4endl;
        //G4cout <<"cry: "<<cry<<G4endl;
        fDetector->fGrapeArray->SetCrystalFlagTrue(det,cry,0);  // The entire crystal
        fDetector->fGrapeArray->AddCrystalEnergy(det,cry,0,edep);

        fDetector->fGrapeArray->SetCrystalFlagTrue(det,cry,seg+1);  // The individual segments
        fDetector->fGrapeArray->AddCrystalEnergy(det,cry,seg+1,edep);
      }  
    }
    //___________________________________________
    if(trackVolume >=300100 && trackVolume < 300500)
    {
      //the copy number is 3xxxyy, with xxx=detector number(0-17),y=coax(0),strip(1-26)
      // for all other material belonging to the detector, y>=30 
      // reminder of the variables:
      // sgtDetectorFlag[2];  //maximum 2 detectors at the moment
      // sgtCoaxFlag[2];
      // sgtPlanarFlag[2][26];
      // sgtCombinedEnergyGamma[2];
      // sgtCoaxEnergyGamma[2];
      // sgtPlanarEnergyGamma[2][26];
      if(trackVolume%100<27)
      {
        if(trackVolume%100==0)
        {
          int det = (trackVolume/100)%100 -1;
          fDetector->fSGTArray->SetDetectorFlagTrue(det);
          fDetector->fSGTArray->SetCoaxFlagTrue(det);
          fDetector->fSGTArray->AddCoaxEnergy(det,edep);
          fDetector->fSGTArray->AddDetectorEnergy(det,edep);
        }
        else
        {
          int strip = (trackVolume%100);  //from 1 to 26
          int det = ((trackVolume-strip)/100)%100 -1;
          fDetector->fSGTArray->SetDetectorFlagTrue(det);
          fDetector->fSGTArray->SetPlanarFlagTrue(det,0); //the total planar detector flag
          fDetector->fSGTArray->AddPlanarEnergy(det,0,edep);
          fDetector->fSGTArray->SetPlanarFlagTrue(det,strip);  // the 25 strips
          fDetector->fSGTArray->AddPlanarEnergy(det,strip,edep);
          fDetector->fSGTArray->AddDetectorEnergy(det,edep);
        }
      } 
    }
    //___________________________________________
    //The sphere:
    if(trackVolume==500000)  {
      fDetector->fSphere->AddEnergy(edep);
      //G4cout<<"Edep: "<<edep<<G4endl;
      if(fDetector->fSphere->GetFlagFirstInteraction()==false)  {
        //G4cout<<"Got hit"<<G4endl;
        fDetector->fSphere->SetFirstInteractionX(theTrack->GetPosition().x()/cm); 
        fDetector->fSphere->SetFirstInteractionY(theTrack->GetPosition().y()/cm);
        fDetector->fSphere->SetFirstInteractionZ(theTrack->GetPosition().z()/cm);
        fDetector->fSphere->SetFlagFirstInteraction(true);
      }
    }
    //___________________________________________
    if(trackVolume == 999997)  {
      if(fTargetHolderFlag==0)  {
        fTargetHolderHit++;
        fTargetHolderFlag = 1;
        //G4cout<<"Hit on TargetHolder number: "<<fTargetHolderHit<<G4endl;
      }
    }
    //___________________________________________
    if(trackVolume == 999999)  {
      if(fBeamPipeFlag==0)  {
        fBeamPipeHit++;
        fBeamPipeFlag = 1;
        //G4cout<<"Hit on BeamPipe number: "<<fBeamPipeHit<<G4endl;
      }
    }
    //___________________________________________
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::StartBeam()  {
  fPosDetFlag       = 0;
  // Heavy ion Energy Detector after the target
  fHIEnDetFlag      = 0;
  fTargetHolderFlag = 0;
  fBeamPipeFlag     = 0;
  
  if(fDetector->GetDali2Include()==1) fDetector->fDali2Array->ResetValues();
  if(fDetector->GetShogunInclude()==1) fDetector->fShogunArray->ResetValues();
  if(fDetector->GetGrapeInclude()==1) fDetector->fGrapeArray->ResetValues();
  if(fDetector->GetSGTInclude()==1)   fDetector->fSGTArray->ResetValues();
  if(fDetector->GetSphereInclude()==1)fDetector->fSphere->ResetValues();
  if(fDetector->GetLaBr3Include()==1) fDetector->fLaBr3Array->ResetValues();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::StopBeam()  {
}
