#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4SteppingManager.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(DetectorConstruction* det,RunAction* RuAct)
:detector(det),runAction(RuAct)
{ 
  time = 0.;
  decayTime = 0.;
  totalEnergy=0.; 
  totalEnergyAtPointOfGammaDecay=0.; 
  totalEnergyAfterTarget=0.; 
  //energyLoss=0.;
  interactionPointX=interactionPointY=interactionPointZ=9999.;
  FlagAfterGammaExcitation = FALSE;
  gammaDecayPointX=gammaDecayPointY=gammaDecayPointZ=9999.;
  FlagDecayPoint = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track * theTrack = aStep->GetTrack(); 
  //totalEnergy = theTrack->GetKineticEnergy();
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();

  //Getting the total energy at the interaction point
  if(theTrack->GetVolume()->GetCopyNo() ==1111  
     && theTrack->GetVolume()->GetName() == "Target"
     && particle->GetParticleType()=="nucleus")
   //G4cout<<"The energy loss is "<<aStep->GetTotalEnergyDeposit()/keV<<G4endl;
  {
    if(FlagAfterGammaExcitation==FALSE && theTrack->GetPosition().z()*mm <= interactionPointZ*mm+0.01)
    //The 0.05 had to be added due to the stepsize in the Physicslist
    {
      FlagTarget=1;
      totalEnergy = theTrack->GetKineticEnergy();
    }
  }
  // Getting the total energy at gamma ray decay point if its in the targert
  // Function has to be passed at least once
  if(particle->GetParticleType()=="nucleus" && FlagAfterGammaExcitation==TRUE && 
    (time <= decayTime || FlagDecayPoint<1) && theTrack->GetVolume()->GetCopyNo() ==1111)
  {
    FlagDecayPoint=1;
    time = time + aStep-> GetDeltaTime()*ns;
    if(time>decayTime) FlagDecayWithinTarget=TRUE;

    movedTrackInTarget =theTrack->GetPosition() - initialPosition;
    pointOfGammaDecay =  theTrack->GetPosition();
    totalEnergyAtPointOfGammaDecay = theTrack->GetKineticEnergy();

    //G4cout<<"The Z position is "<<theTrack->GetPosition().z()*mm<<G4endl;
    //G4cout<<"The time  is "<<time<<G4endl;
    //G4cout<<"The decay time  is "<<decayTime<<G4endl;
  }
  if(particle->GetParticleType()=="nucleus" && theTrack->GetVolume()->GetCopyNo() ==1111) 
  { 
    //The total energy after the target: The last value is the correct one
    totalEnergyAfterTarget = theTrack->GetKineticEnergy();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::StartBeam()
{
  FlagTarget = 0;
  FlagAfterGammaExcitation=FALSE;
  FlagDecayWithinTarget=FALSE;
  FlagDecayPoint = 0;
}
