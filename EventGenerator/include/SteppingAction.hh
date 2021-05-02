#ifndef SteppingAction_h
#define SteppingAction_h 1

//20161026 Pieter
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
class DetectorConstruction;
class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(DetectorConstruction*,RunAction*);
   ~SteppingAction();

    void UserSteppingAction(const G4Step*);
    void StartBeam();

    void SetTime(G4double input){time = input;};
    void SetDecayTime(G4double input){decayTime=input;};
    void SetFlagDecayWithinTarget(bool input){FlagDecayWithinTarget=input;};
    //The interaction point, e.g. Fragmentation
    void SetInteractionPoint(G4double X,G4double Y,G4double Z) {interactionPointX=X*mm;interactionPointY=Y*mm;interactionPointZ=Z*mm;} 
    void SetAfterGammaExcitation(bool State){FlagAfterGammaExcitation = State;}

    void InitialParticlePosition(G4ThreeVector input){initialPosition = input;};
    bool GetFlagDecayWithinTarget(){return FlagDecayWithinTarget;};
    G4int GetFlagTarget(){return FlagTarget;};
    G4double GetTime(){return time;}
    G4double GetTotalEnergy(){return totalEnergy;}; 
    G4double GetTotalEnergyAtPointOfGammaDecay(){return totalEnergyAtPointOfGammaDecay;};
    G4double GetTotalEnergyAfterTarget(){return totalEnergyAfterTarget;}; 
    G4ThreeVector GetPointOfGammaDecay(){ return pointOfGammaDecay;};
    G4ThreeVector GetMovedTrackInTarget(){ return movedTrackInTarget;};
 
  private:
    DetectorConstruction* detector;
    RunAction*            runAction;

    bool      FlagAfterGammaExcitation;
    bool      FlagDecayWithinTarget;
    int       FlagDecayPoint;
    G4int     FlagTarget;
    G4double  time;
    G4double  decayTime;
    G4double  totalEnergy;
    G4double  totalEnergyAtPointOfGammaDecay;
    G4double  totalEnergyAfterTarget;
    //G4double  energyLoss;
    G4double  interactionPointX,interactionPointY,interactionPointZ;
    G4double  targetEndZ;
    G4double  gammaDecayPointX,gammaDecayPointY,gammaDecayPointZ;
    G4ThreeVector pointOfGammaDecay;
    G4ThreeVector initialPosition;
    G4ThreeVector movedTrackInTarget;
};
#endif
