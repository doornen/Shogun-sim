#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SteppingAction : public G4UserSteppingAction  {
 public:

  SteppingAction(DetectorConstruction*,RunAction*);
  ~SteppingAction();

  int      GetHIEnDetFlag() {return fHIEnDetFlag;};
  int      GetPosDetFlag() {return fPosDetFlag;};
  G4double GetPosDetX() {return fPosDetX;};
  G4double GetPosDetY() {return fPosDetY;};
  G4double GetPosDetZ() {return fPosDetZ;};

  void     StartBeam();
  void     SetGammaShot(int dummy){fGammaShot = dummy;};
  void     StopBeam();
  void     UserSteppingAction(const G4Step*);
    
 private:
  
  int                   fBeamPipeFlag;
  int                   fBeamPipeHit;
  DetectorConstruction* fDetector;
  int                   fGammaShot;
  int                   fHIEnDetFlag;
  int                   fPosDetFlag; 
  int                   fTargetHolderFlag;
  int                   fTargetHolderHit;
  G4double              fPosDetX;
  G4double              fPosDetY;
  G4double              fPosDetZ;
  RunAction*            fRunAction;
};
#endif
