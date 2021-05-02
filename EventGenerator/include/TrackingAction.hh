#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction {

  public:  
    TrackingAction(RunAction*);
   ~TrackingAction() {};
   
    void PostUserTrackingAction(const G4Track*);
    
  private:
    RunAction* runAction;    
};
#endif
