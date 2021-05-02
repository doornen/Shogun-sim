
#include "TrackingAction.hh"

#include "RunAction.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackingAction::TrackingAction(RunAction* run)
:runAction(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // extract Projected Range of primary particle
  if (aTrack->GetTrackID() == 1) {
    G4double x = aTrack->GetPosition().x() + runAction->GetOffsetX();
    if(x > runAction->GetLength()) x = runAction->GetLength(); 
    runAction->AddProjRange(x);
  }  
}
