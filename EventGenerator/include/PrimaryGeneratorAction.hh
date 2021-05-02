#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);    
   ~PrimaryGeneratorAction();

  public:  
    void SetRndmBeam(G4double val)  {rndmBeam = val;}   
    void GeneratePrimaries(G4Event*);
    
    void   ResetEbeamCumul() {EbeamCumul = 0.;}
    G4double GetEbeamCumul() {return EbeamCumul;}
     
    G4ParticleGun* GetParticleGun() {return particleGun;}
    
  private:
    G4ParticleGun*             particleGun;
    DetectorConstruction*      detector;
    G4double                   rndmBeam;
    G4double                   EbeamCumul;       
    PrimaryGeneratorMessenger* gunMessenger;     
};
#endif


