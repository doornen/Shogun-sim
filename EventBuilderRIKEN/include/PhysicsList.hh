#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class StepMax;
class PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList();
  virtual ~PhysicsList();

  void ConstructParticle();
    
  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
        
  void AddPhysicsList(const G4String& name);
  void ConstructProcess();
    
  void AddStepMax();       
  StepMax* GetStepMaxProcess() {return stepMaxProcess;};

private:
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4bool helIsRegisted;
  G4bool bicIsRegisted;
  G4bool biciIsRegisted;
    
  G4String                             emName;
  G4VPhysicsConstructor*               emPhysicsList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;
    
  StepMax* stepMaxProcess;
    
  PhysicsListMessenger* pMessenger;
};
#endif

