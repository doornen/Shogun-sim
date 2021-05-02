#include "PhysListEmStandardSS.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4CoulombScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhysListEmStandardSS::PhysListEmStandardSS(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhysListEmStandardSS::~PhysListEmStandardSS()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PhysListEmStandardSS::ConstructProcess()
{
  // Add standard EM Processes

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddDiscreteProcess(new G4CoulombScattering);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
	    
    } else if (particleName == "e+") {
      //positron
      pmanager->AddDiscreteProcess(new G4CoulombScattering);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
            
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddDiscreteProcess(new G4CoulombScattering);
      pmanager->AddProcess(new G4MuIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3,3);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 4,4);
             
    } else if (particleName == "alpha" || particleName == "He3") {
      pmanager->AddDiscreteProcess(new G4CoulombScattering);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);

    } else if (particleName == "GenericIon" ) { 

      G4CoulombScattering* cs = new G4CoulombScattering();
      cs->SetBuildTableFlag(false);
      pmanager->AddDiscreteProcess(cs);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddDiscreteProcess(new G4CoulombScattering);
      pmanager->AddProcess(new G4hIonisation,        -1,2,2);
    }
  }
}
