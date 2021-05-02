#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                                   PrimaryGeneratorAction* Gun)
:Action(Gun)
{ 
  gunDir = new G4UIdirectory("/testem/gun/");
  gunDir->SetGuidance("gun control");

  RndmCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/rndm",this);
  RndmCmd->SetGuidance("random lateral extension on the beam");
  RndmCmd->SetParameterName("rBeam",false);
  RndmCmd->SetRange("rBeam>=0.");
  RndmCmd->SetUnitCategory("Length");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete gunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{ 
  if (command == RndmCmd)
   {Action->SetRndmBeam(RndmCmd->GetNewDoubleValue(newValue));}   
}
