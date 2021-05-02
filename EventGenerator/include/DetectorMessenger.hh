#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class DetectorMessenger: public G4UImessenger
{
  public:
  
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
    DetectorConstruction*   Detector;
    
    G4UIdirectory*             testemDir;
    G4UIdirectory*             detDir;    
    G4UIcmdWithAString*        MaterCmd;
    G4UIcmdWithADoubleAndUnit* SizeXCmd;
    G4UIcmdWithADoubleAndUnit* SizeYZCmd;    
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    
    G4UIcmdWithAString*        TalMateCmd;
    G4UIcmdWith3VectorAndUnit* TalSizeCmd;
    G4UIcmdWith3VectorAndUnit* TalPosiCmd;
            
    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif

