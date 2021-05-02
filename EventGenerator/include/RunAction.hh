#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class PhysicsList;
class PrimaryGeneratorAction;
class G4Run;

namespace AIDA {
 class IAnalysisFactory;
 class ITree;
 class IHistogram1D;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class RunAction : public G4UserRunAction
{
public:
  RunAction(DetectorConstruction*, PhysicsList*,PrimaryGeneratorAction*);
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);
    
  void FillTallyEdep(G4int n, G4double e) {tallyEdep[n] += e;};
       
  G4double GetBinLength() {return binLength;};
  G4double GetLength()    {return length;};
  G4double GetOffsetX()   {return offsetX;} 
  void     FillHisto(G4int id, G4double x, G4double weight = 1.0);
    
  void AddProjRange (G4double x) {projRange += x; projRange2 += x*x;};
  void AddPrimaryStep() {nPrimarySteps++;};
                   
private:  
  void bookHisto();
  void cleanHisto();
    
private:
  DetectorConstruction*   detector;
  PhysicsList*            physics;
  PrimaryGeneratorAction* kinematic;
  G4double*               tallyEdep;   
  G4double                binLength;
  G4double                offsetX;
  G4double                length;
  G4double                projRange, projRange2;
  G4int                   nPrimarySteps;
           
  AIDA::IAnalysisFactory* af;  
  AIDA::ITree*            tree;
  AIDA::IHistogram1D*     histo[1];        
};
#endif

