#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "Dali2.hh"
#include "Shogun.hh"
#include "Grape.hh"
#include "SGT.hh"
#include "Sphere.hh"
#include "LaBr3Array.hh"
#include "STQ.hh"
#include "TargetHolder.hh"
#include "BeamPipe.hh"
#include "MaterialList.hh"
#include "Collimator.hh"

class G4LogicalVolume;
class G4Material;
class G4Box;
class G4Tubs;

class G4VPhysicalVolume;
class DetectorMessenger;

class Dali2;
class Shogun;
class Grape;
class LaBr3Array;
class Sphere;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class DetectorConstruction : public G4VUserDetectorConstruction  {
public:
  
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume* ConstructVolumes(); 

  int         GetDali2Include()     {return fDali2Include;}
  int         GetShogunInclude()    {return fShogunInclude;}
  int         GetGrapeInclude()     {return fGrapeInclude;}
  int         GetSGTInclude()       {return fSGTInclude;}
  int         GetSphereInclude()    {return fSphereInclude;}
  int         GetLaBr3Include()     {return fLaBr3Include;}
  int         GetSTQInclude()       {return fSTQInclude;}
  int         GetCollimatorInclude(){return fCollimatorInclude;}
  G4Material* GetTargetMaterial()   {return fTargetMaterial;} 
  G4double    GetTargetSizeX()      {return fTargetSizeX;}
  G4double    GetTargetSizeY()      {return fTargetSizeY;}
  G4double    GetTargetSizeZ()      {return fTargetSizeZ;}
  G4Material* GetWorldMaterial()    {return fWorldMaterial;} 
  G4double    GetWorldSizeX()       {return fWorldSizeX;}
  G4double    GetWorldSizeY()       {return fWorldSizeY;}
  G4double    GetWorldSizeZ()       {return fWorldSizeZ;}

  void        SetBeamPipeInclude(int a, double b, double c){fBeamPipeInclude = a; bpinsideRadius = b;bpoutsideRadius = c;}
  void        SetDetectorInclude(int a, int b, int c, int d, int e,int f) {
    fDali2Include=a;fShogunInclude=b;fGrapeInclude=c;fSGTInclude=d; fSphereInclude=e;fLaBr3Include=f;}
  void        SetSTQInclude(int a){fSTQInclude = a;}
  void        SetTargetHolderInclude(int a){fTargetHolderInclude = a;}
  void        SetTargetMaterial (G4String);
  void        SetTargetSize(double X,double Y, double Z)  {
    fTargetSizeX=X; fTargetSizeY=Y; fTargetSizeZ=Z;}
  void        SetCollimatorInclude(int a){fCollimatorInclude = a;}
 
  void        UpdateGeometry();
    
  Dali2*               fDali2Array;
  Shogun*              fShogunArray;
  Grape*               fGrapeArray;
  SGT*                 fSGTArray;
  Sphere*              fSphere;
  LaBr3Array*          fLaBr3Array;
  STQ*                 fSTQ;
  TargetHolder*        fTargetHolder;
  BeamPipe*            fBeamPipe;
  Collimator*          fCollimator;

private:

  G4Box*               sWorld;
  G4LogicalVolume*     lWorld;
  G4VPhysicalVolume*   pWorld;
  
  //G4Box*               sTarget; 
  G4Tubs*              sTarget;
  G4LogicalVolume*     lTarget; 
  G4VPhysicalVolume*   pTarget; 
    
  int                  fBeamPipeInclude;
  int                  fDali2Include;
  int                  fShogunInclude;
  DetectorMessenger*   fDetectorMessenger;
  int                  fGrapeInclude;
  MaterialList*        fMaterialList;
  int                  fSGTInclude;
  int                  fSphereInclude;
  int                  fLaBr3Include;
  int                  fSTQInclude;
  int                  fCollimatorInclude;
  int                  fTargetHolderInclude;
  G4Material*          fTargetMaterial;
  G4double             fTargetSizeX;
  G4double             fTargetSizeY;
  G4double             fTargetSizeZ;
  char                 fTempname[200];
  G4Material*          fWorldMaterial;
  G4double             fWorldSizeX;
  G4double             fWorldSizeY;
  G4double             fWorldSizeZ;
  G4double             bpinsideRadius;
  G4double             bpoutsideRadius;
};
#endif

