#ifndef TargetHolder_h
#define TargetHolder_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"

#include "G4UnitsTable.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "Randomize.hh"

#include "G4ios.hh"
#include <iostream>
#include "TMath.h"

#include  "MaterialList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class TargetHolder
{
 public:

  TargetHolder(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld);
  ~TargetHolder();
  void CreateTargetHolder();
  
 private:

  G4LogicalVolume*       lWorld;                  // Need a pointer to the logic world volume
  MaterialList*          fMaterialList;           // Pointer to the material list

  G4Tubs*                sTopFlange;              // The flange on to close the beam pipe
  G4LogicalVolume*       lTopFlange; 
  G4VPhysicalVolume*     pTopFlange;

  G4Box*                 sTargetBoxRaw;           // The plastic holding structure,
  G4Tubs*                sTargetHole;             // The six target holes
  G4SubtractionSolid*    sTargetBox;              // The target box without the target holes
  G4LogicalVolume*       lTargetBox; 
  G4VPhysicalVolume*     pTargetBox[6];           // The physical volume

  G4Tubs*                sTarget;                 // The target
  G4LogicalVolume*       lTarget[6]; 
  G4VPhysicalVolume*     pTarget[6];              // The six target physical volumes

  FILE*                  fFileIn;                 // The input file for the target thicknesses 
  float                  fTargetThickness[6];     // Thickness of the taret in g/cm^2
  char                   fTargetMaterial[6][20];  // The target material 
  char                   fTemp[100];      
  int                    fSetTarget;              // The set target position.        

};

#endif
