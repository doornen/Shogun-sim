#ifndef Collimator_h
#define Collimator_h 1

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
class Collimator
{
 public:

  Collimator(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld);
  ~Collimator();
  void CreateCollimator();
  
 private:

  G4LogicalVolume*       lWorld;                  // Need a pointer to the logic world volume
  MaterialList*          fMaterialList;           // Pointer to the material list

  G4Cons*                sCollHole;               // The hole through the collimator
  G4Box*                 sCollBlock;              // Form of the collimator: 1 -> box
  G4Tubs*                sCollTubs;               // 2 -> tubs
  G4SubtractionSolid*    sCollimator;             // The target box without the target holes
  G4LogicalVolume*       lCollimator; 
  G4VPhysicalVolume*     pCollimator;             // The physical volume

  bool                   fBoxShape;               // option if box or tubs shape 
  FILE*                  fFileIn;                 // The input file for the target thicknesses 
  float                  fCollHoleTarget;         // Thickness of the hole towards at the source 
  float                  fCollHoleDetector;       // Thickness of the hole towards at the detector 

  float                  fCollBlockSize[3];       // The size of the collimator block;
  float                  fCollBlockComposition[3];// The composition in percentage for W,Ni,Fe
  float                  fCollBlockDensity;       // The density of the collimator in g/cm^3
  float                  fCollimatorPos[3];       // The position of the collimator in cm.
};

#endif
