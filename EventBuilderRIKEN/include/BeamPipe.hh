#ifndef BeamPipe_h
#define BeamPipe_h 1

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
class BeamPipe  {
 public:
  BeamPipe(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld, double bpinsideRadius, double OutsideRadius, bool oldPipe);
  ~BeamPipe();
  
  void ConstructShield();
  void SetPbShieldThickness(float a){PbShieldThickness = a;}
  void SetSnShieldThickness(float a){SnShieldThickness = a;}
  void SetShieldRadius(float a){ShieldRadius = a;}

  
 private:

  G4LogicalVolume*       lWorld;                  //need a pointer to the logic world volume
  MaterialList*          fMaterialList;           //Pointer to the material list
  
  G4Tubs*                sBeamPipeTubs[4]; 
  G4LogicalVolume*       lBeamPipeTubs[4]; 
  G4VPhysicalVolume*     pBeamPipeTubs[4];

  G4Tubs*                sBeamPipeOut[2];
  G4Tubs*                sBeamPipeIn[2];

  G4UnionSolid*          sBeamPipeUnion;
  G4SubtractionSolid*    sBeamPipeSubtraction[2];

  G4LogicalVolume*       lBeamPipe; 
  G4VPhysicalVolume*     pBeamPipe;

  G4Tubs*                sPbShield;
  G4LogicalVolume*       lPbShield; 
  G4VPhysicalVolume*     pPbShield;

  G4Tubs*                sSnShield;
  G4LogicalVolume*       lSnShield; 
  G4VPhysicalVolume*     pSnShield;

  float PbShieldThickness;
  float SnShieldThickness;
  float ShieldRadius;

  double bpinsideRadius;
  double bpoutsideRadius;
};

#endif
