#ifndef STQ_h
#define STQ_h 1

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
class STQ
{
  public:
  
   // STQ(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld,float fSizeX,float fSizeY,float fSizeZ,float posX,float posY,float posZ,float fHoleDiameter);
   STQ(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld, float posX, float posY, float posZ, float fHole, float fDiameter,float fLength);
   ~STQ();


  private:

     G4LogicalVolume*    lWorld;                   //need a pointer to the logic world volume
     MaterialList*       fMaterialList;            //Pointer to the material list

     //G4Box*              sSTQBox;       // Basic box for the STQ
     G4Tubs*             sSTQ;
     //G4SubtractionSolid* sSTQ;
     G4LogicalVolume*    lSTQ;
     G4VPhysicalVolume*  pSTQ;
};

#endif
