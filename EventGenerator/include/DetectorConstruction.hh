#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "MaterialList.hh"
//20161026 Pieter
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


class G4LogicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class DetectorConstruction : public G4VUserDetectorConstruction  {
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
     
     G4VPhysicalVolume* Construct();
     void               UpdateGeometry();
     
  public:  
           
     //Information about the world         
     G4double     GetWorldSizeX()     {return fWorldSizeX;};
     G4double     GetWorldSizeY()     {return fWorldSizeY;};
     G4double     GetWorldSizeZ()     {return fWorldSizeZ;};
     G4Material*  GetWorldMaterial()  {return fWorldMaterial;}; 
     //Information about the target
     G4double     GetTargetSizeX()    {return fTargetSizeX;};
     G4double     GetTargetSizeY()    {return fTargetSizeY;};
     G4double     GetTargetSizeZ()    {return fTargetSizeZ;};
     G4Material*  GetTargetMaterial() {return fTargetMaterial;}; 
     inline void SetTargetSize(double X,double Y, double Z)    
     {fTargetSizeX=X;fTargetSizeY=Y;fTargetSizeZ=Z;return;};
     void SetTargetMaterial (G4String);
                 
  private:
     MaterialList*       fMaterialList;
     //The world
     G4Box*              sWorld;
     G4LogicalVolume*    lWorld;
     G4VPhysicalVolume*  pWorld;
     G4double            fWorldSizeX;
     G4double            fWorldSizeY;
     G4double            fWorldSizeZ;
     G4Material*         fWorldMaterial;
     //The target
     G4double            fTargetSizeX;
     G4double            fTargetSizeY;
     G4double            fTargetSizeZ;
     G4Material*         fTargetMaterial;
     G4Box*              sTarget; 
     G4LogicalVolume*    lTarget; 
     G4VPhysicalVolume*  pTarget; 

     DetectorMessenger* detectorMessenger;

  private:
     G4VPhysicalVolume* ConstructVolumes();     
};
#endif

