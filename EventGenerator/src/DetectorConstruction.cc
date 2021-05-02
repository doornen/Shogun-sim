#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
{
  // default parameter values
  //worldsize
  fWorldSizeX = 4.0;
  fWorldSizeY = 4.0;
  fWorldSizeZ = 4.0;
  //targetsize
  fTargetSizeX = fTargetSizeY = 7.;
  fTargetSizeZ = .4;

  // create commands for interactive definition of the detector  
  detectorMessenger = new DetectorMessenger(this);
  fMaterialList     = new MaterialList; 
  
  fTargetMaterial    = fMaterialList->GetMaterial("Fe");
  fWorldMaterial     = fMaterialList->GetMaterial("Galactic");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // World
  sWorld = new G4Box("sWorld",fWorldSizeX*m/2,fWorldSizeY*m/2,fWorldSizeZ*m/2); 
  lWorld = new G4LogicalVolume(sWorld,fWorldMaterial,"lWorld",0, 0, 0);
  pWorld = new G4PVPlacement(0,G4ThreeVector(),"pWorld",lWorld,0,false,0);
 
 // Defining the target:
 sTarget  = new G4Box("Target",fTargetSizeX*cm/2,fTargetSizeY*cm/2,fTargetSizeZ*cm/2);
 lTarget  = new G4LogicalVolume( sTarget, fTargetMaterial, "Target", 0, 0, 0);
 pTarget  = new G4PVPlacement(0,G4ThreeVector(0.0*m,0.0*m,0.0*m),lTarget,"Target",lWorld,false,1111);

  //
  //always return the World volume
  //  
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTargetMaterial(G4String materialName)
{
   // search the material by its name 
  G4cout << "Material set is " << materialName << G4endl;
  fTargetMaterial = fMaterialList->GetMaterial(materialName);
  lTarget->SetMaterial(fTargetMaterial);
  
    G4cout << "\n----> The target is " << fTargetSizeZ << " cm of "
           << materialName << G4endl;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4RunManager.hh" 
 
void DetectorConstruction::UpdateGeometry()
{
G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}
