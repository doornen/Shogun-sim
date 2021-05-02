#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"

#include "G4UnitsTable.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include <iostream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()  {
  // default parameter values
  //worldsize
  fWorldSizeX = 4.0;
  fWorldSizeY = 4.0;
  fWorldSizeZ = 4.0;
  //targetsize
  fTargetSizeX = fTargetSizeY = 10.;
  fTargetSizeZ = .4;
  
  fSTQInclude           = 0;
  fTargetHolderInclude  = 0;
  fBeamPipeInclude      = 0;
  // create commands for interactive definition of the detector  
  // Not used.
  fDetectorMessenger = new DetectorMessenger(this);
  //The materials
  fMaterialList      = new MaterialList;
  fTargetMaterial    = fMaterialList->GetMaterial("Be");
  fWorldMaterial     = fMaterialList->GetMaterial("Galactic");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()  {
  delete fDetectorMessenger;
  delete fMaterialList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()  {
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()  {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  //---------------------------------------------------------------
  // World---------------------------------------------------------
  //---------------------------------------------------------------
  sWorld = new G4Box("sWorld",fWorldSizeX*m/2,fWorldSizeY*m/2,fWorldSizeZ*m/2);   			                      
  lWorld = new G4LogicalVolume(sWorld,fWorldMaterial,"lWorld",0,0,0);                          
  pWorld = new G4PVPlacement(0,G4ThreeVector(),"pWorld",lWorld,0,false,0);				
  //--------------------------------------------------------------- 
  // Target--------------------------------------------------------
  //--------------------------------------------------------------- 
  // The target is only created if there is no target holding structure!!!
  //sTarget = new G4Box("Target",fTargetSizeX*cm/2,fTargetSizeY*cm/2,fTargetSizeZ*cm/2);
  sTarget = new G4Tubs("Target",0.0*mm/2,fTargetSizeX*cm/2,fTargetSizeZ*cm/2,0*deg,360*deg);
  cout<<"fTargetSizeX = "<<fTargetSizeX<<"; fTargetSizeY = " <<fTargetSizeY<<"; fTargetSizeZ = "<<fTargetSizeZ<<endl;
  lTarget = new G4LogicalVolume(sTarget, fTargetMaterial, "Target", 0, 0, 0);
  if(fTargetHolderInclude==0)pTarget = new G4PVPlacement(0,G4ThreeVector(0.0*m,0.0*m,0.0*m),lTarget,"Target",lWorld,false,1111);

  //---------------------------------------------------------------
  // The Dali2 array
  if(fDali2Include==1) fDali2Array = new Dali2(fMaterialList,lWorld);
  // The Shogun array
  if(fShogunInclude==1) fShogunArray = new Shogun(fMaterialList,lWorld);
  // The Grape array
  if(fGrapeInclude==1) fGrapeArray = new Grape(fMaterialList,lWorld);
  // The SGT array
  if(fSGTInclude==1)   fSGTArray   = new SGT(fMaterialList,lWorld);
  // The sphere
  if(fSphereInclude==1)fSphere     = new Sphere(fMaterialList,lWorld);
  // The LaBr3Array
  if(fLaBr3Include==1)fLaBr3Array     = new LaBr3Array(fMaterialList,lWorld);
  //---------------------------------------------------------------

  //---------------------------------------------------------------
  // The STQ
  if(fSTQInclude==1)          fSTQ          = new STQ(fMaterialList,lWorld,0.,0.,202.3,28.,130.8,264.6);
  if(fTargetHolderInclude==1) fTargetHolder = new TargetHolder(fMaterialList,lWorld);
  if(fBeamPipeInclude==1)     fBeamPipe     = new BeamPipe(fMaterialList,lWorld,bpinsideRadius,bpoutsideRadius,false);

  //---------------------------------------------------------------
  // The Collimator
  if(fCollimatorInclude==1)  fCollimator    = new Collimator(fMaterialList,lWorld);


  G4VisAttributes* visAttTarget= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  G4VisAttributes* visAttBox= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  
  lWorld->SetVisAttributes(G4VisAttributes::Invisible); 
  //The target
  lTarget->SetVisAttributes(visAttTarget);
  
  cout<<"Created all detectors"<<endl;
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTargetMaterial(G4String materialName)  {
  // search the material by its name 
  G4cout << "Material set is " << materialName << G4endl;
  fTargetMaterial = fMaterialList->GetMaterial(materialName);
  lTarget->SetMaterial(fTargetMaterial);
  
    G4cout << "\n----> The target is " << fTargetSizeZ << " cm of "
           << materialName << G4endl;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4RunManager.hh" 
 
void DetectorConstruction::UpdateGeometry()  {
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}
