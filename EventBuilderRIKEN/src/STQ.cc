#include "STQ.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//STQ::STQ(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld,float fSizeX,float fSizeY,float fSizeZ,float posX,float posY,float posZ,float fHoleDiameter) 
STQ::STQ(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld, float posX, float posY, float posZ, float fHole, float fDiameter,float fLength)
{
  cout<<"STQ Constructor"<<endl;

  fMaterialList        = ptMaterialList;
  lWorld               = ptWorld;

  G4RotationMatrix Rot3D;
  Rot3D.set(0, 0, 0);

  //sSTQBox = new G4Box("sSTQBox",fSizeX*cm/2.0,fSizeY*cm/2.0,fSizeZ*cm/2.0);
  sSTQ= new G4Tubs("sTQHole",fHole*cm/2,fDiameter*cm/2,fLength*cm/2,0*deg,360*deg);
 // sSTQ    = new G4SubtractionSolid("sSTQ",sSTQBox,sSTQHole);
  lSTQ    = new G4LogicalVolume(sSTQ,fMaterialList->GetMaterial("Fe"),"lSTQ",0,0,0);

  // Setting the vis attributes:
  G4VisAttributes* visAttSTQ = new G4VisAttributes(G4Colour(0.0,1.0,1.0)); 
  lSTQ->SetVisAttributes(visAttSTQ);

  G4ThreeVector Pos(posX*cm,posY*cm,posZ*cm); //Position of the center

  pSTQ = new G4PVPlacement(G4Transform3D(Rot3D,Pos),lSTQ,"pSTQ",lWorld,false,999998);
  cout<<"Exiting STQ Constructor"<<endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
STQ::~STQ()
{
}
