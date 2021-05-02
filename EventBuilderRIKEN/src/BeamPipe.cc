#include "BeamPipe.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
BeamPipe::BeamPipe(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld, double inside, double outside, bool oldPipe)  {
  cout<<"BeamPipe Constructor"<<endl;

  fMaterialList        = ptMaterialList;
  lWorld               = ptWorld;

  PbShieldThickness = SnShieldThickness = 0.0;

  G4RotationMatrix Rot3D;
  Rot3D.set(0, 0, 0);
 
  sBeamPipeTubs[0]= new G4Tubs("sBeamPipeTubs[0]",144.0*mm/2,185.0*mm/2,20.0*mm/2,0*deg,360*deg);
  lBeamPipeTubs[0] = new G4LogicalVolume(sBeamPipeTubs[0],fMaterialList->GetMaterial("Al"),"lBeamPipeTubs[0]",0,0,0);
  
  sBeamPipeTubs[1]= new G4Tubs("sBeamPipeTubs[1]",106.0*mm/2,150.0*mm/2,16.0*mm/2,0*deg,360*deg);
  lBeamPipeTubs[1] = new G4LogicalVolume(sBeamPipeTubs[1],fMaterialList->GetMaterial("Al"),"lBeamPipeTubs[1]",0,0,0);

  sBeamPipeTubs[2]= new G4Tubs("sBeamPipeTubs[2]",144.0*mm/2,400.0*mm/2,16.0*mm/2,0*deg,360*deg);
  lBeamPipeTubs[2] = new G4LogicalVolume(sBeamPipeTubs[2],fMaterialList->GetMaterial("Al"),"lBeamPipeTubs[2]",0,0,0);

  // The bottom of the 90 degree target cylinder
  sBeamPipeTubs[3]= new G4Tubs("sBeamPipeTubs[3]",0.0*mm/2,150.0*mm/2,3.0*mm/2,0*deg,360*deg);
  lBeamPipeTubs[3] = new G4LogicalVolume(sBeamPipeTubs[3],fMaterialList->GetMaterial("Al"),"lBeamPipeTubs[3]",0,0,0);


  sBeamPipeOut[0] = new G4Tubs("sBeamPipeOut[0]",0.0*mm/2,outside*mm/2,1184.0*mm/2,0*deg,360*deg);
  sBeamPipeOut[1] = new G4Tubs("sBeamPipeOut[1]",0.0*mm/2,150.0*mm/2, 222.0*mm/2,0*deg,360*deg);


  if(oldPipe==true)
    sBeamPipeIn[0] = new G4Tubs("sBeamPipeIn[0]",0.0*mm,144.0*mm/2,1200.0*mm/2,0*deg,360*deg);
  else sBeamPipeIn[0]  = new G4Tubs("sBeamPipeIn[0]",0.0*mm,inside*mm/2,1200.0*mm/2,0*deg,360*deg);
  sBeamPipeIn[1]  = new G4Tubs("sBeamPipeIn[1]",0.0*mm,144.0*mm/2, 250.0*mm/2,0*deg,360*deg);
 
 
  G4RotationMatrix rot3DBeamPipe;
  G4RotationMatrix rot3DBeamPipe2;
  rot3DBeamPipe.set(0, 0, 0);
  rot3DBeamPipe.rotateX(90.0*degree); 

  G4VisAttributes* visAttBeamPipe = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  if(oldPipe==true)  {
    // Combining the two outer chambers, second one rotated by 90 degrees, shift back by 224 mm the physical placement has to put it to the target position
    sBeamPipeUnion = new G4UnionSolid("sBeamPipeUnion",sBeamPipeOut[0],sBeamPipeOut[1],&rot3DBeamPipe,G4ThreeVector(0.0*mm,-9.0*mm,-176.0*mm));
    
  // The subtraction geometry:
    sBeamPipeSubtraction[0] = new G4SubtractionSolid("sBeamPipeSubtraction[0]",sBeamPipeUnion,sBeamPipeIn[0],0,G4ThreeVector());

    sBeamPipeSubtraction[1] = new G4SubtractionSolid("sBeamPipeSubtraction[1]",sBeamPipeSubtraction[0],sBeamPipeIn[1],&rot3DBeamPipe,
                                                     G4ThreeVector(0.0*mm,-9.*mm,-176.0*mm));
    
        
    //Making the physical parts:
    pBeamPipeTubs[1] = new G4PVPlacement(G4Transform3D(rot3DBeamPipe,G4ThreeVector(0.0*m,110.0*mm,0.0*mm)),lBeamPipeTubs[1],"pBeamPipeTubs[1]",lWorld,false, 999999);
    pBeamPipeTubs[3] = new G4PVPlacement(G4Transform3D(rot3DBeamPipe,G4ThreeVector(0.0*m,-121.5*mm,0.0*mm)),lBeamPipeTubs[3]," pBeamPipeTubs[3]",lWorld,false, 999999);
    
    pBeamPipeTubs[0] = new G4PVPlacement(G4Transform3D(rot3DBeamPipe2,G4ThreeVector(0.0*m,0.0*m,-426.0*mm)),lBeamPipeTubs[0]," pBeamPipeTubs[0]",lWorld,false, 999999);
    pBeamPipeTubs[2] = new G4PVPlacement(G4Transform3D(rot3DBeamPipe2,G4ThreeVector(0.0*m,0.0*m,776.0*mm)),lBeamPipeTubs[2]," pBeamPipeTubs[2]",lWorld,false, 999999);
    
    // Setting the Vis attributes:
    for(int i=0;i<3;i++) lBeamPipeTubs[i]->SetVisAttributes(visAttBeamPipe);
  }

  // The new beam pipe:
  else {
    sBeamPipeSubtraction[1] = new G4SubtractionSolid("sBeamPipeSubtraction[1]",sBeamPipeOut[0],sBeamPipeIn[0],0,G4ThreeVector());
  }
   
  // Making the logic and physical beam pipe:
  lBeamPipe = new G4LogicalVolume(sBeamPipeSubtraction[1],fMaterialList->GetMaterial("Al"),"lBeamPipe",0,0,0);
  lBeamPipe->SetVisAttributes(visAttBeamPipe);
  pBeamPipe = new G4PVPlacement(G4Transform3D(rot3DBeamPipe2,G4ThreeVector(0.0*m,.0*mm,176.0*mm)),lBeamPipe,"pBeamPipe",lWorld,false, 999999);

  cout<<"Exiting BeamPipe Constructor"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BeamPipe::~BeamPipe()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BeamPipe::ConstructShield()  {

  G4RotationMatrix rot3DShield;
  rot3DShield.set(0, 0, 0);

  if(ShieldRadius>0)  {
    if(PbShieldThickness>0)  {
      sPbShield = new G4Tubs("sPbShield",(ShieldRadius)*10*mm,(ShieldRadius+PbShieldThickness)*10*mm,1184.0*mm/2,0*deg,360*deg);
      lPbShield = new G4LogicalVolume(sPbShield,fMaterialList->GetMaterial("Pb"),"lPbShield",0,0,0);
      pPbShield = new G4PVPlacement(G4Transform3D(rot3DShield,G4ThreeVector(0.0*m,.0*mm,176.0*mm)),lPbShield,"pPbShield",lWorld,false, 999998);
      G4VisAttributes* visAttPbShield = new G4VisAttributes(G4Colour(.0,1.0,1.0));
      lPbShield->SetVisAttributes(visAttPbShield);
    }
    
    if(SnShieldThickness>0)  {
      sSnShield = new G4Tubs("sSnShield",(ShieldRadius+PbShieldThickness)*10*mm,
                             (ShieldRadius+PbShieldThickness+SnShieldThickness)*10*mm,1184.0*mm/2,0*deg,360*deg);
      lSnShield = new G4LogicalVolume(sSnShield,fMaterialList->GetMaterial("Sn"),"lSnShield",0,0,0);
      pSnShield = new G4PVPlacement(G4Transform3D(rot3DShield,G4ThreeVector(0.0*m,.0*mm,176.0*mm)),lSnShield,"pSnShield",lWorld,false, 999997);
      G4VisAttributes* visAttSnShield = new G4VisAttributes(G4Colour(.5,1.0,1.0));
      lSnShield->SetVisAttributes(visAttSnShield);
    }
  }
}




