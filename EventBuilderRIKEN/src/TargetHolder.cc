#include "TargetHolder.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TargetHolder::TargetHolder(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld)  {
  cout<<"TargetHolder Constructor"<<endl;

  fMaterialList        = ptMaterialList;
  lWorld               = ptWorld;

  int idummy;
  char cdummy[20];
  float fdummy;
  fFileIn = fopen("./input/TargetHolder.in","r");
  while(!feof(fFileIn))  {
    fscanf(fFileIn,"%i %s %f",&idummy,&cdummy,&fdummy);
    if(idummy==0) fSetTarget = (int)fdummy;
    else if(idummy>=1 && idummy <=6)  {
      for (int i =0;i<20;i++) fTargetMaterial[idummy-1][i] = cdummy[i];
      fTargetThickness[idummy-1] = fdummy*10.;
    }
    else break; 
  }
  fclose(fFileIn);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TargetHolder::~TargetHolder()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TargetHolder::CreateTargetHolder()  {
  G4cout<<"Creating the target holder"<<endl;

  //______________________________________________
  // The flange:
  G4RotationMatrix rot3DFlange;
  rot3DFlange.set(0, 0, 0);
  rot3DFlange.rotateX(90.0*degree); 

  sTopFlange    = new G4Tubs("sTopFlange",0.0*mm/2,150.0*mm/2,16.0*mm/2,0*deg,360*deg);
  lTopFlange    = new G4LogicalVolume(sTopFlange,fMaterialList->GetMaterial("Al"),"lTopFlange",0,0,0);
  pTopFlange    = new G4PVPlacement(G4Transform3D(rot3DFlange,G4ThreeVector(0.0*mm,126.0*mm,0.0*mm)),lTopFlange,"pTopFlange",lWorld,false, 999998);

  // Creating the plastic holder:
  sTargetBoxRaw = new G4Box("sTargetBoxRaw",40.*mm/2.0,40.*mm/2.0,6.0*mm/2.0);
  sTargetHole   = new G4Tubs("sTargetHole",0.0*mm/2,30.0*mm/2,7.0*mm/2,0*deg,360*deg);

  sTargetBox    = new G4SubtractionSolid("sTargetBox",sTargetBoxRaw,sTargetHole,0,G4ThreeVector());
  lTargetBox    = new G4LogicalVolume(sTargetBox,fMaterialList->GetMaterial("Delrin"),"lTargetBox",0,0,0);


  float holeCenterX,holeCenterY;

  float targetX;
  if(fSetTarget<4) targetX = fSetTarget;
  else targetX = fSetTarget - 3;

  for(int i = 0;i<6;i++)  {
    if(i<3){holeCenterX =   0.0; holeCenterY = (1. - i+(targetX-2))*40.;}
    else   {holeCenterX = -48.0; holeCenterY = (4. - i+(targetX-2))*40.;}

    G4ThreeVector holeCenter(holeCenterX*mm,holeCenterY*mm,0.0);
    sprintf(fTemp,"pTargetBox[%i]",i);
    pTargetBox[i] = new G4PVPlacement(0,holeCenter,lTargetBox,fTemp,lWorld,false, 999998);
 
    //Placing the different targets:
    sTarget   = new G4Tubs("sTarget",0.0*mm/2,30.0*mm/2,fTargetThickness[i]*mm/2,0*deg,360*deg);
    lTarget[i]= new G4LogicalVolume(sTarget,fMaterialList->GetMaterial(fTargetMaterial[i]),"lTarget",0,0,0);

    G4RotationMatrix rot3DTarget;
    
    int zetShift = 1;
    if(fSetTarget<4) {if(i<3) holeCenterX =  0.0;else { holeCenterX = -48.0; zetShift = -1;}}
    else           {if(i<3){holeCenterX =-48.0; zetShift = -1;} else holeCenterX = 0.0;}

    
      G4ThreeVector targetPosition(holeCenterX*mm,holeCenterY*mm,(zetShift*(fTargetThickness[i]/2+3))*mm);
    sprintf(fTemp,"pTarget[%i]",i);
    pTarget[i] = new G4PVPlacement(G4Transform3D(rot3DTarget,targetPosition),lTarget[i],fTemp,lWorld,false, 999997);
  }   
  //______________________________________________
  // Setting the Vis attributes:
  G4VisAttributes* visAttTopFlange = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* visAttTargetBox = new G4VisAttributes(G4Colour(0.8,1.0,1.0));
  G4VisAttributes* visAttTargetC   = new G4VisAttributes(G4Colour(0.,.0,.0));
  G4VisAttributes* visAttTargetPb  = new G4VisAttributes(G4Colour(0.9,1.0,1.0));

  lTopFlange   ->SetVisAttributes(visAttTopFlange);
  lTargetBox   ->SetVisAttributes(visAttTargetBox);
  for(int i=0;i<6;i++)  {
   if(strcmp(fTargetMaterial[i],"Pb")==0) lTarget[i]->SetVisAttributes(visAttTargetPb);
   else if(strcmp(fTargetMaterial[i],"C")==0) lTarget[i]->SetVisAttributes(visAttTargetC);
  }
}
