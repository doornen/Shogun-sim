#include "Collimator.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Collimator::Collimator(MaterialList *ptMaterialList, G4LogicalVolume *ptWorld) {
  cout<<"Collimator Constructor"<<endl;

  fMaterialList        = ptMaterialList;
  lWorld               = ptWorld;

  // Some start values:
  fBoxShape = 0;
  fCollHoleTarget=0.2;
  fCollHoleDetector=0.2;
  fCollBlockSize[0]=10.;
  fCollBlockSize[1]=10.;
  fCollBlockSize[2]=10.;
  fCollBlockComposition[0]=95.;
  fCollBlockComposition[1]=2.;
  fCollBlockComposition[2]=35.;
  fCollBlockDensity=17;
  fCollimatorPos[0]=0;
  fCollimatorPos[1]=0;
  fCollimatorPos[2]=0;
  
  char temp[100];

  FILE *fFileIn = fopen("./input/Collimator.in","r");
  while(!feof(fFileIn)) {
    fscanf(fFileIn,"%s ",temp);
    if(strcmp(temp,"GEOMETRY")==0)  { 
      fscanf(fFileIn,"%s ",temp);
      if(strcmp(temp,"BOX")==0) { 
        fBoxShape = true;
        fscanf(fFileIn,"%f %f %f",&fCollBlockSize[0],&fCollBlockSize[1],&fCollBlockSize[2]);
      }
      else fscanf(fFileIn,"%f %f",&fCollBlockSize[0],&fCollBlockSize[1]);   
    }
    if(strcmp(temp,"HOLESIZE")==0) 
      fscanf(fFileIn,"%f %f",&fCollHoleTarget,&fCollHoleDetector);
    if(strcmp(temp,"COMPOSITION")==0) 
      fscanf(fFileIn,"%f %f %f %f",&fCollBlockComposition[0],&fCollBlockComposition[1],
             &fCollBlockComposition[2],&fCollBlockDensity);
    if(strcmp(temp,"POS")==0) 
      fscanf(fFileIn,"%f %f %f",&fCollimatorPos[0],&fCollimatorPos[1],&fCollimatorPos[2]);
    else if(strcmp(temp,"END")==0) break;
  }
  fclose(fFileIn);
  
  //Set the collimator material:
  fMaterialList->SetCollimatorMaterial(fCollBlockComposition,fCollBlockDensity);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Collimator::~Collimator()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Collimator::CreateCollimator()  {
  G4cout<<"Creating the Collimator"<<endl;

  //______________________________________________
  // The flange:
  G4RotationMatrix rot3DCollimator;
  rot3DCollimator.set(0, 0, 0);
  //rot3DCollimator.rotateX(90.0*degree); 

  if(fBoxShape==true)
    sCollBlock  = new G4Box("sCollBlock",fCollBlockSize[0]*cm/2.0,fCollBlockSize[1]*cm/2.0,
                            fCollBlockSize[2]*cm/2.0);
  else sCollTubs = new G4Tubs("sCollTubs",0.0,fCollBlockSize[0]*cm/2.0,fCollBlockSize[1]*cm/2.0,0.0*deg,360.0*deg);

  sCollHole   = new G4Cons("sCollHole", 0,fCollHoleTarget*cm/2,0,fCollHoleDetector*cm/2,
                           (fCollBlockSize[2]+0.1)*cm/2,0*deg,360*deg);
  if(fBoxShape==true)
    sCollimator = new G4SubtractionSolid("sCollimator",sCollBlock,sCollHole,0,G4ThreeVector());
  else sCollimator = new G4SubtractionSolid("sCollimator",sCollTubs,sCollHole,0,G4ThreeVector());

  lCollimator = new G4LogicalVolume(sCollimator,fMaterialList->GetMaterial("DENSIMET"),
                                    "lCollimator",0,0,0);
 
  pCollimator = new G4PVPlacement(G4Transform3D(rot3DCollimator, G4ThreeVector(fCollimatorPos[0]*cm,
                                                                               fCollimatorPos[1]*cm,
                                                                               fCollimatorPos[2]*cm)),
                                  lCollimator,"pCollimator",lWorld,false, 999990);
  
  //______________________________________________
  // Setting the Vis attributes:
  G4VisAttributes* visAttCollBlock = new G4VisAttributes(G4Colour(0.8,1.0,1.0));
  lCollimator->SetVisAttributes(visAttCollBlock); 
}
