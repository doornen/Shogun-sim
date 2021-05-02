#include "Sphere.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Sphere::Sphere(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld)  {
  cout<<"Sphere Constructor"<<endl;
  fMaterialList = ptMaterialList;
  lWorld = ptWorld;

  //Some initial values:
  fInnerRadius = 20.;
  fOuterRadius = 26.;
  
  //CreateArrayXYZPsiThetaPhi();
  cout<<"Exiting Sphere Constructor"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Sphere::~Sphere() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Sphere::ResetValues() { 
  fFlagFirstInteraction = false;
  fEnergy = 0.;
  fTime = 0.;
  fFirstInteractionX = 0.;
  fFirstInteractionY = 0.;
  fFirstInteractionZ = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Sphere::CreateSphere()  {
  sSphere = new G4Sphere("sSphere",fInnerRadius*cm, fOuterRadius*cm,
                         0.,100.,//if the delta angle is >=2*pi, or >=pi the shape is treated as
                         0.,100.);  // continuous in phi or theta respectively.
  lSphere = new G4LogicalVolume(sSphere,fMaterialList->GetMaterial("LaBr3"),"lSphere",0,0,0);
  pSphere = new G4PVPlacement(0,G4ThreeVector(0.0*m,0.0*m,0.0*m),lSphere,"pSphere",lWorld,false,500000);
 
  // Setting the vis attributes:
  G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); 
  lSphere->SetVisAttributes(visAtt);
}
