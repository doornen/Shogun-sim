#ifndef Sphere_h
#define Sphere_h 1

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"

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

#include "MaterialList.hh"
#include "Globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class Sphere  {
  public:
  
    Sphere(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld);
   ~Sphere();

    void  AddEnergy(float a){fEnergy = fEnergy + a;}
    void  CreateSphere();

    float GetEnergy(){return fEnergy;}
    float GetFirstInteractionX(){return fFirstInteractionX;}
    float GetFirstInteractionY(){return fFirstInteractionY;}
    float GetFirstInteractionZ(){return fFirstInteractionZ;}
    bool  GetFlagFirstInteraction(){return fFlagFirstInteraction;}
    float GetTime(){return fTime;};

    void  ResetValues();

    void  SetDimensions(float a,float b){fOuterRadius=a;fInnerRadius=b;}
    void  SetFirstInteractionX(float a){fFirstInteractionX = a;}
    void  SetFirstInteractionY(float a){fFirstInteractionY = a;}
    void  SetFirstInteractionZ(float a){fFirstInteractionZ = a;}    
    void  SetFlagFirstInteraction(bool a){fFlagFirstInteraction = a;}
   
  private:

     G4LogicalVolume*    lWorld;                   //need a pointer to the logic world volume
     MaterialList*       fMaterialList;            //Pointer to the material list
     
     G4Sphere*          sSphere;
     G4LogicalVolume*   lSphere;
     G4VPhysicalVolume* pSphere;

     float fEnergy;
     float fFirstInteractionX;
     float fFirstInteractionY;
     float fFirstInteractionZ;
     bool  fFlagFirstInteraction;
     float fInnerRadius;
     float fOuterRadius;  
     char  fTempname[200];   
     float fTime;
};
#endif

