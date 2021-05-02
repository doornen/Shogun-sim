#ifndef SGT_h
#define SGT_h 1

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
#include "TMath.h"
#include "Randomize.hh"
#include "G4ios.hh"
#include <iostream>

#include  "MaterialList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SGT  {
  public:
  
    SGT(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld);
   ~SGT();

    void   AddCoaxEnergy(int a,float b){fCoaxEnergy[a] = fCoaxEnergy[a] + b;}
    void   AddDetectorEnergy(int a,float b){fDetectorEnergy[a] = fDetectorEnergy[a] + b;}
    void   AddPlanarEnergy(int a,int b,float c){fPlanarEnergy[a][b] = fPlanarEnergy[a][b] + c;}
    void   CreateArrayXYZPsiThetaPhi();
    void   DetermineDetectorMult(){for(int i=0;i<10;i++)fDetectorMult = fDetectorMult + fDetectorFlag[i];}
    
    float  GetCoaxEnergy(int i){return fCoaxEnergy[i];}
    int    GetCoaxFlag(int i){return fCoaxFlag[i];}
    float  GetCoaxMeasuredEnergy(int a);
    float  GetCoaxMeasuredTime(int a);
    float  GetCoaxTime(int i){return fCoaxTime[i];}
    float  GetDetectorEnergy(int a){return fDetectorEnergy[a];}
    float  GetDetectorMult();
    bool   GetDetectorFlag(int a){return fDetectorFlag[a];}
    float  GetPlanarEnergy(int i,int j){return fPlanarEnergy[i][j];}
    int    GetPlanarFlag(int i,int j){return fPlanarFlag[i][j];}
    float  GetPlanarMeasuredEnergy(int i,int j);
    float  GetPlanarMeasuredTime(int i,int j);
    float  GetPlanarTime(int i,int j){return fPlanarTime[i][j];}
    double GetPosXCrystal(int i,int j){return fPosX[i][j];}
    double GetPosYCrystal(int i,int j){return fPosY[i][j];}
    double GetPosZCrystal(int i,int j){return fPosZ[i][j];}

    void   ResetValues();

  //void   SetAlAbsorberThickness(double thicknessAl){fAbsorberThicknessAl = 10.*thicknessAl*mm;};
    void   SetCoaxFlagTrue(int a){fCoaxFlag[a] = true;};
    void   SetCoaxTime(int a,float b){fCoaxTime[a] = b;};
    void   SetDetectorFlagTrue(int a){fDetectorFlag[a] = true;};
    void   SetEnergyResolution(int a,float b,float c){fTypeOfEnergyResolution=a;fEnergyResolution[0]=b;fEnergyResolution[1]=c;};
    void   SetPlanarFlagTrue(int a, int b){fPlanarFlag[a][b] = true;};
  //void   SetPbAbsorberThickness(double thicknessPb){fAbsorberThicknessPb = 10.*thicknessPb*mm;};
  //void   SetSnAbsorberThickness(double thicknessSn){fAbsorberThicknessSn = 10.*thicknessSn*mm;};
    void   SetTimeResolution(float a, float b){fTimeResolution[0] = a;fTimeResolution[1]=b;};

  private:

     G4LogicalVolume*    lWorld;                   //need a pointer to the logic world volume
     MaterialList*       fMaterialList;            //Pointer to the material list

     G4Box*              sSGTStrip;
     G4LogicalVolume*    lSGTStrip; 
     G4VPhysicalVolume*  pSGTStrip[10][25]; 
     G4Tubs*             sSGTCoaxial;
     G4LogicalVolume*    lSGTCoaxial;
     G4VPhysicalVolume*  pSGTCoaxial[10];

  //double fAbsorberThicknessAl; 
  //double fAbsorberThicknessPb; 
  //double fAbsorberThicknessSn;
     float  fCoaxEnergy[10];
     bool   fCoaxFlag[10];
     float  fCoaxTime[10];
     float  fDetectorEnergy[10];
     bool   fDetectorFlag[10];
     int    fDetectorMult;
     float  fEnergyResolution[2];
     FILE*  fFileIn;
     FILE*  fFileOut;
     double fPosX[10][30];
     double fPosY[10][30];
     double fPosZ[10][30];
     float  fPlanarEnergy[10][26];
     bool   fPlanarFlag[10][26];
     float  fPlanarTime[10][26];
     char   fTempname[200];
     int    fTypeOfEnergyResolution;
     float  fTimeResolution[2];
};
#endif
