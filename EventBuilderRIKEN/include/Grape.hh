#ifndef Grape_h
#define Grape_h 1

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

#include "G4ios.hh"
#include <iostream>
#include "Randomize.hh"
#include "TMath.h"

#include "MaterialList.hh"
#include "Globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// A Crystal is composed out of 9 segments. The members are ordered as follows:
// parameter[x][y][0] for the entire crystal
// parameter[x][y][1-9] for the 9 segments
class Grape  {
  public:
  
    Grape(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld);
   ~Grape();

    void   AddCrystalEnergy(int a,int b,int c,float d){fCrystalEnergy[a][b][c] = fCrystalEnergy[a][b][c] + d;}
    void   CreateArrayXYZPsiThetaPhi();
    void   DetermineCrystalMult(){
      for(int i=0;i<20;i++)  {
        for(int j=0;j<2;j++)  {
          fCrystalMult=fCrystalMult+fCrystalFlag[i][j][0];
        }
      }
    }
    float  GetCrystalEnergy(int a,int b,int c){return fCrystalEnergy[a][b][c];}
    bool   GetCrystalFlag(int a,int b,int c){return fCrystalFlag[a][b][c];}
    float  GetCrystalMeasuredEnergy(int a,int b,int c);
    float  GetCrystalMeasuredTime(int a,int b,int c);
    int    GetCrystalMult(){return fCrystalMult;}
    float  GetCrystalTime(int a,int b,int c){return fCrystalTime[a][b][c];}
    double GetPosXCrystal(int a,int b,int c){return fPosX[a][b][c];}
    double GetPosYCrystal(int a,int b,int c){return fPosY[a][b][c];}
    double GetPosZCrystal(int a,int b,int c){return fPosZ[a][b][c];}

    void   ResetValues();

  //void   SetAlAbsorberThickness(double thicknessAl){fAbsorberThicknessAl = 10.*thicknessAl*mm;}
    void   SetCrystalFlagTrue(int a,int b,int c){fCrystalFlag[a][b][c] = true;}
    void   SetCrystalTime(int a,int b,int c,float d){fCrystalTime[a][b][c] = d;}
    void   SetEnergyResolution(int a,float b,float c)  {
      fTypeOfEnergyResolution = a; fEnergyResolution[0] = b;fEnergyResolution[1] = c;
    }
  //void   SetPbAbsorberThickness(double thicknessPb){fAbsorberThicknessPb = 10.*thicknessPb*mm;}
  //void   SetSnAbsorberThickness(double thicknessSn){fAbsorberThicknessSn = 10.*thicknessSn*mm;}
    void   SetTimeResolution(float a, float b){fTimeResolution[0] = a;fTimeResolution[1]=b;}

  private:

     G4LogicalVolume*    lWorld;                   //need a pointer to the logic world volume
     MaterialList*       fMaterialList;            //Pointer to the material list
    
     G4Tubs*             sGrapeCrystalTube;
     G4Box*              sGrapeCrystalBox;
     G4IntersectionSolid* sGrapeSegment[9];
     G4LogicalVolume*    lGrapeSegment[9];
     G4VPhysicalVolume*  pGrapeSegment[20][2][9];
     // The dewar-----------------
     G4Tubs*             sGrapeDewar;
     G4LogicalVolume*    lGrapeDewar;
     G4VPhysicalVolume*  pGrapeDewar[20];
     //The Al bar
     G4Tubs*             sGrapeAlBar;
     G4LogicalVolume*    lGrapeAlBar;
     G4VPhysicalVolume*  pGrapeAlBar[20];

  //double fAbsorberThicknessAl; 
  //double fAbsorberThicknessPb; 
  //double fAbsorberThicknessSn;
     bool   fCrystalFlag[20][2][10];
     float  fCrystalEnergy[20][2][10];
     int    fCrystalMult;
     float  fCrystalTime[20][2][10];
     float  fEnergyResolution[2];
     FILE*  fFileIn;
     FILE*  fFileOut;
     double fPosX[20][2][10];
     double fPosY[20][2][10];
     double fPosZ[20][2][10];
     char   fTempname[200];
     int    fTypeOfEnergyResolution;
     float  fTimeResolution[2];
};
#endif
