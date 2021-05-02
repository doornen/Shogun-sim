#ifndef Shogun_h
#define Shogun_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
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
#include "Randomize.hh"
#include "G4ios.hh"
#include <iostream>
#include "TMath.h"

#include "MaterialList.hh"
#include "Globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class Shogun  {
public:
  
  Shogun(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld);
  ~Shogun();
  
  void   AddCrystalEnergy(int a,float b){fCrystalEnergy[a] = fCrystalEnergy[a] + b;}
  void   CreateArrayXYZPsiThetaPhi();
  void   DetermineCrystalMult(){for(int i=0;i<NUMBEROFSHOGUNCRYSTALS;i++)fCrystalMult = fCrystalMult + fCrystalFlag[i];}
                                              
  float  GetCrystalEnergy(int a){return fCrystalEnergy[a];}
  bool   GetCrystalFlag(int a){return fCrystalFlag[a];}
  float  GetCrystalMeasuredEnergy(int a);
  float  GetCrystalMeasuredTime(int a);
  int    GetCrystalMult(){return fCrystalMult;};
  float  GetCrystalTime(int a){return fCrystalTime[a];}
  int    GetFIPointOption(){return fFIPointOption;}
  double GetFITime(int a){return fFITime[a];}
  float  GetFIX(int a){return fFIX[a];}
  float  GetFIY(int a){return fFIY[a];}
  float  GetFIZ(int a){return fFIZ[a];}   
  double GetPosXCrystal(int i){return fPosX[i];}
  double GetPosYCrystal(int i){return fPosY[i];}
  double GetPosZCrystal(int i){return fPosZ[i];}

  void   ResetValues();

  //void   SetAlAbsorberThickness(double thicknessAl){fAbsorberThicknessAl = 10.*thicknessAl*mm;}
  void   SetCrystalFlagTrue(int a){fCrystalFlag[a] = true;}
  // void   SetCrystalSize(float X, float Y, float Z)
  // {
  //   fCrystalSizeX = X; fCrystalSizeY = Y; fCrystalSizeZ = Z;
  //   G4cout<<"The Size of the Shogun detectors is : "<<fCrystalSizeX<<" "<<fCrystalSizeY<<" "<<fCrystalSizeZ<<G4endl;
  //   return;
  // } 
  void   SetCrystalTime(int a,float b){fCrystalTime[a] = b;}
  void   SetEnergyResolution(int a,float b,float c){fTypeOfEnergyResolution = a; fEnergyResolution[0] = b;fEnergyResolution[1] = c;}
  void   SetFIPointOption(int a){fFIPointOption = a;}
  void   SetFITime(int a, double b){fFITime[a] = b;}
  void   SetFIX(int a, float b){fFIX[a] = b;}
  void   SetFIY(int a, float b){fFIY[a] = b;}
  void   SetFIZ(int a, float b){fFIZ[a] = b;}
  void   SetHousingThickness(float X, float Y, float Z) {
    fHousingThicknessX = X; fHousingThicknessY = Y; fHousingThicknessZ = Z;
    G4cout<<"The thickness of the Shogun housing is: "<<fHousingThicknessX<<" "<<fHousingThicknessY<<" "<<fHousingThicknessZ<<G4endl;
    return;
  } 
  void   SetMgOThickness(float X, float Y, float Z) {
    fMgOThicknessX = X; fMgOThicknessY = Y; fMgOThicknessZ = Z;
    G4cout<<"The thickness of the Shogun MgO insulation is: "<<fMgOThicknessX<<" "<<fMgOThicknessY<<" "<<fMgOThicknessZ<<G4endl;
    return;
  } 
  void SetZPosShift(float a){
    fZPosShift = a;
    G4cout<<"Zpos shift is: "<<fZPosShift<<G4endl;
    return;}

  //void   SetPbAbsorberThickness(double thicknessPb){fAbsorberThicknessPb = 10.*thicknessPb*mm;}
  //void   SetSnAbsorberThickness(double thicknessSn){fAbsorberThicknessSn = 10.*thicknessSn*mm;}
  void   SetTimeResolution(float a, float b){fTimeResolution[0] = a; fTimeResolution[1] = b;}
  
private:
  
  G4LogicalVolume*    lWorld;                   //need a pointer to the logic world volume
  MaterialList*       fMaterialList;            //Pointer to the material list
     
  G4Trap*             sShogunCrystalTrap; 
  G4Tubs*             sShogunCrystalTubs;
  G4LogicalVolume*    lShogunCrystal;         
  G4VPhysicalVolume*  pShogunCrystal[NUMBEROFSHOGUNCRYSTALS];
  //____________________________________________________________________________________________
  // The housing
  G4Trap*             sShogunAlHouseOut;
  G4Trap*             sShogunAlHouseIn;
  G4SubtractionSolid* sShogunAlHouse;
  G4LogicalVolume*    lShogunAlHouse;
  G4VPhysicalVolume*  pShogunAlHouse[NUMBEROFSHOGUNCRYSTALS];
  //____________________________________________________________________________________________
  // The MgO insulation
  G4Trap*             sShogunMgOOutTrap;
  G4Trap*             sShogunMgOInTrap;
  G4Tubs*             sShogunMgOOutTubs;
  G4Tubs*             sShogunMgOInTubs;
  G4SubtractionSolid* sShogunMgO;
  G4LogicalVolume*    lShogunMgO;
  G4VPhysicalVolume*  pShogunMgO[NUMBEROFSHOGUNCRYSTALS];
  //____________________________________________________________________________________________
  // The lightGuide
  G4Trap*             sShogunLightGuide;
  G4LogicalVolume*    lShogunLightGuide;
  G4VPhysicalVolume*  pShogunLightGuide[1];
  //____________________________________________________________________________________________
  // The PMT
  G4Tubs*             sShogunPMT;
  G4LogicalVolume*    lShogunPMT;
  G4VPhysicalVolume*  pShogunPMT[4];
     
  //double fAbsorberThicknessAl; 
  //double fAbsorberThicknessPb; 
  //double fAbsorberThicknessSn;
  bool   fCrystalFlag[NUMBEROFSHOGUNCRYSTALS];
  float  fCrystalEnergy[NUMBEROFSHOGUNCRYSTALS];
  int    fCrystalMult;

  float  fCrystalTime[NUMBEROFSHOGUNCRYSTALS];
  float  fEnergyResolution[2];
  int    fFIPointOption;
  double fFITime[NUMBEROFSHOGUNCRYSTALS];
  float  fFIX[NUMBEROFSHOGUNCRYSTALS];
  float  fFIY[NUMBEROFSHOGUNCRYSTALS];
  float  fFIZ[NUMBEROFSHOGUNCRYSTALS];
  FILE*  fFileIn;
  FILE*  fFileOut;
  float  fHousingThicknessX;
  float  fHousingThicknessY;
  float  fHousingThicknessZ;
  float  fMgOThicknessX;
  float  fMgOThicknessY;
  float  fMgOThicknessZ;
  int    fNumberOfCrystals;
  double fPosX[NUMBEROFSHOGUNCRYSTALS];
  double fPosY[NUMBEROFSHOGUNCRYSTALS];
  double fPosZ[NUMBEROFSHOGUNCRYSTALS];
  char   fTempname[200];
  int    fTypeOfEnergyResolution;
  float  fTimeResolution[2];
  float  fZPosShift;
};
#endif
