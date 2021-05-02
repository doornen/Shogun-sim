#ifndef LaBr3Array_h
#define LaBr3Array_h 1

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
class LaBr3Array  {
public:
  
  LaBr3Array(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld);
  ~LaBr3Array();
  
  void   AddCrystalEnergy(int a,float b){fCrystalEnergy[a] = fCrystalEnergy[a] + b;}
  void   CreateArrayXYZPsiThetaPhi();
  void   DetermineCrystalMult()  {
    for(int i=0;i<NUMBEROFLABR3ARRAYCRYSTALS;i++)fCrystalMult = fCrystalMult + fCrystalFlag[i];
  }

  float  GetCrystalEnergy(int a){return fCrystalEnergy[a];}
  bool   GetCrystalFlag(int a){return fCrystalFlag[a];}
  float  GetCrystalMeasuredEnergy(int a);
  float  GetCrystalMeasuredTime(int a);
  int    GetCrystalMult(){return fCrystalMult;};
  float  GetCrystalTime(int a){return fCrystalTime[a];}
  int    GetNumberOfCrystals(){return fNumberOfCrystals;}
  double GetPosXCrystal(int i){return fPosX[i];}
  double GetPosYCrystal(int i){return fPosY[i];}
  double GetPosZCrystal(int i){return fPosZ[i];}

  void   ResetValues();

  //void   SetAlAbsorberThickness(double thicknessAl){fAbsorberThicknessAl = 10.*thicknessAl*mm;}
  void   SetCrystalFlagTrue(int a){fCrystalFlag[a] = true;}

  void   SetCrystalTime(int a,float b){fCrystalTime[a] = b;}
  void   SetEnergyResolution(int a,float b,float c){fTypeOfEnergyResolution = a; fEnergyResolution[0] = b;fEnergyResolution[1] = c;}
  void   SetEnergyResolution(int a,float b,float c,float d){fTypeOfEnergyResolution = a; fEnergyResolution[0] = b;fEnergyResolution[1] = c;fEnergyResolution[2] = d;}
  void   SetHousingThickness(float X, float Y) {
    fHousingThicknessFront = X; fHousingThicknessSide = Y;
    G4cout<<"The thickness of the LaBr3Array housing is: "<<fHousingThicknessFront<<" "<<fHousingThicknessSide<<G4endl;
    if(fHousingThicknessFront<0 || fHousingThicknessSide<0)  {
      G4cout<<"At least one of the values is smaller than zero. Aborting programm... "<<G4endl;
      abort();
    }
    else{fHousing=true;}
    return;
  }
  void   SetHousing(bool doHousing, float X, float Y, bool transparent);
  
  void   SetInsulationThickness(float X, float Y) {
    fInsulationThicknessFront = X; fInsulationThicknessSide = Y;
    G4cout<<"The thickness of the LaBr3Array insulation is: "<<fInsulationThicknessFront<<" "<<fInsulationThicknessSide<<G4endl;
    if(fInsulationThicknessFront<0 || fInsulationThicknessSide<0)  {
      G4cout<<"At least one of the values is smaller than zero. Aborting programm... "<<G4endl;
      abort();
    }
    else{fInsulation=true;}
    return;
  }
  void   SetInsulation(bool doInsulation, float X, float Y, bool transparent);
  
  //void   SetPbAbsorberThickness(double thicknessPb){fAbsorberThicknessPb = 10.*thicknessPb*mm;}
  //void   SetSnAbsorberThickness(double thicknessSn){fAbsorberThicknessSn = 10.*thicknessSn*mm;}
  void   SetTimeResolution(float a, float b){fTimeResolution[0] = a; fTimeResolution[1] = b;}
  void   SetZPosShift(float a){fZPosShift = a;}
  
private:
  
  G4LogicalVolume*    lWorld;                   //need a pointer to the logic world volume
  MaterialList*       fMaterialList;            //Pointer to the material list
  
  G4Tubs*             sLaBr3ArrayCrystal; 
  G4LogicalVolume*    lLaBr3ArrayCrystal;         
  G4VPhysicalVolume*  pLaBr3ArrayCrystal[NUMBEROFLABR3ARRAYCRYSTALS];
  //____________________________________________________________________________________________
  // The housing
  G4Tubs*             sLaBr3ArrayAlHouseOut;
  G4Tubs*             sLaBr3ArrayAlHouseIn;
//   G4SubtractionSolid* sLaBr3ArrayAlHouse; //real housing around the crystal
  G4UnionSolid*	      sLaBr3ArrayAlHouse; //real housing around the crystal
  G4LogicalVolume*    lLaBr3ArrayAlHouse;
  G4VPhysicalVolume*  pLaBr3ArrayAlHouse[NUMBEROFLABR3ARRAYCRYSTALS];
  
  G4Tubs*             sLaBr3ArrayAlHouseFront; //front cap
  G4Tubs*             sLaBr3ArrayAlHouseBackFlange;
  G4Tubs*             sLaBr3ArrayAlHouseBackPMT;
  G4UnionSolid*       sLaBr3ArrayAlHouseBack; //flange and pnt at the end
  G4UnionSolid*       sLaBr3ArrayAlHouseOutBack;
  G4SubtractionSolid* sLaBr3ArrayAlHouseOutBackIn;
  
  //____________________________________________________________________________________________
  // The insulation
  G4Tubs*             sLaBr3ArrayInsulationOut;
  G4Tubs*             sLaBr3ArrayInsulationIn;
  G4SubtractionSolid* sLaBr3ArrayInsulation;
  G4LogicalVolume*    lLaBr3ArrayInsulation;
  G4VPhysicalVolume*  pLaBr3ArrayInsulation[NUMBEROFLABR3ARRAYCRYSTALS];
  
  //double fAbsorberThicknessAl; 
  //double fAbsorberThicknessPb; 
  //double fAbsorberThicknessSn;
  bool   fCrystalFlag[NUMBEROFLABR3ARRAYCRYSTALS];
  float  fCrystalEnergy[NUMBEROFLABR3ARRAYCRYSTALS];
  int    fCrystalMult;

  float  fCrystalTime[NUMBEROFLABR3ARRAYCRYSTALS];
//   float  fEnergyResolution[2];
  float  fEnergyResolution[3]; //pschrock: 3 par needed for resolution function from NIM paper
  FILE*  fFileIn;
  FILE*  fFileOut;
  bool   fHousing;
  bool   fInsulation;
  bool   fTransparentHousing;
  bool   fTransparentInsulation;
  float  fHousingThicknessFront;
  float  fHousingThicknessSide;
  char   fInsulationMaterialName[100];
  float  fInsulationThicknessFront;
  float  fInsulationThicknessSide;
  int    fNumberOfCrystals;
  double fPosX[NUMBEROFLABR3ARRAYCRYSTALS];
  double fPosY[NUMBEROFLABR3ARRAYCRYSTALS];
  double fPosZ[NUMBEROFLABR3ARRAYCRYSTALS];
  char   fTempname[200];
  int    fTypeOfEnergyResolution;
  float  fTimeResolution[2];
  float  fZPosShift;
};
#endif
