#include "Shogun.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Shogun::Shogun(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld)   {
  cout<<"Shogun Constructor"<<endl;

  fMaterialList        = ptMaterialList;
  lWorld               = ptWorld;
  //fAbsorberThicknessPb = fAbsorberThicknessSn = fAbsorberThicknessAl = 0.0;
  fNumberOfCrystals    = 0;
  
  fHousingThicknessX   = 0.;
  fHousingThicknessY   = 0.;
  fHousingThicknessZ   = 0.;

  fMgOThicknessX   = 0.;
  fMgOThicknessY   = 0.;
  fMgOThicknessZ   = 0.;

  sShogunPMT = new G4Tubs("sDaliPMT",0.0*mm,15.0*mm/2.0,60.0*mm/2.0,0.0*deg,360.0*deg);
  lShogunPMT = new G4LogicalVolume(sShogunPMT,fMaterialList->GetMaterial("Air"),"lShogunPMT",0,0,0);

  sShogunLightGuide = new G4Trap("sShogunLightGuide",1.5*cm/2.0,1.5*cm/2.0,1.5*cm/2.0,4.*cm/2.0,4*cm/2.0);
  lShogunLightGuide = new G4LogicalVolume(sShogunLightGuide,fMaterialList->GetMaterial("Air"),"lShogunLightGuide",0,0,0);

  G4VisAttributes* visAttPMT = new G4VisAttributes(G4Colour(1.0,1.0,1.0)); 
  lShogunPMT->SetVisAttributes(visAttPMT);

  for(int i=0;i<NUMBEROFSHOGUNCRYSTALS;i++)  {
    fPosX[i] = fPosY[i] = fPosZ[i] = 0.0;
  }
 
  //CreateArrayXYZPsiThetaPhi();
  cout<<"Exiting Shogun Constructor"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Shogun::~Shogun()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Shogun::ResetValues()  {
  fCrystalMult = 0;
  for(int i=0;i<NUMBEROFSHOGUNCRYSTALS;i++)  {
    fCrystalFlag[i]   = false;
    fCrystalEnergy[i] = 0.0;
    fCrystalTime[i]   = 0.0;
    fFIX[i]           = -999.0;
    fFIY[i]           = -999.0;
    fFIZ[i]           = -999.0;
    fFITime[i]        = 1000000;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float Shogun::GetCrystalMeasuredEnergy(int a)  {
  if(fCrystalEnergy[a]==0) return 0.;
  float dummy;
  if(fTypeOfEnergyResolution==1) dummy = (fEnergyResolution[0] + fEnergyResolution[1]*fCrystalEnergy[a]);
  else if(fTypeOfEnergyResolution==2) dummy = fEnergyResolution[0]*TMath::Power(fCrystalEnergy[a],fEnergyResolution[1]);
  else {cout<<"Wrong resolution option for SHOGUN. Aborting program."<<endl; abort();}
  //Observed energy in the detector
  return G4RandGauss::shoot(fCrystalEnergy[a],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float Shogun::GetCrystalMeasuredTime(int a)  {
  if(fCrystalTime[a]==0) return 0.;
  //Observed time in the detector
  float dummy = (fTimeResolution[0] - fTimeResolution[1]*TMath::Power(fCrystalEnergy[a],0.5));
  return G4RandGauss::shoot(fCrystalTime[a],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Shogun::CreateArrayXYZPsiThetaPhi()  {
  G4cout<<"Creating the Shogun array"<<endl;
  float x,y,z,psi,theta,phi,sizeX1,sizeX2,sizeY1,sizeY2,sizeZ;
  x = y = z = psi = theta = phi = sizeX1 = sizeX2 = sizeY1 = sizeY2 = sizeZ = 0.;
  int id;
  int ring;
  int type;
  float radius;
  int ShogunDet;
  int crystalsInDetector, crystalNumber;
  G4RotationMatrix Rot3D, Rot3D2;
    
  fFileIn  = fopen("./geometry/Shogun_geometry_in.txt","r");
  fFileOut = fopen("./geometry/Shogun_geometry_out.txt","w");
  for(int i=0;!feof(fFileIn)&& i<NUMBEROFSHOGUNCRYSTALS;i++)  {
    fscanf(fFileIn,"%i %i %i %f %f %f %f %f %f %f %f %f %f %f %i %i",&id,&ring,&type,&x,&y,&z,&radius,&theta,&phi,&sizeX2,&sizeX1,&sizeY2,&sizeY1,&sizeZ,&crystalsInDetector,&crystalNumber);
    G4cout<<"Read detector "<<i<<": "<<x<<" "<<y<<" "<<z<<" "<<psi<<" "<<theta<<" "<<phi<<endl;
    if(id==-1) break;
    //Type =1-> Trap, Type = 2 -> Tubs
    if(type==1) {
      sShogunCrystalTrap = new G4Trap("ShogunCrystal",sizeX1*cm/2.0,sizeX2*cm/2.0,sizeY1*cm/2.0,sizeY2*cm/2.0,sizeZ*cm/2.0);
      lShogunCrystal = new G4LogicalVolume(sShogunCrystalTrap,fMaterialList->GetMaterial("CEGAGG"),"lShogunCrystal",0,0,0);
      //lShogunCrystal = new G4LogicalVolume(sShogunCrystalTrap,fMaterialList->GetMaterial("CsI"),"lShogunCrystal",0,0,0);
      //lShogunCrystal = new G4LogicalVolume(sShogunCrystalTrap,fMaterialList->GetMaterial("CsI"),"lShogunCrystal",0,0,0);
    }
    else if(type==2) {
      sShogunCrystalTubs = new G4Tubs("LaBr3ArrayCrystal",0.0,sizeX1*cm/2.0,sizeZ*cm/2.0,0.0*deg,360.0*deg);
      lShogunCrystal = new G4LogicalVolume(sShogunCrystalTubs,fMaterialList->GetMaterial("CEGAGG"),"lShogunCrystal",0,0,0);
      //lShogunCrystal = new G4LogicalVolume(sShogunCrystalTubs,fMaterialList->GetMaterial("CsI"),"lShogunCrystal",0,0,0);
      //lShogunCrystal = new G4LogicalVolume(sShogunCrystalTubs,fMaterialList->GetMaterial("CsI"),"lShogunCrystal",0,0,0);
    }
    else break;
   
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5,0.5,1.0)); 
    lShogunCrystal->SetVisAttributes(visAtt);

    //The Housing
    if(fHousingThicknessX>0. && fHousingThicknessY>0. &&fHousingThicknessZ >0. && crystalNumber==1)  {
      sShogunAlHouseOut = new G4Trap("sShogunAlHouseOut",
                                     (sizeX1*crystalsInDetector+2*(fHousingThicknessX+fMgOThicknessX))*cm/2.0, 
                                      (sizeX2*crystalsInDetector+2*(fHousingThicknessX+fMgOThicknessX))*cm/2.0,
                                     (sizeY1+2*(fHousingThicknessY+fMgOThicknessY))*cm/2.0,
                                     (sizeY2+2*(fHousingThicknessY+fMgOThicknessY))*cm/2.0,
                                     (sizeZ+2*(fHousingThicknessZ+fMgOThicknessZ))*cm/2);
      sShogunAlHouseIn = new G4Trap("sShogunAlHouseIn",
                                    (sizeX1*crystalsInDetector+2*fMgOThicknessX)*cm/2.0, 
                                    (sizeX2*crystalsInDetector+2*fMgOThicknessX)*cm/2.0,
                                    (sizeY1+2*fMgOThicknessY)*cm/2.0,
                                    (sizeY2+2*fMgOThicknessY)*cm/2.0,
                                    (sizeZ+2*fMgOThicknessZ)*cm/2);

      sShogunAlHouse   = new G4SubtractionSolid("sShogunAlHouse",sShogunAlHouseOut,sShogunAlHouseIn);
      lShogunAlHouse   = new G4LogicalVolume(sShogunAlHouse,fMaterialList->GetMaterial("Al"),"lShogunAlHouse",0,0,0);
      
      // Setting the vis attributes:
      G4VisAttributes* visAttHouse = new G4VisAttributes(G4Colour(1.0,1.0,1.0)); 
      lShogunAlHouse->SetVisAttributes(visAttHouse);
    }
    // The MgO insulation:
    if(fMgOThicknessX>0. && fMgOThicknessY>0. &&fMgOThicknessZ >0. && crystalNumber==1)  {
      
      if(type==1) {
        sShogunMgOOutTrap = new G4Trap("sShogunMgOOut",
                                   (sizeX1*crystalsInDetector+2*fMgOThicknessX)*cm/2.0, 
                                   (sizeX2*crystalsInDetector+2*fMgOThicknessX)*cm/2.0,
                                   (sizeY1+2*fMgOThicknessY)*cm/2.0,
                                   (sizeY2+2*fMgOThicknessY)*cm/2.0,
                                   (sizeZ+2*fMgOThicknessZ)*cm/2);
        sShogunMgOInTrap = new G4Trap("sShogunMgOIn",
                                  (sizeX1*crystalsInDetector)*cm/2.0, 
                                  (sizeX2*crystalsInDetector)*cm/2.0,
                                  sizeY1*cm/2.0,
                                  sizeY2*cm/2.0,
                                  sizeZ*cm/2);
        sShogunMgO   = new G4SubtractionSolid("sShogunMgO",sShogunMgOOutTrap,sShogunMgOInTrap);
      }
      else if(type==2) {
        sShogunMgOOutTubs = new G4Tubs("sShogunMgOOut",
                                       0.0,(sizeX1+2*fMgOThicknessX)*cm/2.0,
                                       (sizeZ+2*fMgOThicknessX)*cm/2.0,
                                       0*deg,360*deg);
        
        sShogunMgOInTubs = new G4Tubs("sShogunMgOIn",
                                      0.0,sizeX1*cm/2.0,
                                      sizeZ*cm/2.0,
                                      0*deg,360*deg);
        sShogunMgO   = new G4SubtractionSolid("sShogunMgO",sShogunMgOOutTubs,sShogunMgOInTubs);
      }
      else break;
      
      lShogunMgO   = new G4LogicalVolume(sShogunMgO,fMaterialList->GetMaterial("MgO"),"lShogunMgO",0,0,0);
      
      // Setting the vis attributes:
      G4VisAttributes* visAttMgO = new G4VisAttributes(G4Colour(1.0,1.0,.0)); 
      lShogunMgO->SetVisAttributes(visAttMgO);
    }

    fNumberOfCrystals++;

    //The first digit is for the array type, The second digit is for the different
    //materials in the detector. The Last four digits are for the detector number
    ShogunDet = id -1;
    id       =  400000 + ShogunDet;

    z = z + fZPosShift; // Shift along beam axis
    cout<<"z: "<<z<<endl;

    G4ThreeVector Pos(x*cm,y*cm,z*cm); //Position of the center

    Rot3D.set(0, 0, 0);
    Rot3D.rotateX(psi*degree);
    Rot3D.rotateY(theta*degree);  
    Rot3D.rotateZ(phi*degree);  

    sprintf(fTempname,"pShogunCrystal[%i]",ShogunDet);

    pShogunCrystal[ShogunDet] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                                                  lShogunCrystal,fTempname,lWorld,false,id);

    // Putting the housing
    if(fHousingThicknessX>0. && fHousingThicknessY>0. && fHousingThicknessZ>0. && crystalNumber==1)  {
      sprintf(fTempname,"pShogunAlHouse[%i]",ShogunDet);

      G4ThreeVector Pos2((x+(-0.5+0.5*crystalsInDetector)*sizeX1*cos((theta)*degree)*cos(phi*degree))*cm,
                         (y+(-0.5+0.5*crystalsInDetector)*sizeX1*cos((theta)*degree)*sin(phi*degree))*cm,
                         (z+(0.5-0.5*crystalsInDetector)*sizeX1*sin((theta)*degree))*cm); 
      // Pieter
      pShogunAlHouse[ShogunDet] = new G4PVPlacement(G4Transform3D(Rot3D,Pos2),
                                                    lShogunAlHouse,fTempname,lWorld,false,id+10000);
    }
    
    // Putting the MgO insulator:
    if(fMgOThicknessX>0. && fMgOThicknessY>0. && fMgOThicknessZ>0. && crystalNumber==1)  {
      sprintf(fTempname,"pShogunMgO[%i]",ShogunDet);
      
      G4ThreeVector Pos2((x+(-0.5+0.5*crystalsInDetector)*sizeX1*cos((theta)*degree)*cos(phi*degree))*cm,
                         (y+(-0.5+0.5*crystalsInDetector)*sizeX1*cos((theta)*degree)*sin(phi*degree))*cm,
                         (z+(0.5-0.5*crystalsInDetector)*sizeX1*sin((theta)*degree))*cm); 
      
      pShogunMgO[ShogunDet] = new G4PVPlacement(G4Transform3D(Rot3D,Pos2),
                                                lShogunMgO,fTempname,lWorld,false,id+20000);
    }

    // Assigning the positions:
    fPosX[ShogunDet] = Pos.getX();
    fPosY[ShogunDet] = Pos.getY();
    fPosZ[ShogunDet] = Pos.getZ();
    //-----------------------------------------
    //------------------------------------------------
    cout << id << " " << fPosX[ShogunDet]/cm 
         << " " << fPosY[ShogunDet]/cm 
         << " " << fPosZ[ShogunDet]/cm << " " << endl;
    //----------------------------------------------------
    float distance = sqrt(fPosX[ShogunDet]/mm*fPosX[ShogunDet]/mm
                          + fPosY[ShogunDet]/mm*fPosY[ShogunDet]/mm	 
                          + fPosZ[ShogunDet]/mm*fPosZ[ShogunDet]/mm);

    float thetaPlaced = acos(fPosZ[ShogunDet]/mm/distance);
    float phiPlaced = acos(fPosX[ShogunDet]/mm/distance/sin(thetaPlaced));

    if(fPosY[ShogunDet]/mm < 0.0) phiPlaced = 2.0*3.14159-phiPlaced;
    thetaPlaced = thetaPlaced/3.14159*180.0;
    phiPlaced   = phiPlaced/3.14159*180.0;

    if(abs(fPosY[ShogunDet]/mm)<1.0 && fPosX[ShogunDet]/mm>=0.0) phiPlaced = 0.0;
    if(abs(fPosY[ShogunDet]/mm)<1.0 && fPosX[ShogunDet]/mm<0.0)  phiPlaced = 180.0;

    fprintf(fFileOut,"%i %f %f %f \n", ShogunDet, thetaPlaced, phiPlaced, distance);

    //________________________________________________________________________________
    // Placing the PMT:
    /*
      if(i==0)  {
      G4ThreeVector PosPMT(0.*cm,5.*cm,-7.07*cm); //Position of the PMT
      pShogunPMT[0] = new G4PVPlacement(G4Transform3D(Rot3D,PosPMT),
      lShogunPMT,fTempname,lWorld,false,id+10001);
      }
      if(i==1)  {
      G4ThreeVector PosPMT1(1.*cm,0.*cm,-7.07*cm); //Position of the PMT
      pShogunPMT[1] = new G4PVPlacement(G4Transform3D(Rot3D,PosPMT1),
      lShogunPMT,fTempname,lWorld,false,id+10001);

      G4ThreeVector PosPMT2(-1.*cm,0.*cm,-7.07*cm); //Position of the PMT
      pShogunPMT[2] = new G4PVPlacement(G4Transform3D(Rot3D,PosPMT2),
      lShogunPMT,fTempname,lWorld,false,id+10002);

      }
      if(i==2)  {
      G4ThreeVector PosLG(0.*cm,-5.*cm,-6.07*cm); //Position of the PMT
      pShogunLightGuide[0]= new G4PVPlacement(G4Transform3D(Rot3D,PosLG),
      lShogunLightGuide,fTempname,lWorld,false,id+10002);

      G4ThreeVector PosPMT(0.*cm,-5.*cm,-11.07*cm); //Position of the PMT
      pShogunPMT[3] = new G4PVPlacement(G4Transform3D(Rot3D,PosPMT),
      lShogunPMT,fTempname,lWorld,false,id+10001);
      }
    */
    //_________________________________________________________________________________
  }
  cout<<"Created the Shogun array"<<endl;
  fclose(fFileIn);
  fclose(fFileOut);
}
