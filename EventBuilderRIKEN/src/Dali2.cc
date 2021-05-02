#include "Dali2.hh"

using namespace std;

Dali2::Dali2(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld)   {
  fMaterialList = ptMaterialList;
  lWorld        = ptWorld;

  //fAbsorberThicknessPb = fAbsorberThicknessSn=fAbsorberThicknessAl = 0.0;
  for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
    fPosX[i] = fPosY[i] = fPosZ[i] = 0.0;
    fFIX[i] = fFIY[i] = fFIZ[i] = -999.0;
    fEnergyResolutionInd[0][i] = 9999;
    fEnergyResolutionInd[1][i] = 1;
  }

  fIndEnergyResOpt = 0;

  // The Saint-Gobain Housing:
  sDali2AlHouseOutSG = new G4Box("sDali2AlHouseOutSG",17.64*cm/2.0,8.52*cm/2.0,5.014*cm/2.0);
  sDali2AlHouseInSG  = new G4Box("sDali2AlHouseInSG",17.44*cm/2.0,8.32*cm/2.0,4.814*cm/2.0);
  sDali2AlHouseSG    = new G4SubtractionSolid("sDali2AlHouseSG",sDali2AlHouseOutSG,sDali2AlHouseInSG);
  lDali2AlHouseSG    = new G4LogicalVolume(sDali2AlHouseSG,fMaterialList->GetMaterial("Al"),"lDali2AlHouseSG",0,0,0);
  // The Scionix Housing:
  sDali2AlHouseOutSC = new G4Box("sDali2AlHouseOutSC",18.*cm/2.0,8.3*cm/2.0,4.5*cm/2.0);
  sDali2AlHouseInSC  = new G4Box("sDali2AlHouseInSC",17.8*cm/2.0,8.1*cm/2.0,4.3*cm/2.0);
  sDali2AlHouseSC    = new G4SubtractionSolid("sDali2AlHouseSC",sDali2AlHouseOutSC,sDali2AlHouseInSC);
  lDali2AlHouseSC    = new G4LogicalVolume(sDali2AlHouseSC,fMaterialList->GetMaterial("Al"),"lDali2AlHouseSC",0,0,0);
  //------------------------------
  // The Saint-Gobain Coating:
  sDali2MgOCoatOutSG = new G4Box("sDali2MgOCoatOutSG",16.28*cm/2.0, 8.28*cm/2.0,4.78*cm/2.0);
  sDali2MgOCoatInSG  = new G4Box("sDali2MgOCoatInSG",16.*cm/2.0, 8.*cm/2.0,4.5*cm/2.0);
  sDali2MgOCoatSG    = new G4SubtractionSolid("sDali2MgOCoatSG",sDali2MgOCoatOutSG,sDali2MgOCoatInSG);
  lDali2MgOCoatSG    = new G4LogicalVolume(sDali2MgOCoatSG,fMaterialList->GetMaterial("MgO"),"lDali2MgOCoatSG",0,0,0);
  // The Scionix Coating:
  sDali2MgOCoatOutSC = new G4Box("sDali2MgOCoatOutSC",16.4*cm/2.0,8.1*cm/2.0,4.3*cm/2.0);
  sDali2MgOCoatInSC  = new G4Box("sDali2MgOCoatInSC",16.*cm/2.0, 8.*cm/2.0,4.0*cm/2.0);
  sDali2MgOCoatSC    = new G4SubtractionSolid("sDali2MgOCoatSC",sDali2MgOCoatOutSC,sDali2MgOCoatInSC);
  lDali2MgOCoatSC    = new G4LogicalVolume(sDali2MgOCoatSC,fMaterialList->GetMaterial("MgO"),"lDali2MgOCoatSC",0,0,0);
  // -----------------------------
  // the Crystals:
  // The Saint-Gobain:
  sDali2CrystalSG    = new G4Box("sDali2CrystalSG",16.*cm/2.0,8.*cm/2.0,4.5*cm/2.0);
  lDali2CrystalSG    = new G4LogicalVolume(sDali2CrystalSG,fMaterialList->GetMaterial("NaI"),"lDali2CrystalSG",0,0,0);
  //lDali2CrystalSG    = new G4LogicalVolume(sDali2CrystalSG,fMaterialList->GetMaterial("CEGAGG"),"lDali2CrystalSG",0,0,0);
  //The Scionix:
  sDali2CrystalSC    = new G4Box("sDali2CrystalSC",16.*cm/2.0, 8.*cm/2.0,4.*cm/2.0);
  lDali2CrystalSC    = new G4LogicalVolume(sDali2CrystalSC,fMaterialList->GetMaterial("NaI"),"lDali2CrystalSC",0,0,0);
  //lDali2CrystalSC    = new G4LogicalVolume(sDali2CrystalSC,fMaterialList->GetMaterial("CEGAGG"),"lDali2CrystalSC",0,0,0);

  //THE PMT
  sPMT = new G4Tubs("sPMT",19.5*mm,20.*mm,75.0*mm,0*deg,360*deg);
  lPMT = new G4LogicalVolume(sPMT,fMaterialList->GetMaterial("Al"),"lPMT",0,0,0);

  //_____________________________________________
  // Including also the old dali1 crystals...
  sDali1AlHouseOut = new G4Box("sDali1AlHouseOut",12.93*cm/2.0,6.604*cm/2.0,6.604*cm/2.0);
  sDali1AlHouseIn  = new G4Box("sDali1AlHouseIn",12.83*cm/2.0,6.504*cm/2.0,6.504*cm/2.0);
  sDali1AlHouse    = new G4SubtractionSolid("sDali1AlHouse",sDali1AlHouseOut,sDali1AlHouseIn);
  lDali1AlHouse    = new G4LogicalVolume(sDali1AlHouse,fMaterialList->GetMaterial("Al"),"lDali1AlHouse",0,0,0);

  // The DALI1 Coating:
  sDali1MgOCoatOut = new G4Box("sDali1MgOCoatOut",12.8*cm/2.0,6.5*cm/2.0,6.5*cm/2.0);
  sDali1MgOCoatIn  = new G4Box("sDali1MgOCoatIn",12.2*cm/2.0, 6.1*cm/2.0,6.1*cm/2.0);
  sDali1MgOCoat    = new G4SubtractionSolid("sDali2MgOCoat",sDali1MgOCoatOut,sDali1MgOCoatIn);
  lDali1MgOCoat    = new G4LogicalVolume(sDali1MgOCoat,fMaterialList->GetMaterial("MgO"),"lDali1MgOCoat",0,0,0);
  
  sDali1Crystal    = new G4Box("sDali1Crystal",12.192*cm/2.0, 6.096*cm/2.0,6.096*cm/2.0);
  lDali1Crystal    = new G4LogicalVolume(sDali1Crystal,fMaterialList->GetMaterial("NaI"),"lDali1Crystal",0,0,0);
  //lDali1Crystal    = new G4LogicalVolume(sDali1Crystal,fMaterialList->GetMaterial("CEGAGG"),"lDali1Crystal",0,0,0);

  //_____________________________________________
  // New crystals for DALI2 upgrade, same size as old DALI1 crystals
  sDali3AlHouseOut = new G4Box("sDali3AlHouseOut",12.93*cm/2.0,6.604*cm/2.0,6.604*cm/2.0);
  sDali3AlHouseIn  = new G4Box("sDali3AlHouseIn",12.83*cm/2.0,6.504*cm/2.0,6.504*cm/2.0);
  sDali3AlHouse    = new G4SubtractionSolid("sDali3AlHouse",sDali3AlHouseOut,sDali3AlHouseIn);
  lDali3AlHouse    = new G4LogicalVolume(sDali3AlHouse,fMaterialList->GetMaterial("Al"),"lDali3AlHouse",0,0,0);

  sDali3MgOCoatOut = new G4Box("sDali3MgOCoatOut",12.8*cm/2.0,6.5*cm/2.0,6.5*cm/2.0);
  sDali3MgOCoatIn  = new G4Box("sDali3MgOCoatIn",12.2*cm/2.0, 6.1*cm/2.0,6.1*cm/2.0);
  sDali3MgOCoat    = new G4SubtractionSolid("sDali2MgOCoat",sDali1MgOCoatOut,sDali1MgOCoatIn);
  lDali3MgOCoat    = new G4LogicalVolume(sDali3MgOCoat,fMaterialList->GetMaterial("MgO"),"lDali3MgOCoat",0,0,0);

  sDali3Crystal    = new G4Box("sDali3Crystal",12.192*cm/2.0, 6.096*cm/2.0,6.096*cm/2.0);
  lDali3Crystal    = new G4LogicalVolume(sDali3Crystal,fMaterialList->GetMaterial("NaI"),"lDali3Crystal",0,0,0);
  // lDali3Crystal    = new G4LogicalVolume(sDali3Crystal,fMaterialList->GetMaterial("LUAG"),"lDali3Crystal",0,0,0);

  //_____________________________________________

  // Setting the vis attributes:
  G4VisAttributes* visAttRed   = new G4VisAttributes(G4Colour(1.0,0.0,0.0));  
  G4VisAttributes* visAttGreen = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* visAttBlue  = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* visAttHouse = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* visAttPMT   = new G4VisAttributes(G4Colour(0.0,0.0,0.0));

  lDali2MgOCoatSG->SetVisAttributes(visAttGreen);lDali2MgOCoatSC->SetVisAttributes(visAttGreen);
  lDali2AlHouseSG->SetVisAttributes(visAttHouse);lDali2AlHouseSC->SetVisAttributes(visAttHouse);

  lPMT->SetVisAttributes(visAttPMT);

  lDali1AlHouse->SetVisAttributes(visAttHouse);lDali1Crystal->SetVisAttributes(visAttGreen);
  cout<<"Exiting Dali2 Constructor"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Dali2::~Dali2()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dali2::ResetValues()  {
  fCrystalMult = 0;
  for(int i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
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
float Dali2::GetCrystalMeasuredEnergy(int a)  {
  if(fCrystalEnergy[a]==0) return 0.;
  float dummy;
  if(fIndEnergyResOpt == 0)  {
    //Here the Energy resolution is given in sigma
    if(fTypeOfEnergyResolution==1)      dummy = (fEnergyResolution[0] + fEnergyResolution[1]*fCrystalEnergy[a]);
    else if(fTypeOfEnergyResolution==2) dummy = fEnergyResolution[0]*TMath::Power(fCrystalEnergy[a],fEnergyResolution[1]);
    else {cout<<"Wrong global resolution option for DALI2. Aborting program."<<endl; abort();}
  }
  else if(fIndEnergyResOpt == 1)  {
    //Here, the energy resolution is given in sigma
    dummy = fEnergyResolutionInd[0][a]*TMath::Power(fCrystalEnergy[a],fEnergyResolutionInd[1][a]);
  }
  else { cout<<"Individual resolution option for DAlI2 set to: "<<fIndEnergyResOpt<<"."<<endl;
    cout<<"Only the values 0 and 1 are accepted. Aborting program."<<endl; abort();}
  //cout<<"dummy: "<<dummy<<endl;
    //Observed energy in the detector
  return G4RandGauss::shoot(fCrystalEnergy[a],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float Dali2::GetCrystalMeasuredTime(int a)  {
  if(fCrystalTime[a]==0) return 0.;
  //Observed time in the detector
  float dummy = (fTimeResolution[0] -fTimeResolution[1]*TMath::Power(fCrystalEnergy[a],0.5));
  return G4RandGauss::shoot(fCrystalTime[a],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dali2::CreateArrayXYZPsiThetaPhi()  {
  G4cout<<"Creating the Dali2 array"<<endl;
  float x,y,z,psi,theta,phi,rotSign;
  x = y = z = psi = theta = phi = rotSign = 0.;
  int id;
  int detType = 0;
  G4RotationMatrix Rot3D;
    
  fFileIn  = fopen("./geometry/dali2_geometry_in.txt","r");
  fFileOut = fopen("./geometry/dali2_geometry_out.txt","w");
  for(int i=0;!feof(fFileIn)&&i<NUMBEROFDALI2CRYSTALS;i++)  {
    fscanf(fFileIn,"%f %f %f %f %f %f %f %i",&x,&y,&z,&psi,&theta,&phi,&rotSign,&detType);
    G4cout<<"Read detector "<<i<<": "<<x<<" "<<y<<" "<<z<<" "<<psi<<" "<<theta<<" "<<phi<<" "<<rotSign<<" "<<detType<<endl;

    //The first digit is for the array type, digits 2-4 are used for the detector number
    //The last two digits are for the different physical materials of one detector
    id =  100000 +(i+1)*100;

    if(detType==-1) break;
    //if(i<80) detType = 1;
    //else     detType = 0;  // There are two different types of detectors for the DALI 2 array
    //if(i>=170) detType = 2;  //The old dali1 crystals
 
    int dali2Det = (id-100000)/100-1;
    z = z + fZPosShift; // Shift along beam axis because of target thickness and/or offset of the DALI2 frame.
 
    G4ThreeVector Pos(x*cm,y*cm,z*cm); //Position of the center

    Rot3D.set(0, 0, 0);
    Rot3D.rotateX(psi*degree);
    Rot3D.rotateY(theta*degree);  
    Rot3D.rotateZ(phi*degree);  

    G4ThreeVector AlHousePosition;
    G4ThreeVector  PMTPosition;

    //The detector type (SG==0,SC==1):
   if(detType==0)  {
      AlHousePosition = Pos + Rot3D(G4ThreeVector(rotSign*0.535*cm,0.,0)); //cout<<"rotSign: "<<rotSign<<endl;
      sprintf(fTempname,"pDali2AlHouse[%i]",dali2Det);
      pDali2AlHouse[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,AlHousePosition),lDali2AlHouseSG,fTempname,lWorld,false,id+2);
      
      sprintf(fTempname,"pDali2MgOCoat[%i]",dali2Det);
      pDali2MgOCoat[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),lDali2MgOCoatSG,fTempname,lWorld,false,id+1);
      
      sprintf(fTempname,"pDali2Crystal[%i]",dali2Det);
      pDali2Crystal[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),lDali2CrystalSG,fTempname,lWorld,false,id);
/*
      sprintf(fTempname,"pPMT[%i]",dali2Det);
      PMTPosition =  AlHousePosition + Rot3D(G4ThreeVector(15.*cm,0.,0));
      Rot3D.set(0, 0, 0);
      Rot3D.rotateY(90.*degree);  
      Rot3D.rotateZ(phi*degree);
      pPMT[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,PMTPosition),lPMT,fTempname,lWorld,false,id+3);
*/

    }
    else if(detType==1)  {
      AlHousePosition = Pos + Rot3D(G4ThreeVector(rotSign*0.7*cm,0.,0)); 
      sprintf(fTempname,"pDali2AlHouse[%i]",dali2Det);
      pDali2AlHouse[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,AlHousePosition),
                           lDali2AlHouseSC,fTempname,lWorld,false,id+2);
      sprintf(fTempname,"pDali2MgOCoat[%i]",dali2Det);
      pDali2MgOCoat[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali2MgOCoatSC,fTempname,lWorld,false,id+1);
      sprintf(fTempname,"pDali2Crystal[%i]",dali2Det);
      pDali2Crystal[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali2CrystalSC,fTempname,lWorld,false,id);
    }
    else if(detType==2)  {
      sprintf(fTempname,"pDali2AlHouse[%i]",dali2Det);
      pDali2AlHouse[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali1AlHouse,fTempname,lWorld,false,id+2);
      sprintf(fTempname,"pDali2MgOCoat[%i]",dali2Det);
      pDali2MgOCoat[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali1MgOCoat,fTempname,lWorld,false,id+1);
      sprintf(fTempname,"pDali2Crystal[%i]",dali2Det);
      pDali2Crystal[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali1Crystal,fTempname,lWorld,false,id);
    }
    else if(detType==3)  {
      sprintf(fTempname,"pDali2AlHouse[%i]",dali2Det);
      pDali2AlHouse[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali3AlHouse,fTempname,lWorld,false,id+2);
      sprintf(fTempname,"pDali2MgOCoat[%i]",dali2Det);
      pDali2MgOCoat[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali3MgOCoat,fTempname,lWorld,false,id+1);
      sprintf(fTempname,"pDali2Crystal[%i]",dali2Det);
      pDali2Crystal[dali2Det] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                           lDali3Crystal,fTempname,lWorld,false,id);
    }

    else {cout<<"You chose the wrong detector type"<<endl;return;}
           
    // Assigning the positions:
    fPosX[dali2Det] = Pos.getX();
    fPosY[dali2Det] = Pos.getY();
    fPosZ[dali2Det] = Pos.getZ();
    //-----------------------------------------
    //------------------------------------------------
    cout << id << " " << fPosX[dali2Det]/cm 
               << " " << fPosY[dali2Det]/cm 
               << " " << fPosZ[dali2Det]/cm << " " << endl;
    //----------------------------------------------------

    float distance = sqrt(fPosX[dali2Det]/mm*fPosX[dali2Det]/mm  +
                          fPosY[dali2Det]/mm*fPosY[dali2Det]/mm	 +
                          fPosZ[dali2Det]/mm*fPosZ[dali2Det]/mm);

    float thetaPlaced = acos(fPosZ[dali2Det]/mm/distance);
    float phiPlaced   = acos(fPosX[dali2Det]/mm/distance/sin(thetaPlaced));

    if(fPosY[dali2Det]/mm < 0.0) phiPlaced = 2.0*3.14159-phiPlaced;
    thetaPlaced = thetaPlaced/3.14159*180.0;
    phiPlaced   = phiPlaced/3.14159*180.0;

    if(abs(fPosY[dali2Det]/mm)<1.0 && fPosX[dali2Det]/mm>=0.0) phiPlaced = 0.0;
    if(abs(fPosY[dali2Det]/mm)<1.0 && fPosX[dali2Det]/mm<0.0)  phiPlaced = 180.0;

    fprintf(fFileOut,"%i %i %f %f %f \n", dali2Det, detType, thetaPlaced, phiPlaced, distance);
  }
  fclose(fFileIn);
  fclose(fFileOut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dali2::SetEnergyResolutionInd(int a)  {
  fIndEnergyResOpt = a;
  if(fIndEnergyResOpt!=1) return;
  cout<<"Setting the Energy resolution."<<endl;
  fFileIn  = fopen("./input/Dali2Resolution.txt","r");
  int detector;
  float FWHM = 1.;
  char resolution[100];
  int resolutionType;
  cout<<"Opened file."<<endl;
  fscanf(fFileIn,"%s %i",&resolution,&resolutionType); //Check if the given resolution is in SIGMA or FWHM
  if(resolutionType==1) FWHM=2.35;
  for(int i=0;!feof(fFileIn)&&i<NUMBEROFDALI2CRYSTALS;i++)  {
    fscanf(fFileIn,"%i",&detector);
    fscanf(fFileIn,"%f %f",&fEnergyResolutionInd[0][detector-1],&fEnergyResolutionInd[1][detector-1]);
    fEnergyResolutionInd[0][detector-1] = fEnergyResolutionInd[0][detector-1]/FWHM;
    G4cout<<"Energy resolution detector "<<detector<<": "<<fEnergyResolutionInd[0][detector-1]<<" "<<fEnergyResolutionInd[1][detector-1]<<endl;
  }
  fclose(fFileIn);
}

