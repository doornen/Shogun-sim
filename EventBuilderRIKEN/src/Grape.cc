#include "Grape.hh"
#include <string>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Grape::Grape(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld) {
  cout<<"Grape Constructor"<<endl;
  fMaterialList = ptMaterialList;
  lWorld = ptWorld;
  //fAbsorberThicknessPb = fAbsorberThicknessSn = fAbsorberThicknessAl = 0.0;

  for(int i=0;i<20;i++)  {
    for(int j=0;j<2;j++)  {
      for(int k=0;k<10;k++)  {
        fPosX[i][j][k] = fPosY[i][j][k] = fPosZ[i][j][k] = 0.0;
      }
    }
  }
  // So far only the crystal without any housing
  G4ThreeVector segmentPos;
  G4RotationMatrix segmentRot3D;
  segmentRot3D.set(0, 0, 0);
 
  //Starting with a tube:
  sGrapeCrystalTube = new G4Tubs("sGrapeCrystalTube",0.0*mm,60.0*mm/2.0,20.0*mm/2.0,0*deg, 360*deg);
  //Making the boxes, the intersections with the tube will correspond to the segments.
  sGrapeCrystalBox = new G4Box("sGrapeCrystalBox",20*mm/2.0,20*mm/2.0,20*mm/2.0);
  //Making the nine crystals:
  for(int i=0;i<9;i++)  {
    sprintf(fTempname,"sGrapeSegment[%i]",i);
    G4double xPos = ((i%3)-1)*20.0*mm;
    G4double yPos = ((int)(i/3)-1)*20.0*mm;
    
    segmentPos = G4ThreeVector(xPos, yPos, 0.0);
    sGrapeSegment[i] = new G4IntersectionSolid(fTempname,sGrapeCrystalTube,sGrapeCrystalBox,&segmentRot3D,segmentPos);
    // Making a logic volume out og the nine segments:
    sprintf(fTempname,"lGrapeSegment[%i]",i);
    lGrapeSegment[i] = new G4LogicalVolume(sGrapeSegment[i],fMaterialList->GetMaterial("Ge"),fTempname,0,0,0);
  }
  //_______________________________________________________________________________________________
  //------------------------------
  //Grape Dewar  
  sGrapeDewar = new G4Tubs("sGrapeDewar",0.0*mm,210.0*mm/2.0,280.0*mm/2.0,0.0*deg,360.0*deg);
  lGrapeDewar = new G4LogicalVolume(sGrapeDewar,fMaterialList->GetMaterial("Air"),"lGrapeDewar",0,0,0);
  //------------------------------
  //Grape AlBar  
  sGrapeAlBar = new G4Tubs("sGrapeAlBar",0.0*mm,35.0*mm/2.0,130.0*mm/2.0,0.0*deg,360.0*deg);
  lGrapeAlBar = new G4LogicalVolume(sGrapeAlBar,fMaterialList->GetMaterial("Air"),"lGrapeAlBar",0,0,0);

  G4VisAttributes* visAttRed   = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); 
  G4VisAttributes* visAttGreen = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* visAttBlue  = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* visAttDewar = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  G4VisAttributes* visAttAl    = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  G4VisAttributes* visAttGe    = new G4VisAttributes(G4Colour(1.0,0.5,0.6));
  
  int i3=0;
  for(i3=0;i3<3;i3++)  {
    lGrapeSegment[i3]->SetVisAttributes(visAttGe);
    lGrapeSegment[i3+3]->SetVisAttributes(visAttGe);
    lGrapeSegment[i3+6]->SetVisAttributes(visAttGe);
  }
  lGrapeDewar->SetVisAttributes(visAttDewar);
  lGrapeAlBar->SetVisAttributes(visAttAl);

  cout<<"Exiting Grape Constructor"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Grape::~Grape()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Grape::ResetValues()  {
  int i, j, k;
  fCrystalMult = 0;
  for(i=0;i<20;i++)  {
    for(j=0;j<2;j++)  {
      for(k=0;k<10;k++)  {
        fCrystalFlag[i][j][k] = false;
        fCrystalEnergy[i][j][k] = 0.0;
        fCrystalTime[i][j][k] = 0.0;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float Grape::GetCrystalMeasuredEnergy(int a, int b, int c)  {
  float dummy;
  if(fTypeOfEnergyResolution==1) dummy = (fEnergyResolution[0] + fEnergyResolution[1]*fCrystalEnergy[a][b][c]);
  else if(fTypeOfEnergyResolution==2) dummy = fEnergyResolution[0]*TMath::Power(fCrystalEnergy[a][b][c],fEnergyResolution[1]);
  else {cout<<"Wrong resolution option for GRAPE. Aborting program."<<endl; abort();}
  //Observed energy in the detector
  return G4RandGauss::shoot(fCrystalEnergy[a][b][c],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float Grape::GetCrystalMeasuredTime(int a,int b, int c)  {
  if(fCrystalTime[a][b][c]==0) return 0.;
  //Observed time in the detector
  float dummy = (fTimeResolution[0] -fTimeResolution[1]*TMath::Power(fCrystalEnergy[a][b][c],0.5));
  return G4RandGauss::shoot(fCrystalTime[a][b][c],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Grape::CreateArrayXYZPsiThetaPhi()  {
  G4cout<<"Creating the GRAPE array"<<endl;
  float x,y,z,psi,theta,phi;
  int id;
  string nameList = "ABCDEFGHIJKLMNOPQR";

  fFileIn = fopen("./geometry/grape_geometry_in.txt","r");
  fFileOut = fopen("./geometry/grape_geometry_out.txt","w");

  int i = 0;
  while (!feof(fFileIn)&&i<NUMBEROFGRAPEDETECTORS)  {//Maximum 18 Detectors  
    fscanf(fFileIn,"%f %f %f %f %f %f",&x,&y,&z,&psi,&theta,&phi); 
    
    //The first digit is for the array type, digits 2-4 are used for the detector number
    //The last two digits are for the different physical materials of one detector
    id =  200000 +(i+1)*100;    
    //The angle values for phi and theta are 90 degrees less than "expected" for a right-handed system
    //due to the intrinsic orientation of the detector!
    // OK the values are corrected in the Grape creator now!
    // Also the theta angle is defined as 0 for antiparallel to the beam direction (crystal looking upstream) 
    //and 180 parallel to the beam direction (crystal looking downstream)

    int grapeDet = (id-200000)/100-1;

    G4RotationMatrix Rot3D;
    G4ThreeVector centerPos(x*cm,y*cm,z*cm); //Position of the center

    G4ThreeVector crystalPos;  //Position of the cylindars centers 

    Rot3D.set(0, 0, 0);
    //Rot3D.rotateY(90.*degree); 
    Rot3D.rotateX(psi*degree);
    Rot3D.rotateY(theta*degree+90*degree);  
    Rot3D.rotateZ(phi*degree);  

    // The two crystals:
    for(int iii=0;iii<2;iii++)  {
      // All the segments:
      for(int j=0;j<9;j++)  { 
        float  zShift = 0.;

        if(iii==0) zShift = -1.0*cm;
        if(iii==1) zShift = +1.0*cm;

        crystalPos = centerPos + Rot3D(G4ThreeVector(0.,0.,zShift));

        sprintf(fTempname,"pGrapeSegment[%i][%i][%i]",grapeDet,iii,j);
        pGrapeSegment[grapeDet][iii][j] = new G4PVPlacement(G4Transform3D(Rot3D,crystalPos),
                                           lGrapeSegment[j],fTempname,lWorld,false,(id+10*iii+j));

        //cout<<"id_start+10*iii+j = "<< (id_start+10*iii+j)<<endl;
     
        G4double xPos = ((j%3)-1)*20.0*mm;
        G4double yPos = ((int)(j/3)-1)*20.0*mm;

        G4RotationMatrix segmentRotation;
        G4ThreeVector segmentPos;

        segmentPos = crystalPos + Rot3D(G4ThreeVector(xPos, yPos, 0)); 
      
        // Assigning the positions:
        fPosX[grapeDet][iii][j] = segmentPos.getX();
        fPosY[grapeDet][iii][j] = segmentPos.getY();
        fPosZ[grapeDet][iii][j] = segmentPos.getZ();
        //-----------------------------------------
        //cout------------------------------------------------
        cout << id+10*iii+j << " " << fPosX[grapeDet][iii][j]/cm 
                            << " " << fPosY[grapeDet][iii][j]/cm 
                            << " " << fPosZ[grapeDet][iii][j]/cm << " " << endl;
        //----------------------------------------------------

        float dis_this = sqrt(fPosX[grapeDet][iii][j]/mm
                             *fPosX[grapeDet][iii][j]/mm
                            + fPosY[grapeDet][iii][j]/mm
                             *fPosY[grapeDet][iii][j]/mm	 
                            + fPosZ[grapeDet][iii][j]/mm
                             *fPosZ[grapeDet][iii][j]/mm);

        float thetaThis = acos(fPosZ[grapeDet][iii][j]/mm/dis_this);
        float phiThis = acos(fPosX[grapeDet][iii][j]/mm/dis_this/sin(thetaThis));

        if(fPosY[grapeDet][iii][j]/mm < 0.0) phiThis = 2.0*3.14159-phiThis;
        thetaThis = thetaThis/3.14159*180.0;
        phiThis = phiThis/3.14159*180.0;
 
        if(abs(fPosY[grapeDet][iii][j]/mm)<1.0 && fPosX[grapeDet][iii][j]/mm>=0.0) phiThis = 0.0;
        if(abs(fPosY[grapeDet][iii][j]/mm)<1.0 && fPosX[grapeDet][iii][j]/mm<0.0) phiThis = 180.0;

        fprintf(fFileOut,"%i %i %i %f %f %f \n",i , iii, j, thetaThis, phiThis, dis_this);
      }
    }
    //Placing the dewar:

    G4ThreeVector dewarPos;
    G4ThreeVector AlBarPos;
    G4RotationMatrix dummyMatrix;
    dummyMatrix.set(0, 0, 0); 
    dummyMatrix.rotateX(psi*degree);
    dummyMatrix.rotateY(theta*degree);  
    dummyMatrix.rotateZ(phi*degree); 
    dewarPos = centerPos + dummyMatrix(G4ThreeVector(0.,0.,30*cm));
    AlBarPos = centerPos + dummyMatrix(G4ThreeVector(0.,0.,9.5*cm));

    sprintf(fTempname,"pGrapeAlBar[%i]",grapeDet);
    pGrapeAlBar[grapeDet] = new G4PVPlacement(G4Transform3D(dummyMatrix,AlBarPos),
                           lGrapeAlBar,fTempname,lWorld,false,(id+31));

    sprintf(fTempname,"pGrapeDewar[%i]",grapeDet);
    pGrapeDewar[grapeDet] = new G4PVPlacement(G4Transform3D(dummyMatrix,dewarPos),
                           lGrapeDewar,fTempname,lWorld,false,(id+30));
    i++;
  }  
  fclose(fFileOut);
  fclose(fFileIn);
}
