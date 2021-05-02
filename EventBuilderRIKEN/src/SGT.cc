#include "SGT.hh"

using namespace std;

SGT::SGT(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld)  {
  fMaterialList = ptMaterialList;
  lWorld = ptWorld;

  //fAbsorberThicknessPb = fAbsorberThicknessSn=fAbsorberThicknessAl = 0.0;

  for(int i=0;i<10;i++)  {
    for(int j=0;j<30;j++)  {
      fPosX[i][j] = fPosY[i][j] = fPosZ[i][j] = 0.0;
    }
  }
  //Ich weiss nicht mehr, wozu ich die naechsten beiden Objekte gebraucht habe.
  G4ThreeVector stripPos;
  G4RotationMatrix stripRot3D;
  stripRot3D.set(0, 0, 0);

  sSGTStrip = new G4Box("sSGTStrip",20.*mm/2.0,50.*mm/2.0,2.*mm/2.0);
  lSGTStrip = new G4LogicalVolume(sSGTStrip,fMaterialList->GetMaterial("Ge"),"lSGTStrip",0,0,0);
  sSGTCoaxial = new G4Tubs("sSGTCoaxial",0.0*mm,70.0*mm/2.0,70.0*mm/2.0,0*deg,360*deg);
  lSGTCoaxial = new G4LogicalVolume(sSGTCoaxial,fMaterialList->GetMaterial("Ge"),"lSGTCoaxial",0,0,0);

  // Setting the vis attributes:
  G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.2,0.5,1.0)); 
  lSGTStrip->SetVisAttributes(visAtt);
  lSGTCoaxial->SetVisAttributes(visAtt);

  cout<<"Exiting SGT Constructor"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SGT::~SGT()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SGT::ResetValues()  {
  fDetectorMult = 0;
  for(int i=0;i<10;i++)  {
    fCoaxEnergy[i]= 0.;
    fCoaxFlag[i]= 0;
    fCoaxTime[i]= 0.;
    fDetectorEnergy[i]= 0;
    fDetectorFlag[i]= 0;
    for(int j=0;j<30;j++)  {
      fPlanarEnergy[i][j] = 0.;
      fPlanarFlag[i][j] = 0.;
      fPlanarTime[i][j] = 0.0;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float SGT::GetCoaxMeasuredEnergy(int a)  {
  if(fCoaxEnergy[a]==0) return 0.;
  float dummy;
  if(fTypeOfEnergyResolution==1) dummy= (fEnergyResolution[0] + fEnergyResolution[1]*fCoaxEnergy[a]);
  else if(fTypeOfEnergyResolution==2) dummy= fEnergyResolution[0]*TMath::Power(fCoaxEnergy[a],fEnergyResolution[1]);
  else {cout<<"Wrong resolution option for SGT. Aborting program."<<endl; abort();}
  //Observed energy in the detector
  return G4RandGauss::shoot(fCoaxEnergy[a],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float SGT::GetCoaxMeasuredTime(int a)  {
  if(fCoaxTime[a]==0) return 0.;
  //Observed time in the detector
  float dummy = (fTimeResolution[0] - fTimeResolution[1]*TMath::Power(fCoaxEnergy[a],0.5));
  return G4RandGauss::shoot(fCoaxTime[a],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float SGT::GetPlanarMeasuredEnergy(int a,int b)  {
  if(fPlanarEnergy[a][b]==0) return 0.;
  float dummy;
  if(fTypeOfEnergyResolution==1) dummy = (fEnergyResolution[0] + fEnergyResolution[1]*fPlanarEnergy[a][b])/2.35;
  else if(fTypeOfEnergyResolution==2) dummy = fEnergyResolution[0]*TMath::Power(fPlanarEnergy[a][b],fEnergyResolution[1])/2.35;
  else {cout<<"Wrong resolution option for SGT. Aborting program."<<endl; abort();}
  //Observed energy in the detector
  return G4RandGauss::shoot(fPlanarEnergy[a][b],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float SGT::GetPlanarMeasuredTime(int a,int b)  {
  if(fPlanarTime[a][b]==0) return 0.;
  //Observed time in the detector
  float dummy = (fTimeResolution[0] - fTimeResolution[1]*TMath::Power(fPlanarEnergy[a][b],0.5))/2.35;
  return G4RandGauss::shoot(fPlanarTime[a][b],dummy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float SGT::GetDetectorMult()  {
  for(int i=0;i<10;i++)  {
    if(fCoaxFlag[i]==1) {fDetectorMult++; continue;}
    if(fPlanarFlag[i][0]==1)fDetectorMult++;
  }
  return fDetectorMult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SGT::CreateArrayXYZPsiThetaPhi()  {
  G4cout<<"Creating the SGT array"<<endl;
  float x,y,z,phi;
  int id;
  char sgt_name[10];
  
  //The position of the sgt is always given from the center of the 25 strip detectors
  fFileIn = fopen("./geometry/sgt_geometry_in.txt","r");
  fFileOut = fopen("./geometry/sgt_geometry_out.txt","w");

  for(int i=0;!feof(fFileIn)&&i<4;i++)  {
    fscanf(fFileIn,"%f %f %f %f",&x,&y,&z,&phi);
    id =  300000 +(i+1)*100;
    cout<<x<<" "<<y<<" "<<z<<" "<<phi<<endl;

    if(id==300100) sprintf(sgt_name,"A");
    if(id==300200) sprintf(sgt_name,"B");
    if(id==300300) sprintf(sgt_name,"C");
    if(id==300400) sprintf(sgt_name,"D");

    int sgtDet = (id-300000)/100-1;
 
    //Positioning the 25 Strips
    G4RotationMatrix Rot3D;
    G4ThreeVector Pos(x*cm,y*cm,z*cm);

    Rot3D.set(0, 0, 0);
    Rot3D.rotateZ(phi*degree);  

    G4RotationMatrix dummyMatrix;
    G4ThreeVector dummyVector;

    dummyMatrix = Rot3D;

    //Assigning the position for the center of the planar detector
    fPosX[sgtDet][1] = Pos.getX();
    fPosY[sgtDet][1] = Pos.getY();
    fPosZ[sgtDet][1] = Pos.getZ();

    // The 25 segments:
    for(int iii=0;iii<25;iii++)  {
      dummyVector = Pos+ G4ThreeVector(0,0,-24*mm+iii*2.*mm);
      sprintf(fTempname,"pSGTStrip[%i][%i]",sgtDet,iii);
      pSGTStrip[sgtDet][iii] = new G4PVPlacement(G4Transform3D(dummyMatrix,dummyVector),
                                                 lSGTStrip,fTempname,lWorld,false,(id+iii+1));

      // Assigning the positions:
      fPosX[sgtDet][iii+2]= dummyVector.getX();
      fPosY[sgtDet][iii+2]= dummyVector.getY();
      fPosZ[sgtDet][iii+2]= dummyVector.getZ();
      //------------------------------------
      cout << id+iii+1 << " " << fPosX[sgtDet][iii+2]/cm 
           << " " << fPosY[sgtDet][iii+2]/cm 
           << " " << fPosZ[sgtDet][iii+2]/cm << " " << endl;
      //---------------------------------------------------
      float dis_this = sqrt(fPosX[sgtDet][iii+2]/mm
                            *fPosX[sgtDet][iii+2]/mm
                            + fPosY[sgtDet][iii+2]/mm
                            *fPosY[sgtDet][iii+2]/mm	 
                            + fPosZ[sgtDet][iii+2]/mm
                            *fPosZ[sgtDet][iii+2]/mm);

      float thetaThis = acos(fPosZ[sgtDet][iii+2]/mm/dis_this);
      float phiThis = acos(fPosX[sgtDet][iii+2]/mm/dis_this/sin(thetaThis));

      if(fPosY[sgtDet][iii+2]/mm < 0.0) phiThis = 2.0*3.14159-phiThis;
      thetaThis = thetaThis/3.14159*180.0;
      phiThis = phiThis/3.14159*180.0;

      if(abs(fPosY[sgtDet][iii+2]/mm)<1.0 
         && fPosX[sgtDet][iii+2]/mm>=0.0) phiThis = 0.0;
      if(abs(fPosY[sgtDet][iii+2]/mm)<1.0 
         && fPosX[sgtDet][iii+2]/mm<0.0) phiThis = 180.0;

      fprintf(fFileOut,"%s %i %f %f %f \n", sgt_name, (iii+1), thetaThis, phiThis, dis_this);
    } 
    // --------------------------------------------------------------
    // The coaxial:
    dummyMatrix.set(0, 0, 0);
    dummyMatrix.rotateY(90.0*degree);
    dummyMatrix.rotateZ(phi*degree); 
    dummyVector = Pos + dummyMatrix(G4ThreeVector(0.,0.,50*mm));

    sprintf(fTempname,"pSGTCoaxial[%i]",sgtDet);
    pSGTCoaxial[sgtDet] = new G4PVPlacement(G4Transform3D(dummyMatrix,dummyVector),
                                            lSGTCoaxial,fTempname,lWorld,false,(id));
    // Assigning the positions:
    fPosX[sgtDet][0] = dummyVector.getX();
    fPosY[sgtDet][0] = dummyVector.getY();
    fPosZ[sgtDet][0] = dummyVector.getZ();
    //-----------------------------------------
    cout << id << " " << fPosX[sgtDet][0]/cm 
         << " " << fPosY[sgtDet][0]/cm 
         << " " << fPosZ[sgtDet][0]/cm << " " << endl;
    //----------------------------------------------------
    float dis_this = sqrt(fPosX[sgtDet][0]/mm
                          *fPosX[sgtDet][0]/mm
                          + fPosY[sgtDet][0]/mm 
                          *fPosY[sgtDet][0]/mm	 
                          + fPosZ[sgtDet][0]/mm
                          *fPosZ[sgtDet][0]/mm);

    float thetaThis = acos(fPosZ[sgtDet][0]/mm/dis_this);
    float phiThis = acos(fPosX[sgtDet][0]/mm/dis_this/sin(thetaThis));

    if(fPosY[sgtDet][0]/mm < 0.0) phiThis = 2.0*3.14159-phiThis;
    thetaThis = thetaThis/3.14159*180.0;
    phiThis = phiThis/3.14159*180.0;

    if(abs(fPosY[sgtDet][0]/mm)<1.0 
       && fPosX[sgtDet][0]/mm>=0.0) phiThis = 0.0;
    if(abs(fPosY[sgtDet][0]/mm)<1.0 
       && fPosX[sgtDet][0]/mm<0.0) phiThis = 180.0;

    fprintf(fFileOut,"%s 0 %f %f %f \n", sgt_name, thetaThis, phiThis, dis_this);
  }
}
