#include "LaBr3Array.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LaBr3Array::LaBr3Array(MaterialList *ptMaterialList,G4LogicalVolume *ptWorld)   {
  cout<<"LaBr3Array Constructor"<<endl;

  fMaterialList        = ptMaterialList;
  lWorld               = ptWorld;
  //fAbsorberThicknessPb = fAbsorberThicknessSn = fAbsorberThicknessAl = 0.0;
  fNumberOfCrystals    = 0;
  
  fHousing = false;
  fHousingThicknessFront  = 0.;
  fHousingThicknessSide   = 0.;
  
  fInsulation = false;
  fInsulationThicknessFront = 0.;
  fInsulationThicknessSide = 0.;
  
  fTransparentHousing = false;
  fTransparentInsulation = false;

  for(int i=0;i<NUMBEROFLABR3ARRAYCRYSTALS;i++)  {
    fPosX[i] = fPosY[i] = fPosZ[i] = 0.0;
  }
 
  //CreateArrayXYZPsiThetaPhi();
  cout<<"Exiting LaBr3Array Constructor"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LaBr3Array::~LaBr3Array()  {
}


void LaBr3Array::SetHousing(bool doHousing, float X, float Y, bool transparent){
  fHousingThicknessFront = X; fHousingThicknessSide = Y;
  fTransparentHousing = transparent;
  G4cout<<"The thickness of the LaBr3Array housing is: "<<fHousingThicknessFront<<" "<<fHousingThicknessSide<<G4endl;
  if(fHousingThicknessFront<0 || fHousingThicknessSide<0)  {
    G4cout<<"At least one of the values is smaller than zero. HOUSING WILL BE OMITTED!!!!! "<<G4endl;
    fHousing=false;
  }
  else{fHousing=doHousing;}
  return;
}

void LaBr3Array::SetInsulation(bool doInsulation, float X, float Y, bool transparent){
  fInsulationThicknessFront = X; fInsulationThicknessSide = Y;
  fTransparentInsulation = transparent;
  G4cout<<"The thickness of the LaBr3Array insulation is: "<<fInsulationThicknessFront<<" "<<fInsulationThicknessSide<<G4endl;
  if(fInsulationThicknessFront<0 || fInsulationThicknessSide<0)  {
    G4cout<<"At least one of the values is smaller than zero. INSULATION WILL BE OMITTED!!!!! "<<G4endl;
    fInsulation = false;
  }
  else{fInsulation = doInsulation;}
  return;
} 








//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Array::ResetValues()  {
  fCrystalMult = 0;
  for(int i=0;i<fNumberOfCrystals;i++)  {
    fCrystalFlag[i]   = false;
    fCrystalEnergy[i] = 0.0;
    fCrystalTime[i]   = 0.0;
  }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float LaBr3Array::GetCrystalMeasuredEnergy(int a)  {
  if(fCrystalEnergy[a]==0) return 0.;
  if(fTypeOfEnergyResolution==0){return fCrystalEnergy[a];} //perfect resolution
  float dummy;
  if(fTypeOfEnergyResolution==1) dummy = (fEnergyResolution[0] + fEnergyResolution[1]*fCrystalEnergy[a]);
  else if(fTypeOfEnergyResolution==2) dummy = fEnergyResolution[0]*TMath::Power(fCrystalEnergy[a],fEnergyResolution[1]);
  else {
    if(fTypeOfEnergyResolution==3){
      dummy = TMath::Sqrt( fEnergyResolution[0] + fEnergyResolution[1]*fCrystalEnergy[a] + fEnergyResolution[2]*TMath::Power(fCrystalEnergy[a],2) );
      dummy/=2.35; //get sigma from fwhm
    }else{
    cout<<"Wrong resolution option '" << fTypeOfEnergyResolution << "' for LaBr3Array. Aborting program."<<endl; 
    abort();
    }
  }
  //Observed energy in the detector
  return G4RandGauss::shoot(fCrystalEnergy[a],dummy);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float LaBr3Array::GetCrystalMeasuredTime(int a)  {
  if(fCrystalTime[a]==0) return 0.;
  //Observed time in the detector
  float dummy = (fTimeResolution[0] - fTimeResolution[1]*TMath::Power(fCrystalEnergy[a],0.5));
  return G4RandGauss::shoot(fCrystalTime[a],dummy);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Array::CreateArrayXYZPsiThetaPhi()  {
  G4cout<<"Creating the LaBr3Array array"<<endl;
  float x,y,z,psi,theta,phi,crystalDiameter,crystalLength;
  x = y = z = psi = theta = phi = crystalDiameter = crystalLength = 0.;
  int id;
  G4RotationMatrix Rot3D, Rot3D2;
  
  fFileIn  = fopen("./geometry/LaBr3Array_geometry_in.txt","r");
  fFileOut = fopen("./geometry/LaBr3Array_geometry_out.txt","w");
  for(int i=0;!feof(fFileIn)&& i<NUMBEROFLABR3ARRAYCRYSTALS;i++)  {
    fscanf(fFileIn,"%i %f %f %f %f %f %f %f %f %s",&id,&x,&y,&z,&psi,&theta,&phi,&crystalDiameter,&crystalLength,fInsulationMaterialName);
    G4cout<<"Read detector "<<i<<": "<<x<<" "<<y<<" "<<z<<" "<<psi<<" "<<theta<<" "<<phi<<endl;
    if(id==-1) break;
    if(id<1) continue;

    sLaBr3ArrayCrystal = new G4Tubs("LaBr3ArrayCrystal",0.0,crystalDiameter*cm/2.0,crystalLength*cm/2.0,0.0*deg,360.0*deg);
    lLaBr3ArrayCrystal = new G4LogicalVolume(sLaBr3ArrayCrystal,fMaterialList->GetMaterial("LaBr3"),"lLaBr3ArrayCrystal",0,0,0);
    
    G4VisAttributes* visAtt = new G4VisAttributes(G4Colour(0.5,0.5,1.0)); //blue
    lLaBr3ArrayCrystal->SetVisAttributes(visAtt);

    //The Housing
    
    if(fHousing)  {
      
      //pschrock: there is excess material at the front edge, 5 mm longer than front cap
      float excessMaterialAtFront=0.5;
      
      float flangeLength=2.0, pmtLength=19.0; //pschrock: these are full-lengths in cm //flange: 2 peaces, each 1 cm thick
      
      
      ///housing of LaBr crystal
      sLaBr3ArrayAlHouseOut = new G4Tubs("sLaBr3ArrayAlHouseOut",
                                         0.0,(crystalDiameter+2*fInsulationThicknessSide+2*fHousingThicknessSide)*cm/2.0,
//                                         (crystalLength+2*fInsulationThicknessFront+2*fHousingThicknessFront)*cm/2.0,
                                         (crystalLength+2*fInsulationThicknessFront+2*fHousingThicknessFront+2.0*excessMaterialAtFront)*cm/2.0, 
                                         0*deg,360*deg);
//      sLaBr3ArrayAlHouseIn = new G4Tubs("sLaBr3ArrayAlHouseIn",
//                                        0.0,(crystalDiameter+2*fInsulationThicknessSide)*cm/2.0,
//                                        (crystalLength+2*fInsulationThicknessFront)*cm/2.0,
//                                        0*deg,360*deg);
      //pschrock
      sLaBr3ArrayAlHouseIn = new G4Tubs("sLaBr3ArrayAlHouseIn",
                                        0.0,(crystalDiameter+2*fInsulationThicknessSide)*cm/2.0,
//                                        (crystalLength+2*fInsulationThicknessFront+2*fHousingThicknessFront+0.1)*cm/2.0,
//                                         (crystalLength+2*fInsulationThicknessFront+2*fHousingThicknessFront+2.0*excessMaterialAtFront+0.1)*cm/2.0, //cut a bit more
                                        (crystalLength+2*fInsulationThicknessFront+2*fHousingThicknessFront+2.0*excessMaterialAtFront+pmtLength)*cm/2.0, //
                                        0*deg,360*deg);
      cout << "LaBr crystal Diameter is " << crystalDiameter << ", lenght is " << crystalLength << endl;
      cout << "LaBr outer housing radius is: " << (crystalDiameter+2*fInsulationThicknessSide+2*fHousingThicknessSide)*cm/2.0 << " mm (should be 60.0 mm for Milano LaBr)" << endl;
      
       
      
      
      ///create back side
      //simplified structures behind LaBr
      //THIS IS NOT THE EXACT GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //It's just for optical reasons
      //please see technical drawings for correct values
      
      sLaBr3ArrayAlHouseBackFlange = new G4Tubs("sLaBr3ArrayAlHouseBackFlange",
                                        0.0,19.028*cm/2.0,flangeLength*cm/2.0, 
                                        0*deg,360*deg);
      sLaBr3ArrayAlHouseBackPMT = new G4Tubs("sLaBr3ArrayAlHouseBackPMT",
                                        0.0,12.0*cm/2.0,pmtLength*cm/2.0,
                                        0*deg,360*deg);
      
      
      
      ///Create front cap
      //should be 5 mm Al
      sLaBr3ArrayAlHouseFront = new G4Tubs("sLaBr3ArrayAlHouseFront",
                                          0.0,(crystalDiameter+2*fInsulationThicknessSide)*cm/2.0,fHousingThicknessFront*cm/2.0, //half thickness for creation
                                          0*deg,360*deg);
      
      
      
      ///merge the parts
      
      //create the back side:
      //pmt and flange
      //new center is center of flange
      sLaBr3ArrayAlHouseBack = new G4UnionSolid("sLaBr3ArrayAlHouseBack",sLaBr3ArrayAlHouseBackFlange,sLaBr3ArrayAlHouseBackPMT,0,G4ThreeVector(0.0*cm,0.0*cm,(flangeLength+pmtLength)*cm/2.0));
      
      
      
      //merge outer house and backside
      //new center is center of outer house
      sLaBr3ArrayAlHouseOutBack = new G4UnionSolid("sLaBr3ArrayAlHouseOutBack",sLaBr3ArrayAlHouseOut,sLaBr3ArrayAlHouseBack,0,
      						G4ThreeVector(0.0*cm,0.0*cm,(crystalLength+flangeLength)*cm/2.0));
      
      //substract inner part of house
      sLaBr3ArrayAlHouseOutBackIn   = new G4SubtractionSolid("sLaBr3ArrayAlHouseOutBackIn",sLaBr3ArrayAlHouseOutBack,sLaBr3ArrayAlHouseIn);
      
      
      //add front cap
      sLaBr3ArrayAlHouse = new G4UnionSolid("sLaBr3ArrayAlHouse",sLaBr3ArrayAlHouseOutBackIn,sLaBr3ArrayAlHouseFront,0,
      						G4ThreeVector(0.0*cm,0.0*cm,-(crystalLength+2.0*fInsulationThicknessFront+fHousingThicknessFront)*cm/2.0));
      
      //housing of LaBr crystal
//       sLaBr3ArrayAlHouse   = new G4SubtractionSolid("sLaBr3ArrayAlHouse",sLaBr3ArrayAlHouseOut,sLaBr3ArrayAlHouseIn);
      lLaBr3ArrayAlHouse   = new G4LogicalVolume(sLaBr3ArrayAlHouse,fMaterialList->GetMaterial("Al"),"lLaBr3ArrayAlHouse",0,0,0);
      
      
      
      
      /// Setting the vis attributes:
//      G4VisAttributes* visAttHouse = new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //gray
      G4VisAttributes* visAttHouse = new G4VisAttributes(G4Colour(1.0,0.3,0.3)); //red
      
      if(fTransparentHousing) {visAttHouse->SetForceWireframe(true);} //set housing transparent
      lLaBr3ArrayAlHouse->SetVisAttributes(visAttHouse);
      
      
    } //end of housing
    
    
    
    
    // The insulation:
//    if(fInsulationThicknessFront>0. && fInsulationThicknessSide>0.)  {
    if(fInsulation)  {
      sLaBr3ArrayInsulationOut = new G4Tubs("sLaBr3ArrayInsulationOut",
                                            0.0,(crystalDiameter+2*fInsulationThicknessSide)*cm/2.0,
                                            (crystalLength+2*fInsulationThicknessFront)*cm/2.0,
                                            0*deg,360*deg);

      sLaBr3ArrayInsulationIn = new G4Tubs("sLaBr3ArrayInsulationIn",
                                           0.0,crystalDiameter*cm/2.0,
                                           crystalLength*cm/2.0,
                                           0*deg,360*deg);

      sLaBr3ArrayInsulation   = new G4SubtractionSolid("sLaBr3ArrayInsulation",sLaBr3ArrayInsulationOut,sLaBr3ArrayInsulationIn);
      lLaBr3ArrayInsulation   = new G4LogicalVolume(sLaBr3ArrayInsulation,fMaterialList->GetMaterial(fInsulationMaterialName),
                                                    "InsulationMaterial",0,0,0);
      cout << "	Using '" << fInsulationMaterialName << "' as insulation material. " << endl;
      
      // Setting the vis attributes:
      G4VisAttributes* visAttInsulation = new G4VisAttributes(G4Colour(1.0,1.0,.0)); //yellow
      
      if(fTransparentInsulation) {visAttInsulation->SetForceWireframe(true);} //set insulation transparent
      lLaBr3ArrayInsulation->SetVisAttributes(visAttInsulation);
    } //insulation

    fNumberOfCrystals++;

    //The first digit is for the array type, the second for the different materials inserted
    //The last digit is for the different detector numbers
    id       =  600000 + id -1;
    
    //z posiiton shift along the beam axis like Dali2
    z = z + fZPosShift;
    
    G4ThreeVector Pos(x*cm,y*cm,z*cm); //Position of the center
    
    Rot3D.set(0, 0, 0);
    Rot3D.rotateX(psi*degree);
    Rot3D.rotateY(theta*degree);  
    Rot3D.rotateZ(phi*degree);  
    


    sprintf(fTempname,"pLaBr3ArrayCrystal[%i]",i);

    pLaBr3ArrayCrystal[i] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                                              lLaBr3ArrayCrystal,fTempname,lWorld,false,id);

    // Putting the housing
//    if(fHousingThicknessSide>0. && fHousingThicknessFront>0.)  {
    if(fHousing)  {
      sprintf(fTempname,"pLaBr3ArrayAlHouse[%i]",i);
      pLaBr3ArrayAlHouse[i] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                                                lLaBr3ArrayAlHouse,fTempname,lWorld,false,id+10000);
      
    }
    
    // Putting the insulation:
//    if(fInsulationThicknessSide>0. && fInsulationThicknessFront>0.)  {
    if(fInsulation)  {
      sprintf(fTempname,"pLaBr3ArrayInsulation[%i]",i);
      
      pLaBr3ArrayInsulation[i] = new G4PVPlacement(G4Transform3D(Rot3D,Pos),
                                                   lLaBr3ArrayInsulation,fTempname,lWorld,false,id+20000);
    }

    // Assigning the positions:
    fPosX[i] = Pos.getX();
    fPosY[i] = Pos.getY();
    fPosZ[i] = Pos.getZ();
    //-----------------------------------------
    //------------------------------------------------
    cout << id << " " << fPosX[i]/cm 
         << " " << fPosY[i]/cm 
         << " " << fPosZ[i]/cm << " " << endl;
    //----------------------------------------------------
    float distance = sqrt(fPosX[i]/mm*fPosX[i]/mm
                          + fPosY[i]/mm*fPosY[i]/mm	 
                          + fPosZ[i]/mm*fPosZ[i]/mm);

    float thetaPlaced = acos(fPosZ[i]/mm/distance);
    float phiPlaced = acos(fPosX[i]/mm/distance/sin(thetaPlaced));

    if(fPosY[i]/mm < 0.0) phiPlaced = 2.0*3.14159-phiPlaced;
    thetaPlaced = thetaPlaced/3.14159*180.0;
    phiPlaced   = phiPlaced/3.14159*180.0;

    if(abs(fPosY[i]/mm)<1.0 && fPosX[i]/mm>=0.0) phiPlaced = 0.0;
    if(abs(fPosY[i]/mm)<1.0 && fPosX[i]/mm<0.0)  phiPlaced = 180.0;

    fprintf(fFileOut,"%i %f %f %f \n", i, thetaPlaced, phiPlaced, distance);
  } // end of loop over LaBr crystals
  fclose(fFileIn);
  fclose(fFileOut);
} // end of CreateArrayXYZPsiThetaPhi



