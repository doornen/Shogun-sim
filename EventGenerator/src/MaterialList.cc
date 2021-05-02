#include "MaterialList.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

using namespace std;

MaterialList::MaterialList()  {
  // define Elements
  //
  G4double z, a;

  elH = new G4Element("Hydrogen" , "elH", z= 1., a= 1.008*g/mole);
  elBe= new G4Element("Berillium","elBe", z= 4., a= 9.012182*g/mole);
  elC = new G4Element("Carbon"   , "elC", z= 6., a= 12.01*g/mole);
  elN = new G4Element("Nitrogen" , "elN", z= 7., a= 14.01*g/mole);
  elO = new G4Element("Oxygen"   , "elO", z= 8., a= 16.00*g/mole);
  elF = new G4Element("Fluorine" , "elF", z= 9., a= 18.998*g/mole);

  elNa= new G4Element("Sodium"   ,"elNa", z=11., a= 22.98977*g/mole);
  elMg= new G4Element("Magnesium","elMg", z=12., a= 24.305*g/mole);
  elAl= new G4Element("Aluminum" ,"elAl", z=13., a= 26.982*g/mole);

  elFe= new G4Element("Iron"     ,"elFe", z=26., a= 55.845*g/mole);
  elNi= new G4Element("Nickel"   ,"elNi", z=28., a= 58.693*g/mole);

  elGa= new G4Element("Gallium"  ,"elGa", z=31., a= 69.723*g/mole);
  elBr= new G4Element("Bromine"  ,"elBr", z=35., a= 79.904*g/mole);

  elI = new G4Element("Iodine"   , "elI", z=53., a= 126.90447*g/mole);
  elCs= new G4Element("Cesium"   ,"elCs", z=55., a= 132.90545*g/mole);
  elBa= new G4Element("Barium"   ,"elBa", z=56., a= 137.327*g/mole);
  elLa= new G4Element("Lanthanum","elLa", z=57., a= 138.9055*g/mole);
  elCe= new G4Element("Cerium"   ,"elCe", z=58., a= 140.116*g/mole);
  elPr= new G4Element("Praeseodymium","elPr",z=59.,a=104.908*g/mole);

  elGd= new G4Element("Gadolinium","elGd",z=64., a= 157.25*g/mole);

  elLu = new G4Element("Lutetium" ,"elLu",z=71., a= 174.9668*g/mole);
  elW = new G4Element("Tungsten" ,"elW" , z=74., a= 183.84*g/mole);
  //
  // define materials.
  //
  G4double density, temperature, pressure;
  G4int    ncomponents, natoms;
  G4double fractionmass;
  G4String name, symbol;

  //Simple material
  // Be
  a = 9.012182*g/mole;
  density = 1.848*g/cm3;
  Be = new G4Material(name="Be", z=4., a, density);

  // C
  a = 12.011*g/mole;
  density = 1.8*g/cm3;
  C =  new G4Material(name="C",  z=6., a, density);
  
  // Al
  a = 26.98154*g/mole;
  density = 2.7*g/cm3;
  Al = new G4Material(name="Al", z=13., a, density);
 
  // Si
  a = 28.0855*g/mole;
  density = 2.33*g/cm3;
  Si = new G4Material(name="Si", z=14., a, density);

  // Fe
  a = 55.845*g/mole;
  density = 7.874*g/cm3;
  Fe = new G4Material(name="Fe", z=26., a, density);

  // Ge
  a = 72.61*g/mole;
  density = 5.323*g/cm3;
  Ge = new G4Material(name="Ge", z=32., a, density);

  // Zr
  a = 91.224*g/mole;
  density = 6.52*g/cm3;
  Ge = new G4Material(name="Zr", z=40., a, density);

  // Sn
  a = 118.71*g/mole;
  density = 7.31*g/cm3;
  Sn = new G4Material(name="Sn", z=50., a, density);

  // Au
  a = 196.97*g/mole;
  density = 19.3*g/cm3;
  Au = new G4Material(name="Au", z=79., a, density);

  // Pb
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  Pb = new G4Material(name="Pb", z=82., a, density);

  // Compound material

  // Air
  Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  // Liquid Hydrogen
  a = 1.008*g/mole;
  density = 70.99*mg/cm3;
  LH2 = new G4Material(name="LH2",z=1.,a,density);

  // Vacuum
  density     = universe_mean_density; 
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  vacuum = new G4Material("Galactic",z= 1,a= 1.008*g/mole,density,
		   kStateGas,temperature,pressure);

  // Water
  H2O = new G4Material("Water", density= 1.0*g/cm3, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);
  
  // Scintillator
  density = 1.032*g/cm3;
  Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(elC, natoms=8);
  Sci->AddElement(elH, natoms=8);
 
  // BaF2
  density = 4.89*g/cm3;
  BaF2 = new G4Material(name="BaF2", density, ncomponents=2);
  BaF2->AddElement(elBa, natoms=1);
  BaF2->AddElement(elF, natoms=2);

  // MgO
  density = 3.58*g/cm3;
  MgO = new G4Material(name="MgO", density, ncomponents=2);
  MgO->AddElement(elMg, natoms=1);
  MgO->AddElement(elO, natoms=1);

  // NaI
  density = 3.67*g/cm3;
  NaI = new G4Material(name="NaI", density, ncomponents=2);
  NaI->AddElement(elNa, natoms=1);
  NaI->AddElement(elI, natoms=1);

  // CsI
  density = 4.51*g/cm3;
  CsI = new G4Material(name="CsI", density, ncomponents=2);
  CsI->AddElement(elCs, natoms=1);
  CsI->AddElement(elI, natoms=1);

  // LaBr3
  density = 5.08*g/cm3;
  LaBr = new G4Material(name="LaBr3", density, ncomponents=2);
  LaBr->AddElement(elLa,natoms= 1);
  LaBr->AddElement(elBr,natoms= 3);

  //Delrin
  density = 0.88*g/cm3;
  Delrin = new G4Material(name="Delrin", density, ncomponents=3);
  Delrin->AddElement(elH,natoms= 16);
  Delrin->AddElement(elC,natoms= 8);
  Delrin->AddElement(elO,natoms= 8);


  //The standard for DENSIMET:
  density = 17.0*g/cm3;
  DENSIMET = new G4Material(name="DENSIMET",density,ncomponents=3);
  DENSIMET->AddElement(elFe,natoms=5);
  DENSIMET->AddElement(elNi,natoms=5);
  DENSIMET->AddElement(elW,natoms=90);

  //Polyethylene - Marlex:
  density = 0.93*g/cm3;
  CH2 = new G4Material(name="CH2",density,ncomponents=2);
  CH2->AddElement(elH,natoms=4);
  CH2->AddElement(elC,natoms=2);
  
  //cegagg
  density = 6.63*g/cm3;
  CEGAGG = new G4Material(name="CEGAGG",density,ncomponents=4);
  CEGAGG->AddElement(elGd,natoms=3);
  CEGAGG->AddElement(elAl,natoms=2);
  CEGAGG->AddElement(elGa,natoms=3);
  CEGAGG->AddElement(elO,natoms=12);

  //LUAG
  density = 6.73*g/cm3;
  LUAG = new G4Material(name="LUAG",density,ncomponents=3);
  LUAG->AddElement(elLu,natoms=3);
  LUAG->AddElement(elAl,natoms=5);
  LUAG->AddElement(elO,natoms=12);


  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MaterialList::~MaterialList()  {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material *MaterialList::GetMaterial(G4String materialName)  {
  G4Material* material = G4Material::GetMaterial(materialName); 
  return material; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MaterialList::SetCollimatorMaterial(float composition[],float density)  {
  G4String name;
  G4int    ncomponents, natoms;
  //delete DENSIMET;
  
  for(int i=0;i<3;i++)  {
    cout<<"Composition["<<i<<"]: "<<composition[i]<<endl;
  }

  DENSIMET = new G4Material(name="DENSIMET", density*g/cm3, ncomponents=3);
  DENSIMET->AddElement(elFe,natoms= composition[0]);
  DENSIMET->AddElement(elNi,natoms= composition[1]);
  DENSIMET->AddElement(elW,natoms= composition[2]);
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}
