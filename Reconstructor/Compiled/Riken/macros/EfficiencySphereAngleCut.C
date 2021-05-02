#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TH1F.h"
#include "TFile.h"

void EfficiencySphereAngleCut(float energy)  {
    char dummytext[1000];
  char root_file[200];
  TH1F *h[18];
    
  sprintf(dummytext,"./output/FittingResults/SphereLaBr3AngleCut.txt");
  FILE *fFitOutput = fopen(dummytext,"w");
  fprintf(fFitOutput,"Thickness_BeamEnergy_GammaEnergy_Counts \n");
          
  sprintf(root_file,"../../../SimulationResults/Reconstructor/Sphere/SphereLaBr3_8cm_100MeV1000keV.root");
  
  TFile *infile = new TFile(root_file,"READ");
  
  for(int i = 0;i<18;i++){
    sprintf(dummytext,"h_sphere_doppler_cor_angle_cut[%i];1",i);
    infile->GetObject(dummytext,h[i]);
    double dummy = h[i]->Integral((energy-20)/5,(energy+20)/5);
    fprintf(fFitOutput,"%i %f \n",i,dummy);
    
  }
  fclose(fFitOutput); 
}
