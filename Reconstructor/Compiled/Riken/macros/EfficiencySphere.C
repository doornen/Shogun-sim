#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "TH1F.h"
#include "TFile.h"

void EfficiencySphere()  {
    char dummytext[1000];
  char root_file[200];
  TH1F *h[200];
    
  sprintf(dummytext,"./output/FittingResults/SphereLaBr3.txt");
  FILE *fFitOutput = fopen(dummytext,"w");
  fprintf(fFitOutput,"Thickness_BeamEnergy_GammaEnergy_Counts \n");
  
  int l = 0;
   
  for(int i=8;i<=8;i+=2)  {
    sprintf(dummytext,"Sphere");
    for(int j=100;j<=200;j+=50)  {
      for(int k=250;k<=3000;k+=250)  {
        
        cout<<"Reading Root-tree "<<i<<endl;
        
        sprintf(root_file,"../../../SimulationResults/Reconstructor/Sphere/SphereLaBr3_%icm_%iMeV%ikeV.root",i,j,k);
        
        
        TFile *infile = new TFile(root_file,"READ");
        sprintf(dummytext,"h_sphere_doppler_cor;1");
        infile->GetObject(dummytext,h[l]);
                
        //h[l]->GetXaxis()->SetRangeUser((k-20),(k+20));
        //fcanvas->cd(l+1);
        //h[l]->Draw();
        
        
        double dummy = h[l]->Integral((k-20)/5,(k+20)/5);
        fprintf(fFitOutput,"%i %i %i %f \n",i,j,k,dummy);
   
        l++;
       
      }
    }
  }
  fclose(fFitOutput); 
}
