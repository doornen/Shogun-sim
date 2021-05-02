#!/bin/bash
INPUTFILE="INPUTFILE ../SimulationResults/Generator/32mg/32mg32mg200mev"
FILEEXT0="mg"
FILEEXT1="ps.root"
FILEEXT2="_32mg32mg200mev"
OUTPUTFILE="OUTPUTFILE ../SimulationResults/Builder/Shogun/32mg/Case"
SHOGUNINCLUDE="SHOGUNINCLUDE 1"
SHOGUNENERGYRESOLUTION="SHOGUNENERGYRESOLUTION 2 0.771 0.5"
SHOGUNTIMERESOLUTION="SHOGUNTIMERESOLUTION 0.5 0"
SHOGUNHOUSINGTHICKNESSXYZ="SHOGUNHOUSINGTHICKNESSXYZ 1 0.1 0.1 0.1"
SHOGUNMGOTHICKNESSXYZ="SHOGUNMGOTHICKNESSXYZ 1 0.1 0.1 0.1"
POSDETECTORONTARGETRESOLUTION="POSDETECTORONTARGETRESOLUTION .3"
ENERGYDETECTORAFTERTARGETINCLUDE="ENERGYDETECTORAFTERTARGETINCLUDE 0"
POSDETECTORAFTERTARGETDISTANCE="POSDETECTORAFTERTARGETDISTANCE 100."
POSDETECTORAFTERTARGETRESOLUTION="POSDETECTORAFTERTARGETRESOLUTION .3"
BETARESOLUTION="BETARESOLUTION 0.001"
BEAMPIPEINCLUDE="BEAMPIPEINCLUDE 0"
TARGETHOLDERINCLUDE="TARGETHOLDERINCLUDE 0"
END="END"

FILE="./input/EventBuilder.in"


for arg1 in 1 500 1000 1500 2000 2500 
do
for arg2 in 0 10 20 30 40 50 60 70 80 90 100
do

echo ""$INPUTFILE$arg1$FILEEXT0$arg2$FILEEXT1"
"$OUTPUTFILE$1$FILEEXT2$arg1$FILEEXT0$arg2$FILEEXT1"
"$SHOGUNINCLUDE"
"$SHOGUNENERGYRESOLUTION"
"$SHOGUNTIMERESOLUTION"
"$SHOGUNHOUSINGTHICKNESSXYZ"
"$SHOGUNMGOTHICKNESSXYZ"
"$POSDETECTORONTARGETRESOLUTION"
"$ENERGYDETECTORAFTERTARGETINCLUDE"
"$POSDETECTORAFTERTARGETDISTANCE"
"$POSDETECTORAFTERTARGETRESOLUTION"
"$BETARESOLUTION"
"$BEAMPIPEINCLUDE"
"$TARGETHOLDERINCLUDE"
"$END" 
" > $FILE

time EventBuilder run_nothing.mac
done
done

exit 0