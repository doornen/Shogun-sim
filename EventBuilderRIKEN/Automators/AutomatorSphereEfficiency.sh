#!/bin/bash
INPUTFILE="INPUTFILE ../SimulationResults/Generator/Efficiency/"
FILEEXT0="cm_"
FILEEXT1="MeV"
FILEEXT2="keV.root"
OUTPUTFILE="OUTPUTFILE ../SimulationResults/Builder/Sphere/SphereLaBr3_"
SHOGUNINCLUDE="SPHEREINCLUDE 1 18 10"
POSDETECTORONTARGETRESOLUTION="POSDETECTORONTARGETRESOLUTION .3"
ENERGYDETECTORAFTERTARGETINCLUDE="ENERGYDETECTORAFTERTARGETINCLUDE 0"
POSDETECTORAFTERTARGETDISTANCE="POSDETECTORAFTERTARGETDISTANCE 100."
POSDETECTORAFTERTARGETRESOLUTION="POSDETECTORAFTERTARGETRESOLUTION .3"
BETARESOLUTION="BETARESOLUTION 0.001"
BEAMPIPEINCLUDE="BEAMPIPEINCLUDE 0"
TARGETHOLDERINCLUDE="TARGETHOLDERINCLUDE 0"
END="END"

FILE="./input/EventBuilder.in"

for arg1 in  100 150 200
do
for arg2 in 250 500 750 1000 1250 1500 1750 2000 2250 2500 2750 3000
do

echo ""$INPUTFILE$arg1$FILEEXT1$arg2$FILEEXT2"
"$OUTPUTFILE$1$FILEEXT0$arg1$FILEEXT1$arg2$FILEEXT2"
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
