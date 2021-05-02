#!/bin/bash

INPUTFILE="INPUTFILE ../../../SimulationResults/Builder/Shogun/112sn/Case"
FILEEXT0="mg"
FILEEXT1="fs.root"
FILEEXT2="_112sn112sn180mev"
OUTPUTFILE="OUTPUTFILE ../../../SimulationResults/Reconstructor/Shogun/112sn/Case"
SPECTRABINANDRANGE="SPECTRABINANDRANGE 250 0. 2000."
BETADOPPLERAVERAGE="BETADOPPLERAVERAGE "
BETATOFAVERAGE="BETATOFAVERAGE 0.5456 "
DECAYPOSITIONZ="DECAYPOSITIONZ 0.0"
END="END"

FILE="./input/Reconstructor.in"

beta="0"
BETA=(0.5443
    0.5431
    0.5431
    0.5431
    0.5431
    0.5431
    0.5431
    0.5431
    0.5431
    0.5431
    0.5431
    0.5416
    0.5390
    0.5384
    0.5382
    0.5381
    0.5380
    0.5379
    0.5379
    0.5379
    0.5378
    0.5378
    0.5383
    0.5348
    0.5335
    0.5328
    0.5324
    0.5321
    0.5319
    0.5317
    0.5316
    0.5315
    0.5315
    0.5348
    0.5310
    0.5290
    0.5278
    0.5271
    0.5265
    0.5261
    0.5258
    0.5255
    0.5253
    0.5251
    0.5311
    0.5271
    0.5247
    0.5231
    0.5219
    0.5210
    0.5204
    0.5199
    0.5195
    0.5191
    0.5188
    0.5274
    0.5232
    0.5205
    0.5184
    0.5169
    0.5157
    0.5148
    0.5141
    0.5134
    0.5129
    0.5125
    0.5235
    0.5191
    0.5162
    0.5138
    0.5119
    0.5104
    0.5092
    0.5082
    0.5074
    0.5067
    0.5062
    0.5195
    0.5150
    0.5118
    0.5091
    0.5069
    0.5051
    0.5036
    0.5024
    0.5014
    0.5005
    0.4997
    0.5153
    0.5106
    0.5073
    0.5044
    0.5019
    0.4998
    0.4980
    0.4965
    0.4952
    0.4941
    0.4932
    0.5109
    0.5061
    0.5026
    0.4994
    0.4967
    0.4943
    0.4923
    0.4905
    0.4890
    0.4877
    0.4865
    0.5064
    0.5015
    0.4978
    0.4944
    0.4914
    0.4888
    0.4865
    0.4845
    0.4827
    0.4811
    0.4797)
a=0
for arg1 in 1 100 200 300 400 500 600 700 800 900 999
do
    
    for arg2 in 0 100 200 300 400 500 600 700 800 900 999
    
    do
        echo $beta
        echo $arg1
        echo $arg2
        beta=${BETA[$a]}
        echo $beta
        echo ""$INPUTFILE$1$FILEEXT2$arg1$FILEEXT0$arg2$FILEEXT1"
"$OUTPUTFILE$1$FILEEXT2$arg1$FILEEXT0$arg2$FILEEXT1"
"$SPECTRABINANDRANGE"
"$BETADOPPLERAVERAGE$beta"
"$BETATOFAVERAGE"
"$DECAYPOSITIONZ"
"$END" 
" > $FILE
        
        time ./ShogunReconstructor
    let a=a+1
    done

done

exit 0




