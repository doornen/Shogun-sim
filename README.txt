This is a simulation package for gamma ray detectors used at the RIKEN
Nishina Center, in particular the SHOGUN spectrometer.

Details of the package can be found in the file Manual.pdf

Version changes:
1.0.4: Added the option to shift the DALI2 spectrometer along the beam line
       with the command 
       
       ZPOSSHIFT Value

       in the EventBuilder input file. Value is the shift in cm along the beam
       axis relative to the target. The center of the target has always the
       position (0,0,0), Therefore, for off center placed targets, the thickness
       has to be taken into account. 

1.0.5: Added the polyethylene target for Generator and Builder. This target has the 
       index number 8 for the input option TARGET

1.1.2: Fixed some bugs in the Dali2Reconstructor for the calculation of the average
       position
1.1.3: Fixed bugs for GRAPE. Only position information from the segments is included
       Possibility to change the beam pipe thickness is included
1.2.0: Changed for GEANT4.10.2
       Included GAGG and Pr:LuAG as detector materials for SHOGUN
       Changes in the reconstructor of DALI2 and SHOGUN
1.2.1: Included the proper geometry for DALI2+ used in psp17       
1.2.2: The Dali2Reconstructor was not included due to linking problem
