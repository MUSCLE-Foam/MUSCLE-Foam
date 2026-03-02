#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions
source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase
rm -rf 0
rm -rf postProcessing/*
rm -f plots/*.png
rm -rf plots/__pycache__
rm -f *.pdf
rm -f *.png
rm -f *.txt
rm -f *.obj
rm -f log*
rm -f *.msh
rm -rf processor*
rm -fr __pycache__
rm -f slurm*
rm -f constant/neutronFlightDirectionsDict