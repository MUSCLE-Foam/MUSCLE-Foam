#!/bin/bash

# This case corresponds to the CNRS benchmark on multi-physics codes
# as described in Tiberga et al. (2020)
# https://doi.org/10.1016/j.anucene.2020.107428

source $FOAM_SRC/../bin/tools/RunFunctions

cp -r 0.orig 0

runApplication blockMesh