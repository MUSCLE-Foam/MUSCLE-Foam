#!/bin/bash

source $FOAM_SRC/../bin/tools/RunFunctions

wmake code

cp -r 0.orig 0

runApplication blockMesh
runApplication createZones