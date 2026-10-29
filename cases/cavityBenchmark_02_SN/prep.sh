#!/bin/bash

# This case corresponds to the CNRS benchmark on multi-physics codes
# as described in Tiberga et al. (2020)
# https://doi.org/10.1016/j.anucene.2020.107428

ORD="${1:-S6}"

# Highest scattering order a level-symmetric S_N set integrates exactly is N/2
case "$ORD" in
    S2) SCATORD=2 ;;
    S6) SCATORD=4 ;;
    *)
        echo "Error: Invalid SN discretization. Only S2 and S6 are valid" >&2
        exit 1
        ;;
esac

source $FOAM_SRC/../bin/tools/RunFunctions

cp -r 0.orig 0

cat "constant/neutronFlightDirectionsDict_${ORD}" > constant/neutronFlightDirectionsDict

sed -i "s/^\(scatteringOrders[[:space:]]\+\)[0-9]\+;/\1${SCATORD};/" system/controlDict

if ! grep -q "^scatteringOrders[[:space:]]\+${SCATORD};" system/controlDict; then
    echo "Error: Could not set scatteringOrders in system/controlDict" >&2
    exit 1
fi

runApplication blockMesh