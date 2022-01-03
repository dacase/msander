#!/bin/bash

set -e

BUILD_DIR="../../../../../../cmake-build-debug"
PMEMD="${BUILD_DIR}/src/pmemd/src/pmemd.cuda_SPFP"

echo "Using pmemd=\"$(readlink -m ${PMEMD})\""

set -x

echo "Minimization"
"${PMEMD}" -O \
         -i min.in \
         -o min.out \
         -p wbox.prmtop \
         -c wbox.inpcrd \
         -r min.rst

echo "Heating"
"${PMEMD}" -O \
         -i heat.in \
         -o heat.out \
         -p wbox.prmtop \
         -c min.rst \
         -x heat.nc \
         -r heat.rst \
         -ref min.rst \
         -AllowSmallBox

echo "Refinement start"
"${PMEMD}" -O \
         -i run.in \
         -o run.out \
         -p wbox_xray.parm7 \
         -c heat.rst \
         -x run.nc \
         -r run.rst

"${PMEMD}" -O \
         -i cool.in \
         -o run_cool.out \
         -p wbox_xray.parm7 \
         -c run.rst \
         -x run_cool.nc \
         -r run_cool.rst

echo "Refinement end"