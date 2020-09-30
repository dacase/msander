#!/bin/sh

export MSANDERHOME=`pwd`
./configure_conda --no-netcdf --verbose

cd src
make -f Makefile.ap install
cd ..

rsync -av bin dat lib $PREFIX
