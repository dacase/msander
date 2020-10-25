#!/bin/sh

export MSANDERHOME=`pwd`
./configure_conda --verbose --no-netcdf

cd src
make -f Makefile.ap install
cd ..

rsync -av bin dat lib $PREFIX
