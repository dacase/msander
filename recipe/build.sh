#!/bin/sh

export MSANDERHOME=`pwd`
./configure --conda --no-netcdf --openmp

cd src
make -f Makefile.ap install
cd ..

rsync -av bin dat lib $PREFIX
