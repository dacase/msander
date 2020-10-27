#!/bin/sh

export MSANDERHOME=`pwd`
./configure --conda --no-netcdf --no-rism --no-boost

cd src
make install
cd ..

rsync -av bin dat lib $PREFIX
