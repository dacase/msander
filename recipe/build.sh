#!/bin/sh

export MSANDERHOME=`pwd`
./configure --conda --openmp

cd src
make install
cd ..

rsync -av bin lib $PREFIX
