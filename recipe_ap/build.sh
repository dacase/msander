#!/bin/sh

export AMBERHOME=`pwd`
./configure_ambermini --conda --no-netcdf --no-rism

cd AmberTools/src 
make -f Makefile.ambermini install
make -f Makefile.ambermini clean
cd ../..

if [ `uname` == "Darwin" ]; then
    mkdir -p lib/amber_3rd_party # to avoid overwritten by other programs
    python ./recipe_ap/copy_and_fix_gfortran.py \
        /usr/local/gfortran/lib lib/amber_3rd_party
    mkdir -p $PREFIX/lib
    rsync -av lib/amber_3rd_party $PREFIX/lib
fi

rsync -av bin dat lib AmberTools $PREFIX
