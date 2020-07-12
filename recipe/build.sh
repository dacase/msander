#!/bin/sh

export MSANDERHOME=`pwd`
./configure 

make install.ap
make clean

#if [ `uname` == "Darwin" ]; then
#    mkdir -p lib/amber_3rd_party # to avoid overwritten by other programs
#    python ./recipe/copy_and_fix_gfortran.py \
#        /usr/local/gfortran/lib lib/amber_3rd_party
#    mkdir -p $PREFIX/lib
#    rsync -av lib/amber_3rd_party $PREFIX/lib
#fi

rsync -av bin dat lib $PREFIX
