#!/bin/sh

export MSANDERHOME=`pwd`
./configure  --conda --verbose

make install.ap
make clean

if false; then
if [ `uname` == "Darwin" ]; then
    mkdir -p lib/amber_3rd_party # to avoid overwritten by other programs
    /Users/case/miniconda/bin/python ./recipe/copy_and_fix_gfortran.py \
        /usr/local/gfortran/lib lib/amber_3rd_party
    mkdir -p $PREFIX/lib
    rsync -av lib/amber_3rd_party $PREFIX/lib
fi
fi

rsync -av bin dat lib $PREFIX
