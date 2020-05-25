#!/bin/sh

export MSANDERHOME=`pwd`
./configure 

make install.ap
make clean

rsync -av bin dat lib $PREFIX
