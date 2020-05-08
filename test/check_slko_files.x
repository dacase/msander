#!/bin/csh -f
# checks if cc.spl file exists - if it does then it
# assumes that all slko files are installed.
# by default assumes DFTB2 with mio-1-1 parameters

set paramset='mio-1-1'
set method='DFTB2'

if ($1 == 'DFTB3') then
    set paramset='3ob-3-1'
    set method='DFTB3'
endif

if( -r $MSANDERHOME/dat/slko/$paramset/C-C.skf ) then
    #Exists
    exit(0)
else
   echo "$method $paramset SLKO files not found - Skipping DFTB tests..."
   echo ""
   exit(1)
endif


