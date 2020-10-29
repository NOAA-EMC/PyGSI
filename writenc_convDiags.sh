#!/bin/bash

# PyGSIdir=path/to/PyGSI_Directory
PyGSIdir=/home/$LOGNAME/PyGSI
OUTDIR=/home/$LOGNAME/PyGSI/ncfiles

# load environment needed to run python scripts
source $PyGSIdir/modulefiles/modulefile.PyGSI.hera.bash

DATES=(20200920 20200921)
HOURS=(00 06 12 18)

for DATE in ${DATES[@]}
do
    for HOUR in ${HOURS[@]}
    do
        python create_convNetCDF.py -y /home/$LOGNAME/PyGSI/conv_${DATE}${HOUR}_yaml.yaml -o $OUTDIR
    done
done