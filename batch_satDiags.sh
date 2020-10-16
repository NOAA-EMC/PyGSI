#!/bin/bash
#SBATCH -J plot_gsi_diags
#SBATCH -A da-cpu
#SBATCH -q debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH -t 00:30:00
#SBATCH --mail-user=$LOGNAME@noaa.gov

# PyGSIdir=path/to/PyGSI_Directory
PyGSIdir=/home/$LOGNAME/PyGSI
OUTDIR=GSI_diag_figures

# load environment needed to run python scripts
source $PyGSIdir/modulefiles/modulefile.PyGSI.hera.bash

python mp_plot_satDiags.py -n 20 -y sat_test_yaml.yaml

mkdir $OUTDIR
mv *.png ./$OUTDIR