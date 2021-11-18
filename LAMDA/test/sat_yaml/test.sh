#!/bin/bash
set -x
CYCLE=2021100400
exp=test
DIAGDIR=/scratch1/NCEPDEV/da/Kevin.Dougherty/PyGSI/LAMDA/data/diags
SATINFO=/scratch2/NCEPDEV/fv3-cam/Xiaoyan.Zhang/noscrub/gsi_code/gsi_master/fix/global_satinfo.txt
yaml_file=${exp}.yaml
echo $DIAGDIR
python ush/satinfo2yaml.py -d $DIAGDIR -c $CYCLE -i $SATINFO -y $yaml_file

