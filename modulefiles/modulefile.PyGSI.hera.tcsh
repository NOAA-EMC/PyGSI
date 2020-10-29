setenv PATH /scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/bin\:$PATH

set srccmd=($_)
set pwddir=`pwd`
set mydir=`dirname $srccmd[2]`
set PyGSIdir=`cd $mydir/../; pwd`
setenv PYTHONPATH $PyGSIdir\:$PYTHONPATH

source /scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/etc/profile.d/conda.csh

conda init tcsh
source ~/.tcshrc
conda activate PyGSI
