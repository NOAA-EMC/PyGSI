setenv PATH $PATH\:/scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/bin

set PyGSIdir="$(cd "$(dirname "$BASH_SOURCE").//"; pwd)"
setenv PYTHONPATH $PYTHONPATH\:$PyGSIdir

source /scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/etc/profile.d/conda.sh

conda init tcsh
source ~/.tcshrc
conda activate PyGSI
