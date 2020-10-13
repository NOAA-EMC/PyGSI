setenv PATH $PATH\:/scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/bin
setenv PYTHONPATH $PYTHONPATH\:/home/Kevin.Dougherty/PyGSI/

source /scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/etc/profile.d/conda.sh

conda init tcsh
source ~/.tcshrc
conda activate PyGSI
