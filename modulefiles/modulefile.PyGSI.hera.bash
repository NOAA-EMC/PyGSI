export PATH="/scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/bin:$PATH"

PyGSIdir="$(cd "$(dirname "$BASH_SOURCE")/../"; pwd)"
export PYTHONPATH="${PYTHONPATH}:${PyGSIdir}"

source /scratch1/NCEPDEV/da/Kevin.Dougherty/anaconda3/etc/profile.d/conda.sh

conda init bash
source ~/.bashrc
conda activate PyGSI
