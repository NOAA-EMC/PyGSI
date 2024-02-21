#!/bin/sh -xvf
# run_error_model_estimate.sh
# This script is intended to simplify/automate obs error parameters estimations.
# This script will call the python script to:
# - Compute error parameters and produce a csv file
# - Plot mean cloud amount vs FG departure std.dev along with original error parameters

#--------------- User modified options below -----------------

# specify the directory for input config files
config_path=/scratch1/NCEPDEV/stmp4/Azadeh.Gholoubi/GDAS-ops/PyGSI/ush/obsErrorAssignmentTool/config_files/ 

# specify the global_satinfo.txt and cloudy_radiance_info.txt path in GSI-fix directory
global_satinPath=../../fix/global_satinfo.txt  # Path to global_satinfo.txt
cloudy_path=../../fix/cloudy_radiance_info.txt # Path to cloudy_radiance_info.txt

#read arguments : specify sensor name, list of channels, bin size, output file name, and qc flag for filtering
sensor=amsua_n19
Channels="1,2"
bin_size=0.005
bindir=bin005
qc_flag=0

# specify the directory for saving output plots and csv file
output=/scratch1/NCEPDEV/stmp4/Azadeh.Gholoubi/GDAS-ops/PyGSI/ush/obsErrorAssignmentTool/output/$bindir/

mkdir -p $output

machine=${machine:-hera} 

if [ $machine = orion ]; then
   PyGSI=${PyGSI:-/Path/to/PyGSI/} # Change this to your own branch
elif [ $machine = hera ]; then
   PyGSI=${PyGSI:-/scratch1/NCEPDEV/stmp4/Azadeh.Gholoubi/GDAS-ops/PyGSI/} # Change this to your own branch
else
   echo "Machine " $machine "not found"
   exit 1
fi


#-------------- Do not modify below this line ----------------
# Load Modules 
module use $PyGSI/modulefiles
module load PyGSI/$machine

# call the Python script
 
python ObsErrorAssignment.py -n $config_path -o $output -g $global_satinPath -l $cloudy_path -s $sensor -c $Channels -b $bin_size -q $qc_flag 
