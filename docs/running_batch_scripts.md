# Running Batch Scripts

If a user is interested in producing results for a wide range of inputs or analyzing the YAMLs produced by the global fix files, PyGSI has scripts that will submit a batch job. To do so, the user will need three things:

* Path to user's PyGSI directory
* Path to output directory
* Path to input YAML file

The shell script then executes a python script that will create multiple output diagnostic figures. 

### Conventional

Navigate to the shell script [`PyGSI/scripts/batch_conv_diags.sh`](https://github.com/NOAA-EMC/PyGSI/blob/develop/scripts/batch_conv_diags.sh). The example code will look as follows:

```
#!/bin/bash
#SBATCH -J plot_gsi_diags
#SBATCH -A da-cpu
#SBATCH -q debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH -t 00:30:00
#SBATCH --mail-user=$LOGNAME@noaa.gov

# PyGSIdir=path/to/PyGSI_Directory
PyGSIdir=../
OUTDIR=../
YAML=$PyGSIdir/yamls/diag_conv_t_ges.yaml

# load environment needed to run python scripts
# source $PyGSIdir/modulefiles/modulefile.PyGSI.hera.bash

python $PyGSIdir/scripts/mp_plot_conv_diags.py -n 20 -y $YAML -o $OUTDIR
```

Change the `PyGSIdir`, `OUTDIR`, and `YAML` variables to point to the users appropriate paths. The Python script being called uses multi-processing and the number of nodes can be changed with the `-n` input. This example uses `-n 20`.

Once the changes have been made to reflect the appropriate inputs, save the `batch_conv_diags.sh` file and run the following on the command line:

```
sbatch batch_conv_diags.sh
```

### Radiance

The instructions are the same for radiance files. The script can be found at [`PyGSI/scripts/batch_sat_diags.sh`](https://github.com/NOAA-EMC/PyGSI/blob/develop/scripts/batch_sat_diags.sh).

### Ozone

The instructions are the same for ozone files. The script can be found at [`PyGSI/scripts/batch_ozone_diags.sh`](https://github.com/NOAA-EMC/PyGSI/blob/develop/scripts/batch_ozone_diags.sh).