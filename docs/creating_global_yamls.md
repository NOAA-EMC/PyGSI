# Creating Global Input YAMLs

In [`PyGSI/ush`](https://github.com/NOAA-EMC/PyGSI/tree/develop/ush), [`convinfo2yaml.py`](https://github.com/NOAA-EMC/PyGSI/blob/develop/ush/convinfo2yaml.py), [`satinfo2yaml.py`](https://github.com/NOAA-EMC/PyGSI/blob/develop/ush/satinfo2yaml.py), and [`ozinfo2yaml.py`](https://github.com/NOAA-EMC/PyGSI/blob/develop/ush/ozinfo2yaml.py) are scripts used to create all the necessary input yamls for all the input data from the global input files that are found in [`PyGSI/fix`](https://github.com/NOAA-EMC/PyGSI/tree/develop/fix).

## How to run:

Each python file is set up to receive similar inputs to create the desired yaml files. This example will use the `convinfo2yaml.py` script, however the other scripts generate yamls in the exact same way. To generate the yaml files, run the following from the command line:

```
python convinfo2yaml.py -d <> -c <> -i <> -y <> -l <> -v <>
```
The inputs are described below:

    Required:
        -d: (--diagdir) local path to GSI diagnostic netCDF4 files
        -c: (--cycle) cycle being used for GSI diagnostic file in YYYYMMDDHH format
        -i: (--infofile) path to global information file. These files are located in PyGSI/fix
        -y: (--yaml) path to where the users output yamls will be stored
    
    Optional:
        -l: (--loop) whether the user is interested in guess (ges) or analysis (anl) GSI diagnostic files. Default is ‘ges’.
        -v: (--variable) whether the user is interested in reading departures (omf), observations (obs), or H(x) (hofx). Default is ‘omf’. 
     
        ** Note, if -l anl, omf will result in oma values.

Example:
```
python convinfo2yaml.py -d /path/to/gdasfile/gdas.20230223/00/atmos/ -c 2023022300 -i /path/to/PyGSI/fix/global_convinfo.txt -y path/to/output_yamls -l ges -v omf
```
