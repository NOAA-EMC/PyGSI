# Getting Started

## Cloning Repository

The Github repository can be found here: https://github.com/NOAA-EMC/PyGSI
and can be cloned by using the follow command:

```
$ git clone https://github.com/NOAA-EMC/PyGSI.git
```
This will clone the develop branch into a directory called PyGSI.

## Setting Environment

### On Supported Platforms

To load the proper environment when working on Hera, use the following commands:

    cd PyGSI
    module use modulefiles
    module load PyGSI/hera

To load the environment on Orion, use the following commands:

    cd PyGSI
    module use modulefiles
    module load PyGSI/orion

### On Local Machine

If working on a local machine, you can install PyGSI using pip.

```
cd PyGSI
pip install .
```
