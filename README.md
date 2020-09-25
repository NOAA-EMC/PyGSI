__Kevin Dougherty__  
__September 2020__  


# PyGSI
Scripts used to validate GSI diagnostic files for JEDI


# Anaconda Environment
> Please read the documents for managing environments.  
> Reference: https://github.com/Unidata/unidata-users-workshop  
> Reference: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
> Reference: https://www.youtube.com/watch?v=15DNH25UCi0

I installed the Anaconda python distribution and created a new environment using the 'environment.yml' file provided in this directory. The name of the environment is DA_Diags. You can copy the the .yml file and run the following command in your terminal

    conda env create -f ./environment.yml
    
On a windows computer, to activate the `DA_Diags` environment, do this in the Windows Command Prompt:

    activate DA_Diags

Or, if you are in the PowerShell

    cmd
    activate DA_Diags

If you are using a `bash` shell in Linux:

    conda init bash  # Only need to do this once to initialize the correct shell
    conda activate DA_Diags

If you are using a `tcsh` shell in Linux:

    conda init tsch  # Only need to do this once to initialize the correct shell
    conda activate DA_Diags


## Update environment
Deactivate the environment

    conda deactivate DA_Diags

Update the environment.yml file, and update the conda environment

    conda env update -f environment.yml

List all the available environments

    conda info --envs