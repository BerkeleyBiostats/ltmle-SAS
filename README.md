# LTMLE SAS macro

[ltmle-SAS](https://github.com/BerkeleyBiostats/ltmle-SAS) is a SAS implementation Longitudinal Targeted Maximum Likelihood Estimation (LTMLE).

Author: Jordan Brooks


## Run the macro

This project depends on [SuperLearner-SAS](https://github.com/BerkeleyBiostats/SuperLearner-SAS).

The `./data` folder is where the simulated data is deposited as a permanent SAS dataset. For data analyses, this is were the data sample would be deposited prior to analysis. It is empty and could be deleted, but anyone actually using the SAS code would need to create it again. Alternatively, SAS can be coded to check for the directory and then create it if does not already exist but I have not implemented that.