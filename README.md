# Federated Principal Component Analysis

This repository contains additional information regarding federated principal component analysis. Library code for the simulation of federated PCA is available in [vertical](/python/PCA/vertical). The folder [misc_scripts](./misc_scripts) contains examples on how to run the code and a basic workflow. Dependencies are managed using a conda environment and can be installed using [environment.yaml](./environmnent.yaml) 

The specific use case for principal component analysis handeled by this version of principal compoenent analysis is patient stratification for medical data. This is a very special use case, where the patients are distributed over serveral sites, which cannot share their data. 

If you would like to try federated principal component analysis in practice, please visit https://federated.compbio.sdu.dk/. The web service allows you to set up a federated study, which includes using the web interface to set the parameters, and downloading a user client, which handles the communication on the client site. More detailed information can be found on the web site. Please note that it is experimental software.
