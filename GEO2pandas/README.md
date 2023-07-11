# GEO2pandas
GEO2pandas is a simple python module that makes it easy to import transcriptome data from the GEO repository as pandas DataFrame object.

Written by Francesc Font-Clos
Please report bugs at
francesc.font@gmail.com

### Main features
+ Downloads .soft.gz files from GEO server only if not found locally
+ Creates two dataframes: **expr** for the gene expression data, and **meta** for the samples metadata
+ Tries to automatically guess which field holds the EntrezID gene name
+ Tries to autoamtically guess if the data is log2 transformed or not
+ Can average out duplicate columns

### TO DO list
+ Create unit tests
