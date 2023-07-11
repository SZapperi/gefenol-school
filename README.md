# Hands-on activity for the XI GEFENOL summer school
Course by Caterina A. M. La Porta and Stefano Zapperi

In this hands-on activity is based on a tutorial showing the basic functionality of `GEO2pandas`.

`GEO2pandas` is a python package written by Francesc Font-Clos
to retrieve data and metadata from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
into `pandas.DataFrame` objects.

The tutorial is part of the XI GEFENOL summer school in Barcelona and is adapted from a
tutorial set up by Francesc Font-Clos at the University of Milan.

### Before starting
To follow this lecture, you need a modern installation of python, together with jupyter, numpy, matplotlib and some other standard python libraries. The simplest way to install all these packages without interfeering with your current python installation is the [Anaconda distribution](https://www.anaconda.com/download/). 
Choose python 3.x and your OS, download, install, and you should be good to go.

Once you have install python you should download the entire repository and local install GEO2pandas.

Use a terminal window, cd to the GEO2pandas directory and install the following:

- pip install geoparse == 2.0.1
- pip install pandas == 1.0.3
- pip install numpy == 1.18.1
- pip install -e ./

the last command will install GEO2pandas 



### Hands-on structure

*First part*  
Starting just with from a keyword of interest (e. g. **obesity**),
I will show how to do some basic exploratory analysis of public transcriptomic data related
to our keyword.

1. Searching for datasets: GEO and ArrayExpress
2. Getting the data. Possible issues.
3. Basic analysis
  3.1 Differentially expressed genes
4. Basic visualization
 4.1 Heatmap of DE genes

*Second part*  
Starting with a keyword of your choice, you will try
to reproduce the steps 1-5 above. Problems might arise. 

