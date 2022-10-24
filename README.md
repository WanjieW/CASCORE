# CASCORE
Data and Rcode for Covariate Assisted SCORE algorithm. 


Simulation folder: the simulation code, including

  -ExpSets.R: the main file to run, which gives the result and plot for two experimental settings.
  
  -TuningParameter.R: the main file to explore the effects of tuning parameter alpha.
  
  -Simulation.R: main simulation file that generate data.
  
  -VariousFunctions.R: the files containing several community detection methods
  
  -Mtable.R: the file to calculate the error rate
  
  -All other files: required in VariousFunctions.R. 
  
  
  
LastFM folder: the data and code about LastFM app user data
  
   -lastfm_asia_edges.csv and lastfm_asia_target.csv: the edge information of LastFM network and the unspecified country labels.
   
   -CovariatesMatrix.Rdata: the list of liked artists in LastFM data.
   
   -LastFM_Main.R: the main code file to process the data. Require file: DataSubsets.R and Dataplot.R.
   
   -VariousFunctions.R: the files containing several community detection methods
   
   -Mtable.R: the file to calculate the error rate
  
   -All other files: required in VariousFunctions.R. 
