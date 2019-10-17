# DGM_parcomp software
Matlab code for Bayesian inference of discrete graphical models using parallelization. 
The calling function is par_discretegraph

Requires five input values
data : the mat file consisting the data without headers (data should consist of numerical values)
function cannot handle missing values. The minimum value in the data should be either one or zero.
Niter: the total number of iterations that the MCMC chain will run
burn: the number of iterations allowed for burn in of the chain.
thin: the thinning value of the markov chain
nr_of_cores: the number of cores to be utilized for parallel operations
PG: if set at 0 will utilize PG algorithm for Ising model, else will use MHwG algorithm
default value set at 1.

Algorithm automatically uses Ising model if state space is {0,1}.

Output file has six colums.
column 1 and 2 provides the row column index
column 3 provides the estimate
column 4 and 5 provides lower and upper bounds of credible intervals respectively.
column 6 provides probability of the parameter being active 

the values are arranged in increasing order of estimated strength. 

# data.csv
PTSD symptoms data used by Epskamp et al(2018).

# distmiss.csv
16PF data filtered and adjusted by missing values
