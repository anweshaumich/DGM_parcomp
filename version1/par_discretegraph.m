%%% Automatically specifies most tuning parameters
%%% If required one can change some of the parameters involved


function par_discretegraph(inputfile,Niter,burn,thin,nr_of_cores)
X = csvread(inputfile,1,0);
if(min(unique(X))>1)
    fprintf('error in input format!')
elseif(min(unique(X))==1)
    X = X-1;
end

[n,p] = size(X);
gam1 = 300; %% penalty parameter for optimization algo to initiate mcmc
gam = 0.1/(n*p); rho_in = .5*sqrt(log(p)/n); rho_fixed = 0;
%%% prior parameters (gam for inactives) and (rho_in for initialization/fixed laplace prior) 

c = 10; %% stabilizing gradient for the langevin algorithm

step = .5; %%% step size for gradient search in langevin algorithm

tow_in = 1;compact_in = 2; mc_rate = .6; inc = .05;%%% parameters for tuning the independent sampler
 
sam = 100; block = 20; %%% parameters for approximating double derivative and inactive theta update
  

if(length(unique(X))>2)
    cd ./potts
else
    cd ./ising
end
[~,~, delta, theta] = langevin_dat_gaussian(X,gam,sam,step,c,block,...
                                                    Niter,burn,gam1,rho_in,rho_fixed, tow_in,...
                                                         compact_in,mc_rate,inc,nr_of_cores);
 cd ../
 
 [confidence] = confidence_edge(theta,delta,thin);
 
 %save('result_gausscor','confidence','edge','theta_init','delta_init','dumbvar','-v7.3') ;
 csvwrite('confidence.txt',confidence);
 
 

end