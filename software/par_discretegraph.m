%%% Automatically specifies most tuning parameters
%%% If required one can change some of the parameters involved


%function par_discretegraph(data,Niter,burn,thin,nr_of_cores,PG)
if nargin > 6
    error('TooManyInputs', ...
        'requires at most 1 optional inputs');
end

% Fill in unset optional values.
switch nargin
    case 5
        PG = 1;
end

if(min(unique(data))>1)
    fprintf('error in input format!')
elseif(min(unique(data))==1)
    data = data-1;
end

[n,p] = size(data);

%% penalty parameter initialization for initialization algo
options.gam1 = 300; 

%% prior parameters (gam for inactives) and (rho_in for initialization of prior for active) 
options.gam = 0.1/max(p,n); options.rho = sqrt(n/log(p));options.u = 2;

% %% langevin parameters
% options.c = 10; %% stabilizing gradient for the langevin algorithm
% options.step_in = .5; %%% step size for gradient search in langevin algorithm


 

if(length(unique(data))>2)
    cd ./potts
    [ delta, theta] = langevin_dat(data,Niter,burn,thin,nr_of_cores,options);
else
    cd ./ising
    if PG~=0
    [delta,theta] = PG_dat(data,Niter,burn,thin,nr_of_cores,options);    
    else
    [delta, theta] = langevin_dat(data,Niter,burn,thin,nr_of_cores,options);
    end
end

 cd ../
 
 [confidence] = confidence_edge(theta,delta);
 
 %save('result_gausscor','confidence','edge','theta_init','delta_init','dumbvar','-v7.3') ;
 %csvwrite('confidence.txt',confidence);
 
 

%end