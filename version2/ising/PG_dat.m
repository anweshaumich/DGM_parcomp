function [ del_est, thet_est] = PG_dat(Z,Niter,burn,thin,nr_of_cores,options)
diary my_data.out;
fprintf('start')
%clear all;
set (parcluster('local'), 'NumWorkers', 40)
myPool = parpool('local', nr_of_cores);
if isempty(myPool)
    error('pctexample:backslashbench:poolClosed', ...
         ['This example requires a parallel pool. ' ...
          'Manually start a pool using the parpool command or set ' ...
          'your parallel preferences to automatically start a pool.']);
    exit
end
fprintf('pool started')
m = length(unique(Z)); [n,p] = size(Z);
Rg = 0:(m-1); % gam is approximation parameter
u=options.u;
rho = options.rho;
gamma = options.gam;
gam1 = options.gam1;
prior_q = (1/p)^(1+u);
%Res1 = zeros(p,p);Res2 = zeros(p,p);
ln = length(1:thin:(Niter-burn));
thet_est = zeros(p,p,ln);
del_est = zeros(p,p,ln);
%%%loading and declaration ends%%%
fprintf('initialization complete')
%%%initialization%%%
parfor j = 1:p
    fprintf('\n starting %i. \n',j)
    
%%%initialization%%%

    
    gm1 = gam1/(n*p);
    [ theta_old,~] = initial(Z,gm1,j,Rg, 500,.0001);
    fprintf('\n initialized %i. \n',j)
   theta = theta_old';
    delta = abs(theta)>0;
%     theta_init = theta_old;
%     d_init = dd;
    
    res1 = zeros(Niter - burn,p);
    res2= zeros(Niter - burn,p);
    y = Z(:,j); X = Z; X(:,j) = ones(n,1);

    jj= 1; 
    while  jj  <= Niter

    %Update W
    W = pgdraw(abs(X(:,delta)*theta(delta)));
    ybar = (diag(W.^(-0.5)))*(y-0.5);
    Xbar = diag(W.^(0.5))*X;
    Xpy=(Xbar')*ybar; 
    XpX = (Xbar')*Xbar;

    %Update theta
    ind = find(delta==0);
    if isempty(ind)
    else
        theta(ind) = sqrt(gamma)*randn(length(ind),1);
    end
    ind = find(delta==1);
    if isempty(ind)
    else
        %v_vec = (1/rho)*ones(length(ind),1);
        v_vec = rho*ones(length(ind),1);
        Sigma_inv  = (diag(1./v_vec) + XpX(ind,ind));
        C1 = chol(Sigma_inv);
        mu = (C1 \( (C1') \ Xpy(ind)  ));
        theta(ind) = mu + C1\randn(length(ind),1);
    end

    
    %Update delta
    for jj = 1:p
        prop = (rand(1)<=0.5);
        if ( delta(jj) == 0 ) && (prop == 1)
            r = log(1-prior_q) - log(prior_q) - 0.5*log(gamma/rho) ...
                + 0.5*(1/rho - 1/gamma)*(theta(jj)^2)...
                +0.5*(theta(jj)^2)*XpX(jj,jj)- theta(jj)*(Xpy(jj)...
                - sum(theta(delta)'.*XpX(jj,delta)));
            Acc  = exp(-r);
            if rand(1)<= Acc
                delta(jj) = prop;
            end
        elseif ( delta(jj) == 1 ) && (prop == 0)
            r = log(1-prior_q) - log(prior_q) - 0.5*log(gamma/rho) ...
                + 0.5*(1/rho - 1/gamma)*(theta(jj)^2)...
                +0.5*(theta(jj)^2)*XpX(jj,jj)- theta(jj)*(Xpy(jj)...
                - sum(theta(delta)'.*XpX(jj,delta)) + theta(jj)*XpX(jj,jj));
            Acc = exp(r);
            if rand(1)<= Acc
                delta(jj) = prop;
            end
        end
        
    end

        
    

   
    
          
        
        %%Save result
        if jj>burn
            res1(jj-burn,:) = theta';
            res2(jj-burn,:) = delta';
         
        end
        %[jj/1000,sum(dd),  rho,count_reinit]
       
    
        jj = jj+1;
    
    
    end
fprintf('\n ending %i. \n',j)

            
%              Res1(j,:) = theta_init;
%              
%              Res2(j,:) = d_init;
        thet_est(j,:,:) = res1(1:thin:end,:)';
        del_est(j,:,:) = res2(1:thin:end,:)';
    



end
%clearvars -except block burn ct Niter p m n Res1 Res2 Rg step thet_est del_est X q sam

%save('result_ispot','del_est','thet_est','Res1','Res2');
diary off;
delete(myPool);



end
