function [ del_est, thet_est] = langevin_dat(X,Niter,burn,thin,nr_of_cores,options)
diary my_data.out;
fprintf('start')
%clear all;

myPool = parpool('local', nr_of_cores);
if isempty(myPool)
    error('pctexample:backslashbench:poolClosed', ...
         ['This example requires a parallel pool. ' ...
          'Manually start a pool using the parpool command or set ' ...
          'your parallel preferences to automatically start a pool.']);
    exit
end
fprintf('pool started')
m = length(unique(X)); [n,p] = size(X);
Rg = 0:(m-1); % gam is approximation parameter
u=options.u;
rho = options.rho;
gam = options.gam;
gam1 = options.gam1;
q = (1/p)^(1+u);
c=10; step = .5;
%Res1 = zeros(p,p);Res2 = zeros(p,p);
ln = length(1:thin:(Niter-burn));
thet_est = zeros(p,p,ln);
del_est = zeros(p,p,ln);
%%%loading and declaration ends%%%
fprintf('initialization phaze 1 complete')
%%%initialization%%%
parfor j = 1:p
    fprintf('\n starting %i. \n',j)
    
%%%initialization of mean and coupling functions%%%

    h_matj = zeros(n,p); 

    for i = 1:n
        for k = 1:p
            
             h_matj(i,k) = X(i,j)*X(i,k)/(m-1)^2;
            
        end
    end
    H_mat = zeros(n*m,p);
    for i = 1:n
        for s = 1:m
            for k = 1:p
               H_mat((i-1)*m + s,k) = (s-1)*X(i,k)/(m-1)^2;
            end
            H_mat((i-1)*m + s,j) = ((s-1)/(m-1))^2;
        end
    end
    gm1 = gam1/(n*p);
    [ theta_old,~] = initial(X,gm1,h_matj,H_mat,Rg, 500,.0001);
   
    dd = abs(theta_old)>0;
%     theta_init = theta_old;
%     d_init = dd;
    res1 = zeros(Niter - burn,p);
    res2= zeros(Niter - burn,p);

    jj= 1; 
    while  jj  <= Niter

    
    
    %%Update theta
        %Langevin
        s = find(dd>0);
        if(isempty(s)==0)
            lo = likehdder(theta_old,dd,Rg,h_matj,H_mat);
            g_old = lo + (theta_old.*dd)/rho + (theta_old.*(1-dd))/gam;
            G_old = g_old*(c/max(c,norm(g_old)));
            h_old = likhd(theta_old,dd,Rg,h_matj,H_mat) + (.5/rho)*norm(theta_old.*dd)^2 + (.5/gam)*norm(theta_old.*(1-dd))^2;
          
            
        for i = s
            theta_new = theta_old;
             
            theta_new(i) = normrnd(theta_old(i) - .5*step*G_old(i),(step));
            %theta_new(i) = normrnd(theta_old(i) -.5*step(i)*G_old(i),(step(i))); %%%%adaptive case
            ln = likehdder(theta_new,dd,Rg,h_matj,H_mat);
            g_new = ln + (theta_new.*dd)/rho  + (theta_new.*(1-dd))/gam;
            G_new = g_new*(c/max(c,norm(g_new)));
            h_new = likhd(theta_new,dd,Rg,h_matj,H_mat) + (.5/rho)*norm(theta_new.*dd)^2 + (.5/gam)*norm(theta_new.*(1-dd))^2;
          
         %save old value
            k1 = - h_new + h_old;
                
            k2 = -((theta_old(i)-(theta_new(i) - step*G_new(i)/2))^2)/(2*step^2)...
                + ((theta_new(i)-(theta_old(i) - step*G_old(i)/2))^2)/(2*step^2) ;
        %%%% k2 = -((theta_old(i)-(theta_new(i) - step(i)*G_new(i)/2))^2)/(2*step(i)^2)...
           %%%     + ((theta_new(i)-(theta_old(i) - step(i)*G_old(i)/2))^2)/(2*step(i)^2) ; %%% adaptive step
            Acc = min(1,exp(k1+k2));
            b = binornd(1,Acc);
            b(isnan(b)) = 1;
            theta_old(i) = theta_new(i)*b + (1-b)*theta_old(i);
            G_old = b*G_new +(1-b)*G_old; h_old = h_new*b + (1-b)*h_old;
        end
        %step(i) = exp(log(step(i)) + (1/jj^0.6)*(Acc - log(step(i))); %%%% adaptive step
        end
  %normal
        do = p - sum(dd);
        if do >0
            
          theta_old(logical(1-dd)) = sqrt(gam)*randn(1,do);
        end
        %%Update delta
    h_old = likhd(theta_old,dd,Rg,h_matj,H_mat) + sum(dd)*log(sqrt(rho)) + sum(1-dd)*log(sqrt(gam)) ...
         + (.5/rho)*norm(theta_old.*dd)^2 + (.5/gam)*norm(theta_old.*(1-dd))^2 - sum(dd)*log(q/(1-q));
    ddprop = double(rand(1,p)<.5);
    id = find(ddprop~=dd);
     for i = id
        
        
            ddn = dd;
            ddn(i) = ddprop(i);
            h_new = likhd(theta_old,ddn,Rg,h_matj,H_mat)+ sum(ddn)*log(sqrt(rho))+ sum(1-ddn)*log(sqrt(gam)) + ...
                (.5/rho)*norm(theta_old.*ddn)^2 + (.5/gam)*norm(theta_old.*(1-ddn))^2 - sum(ddn)*log(q/(1-q));
            k1 = - h_new + h_old;
            Acc = min(1,exp(k1));
            b = binornd(1,Acc);
            b(isnan(b)) = 1;
            dd = ddn*b + (1-b)*dd;
            h_old = h_new*b + (1-b)*h_old;
        
    end
   
   
    

   
   
          
        
        %%Save result
        if jj>burn
            res1(jj-burn,:) = theta_old;
            res2(jj-burn,:) = dd;
         
        end
       % [jj/1000,sum(dd),  rho,count_reinit]
       
    
        jj = jj+1;
    end
    
    
fprintf('\n ending %i. \n',j)

            
%              Res1(j,:) = theta_init;
%              
%              Res2(j,:) = d_init;
        thet_est(j,:,:) = res1(1:thin:end,:)';
        del_est(j,:,:) = res2(1:thin:end,:)';
    



end

diary off;
delete(myPool);



end
