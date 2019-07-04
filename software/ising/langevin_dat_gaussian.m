function [Res1,Res2, del_est, thet_est] = langevin_dat_gaussian(X,gam,~,step_in,c,block,Niter,burn,gam1,rho_in,rho_fixed, tow_in,compact_in,mc_rate,inc,nr_of_cores)
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
u=log(p);
q = (1/p)^(1+u);
Res1 = zeros(p,p);Res2 = zeros(p,p);
thet_est = zeros(p,p,Niter-burn);
del_est = zeros(p,p,Niter-burn);
%%%loading and declaration ends%%%
fprintf('initialization complete')
%%%initialization%%%
parfor j = 1:p
    fprintf('\n starting %i. \n',j)
    
%%%initialization%%%

    
    gm1 = gam1/(n*p);
    [ theta_old,~] = initial(X,gm1,j,Rg, 500,.0001);
   
    dd = abs(theta_old)>0;
    theta_init = theta_old;
    d_init = dd;
    rho = rho_in; tow_o = tow_in; step = step_in;
    res1 = zeros(Niter - burn,p);
    res2= zeros(Niter - burn,p);
%parameter for Update tow
    compact = compact_in; count_reinit = 0;
    jj= 1; ff = 1;
    while  jj  <= Niter

    
    
    %%Update theta
        %Langevin
        s = find(dd>0);
        if(isempty(s)==0)
             lo = likehdder(theta_old,Rg,X,j);
             g_old = lo + rho*sign(theta_old.*dd)  + (theta_old.*(1-dd))/gam;
             G_old = g_old*(c/max(c,norm(g_old)));
             h_old = likhd(theta_old,Rg,X,j) + rho*sum(abs(theta_old.*dd)) + (.5/gam)*norm(theta_old.*(1-dd))^2;
          
            
        for i = s
            theta_new = theta_old;
             
            theta_new(i) = normrnd(theta_old(i) - .5*step*G_old(i),(step));
            %theta_new(i) = normrnd(theta_old(i) -.5*step(i)*G_old(i),(step(i))); %%%%adaptive case
            ln = likehdder(theta_new,Rg,X,j);
            g_new = ln + rho*sign(theta_new.*dd)  + (theta_new.*(1-dd))/gam;
            G_new = g_new*(c/max(c,norm(g_new)));
            h_new = likhd(theta_new,Rg,X,j) + rho*sum(abs(theta_new.*dd)) + (.5/gam)*norm(theta_new.*(1-dd))^2;
          
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
            G_old = G_new*b + G_old*(1-b); h_old = h_new*b + h_old*(1-b);
        end
        %step(i) = exp(log(step(i)) + (1/jj^0.6)*(Acc - log(step(i))); %%%% adaptive step
        end
  %normal
        do = p - sum(dd);
        if do >0
            
           %form blocks
            x = X; x(:,j) = 1;
            t = theta_old.*dd;
            X_t = x(:,dd==1)*t(dd==1)';
            mat1 = exp(X_t);
            w = mat1./((1+mat1).^2);
             
            w(isnan(w))=1;
            W  =  sparse(1:n,1:n,w,n,n);
            
            set = find(dd==0); cl = ceil(do/block);
            y = repmat((1:cl)',1,block); y = reshape(y',1,block*cl);y = y(1:length(set));
            
            
            %block update
            for s = 1:cl
                act = set(y==s); inact = 1:p;  inact(act) = []; sact = length(act);
                cact = x(:,act); %cinact = x(:,inact);
                ind = ones(1,p);ind(act) = 0;
                ldr = likehdder(theta_old.*ind,Rg,X,j);
                Sig = eye(sact) + (gam)*(cact'*W*cact);
                chl = chol(Sig); dim = size(Sig,1); 
                chl_in = chl\eye(dim); Sig_in = chl_in*chl_in';
                mean = gam*(Sig_in*ldr(act)');
                un  = mean + sqrt(gam)*(chl_in*randn(sact,1));
                uo = theta_old(act);  
                theta_new = theta_old;
                theta_new(act) = un;
                
                k3 = -likhd(theta_new,Rg,X,j) - (.5/gam)*norm(theta_new.*(1-dd))^2 + likhd(theta_old,Rg,X,j) + (.5/gam)*norm(theta_old.*(1-dd))^2;
                k4 = (1/(2*gam))*(uo - mean')*(Sig_in*(uo' - mean)) - (1/(2*gam))*(un' - mean')*(Sig_in*(un - mean));
                k5 = min(1,exp(k3 + k4));
                b = binornd(1,k5);
                theta_old = b*theta_new + (1-b)*theta_old;
               
               
            end
        end
        %%Update delta
    
     
   
        k = rho*abs(theta_old) - (.5/gam)*(theta_old.^2);
        r = log(q/(1-q)) + .5*log(2*pi*gam) - log(2/rho) - k;     
        r = exp(r)./(1 + exp(r));
        r(isnan(r)) = 1;
        dd = binornd(1,r);
    
        if rho_fixed ==0
   
            rho_n = normrnd(rho,exp(tow_o));
            if rho_n>0
            k3 = (rho-rho_n)*sum(abs(theta_old)) + sum(dd)*(log(rho_n/2) - log(rho/2));
            k4 = min(1,exp(k3));
            br = binornd(1,k4);
            rho = rho_n*br + (1-br)*rho;
            tow_o = tow_o + ff*(k4 - .3)/(jj^mc_rate);
            end

            if abs(tow_o)>compact
                compact = compact + inc;
                count_reinit = count_reinit +1;
                jj = 1; ff = ff/sqrt(2);
                %Initialization
                rho = rho_in; tow_o = tow_in;% gam is approximation parameter, rho is prior regularizing parameter, 
                theta_old = theta_init; %initializing mcmc for jth set of parameters
                dd = abs(theta_old)>0; 
                res1 = zeros(Niter - burn,p);
                res2= zeros(Niter - burn,p);

            else


                %%Save result
                if jj>burn
                    res1(jj-burn,:) = theta_old;
                    res2(jj-burn,:) = dd;

                end
                %[jj/1000,sum(dd),  rho,count_reinit]


                jj = jj+1;
            end
            if isnan(rho)
                break;
            end
    
        else
            if jj>burn
                res1(jj-burn,:) = theta_old;
                res2(jj-burn,:) = dd;

            end
            jj = jj+1;
        end
 end
fprintf('\n ending %i. \n',j)

            
             Res1(j,:) = theta_init;
             
             Res2(j,:) = d_init;
        thet_est(j,:,:) = res1';
        del_est(j,:,:) = res2';
    



end
%clearvars -except block burn ct Niter p m n Res1 Res2 Rg step thet_est del_est X q sam

%save('result_ispot','del_est','thet_est','Res1','Res2');
diary off;
delete(myPool);



end
