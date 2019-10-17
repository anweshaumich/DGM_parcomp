
function [theta_now,mse] = initial(X,gamn,h_matj,H_mat,Rg, MAX_ITER,TOL)
% load('data.mat')
[n,p] = size(X);
m = length(Rg);

s = .15:.02:.5;

rho = s*sqrt(n*log(p));lrho = length(rho);
thet = zeros(lrho,p);mse = zeros(1,lrho);

% MAX_ITER = 500;  TOL = .0001;
%beta = .5;
for i = 1:lrho
    
    theta_now = zeros(1,p) ; hprox = zeros(1,MAX_ITER);dd = ones(1,p);
    for k = 1:MAX_ITER
        gam = gamn;
        th = theta_now - gam * likehdder(theta_now,dd, Rg, h_matj, H_mat);
        th1 = th + gam * rho(i);
        th2 = th - gam * rho(i);
        theta_now =  (th1<0).*th1 + (th2>0).*th2;
        hprox(k) = norm(theta_now, 2);
        if k > 1 && abs(hprox(k) - hprox(k-1)) < TOL
            
            break;
        end
    end
    mse(i) = likhd(theta_now,dd, Rg,h_matj,H_mat) + (log(n) + 2*20*log(p))*nnz(theta_now);
    thet(i,:) = theta_now;
    
end
ind = find(mse==min(mse));
theta_now = thet(ind(1),:);
mse = mse(ind(1));
end


