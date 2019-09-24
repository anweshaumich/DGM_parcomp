
function [theta_now,mse2] = initial(X,gamn,j,Rg, MAX_ITER,TOL)

[n,p] = size(X);
m = length(Rg);

s = .05:.05:0.5;
rho = s*sqrt(n*log(p));lrho = length(rho);
thet = zeros(lrho,p);
mse2 = zeros(1,lrho);


for i = 1:lrho
    
    theta_now = zeros(1,p) ; hprox = zeros(1,MAX_ITER);
    for k = 1:MAX_ITER
        gam = gamn;

        th = theta_now - gam * likehdder(theta_now, Rg, X,j);
        th1 = th + gam * rho(i);
        th2 = th - gam * rho(i);
        theta_now =  (th1<0).*th1 + (th2>0).*th2;
        hprox(k) = norm(theta_now, 2);
        if k > 1 && abs(hprox(k) - hprox(k-1)) < TOL
            
            break;
        end
    end

    mse2(i) = likhd(theta_now, Rg,X,j) + (log(n) + 2*20*log(p))*nnz(theta_now);

    thet(i,:) = theta_now;
    
end
ind = find(mse2 == min(mse2));
theta_now = thet(ind(1),:);
mse2 = mse2(ind(1));
end




