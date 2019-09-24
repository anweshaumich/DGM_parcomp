function[theta_conf] = confidence_edge(theta,delta,thin)
theta = theta(:,:,1:thin:end);
delta = delta(:,:,1:thin:end);
prob = mean(delta,3); edge = (prob +prob')/2;
pars_p = size(delta,1);
psq = pars_p^2;
est2 = zeros(psq,5); ind = 1;
    for i = 1:pars_p
        for j = 1:i
            if i==j
            val = theta(i,j,:);
            est2(ind,1) = i; est2(ind,2) = j;
            est2(ind,3) = mean(val);
            est2(ind,4:5) = quantile(val,[.025,.975]);
            ind = ind +1;
            else
            val = theta(i,j,:);
            est2(ind,1) = i; est2(ind,2) = j;
            est2(ind,3) = mean(val);
            est2(ind,4:5) = quantile(val,[.025,.975]);
            val = theta(j,i,:);
            est2(ind+1,1) = j; est2(ind+1,2) = i;
            est2(ind+1,3) = mean(val);
            est2(ind+1,4:5) = quantile(val,[.025,.975]);
            ind = ind+2;
            end
        end
    end
    ind = 1; est3 = zeros((psq+pars_p)/2,6); ind2 = 1;
    while ind<=psq
        if  est2(ind,1)==est2(ind,2)
            est3(ind2,1:5) = est2(ind,:);
            est3(ind2,6) = edge(est2(ind,1),est2(ind,2));
            ind = ind+1; ind2 = ind2+1;
        else
            est3(ind2,1:2) = est2(ind,1:2);
            est3(ind2,3) = (est2(ind,3)+est2(ind+1,3))/2;
            est3(ind2,4) = min(est2(ind,4),est2(ind+1,4));
            est3(ind2,5) = max(est2(ind,5),est2(ind+1,5));
            est3(ind2,6) = edge(est2(ind,1),est2(ind,2));
            ind = ind+2; ind2 = ind2+1;
        end
    end 
    theta_conf = sortrows(est3,3);
    
    
end


