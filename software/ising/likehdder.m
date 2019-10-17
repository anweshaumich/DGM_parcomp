function  l  = likehdder(th,del, Rg,Z,j)
%derivative of negative of loglikelihood
%   Z is the data matrix, th is the jth row of the parameter matrix
%   associated with jth conditional likelihood, j is the index of the
%   conditional likelihood to be used, m is the number of state spaces

z = Z;  
x = Z(:,j); % identifies the jth variable in n observations
z(:,j) = 1;
z = z(:,logical(del));
l = zeros(1,length(th));
p = size(z,2);
X = repmat(x,1,p);
c = X.*z; % creates the derivative of the numerator of the likelihood for each observation
k1 = sum(c);

t = th(logical(del));
% creating derivative of the proportional constant.
exp_mat = exp((z*t')*Rg);
k2 = sum(exp_mat,2);
k4 = exp_mat(:,2);
k6 = ( (k4./k2)' )*z;


l(logical(del)) = -k1 + k6; % returns row vector

end

