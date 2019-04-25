function y = likhd(th, Rg,Z, j)
% jth negative conditional loglikelihood
% Z is the data matrix, th is the jth row of the parameter matrix
%   associated with jth conditional likelihood, j is the index of the
%   conditional likelihood to be used, m is the number of state spaces
z = Z;
x = Z(:,j); % identifies the jth variable in n observations
z(:,j) = 1;
mat_1 = z * th';
k1 = sum(mat_1.*x);
k2 = sum(log( sum( exp(mat_1 * Rg), 2) ) );
y = -k1 + k2;


end

