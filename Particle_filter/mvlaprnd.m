function R = mvlaprnd(d,MU,SIGMA)
% Purpose:  
%       Generate multivariate Laplace random numbers

% Inputs: 
%       d       - dimension of the random vector 
%       MU      - mean vector
%       SIGMA   - covariance matrix

% Output:
%       R - Multivariate Laplacian random number
 
% Reference:
%       T. Eltoft et. al.
%       "On the Multivariate Laplace Distribution"
%       IEEE Signal Processing Letters, Vol. 13, No. 5, May 2006


if nargin ~= 3
    error('Input must have 3 arguments, neither more nor less!');
end
if ~all(eig(SIGMA)> 0)
    error('SIGMA is not positive definite')
end

lambda  = det(SIGMA);
GAMMA   = SIGMA/lambda;
z       = exprnd(lambda);
MU_G    = zeros(d,1); 
SIGMA_G = eye(d);
X       = mvnrnd(MU_G,SIGMA_G)';
R       = MU + sqrt(z)*sqrtm(GAMMA)*X;

end