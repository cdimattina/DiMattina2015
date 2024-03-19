% GaussEval.m
% by
% Christopher DiMattina
% Florida Gulf Coast Univeristy
%
% Description: This program evaluates the multivariate Gaussian density at
%              the points contained in theta
%
%              thetaMat has observations in rows, variables in columns
%              mu is a row vector

function y = GaussEval(thetaMat,mu,SigmaInv)
    
    D = length(mu);                     % Dimensionality
    N = size(thetaMat,1);               % Number of samples
    M = repmat(mu,N,1);
    X = (thetaMat - M) ;
    Z = ((2*pi)^(D/2))*((1/det(SigmaInv))^(0.5));
    
    y = (1/Z)*diag(exp(-0.5*X*SigmaInv*(X')));  
end

