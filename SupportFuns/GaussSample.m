% GaussSample.m
% by
% Christopher DiMattina
% Florida Gulf Coast Univeristy
%
% Description: This program simulates nSamples points from a Gaussian density 
%              having mean mu and covariance matrix Sigma
%
%              Samples are in rows, variables in columns
%
%              mu - row vector
%

function Y = GaussSample(nSamples,mu,Sigma)
    ndim = length(mu);
    X    = randn(ndim,nSamples);
    Y    = repmat(mu,nSamples,1) + (sqrtm(Sigma)*X)';
end