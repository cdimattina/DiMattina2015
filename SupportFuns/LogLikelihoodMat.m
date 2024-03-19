% function [L1,L0] = LogLikelihoodMat(theta,x,r)
% LogLikelihoodMat.m
% by
% Christopher DiMattina
% Florida Gulf Coast University
%
% Inputs:
% 			thetaMat 	= Mx2 matrix of logistic function parameters
% 						  [logb, lambda] (log slope, threshold)
% 			x     	 	= Nx1 vector of inputs
% 			r     	 	= Nx1 vector of binary responses (1,0)
% Outputs:
%			L     		= NxM log-likelihood of data (in rows) given parameters (columns)
%

function [L1,L0] = LogLikelihoodMat(thetaMat,x,r)
	ind1 = find(r==1);
	ind0 = find(r==0); 

	L1   =  log(    gMat( x(ind1), thetaMat ) ) ; 
    L0   =  log(1 - gMat( x(ind0), thetaMat ) ) ;
end
