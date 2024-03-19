% function L = NegLogPosteriorWithGradFast(xAug, R, thetaHatCur, S2HatCur, theta)
% NegLogPosteriorWithGrad.m
% by
% Christopher DiMattina
% Florida Gulf Coast University
%
% Inputs:   xAug     : N x (Dx + 1) matrix of augmented x values
%           R        : N x 1        vector of responses
%
%
%           theta    : 1 x (Dx + 1 = Dt) matrix of theta values to evaluate
%                      each x value at (row vector!)
%          
% Outputs:  L        : log-posterior of data        (scalar)
%			g        : gradient                     (row vector!)
%           
%           
% This version has been modified to be faster

function [L,g] = NegLogPosteriorWithHessianFast(X1,X0,D, thetaMuPrior, S2InvPrior, theta)
  
    L1   =  log(     gMat( X1, theta ) ); 
    L0   =  log( 1 - gMat( X0, theta ) );
    L    =  -(sum(L1) + sum(L0)) + (D/2)*log(2*pi) + 0.5*log(1./det(S2InvPrior)) + ...
             0.5*(theta-thetaMuPrior)*S2InvPrior*((theta-thetaMuPrior)');
   
    D1   = diag(gPrimeMat(X1,theta)./gMat(X1,theta));
    D0   = diag(-gPrimeMat(X0,theta)./(1-gMat(X0,theta)));
    g    = -(sum(D1*X1,1) + sum(D0*X0,1)) + (theta-thetaMuPrior)*S2InvPrior;
    
end


