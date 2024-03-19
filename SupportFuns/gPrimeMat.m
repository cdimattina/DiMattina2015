% function y = gPrimeMat(xAug,thetaMat)
% gPrimeMat.m 
% by
% Christopher DiMattina
% Florida Gulf Coast University
%
% Created: 	3-06-2014
%           4-28-2014
%
% Inputs: 
% 			xAug     : N x (Dx + 1) matrix of augmented x values
%           thetaMat : M x (Dx + 1 = Dt) matrix of theta values to evaluate
%                      each x value at
%
% Outputs   y        : gPrime evaluated at every (x,theta) pair
%

function Y = gPrimeMat(xAug,thetaMat)
    U       = xAug*thetaMat';
    Y       = exp(-U)./((1 + exp(-U)).^2); 
end

