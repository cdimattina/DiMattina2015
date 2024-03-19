% function y = gMat(xAug,thetaMat)
% gMat.m 
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
% Outputs:
%			Y	  	 : N x M matrix of P(R = 1 | x, theta )
%                      with
%
%
%


function Y = gMat(xAug,thetaMat)
    U       = xAug*thetaMat';
    Y  = 1./(1 + exp(-U)); 
end
