% function L = LogLikelihoodWithGrad(xAug,R,theta)
% LogLikelihoodWithGrad.m
% by
% Christopher DiMattina
% Florida Gulf Coast University
%
% Inputs:   xAug     : N x (Dx + 1) matrix of augmented x values
%           R        : N x 1        vector of responses
%           theta    : 1 x (Dx + 1 = Dt) matrix of theta values to evaluate
%                      each x value at (row vector!)
%          
% Outputs:  L        : log-likelihood of data     
%			g        : log-likelihood gradient
%           
%           

function [L,g] = LogLikelihoodWithGrad(xAug,R,theta)

    if(size(theta,1) > size(theta,2))
       error('Input argument THETA must be a row vector!'); 
    end

    ind1 = find(R==1);
	ind0 = find(R==0); 
	
    X1   = xAug(ind1,:); 
    X0   = xAug(ind0,:); 

    L1   =  log(     gMat( X1, theta ) ); 
    L0   =  log( 1 - gMat( X0, theta ) );
    L    =  sum(L1) + sum(L0);
   
    D1   = diag(gPrimeMat(X1,theta)./gMat(X1,theta));
    D0   = diag(-gPrimeMat(X0,theta)./(1-gMat(X0,theta)));
    g    = sum(D1*X1,1) + sum(D0*X0,1);
    
end


