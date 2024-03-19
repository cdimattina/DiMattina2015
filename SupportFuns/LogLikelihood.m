% function L = LogLikelihood(xAug,R,theta)
% LogLikelihood.m
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
%			
%           
%           

function L = LogLikelihood(xAug,R,theta)

    if(size(theta,1) > size(theta,2))
       error('Input argument THETA must be a row vector!'); 
    end

    ind1 = find(R==1);
	ind0 = find(R==0); 
	L1   =  log(     gMat( xAug(ind1,:), theta ) ); 
    L0   =  log( 1 - gMat( xAug(ind0,:), theta ) );
    L    = sum(L1) + sum(L0);
    
end


