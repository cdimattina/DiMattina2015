% CalcEnt.m
% by
% Christopher DiMattina
% Florida Gulf Coast University
%
% Description:  This program calculates the differential entropy of each
%               distribution in the M rows of the matrix P and returns an Mx1 vector
%               of entropies
%
function H = CalcEnt(P)
    
    logP                     = log(P); 
    badind                   = find( (logP==NaN) | (logP==(-Inf)) );
    logP(badind)             = 0; 
    H                        = sum(-P.*logP,2);
end
