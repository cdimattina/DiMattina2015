function q=getprefix(s)
%----------------------
% getprefix.m
% by 
% Christopher DiMattina
% Case Western Reserve University
%
% Description:  This function extracts the prefix of the input string
%
% Inputs:       s - input string
%
% Outputs:      q - output string
%

    ind = find(s=='.'); 
    q   = s(1:(ind-1));
end