% FactorMat.m
% by
% Christopher DiMattina
% Florida Gulf Coast University
%
%
%

function M = FactorMat(varargin)
    nVec    = nargin;
    dimVec  = zeros(nVec,1);
    
    for i=1:nVec
       dimVec(i) = length(varargin{i});  
    end
        
    nRows   = prod(dimVec);
    nCols   = nVec;
    
    M       = zeros(nRows,nCols);
    
    outstr  = MakeOutString(nVec);
    instr   = MakeInString(nVec);
    evalstr = strcat(outstr,'=',instr);
    
    eval(evalstr)
    
    for i=1:nCols
        outstr = sprintf('M(:,%d)',i);
        instr  = sprintf('reshape(X%d,nRows,1);',i);
        eval(strcat(outstr,'=',instr));
    end
    
end

function s = MakeOutString(n)

    s = '[';
    for i=1:(n-1)
       s = strcat(s,sprintf('X%d',i),','); 
    end
    s = strcat(s,sprintf('X%d',n),']');
    
end

function s = MakeInString(n)
    
    s='ndgrid(';
    for i=1:(n-1)
       s = strcat(s,sprintf('varargin{%d}',i),','); 
    end
     s = strcat(s,sprintf('varargin{%d}',n),');');

end