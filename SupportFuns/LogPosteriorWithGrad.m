% LogPosteriorWithGrad.m
% by
% Christopher DiMattina
% Florida Gulf Coast University
%
%
%
%

function [L,g] = LogPosteriorWithGrad(xAug,R,theta, muPrior, SigmaInvPrior)
    [LL,gL] = LogLikelihoodWithGrad(xAug,R,theta);
    L       = LL + -0.5*(theta-muPrior)*SigmaInvPrior*(theta-muPrior)';
    g       = gL + -(theta-muPrior)*SigmaInvPrior;
end

