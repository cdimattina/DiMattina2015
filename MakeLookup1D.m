% MakeLookup1D.m
% by
% Christopher DiMattina
% Florida Gulf Coast University
% 

function MakeLookup1D

clear all; 
addpath ./SupportFuns/
addpath ./MatFiles/

NThetaGrid      = 51;

% Grid limits
LogBetaMin      = -0.5;
LogBetaMax      = 0.5;

LambdaMin       = -2; 
LambdaMax       = 2;

LogBetaVec      = linspace(LogBetaMin,LogBetaMax,NThetaGrid);
BetaVec         = 10.^LogBetaVec; 
LambdaVec       = linspace(LambdaMin,LambdaMax,NThetaGrid); 

theta0vec       = -(LambdaVec.*BetaVec);
theta1vec       = BetaVec;

thetaMat        = FactorMat(theta0vec, theta1vec);
nTheta          = size(thetaMat,1); 

subSampFact     = 2; 
sampInd         = (1:subSampFact:nTheta)';

thetaDesign     = thetaMat(sampInd,:);
nThetaDesign    = size(thetaDesign,1); 

xMat            = zeros(nThetaDesign,2);

% Compute the optimal design for each stimulus
for i=1:nThetaDesign  
    options     = optimset('GradObj','on','Display','off' ); warning off;
    xopt        = fminunc(@(x) NegDetIF1D(x,thetaDesign(i,:)),randn(1,2),options);
    if(xopt(1) > xopt(2))
        xMat(i,:) = xopt;
    else
        xMat(i,:) = fliplr(xopt);  
    end    
end

disp('done computing optimal designs...');

U = zeros(nThetaDesign,nTheta); 

for i=1:nThetaDesign
    for j=1:nTheta
        U(i,j) = DetIF1D(xMat(i,:),thetaMat(j,:));
    end
end
xDes = xMat; 
fnameOut = '1DLookupOut.mat';
save(fnameOut,'xDes','thetaMat','U');

disp('done computing utility table...');

clear all; 

end

function [f,g] = NegDetIF1D(x,theta)
        [f,g]  = DetIF1D(x,theta);
        f      = -f; 
        g      = -g;
end

function [f,g] = DetIF1D(x,theta)
    f          = gPrimeMat([1,x(1)],theta).*gPrimeMat([1,x(2)],theta).*((x(1) - x(2))^2);
    g          = zeros(1,2); % row vector!
    g(1)       = gDoublePrimeMat([1 x(1)],theta)*theta(2)*gPrimeMat([1 x(2)],theta)*(x(1)-x(2))^2 + ... 
                 gPrimeMat([1 x(1)],theta)*gPrimeMat([1 x(2)],theta)*2*(x(1)-x(2));
    g(2)       =  gDoublePrimeMat([1 x(2)],theta)*theta(2)*gPrimeMat([1 x(1)],theta)*(x(1)-x(2))^2 - ... 
                 gPrimeMat([1 x(1)],theta)*gPrimeMat([1 x(2)],theta)*2*(x(1)-x(2));
end

