% MakeLookup2DNew.m
% by
% Christopher DiMattina, PhD 
% Florida Gulf Coast University
%
% Brief Description: This program generates a Lookup table of D-optimal  
%                    designs for the 2D psychometric function. Particles
%                    are sampled from a Gaussian Prior
%
% Copyright (C) 2014 Christopher DiMattina
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 

function MakeLookup2DNew

%% Number of particles
NParticles          = 5000; 

%% Establish parameters of true underlying observer
% LogBetaTrue         = 0; 
% BetaTrue            = 10^LogBetaTrue; 
% LambdaTrue          = 3; 
% ThetaTrue           = [ -LambdaTrue*BetaTrue , BetaTrue, BetaTrue, BetaTrue ];

% ThetaTrue = [-3,1,1,0.4];
% options             = optimset('Display','iter','GradObj','off');
% warning off;
% 
% d0                  = 2*rand(8,1);  % random starting point
% dopt                = fmincon(@(d) NegLogDetIF(d,ThetaTrue),d0,[],[],[],[],zeros(8,1),5*ones(8,1),[],options);
% 
% x1vec = dopt(1:2:end);
% x2vec = dopt(2:2:end);
% 
% figure(2); 
% subplot(2,2,4);
% Make2DContourPlot; hold on;
% plot(x1vec,x2vec,'k.','MarkerSize',20);


%% Define priors on LogBeta, Lambda

% Draw NParticles from prior
% Draw LogBeta from a prior with mean 0, SD 0.5
LogBetaVec = 0.5*randn(NParticles,1);
BetaVec    = 10.^LogBetaVec; BetaVec(find(BetaVec > 1) ) = 1;

% Draw Lambda  from a prior with mean 2, SD 1
LambdaVec  = 2 + randn(NParticles,1);

theta0vec  = -(LambdaVec.*BetaVec);
theta1vec  = BetaVec;

LogBetaVec = 0.5*randn(NParticles,1);
BetaVec    = 10.^LogBetaVec; BetaVec(find(BetaVec > 1) ) = 1;
theta2vec  = BetaVec;

LogBetaVec = 0.5*randn(NParticles,1);
BetaVec    = 10.^LogBetaVec; BetaVec(find(BetaVec > 1) ) = 1;
theta3vec  = BetaVec;

thetaMat   = [theta0vec, theta1vec, theta2vec, theta3vec];

%% Define design matrix
xDim       = 2; 
thetaDim   = 4;

desSamp    = 4; % only compute designs for 1/4 of particles
desInd     = 1:desSamp:NParticles;
nDes       = length(desInd); 

xDes       = zeros(nDes,xDim*thetaDim,'single');

options    = optimset('Display','off','GradObj','off','MaxFunEvals',200);
warning off; 

%% Compute D-optimal deisgns
for i=1:nDes
 
   % disp(thetaMat(desInd(i),:));
    
    d0                  = 2*rand(8,1);  % random starting point
    dopt                = fmincon(@(d) NegLogDetIF(d,thetaMat(desInd(i),:)),d0,[],[],[],[],zeros(8,1),5*ones(8,1),[],options);
    xDes(i,:)           = dopt; 
    % Sequentially construct a D-optimal design using a greedy method
    % x1, x2 in [0 5] range only
    
    % Start with identity matrix and optimize the determinant of 
    % 4 rank-1 perturbations
    
%     F = eye(thetaDim); 
%     for j = 1:thetaDim
%         x0   = 2*rand(2,1);  % random starting point
%         xopt = fmincon(@(x) NegQuadFormOpt(x,thetaMat(desInd(i),:),inv(F)),x0,[],[],[],[],[0;0],[5;5],[],options);
%         F    = F + Sens(xopt,thetaMat(desInd(i),:))*Sens(xopt,thetaMat(desInd(i),:))';  % rank-1 correction
%         xDes(i, (2*(j-1) + 1):(2*(j-1) + 2)) = xopt';
%     end
    
    if(rem(i,10)==0)
        disp(sprintf('Constructed design %d of %d\n',i,nDes));
    end
    % disp(i)
end

%% Compute utility matrix
U = zeros(nDes,NParticles,'single');

for i=1:nDes
   tic,
   for j=1:NParticles
       U(i,j) = det(IF(xDes(i,:),thetaMat(j,:))); 
   end
   t=toc;
   disp(sprintf('Filled out matrix row %d of %d. Time = %.2f\n seconds',i,nDes,t));
end

%% Save outputs
fNameOut = sprintf('2DLookupOut_%d.mat',NParticles);
save(fNameOut,'thetaMat','xDes','U');

end


function v = Sens(x,theta)
    yAug = [ 1 ; x(1) ; x(2) ; x(1)*x(2) ];
    v    = sqrt(gPrimeMat(yAug',theta))*yAug;
end

function f = NegQuadFormOpt(x,theta,A)
    f    = -Sens(x,theta)'*A*Sens(x,theta);
end

function F = IF(d,theta)
    F = zeros(length(theta),length(theta)); 
    for i = 1:4
       x = d((2*(i-1) + 1):(2*(i-1) + 2));
       F = F + Sens(x,theta)*Sens(x,theta)'; 
    end
end

function L = NegLogDetIF(d,theta)
    L = -1*log(det(IF(d,theta)));
end

