% PSI_Laplace_1D.m
% by
% Christopher DiMattina, PhD 
% Florida Gulf Coast University
%
% Brief Description: This program implements the Laplace approximation
%                    method for a 1-D psychometric function
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

clear all;
addpath ./SupportFuns/
addpath ./MatFiles/

%% Set simulation inputs
nSamp                   = 100; 

nmc                     = 1; 
Niter                   = 100;
sampMethod              = 'OED';    % set to 'IID' for method of constant stimuli
sysID                   = computer; if(strcmp(sysID,'PCWIN')) sysID = 'PCWIN32'; end
showPlots               = 1;        % set to 0 to supress graphical output
  
%% Establish parameters of true underlying observer
LogBetaTrue         	= 0; 
BetaTrue            	= 10^LogBetaTrue; 
LambdaTrue          	= 0; 
ThetaTrue           	= [ -LambdaTrue*BetaTrue , BetaTrue ];
D                       = length(ThetaTrue);

thetaMuPrior            = [-1, 1.5];
S2InvPrior              = inv(eye(2)); 

%% Set simulation constants
% Grid density
NXGrid                  = 51;
NThetaGrid              = 51; 

% Grid limits
LogBetaMin              = -1;
LogBetaMax              = 1;

LambdaMin               = -2; 
LambdaMax               = 2;

xMin                    = -7;
xMax                    = 7;
xVec                    = linspace(xMin,xMax,NXGrid)';
xAug                    = [ ones(size(xVec)), xVec ];
nx                      = size(xVec,1); 

LogBetaVec              = linspace(LogBetaMin,LogBetaMax,NThetaGrid);
BetaVec                 = 10.^LogBetaVec; 
LambdaVec               = linspace(LambdaMin,LambdaMax,NThetaGrid); 

theta0vec               = -(LambdaVec.*BetaVec);
theta1vec               = BetaVec;

thetaMat                = FactorMat(theta0vec,theta1vec);
 
%% Declare Matrices which will store the final outputs of the simulation
Hmat  = zeros(nmc,Niter+1);         % Entropies
Emat  = zeros(nmc,Niter+1);         % EMSE
Tmat  = zeros(nmc,Niter);           % Time
Xmat  = zeros(nmc,Niter);           % Stimuli

EVmat = zeros(nmc,2);               % Vector of final parameter estimates - Expectation Value
MLmat = zeros(nmc,2);               % Vector of final parameter estimates - Max Likelihood

%% Monte Carlo loop
for j=1:nmc

    DPlot                   = NaN*ones(Niter+1,1);
    HPlot                   = NaN*ones(Niter+1,1);
    EPlot                   = NaN*ones(Niter+1,1);
    
    % Pick random starting point from prior
    thetaHatCur             = GaussSample(1,thetaMuPrior,inv(S2InvPrior));        % sample from prior
    S2HatInvCur             = S2InvPrior;
    
    % Initilize plot variables
    DPlot(1)                = 1./det(S2HatInvCur);
    HPlot(1)                = (1/2)*log(DPlot(1));
    EPlot(1)                = sum((thetaHatCur - ThetaTrue).^2);
    
    % Plot current posterior density on the grid thetaMat
    PTheta                  = GaussEval(thetaMat,thetaHatCur,S2HatInvCur);
    
    if(showPlots)
        close all; 
        figure(1);
        subplot(2,2,1); plot(0:Niter,HPlot); axis square; ylabel('Entropy'); xlabel('Iteration');
        subplot(2,2,2); plot(0:Niter,EPlot); axis square; ylabel('Error'); xlabel('Iteration');
        subplot(2,2,3); plot(ThetaTrue(1),ThetaTrue(2),'bo','MarkerSize',10,'LineWidth',3); hold on; axis square;
                        plot(thetaHatCur(1),thetaHatCur(2),'k.');
        subplot(2,2,4); imagesc(reshape(PTheta,NThetaGrid,NThetaGrid)); axis square;
        pause(0.1);
    end
    
    tElapsedVec = zeros(Niter,1);
    Xall        = zeros(Niter,2);
    Rall        = zeros(Niter,1);
    
    %% Experiment loop
    for i = 1:Niter
    
        Hexp        = zeros(nx,1);
        thetaSamp   = GaussSample(nSamp,thetaHatCur,inv(S2HatInvCur));
        
        % for every stimulus x in xAug
        if(strcmp(sampMethod,'OED') )
            tic,
            for k = 1:nx
                xNext = xAug(k,:);
                
                pR1x  = mean(gMat(xNext,thetaSamp));
                pR0x  = 1-pR1x;
                
              
                %C1_old    = exp(-thetaHatCur*xNext')*((gMat(xNext,thetaHatCur)).^2)*(xNext'*xNext);
                C1    = gPrimeMat(xNext,thetaHatCur)*(xNext'*xNext);
                C0    = gPrimeMat(xNext,thetaHatCur)*(xNext'*xNext);
                
                
                S2HatInvNext1 = C1 + S2HatInvCur;
                S2HatInvNext0 = C0 + S2HatInvCur;
                
                HR1x  = -0.5*log(det(S2HatInvNext1));
                HR0x  = -0.5*log(det(S2HatInvNext0));
                
                Hexp(k) = HR1x*pR1x + HR0x*pR0x;
               
                
            end % end for k
            
            % Present chosen stimulus + record response
            [dum,ind] = min(Hexp);
            xNext       = xAug(ind,:);
            
            tElapsedVec(i) = toc;
    
        elseif(strcmp(sampMethod,'IID'))
            ind         = randi(nx);
            xNext       = xAug(ind,:);
        else
            error;
        end
        
        if( rand <= gMat(xNext,ThetaTrue))
            R = 1;
        else
            R = 0;
        end
        Rall(i)     = R;
        Xall(i,:)   = xNext;
    
        XsoFar  = Xall((1:i),:);
        RsoFar  = Rall(1:i);
    
        X1      = XsoFar(find(RsoFar==1),:);
        X0      = XsoFar(find(RsoFar==0),:);
    
    
        % Update posterior density - find likelihood maximum thetaHatCur and
        % perform a rank-1 correction of the S2HatCur at that point
    
        % Maximize likelihood
        options      = optimset('GradObj','on','Display','off','TolFun',10e-5,'TolX',10e-5);
        thetaHatCur  = fminunc(@(theta) NegLogPosteriorWithGradFast(X1,X0,D,thetaMuPrior,S2InvPrior,theta),thetaHatCur + 0.1*randn(size(thetaHatCur)),options);
    
        % Rank-1 correction of S2HatCur
        if(R==1)
            C = exp(-thetaHatCur*xNext')*((gMat(xNext,thetaHatCur)).^2)*(xNext'*xNext);
        else
            C = gPrimeMat(xNext,thetaHatCur)*(xNext'*xNext);
        end
    
        S2HatInvCur = C + S2HatInvCur;
    
        thetaSamp   = GaussSample(nSamp,thetaHatCur,inv(S2HatInvCur));
        thetaHatCur = mean(thetaSamp);
    
        DPlot(i+1)               = 1./det(S2HatInvCur);
        HPlot(i+1)               = (1/2)*log(DPlot(i+1));
        EPlot(i+1)               = sum((thetaHatCur - ThetaTrue).^2);
    
        % Plot current posterior density on the grid thetaMat
        PTheta      = GaussEval(thetaMat,thetaHatCur,S2HatInvCur);
    
        % Plot
        if(showPlots)
            figure(1);
            subplot(2,2,1); plot(0:Niter,HPlot); axis square; ylabel('Entropy'); xlabel('Iteration');
            subplot(2,2,2); plot(0:Niter,EPlot); axis square; ylabel('Error'); xlabel('Iteration');
            subplot(2,2,3); plot(ThetaTrue(1),ThetaTrue(2),'bo','MarkerSize',10,'LineWidth',3); hold on; axis square;
            plot(thetaHatCur(1),thetaHatCur(2),'k.');
            subplot(2,2,4); imagesc(reshape(PTheta,NThetaGrid,NThetaGrid)); axis square; axis off; 
            pause(0.1);
        end
    
    end % for i = 1:Niter
      
    Tmat(j,:) = tElapsedVec';
    Hmat(j,:) = HPlot';
    Emat(j,:) = EPlot';
    Xmat(j,:) = Xall(:,2)';

    % obtain final estimate by maximizing likelihood function starting at
    % thetaHatCur
    options     = optimset('GradObj','on','Display','off');
    thetaHatML  = fminunc(@(theta) NegLogLikelihoodWithGrad(Xall,Rall,theta),thetaHatCur,options);
    
    EVmat(j,:) = thetaHatCur;
    MLmat(j,:) = thetaHatML;
    
    if(showPlots)
        figure(1); subplot(2,2,3);  plot(thetaHatML(1),thetaHatML(2),'rd','MarkerSize',10,'LineWidth',3);
    end
    
    disp(sprintf('Iteration: %d',j));
    
 end % for j=1:nmc

%% Save outputs
% save output variables
savstr = strcat('PSI_Laplace_1D_Out_',sampMethod,'_',sysID,'.mat'); 
save(savstr,'Hmat','Emat','Tmat','Xmat','EVmat','MLmat','nmc','Niter','sampMethod','sysID','ThetaTrue'); 

clear all;
