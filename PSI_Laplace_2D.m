% PSI_Laplace_2D.m
% by
% Christopher DiMattina, PhD 
% Florida Gulf Coast University
%
% Brief Description: This program implements the PSI method for a 
%                    2-D psychophysical experiment using the Laplace
%                    implementation
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
nSamp                   = 50; % Number of samples to take when simulating from Gaussian
                              % reducing this number increases speed but decreases accuracy
nmc                     = 1; 
Niter                   = 100;
sampMethod              = 'OED';    % set to 'IID' for method of constant stimuli
sysID                   = computer; if(strcmp(sysID,'PCWIN')) sysID = 'PCWIN32'; end
showPlots               = 1;        % set to 0 to supress graphical output


  
%% Establish parameters of true underlying observer
LogBetaTrue             = 0; 
BetaTrue                = 10^LogBetaTrue; 
LambdaTrue              = 3; 
ThetaTrue               = [ -LambdaTrue*BetaTrue , BetaTrue, BetaTrue, BetaTrue ];
D                       = length(ThetaTrue);

thetaMuPrior            = [ -1 , 0.5, 0.5, 0.5 ];
S2InvPrior              = inv(2*eye(4)); 
 
%% Set simulation constants 
% Grid density
NXGrid                  = 21;
xMin                    = 0;
xMax                    = 5; 

xVec                    = linspace(xMin,xMax,NXGrid)';
xMat                    = FactorMat(xVec,xVec);

xAug                    = ones(size(xMat,1),4);
xAug(:,2:3)             = xMat;
xAug(:,4)               = xMat(:,1).*xMat(:,2);

nx                      = size(xMat,1); 

 
%% Declare Matrices which will store the final outputs of the simulation
Hmat  = zeros(nmc,Niter+1);         % Entropies
Emat  = zeros(nmc,Niter+1);         % EMSE
Tmat  = zeros(nmc,Niter);           % Time
X1mat = zeros(nmc,Niter);           % Stimuli - x1
X2mat = zeros(nmc,Niter);           % Stimuli - x2
EVmat = zeros(nmc,4);               % Vector of final parameter estimates - Expectation Value
MLmat = zeros(nmc,4);               % Vector of final parameter estimates - Max Likelihood
 
 
%% Monte Carlo loop 
for j=1:nmc
    
    % Initialize posterior mode and inverse covariance
    % Compute current entropy and current error (EMS error)
    DPlot                  = NaN*ones(Niter+1,1);
    HPlot                  = NaN*ones(Niter+1,1);
    EPlot                  = NaN*ones(Niter+1,1);
    
    thetaHatCur             = GaussSample(1,thetaMuPrior,inv(S2InvPrior));        % sample from prior
    S2HatInvCur             = inv(2*eye(4));
    
    DPlot(1)               = 1./det(S2HatInvCur);
    HPlot(1)               = (1/2)*log(DPlot(1));
    EPlot(1)               = sum((thetaHatCur - ThetaTrue).^2);
    
    if(showPlots)
        close all;
        figure(1); 
        subplot(2,2,1); plot(0:Niter,HPlot); axis square; ylabel('Entropy'); xlabel('Iteration');
        subplot(2,2,2); plot(0:Niter,EPlot); axis square; ylabel('Error'); xlabel('Iteration');
        pause(0.1);
    end
    
    tElapsedVec = zeros(Niter,1);
    Xall        = zeros(Niter,D);
    Rall        = zeros(Niter,1);
    
    %% Experiment loop
    for i = 1:Niter
    
        Hexp        = zeros(nx,1);
        thetaSamp   = GaussSample(nSamp,thetaHatCur,inv(S2HatInvCur));
    
        % for every stimulus x in xAug
    
        if(strcmp(sampMethod,'OED'))
    
           % tic,
            for k = 1:nx
                xNext = xAug(k,:);
                
                tic,
                pR1x  = mean(gMat(xNext,thetaSamp));
                pR0x  = 1-pR1x;
                
%               C1    = exp(-thetaHatCur*xNext')*((gMat(xNext,thetaHatCur)).^2)*(xNext'*xNext);
                C0    = gPrimeMat(xNext,thetaHatCur)*(xNext'*xNext);
                C1    = C0; 
                
                S2HatInvNext1 = C1 + S2HatInvCur;
                S2HatInvNext0 = C0 + S2HatInvCur;
                
                HR1x  = -0.5*log(det(S2HatInvNext1));
                HR0x  = -0.5*log(det(S2HatInvNext0));
                
                Hexp(k) = HR1x*pR1x + HR0x*pR0x;

%                 et=toc;
%                 disp(sprintf('time = %.2f msec',1000*et));
                
            end % for k
   
                % Present chosen stimulus + record response
                [dum,ind] = min(Hexp);
                xNext     = xAug(ind,:);
                
                tElapsedVec(i) = toc;
                
        elseif(strcmp(sampMethod,'IID'))
                ind         = randi(nx);
                xNext       = xAug(ind,:);
        else
                error('Method not supported!');
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
    
   
        % Plot
        if(showPlots)
            figure(1);
            subplot(2,2,1); plot(0:Niter,HPlot); axis square; ylabel('Entropy'); xlabel('Iteration');
            subplot(2,2,2); plot(0:Niter,EPlot); axis square; ylabel('Error'); xlabel('Iteration');
            subplot(2,2,4); hold on; plot(xNext(2),xNext(3),'k.'); axis square; axis([xMin xMax xMin xMax]);
            pause(0.1);
        end
    
    end % for i = 1:Niter
    
    
    % optimize posterior to get final estimate
    options     = optimset('GradObj','on','Display','off');
    thetaHatML  = fminunc(@(theta) NegLogLikelihoodWithGrad(Xall,Rall,theta),double(thetaHatCur),options);
    
    EVmat(j,:)  = thetaHatCur;
    MLmat(j,:)  = thetaHatML;
    
    disp(sprintf('Iteration: %d',j));
    
    % save everything to output matrices
    Tmat(j,:) = tElapsedVec';
    Hmat(j,:) = HPlot';
    Emat(j,:) = EPlot';
    X1mat(j,:) = Xall(:,2)';
    X2mat(j,:) = Xall(:,3)';
    
    
end % for j=1:nmc
 
%% Save output variables
savstr = strcat('PSI_Laplace_2D_Out_',sampMethod,'_',sysID,'.mat'); 
save(savstr,'Hmat','Emat','Tmat','X1mat','X2mat','EVmat','MLmat','nmc','Niter','sampMethod','sysID','ThetaTrue'); 

clear all;
