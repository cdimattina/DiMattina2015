% PSI_Lookup_1D.m
% by
% Christopher DiMattina, PhD 
% Florida Gulf Coast University
%
% Brief Description: This program implements the look-up method for a 1-D
%                    psychometric function
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
nmc                 = 1; 
Niter               = 25;
sampMethod          = 'OED';    % set to 'IID' for method of constant stimuli
sysID               = computer; if(strcmp(sysID,'PCWIN')) sysID = 'PCWIN32'; end
showPlots           = 1;        % set to 0 to supress graphical output

%% Establish parameters of true underlying observer
LogBetaTrue         	= 0; 
BetaTrue            	= 10^LogBetaTrue; 
LambdaTrue          	= 0; 
ThetaTrue           	= [ -LambdaTrue*BetaTrue , BetaTrue ];

%% Load look-up table
if(~exist('1DLookupOut.mat'))
    MakeLookup1D;
else
    load('1DLookupOut.mat') 
end

%% Declare Matrices which will store the final outputs of the simulation
Hmat    = zeros(nmc,2*Niter+1);       % Entropies - 2 stimuli are presented each trial
Emat    = zeros(nmc,2*Niter+1);       % EMSE      - 2 stimuli are presented each trial
Tmat    = zeros(nmc,Niter);           % Time
Xmat   = zeros(nmc,2*Niter);          % Stimuli

EVmat   = zeros(nmc,2);               % Vector of final parameter estimates - Expectation Value
MLmat   = zeros(nmc,2);               % Vector of final parameter estimates - Max Likelihood

%% Monte Carlo loop

for j=1:nmc
    %% Initialize experiment    
    % Initialize uniform prior over theta
    [NParticles,NThetaDim]  = size(thetaMat);
    PTheta                  = (1/NParticles)*ones(NParticles,1);
    nDes                    = length(xDes);
    
    NThetaGrid              = sqrt(NParticles);
    
    % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
    % (i.e. given the prior)
    thetaHatCur         = ((thetaMat')*PTheta)';
    HCur                = CalcEnt(PTheta');
    ECur                = norm(thetaHatCur-ThetaTrue);
    
    % Initialize plot variables
    HPlot           = NaN*ones(2*Niter+1,1);  % on each iteration you present 2 stimuli
    EPlot           = NaN*ones(2*Niter+1,1);  % on each iteration you present 2 stimuli
    PosteriorPlot   = reshape(PTheta,NThetaGrid,NThetaGrid);
    HPlot(1)        = HCur;
    EPlot(1)        = ECur;
    
    if(showPlots)
        close all;
        figure(1);
        subplot(2,2,1); plot(0:(2*Niter),HPlot); axis square; ylabel('Entropy'); xlabel('Iteration');
        subplot(2,2,2); plot(0:(2*Niter),EPlot); axis square; ylabel('Error'); xlabel('Iteration');
        subplot(2,2,3); plot(ThetaTrue(1),ThetaTrue(2),'bo','MarkerSize',10,'LineWidth',3); hold on; axis square;
        plot(thetaHatCur(1),thetaHatCur(2),'k.');
        subplot(2,2,4); imagesc(reshape(PTheta,NThetaGrid,NThetaGrid)); axis square;
    end
    
    tElapsedVec = zeros(Niter,1);
    Xall        = zeros(2*Niter,2); % on each iteration you present 2 stimuli
    Rall        = zeros(2*Niter,1); % on each iteration you get 2 responses
    
    %% Experiment Loop
    for i=1:Niter
        % Find and present next design
        if(strcmp(sampMethod,'OED'))
            tic,
            % Find best design
            uExp        = U*PTheta;
            [dum,mxind] = max(uExp);
            xAug1       = [1,xDes(mxind,1)];
            xAug2       = [1,xDes(mxind,2)];
            tElapsedVec(i)    = toc;
        elseif(strcmp(sampMethod,'IID'))
            xAug1       = [1, -7 + 14*rand];
            xAug2       = [1, -7 + 14*rand];
        elseif(strcmp(sampMethod,'OPT'))    
            xAug1       = [1,1.54];
            xAug2       = [1,-1.54];
        end
                
        % Present stimuli + record responses
        if( rand <= gMat(xAug1,ThetaTrue))
            R1 = 1;
        else
            R1 = 0;
        end
        Rall(2*(i-1)+1,:) = R1;
        Xall(2*(i-1)+1,:) = xAug1;
        
        % Update posterior density
        if(R1==1)
            PTheta = (gMat(xAug1,thetaMat)').*PTheta;
        else
            PTheta = (1-gMat(xAug1,thetaMat)').*PTheta;
        end
        PTheta = (1/sum(PTheta))*PTheta;
        
        % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
        thetaHatCur = ((thetaMat')*PTheta)';
        HCur        = CalcEnt(PTheta');
        ECur        = norm(thetaHatCur-ThetaTrue);
        
        HPlot(2*(i-1)+2)          = HCur;
        EPlot(2*(i-1)+2)          = ECur;
        
        if( rand <= gMat(xAug2,ThetaTrue))
            R2 = 1;
        else
            R2 = 0;
        end
        Rall(2*(i-1)+2,:) = R2;
        Xall(2*(i-1)+2,:) = xAug2;
        
        % Update posterior density
        if(R2==1)
            PTheta = (gMat(xAug2,thetaMat)').*PTheta;
        else
            PTheta = ((1-gMat(xAug2,thetaMat))').*PTheta;
        end
        PTheta = (1/sum(PTheta))*PTheta;
        
        % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
        thetaHatCur = ((thetaMat')*PTheta)';
        HCur        = CalcEnt(PTheta');
        ECur        = norm(thetaHatCur-ThetaTrue);
        
        HPlot(2*(i-1)+3)          = HCur;
        EPlot(2*(i-1)+3)          = ECur;
        
        % Plot
        if(showPlots)
            figure(1);
            subplot(2,2,1); plot(0:(2*Niter),HPlot); axis square;
            subplot(2,2,2); plot(0:(2*Niter),EPlot); axis square; ylabel('Error'); xlabel('Iteration');
            subplot(2,2,3);  plot(ThetaTrue(1),ThetaTrue(2),'bo','MarkerSize',10,'LineWidth',3);  hold on; axis square;
            plot(thetaHatCur(1),thetaHatCur(2),'k.');
            subplot(2,2,4); imagesc(reshape(PTheta,NThetaGrid,NThetaGrid)); axis square; axis off;
            pause(0.1);
        end
        
        
    end % for i=1:Niter
    
    Tmat(j,:) = tElapsedVec';
    Hmat(j,:) = HPlot';
    Emat(j,:) = EPlot';
    Xmat(j,:) = Xall(:,2)';

    % end for i = 1:Niter
    % optimize posterior to get final estimate
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
savstr = strcat('PSI_Lookup_1D_Out_',sampMethod,'_',sysID,'.mat'); 
save(savstr,'Hmat','Emat','Tmat','Xmat','EVmat','MLmat','nmc','Niter','sampMethod','sysID','ThetaTrue','NParticles'); 

clear all; 



