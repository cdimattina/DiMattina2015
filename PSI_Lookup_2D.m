% PSI_Lookup_2D.m
% by
% Christopher DiMattina, PhD 
% Florida Gulf Coast University
%
% Brief Description: This program implements the PSI method for a 
%                    2-D psychophysical experiment using the Lookup
%                    table implementation
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
NParticles          = 5000;     % Number of supports for evolving posterior
nmc                 = 1; 
Niter               = 25;       % 4 stimuli are presented each iteration
sampMethod          = 'OED';    % set to 'IID' for method of constant stimuli
sysID               = computer; if(strcmp(sysID,'PCWIN')) sysID = 'PCWIN32'; end
showPlots           = 1;        % set to 0 to supress graphical output

Nstim               = 4*Niter;  % 4 stimuli are presented each iteration
xMax                = 5;

%% Establish parameters of true underlying observer
LogBetaTrue         = 0; 
BetaTrue            = 10^LogBetaTrue; 
LambdaTrue          = 3; 
ThetaTrue           = [ -LambdaTrue*BetaTrue , BetaTrue, BetaTrue, BetaTrue ];

%% Load pre-computed lookup-table
load(sprintf('2DLookupOut_%d.mat',NParticles));

%% Clean lookup table to remove NaN
[nrow,ncol] = size(xDes);
gdind = [];
for i=1:nrow
    if(sum(isnan(xDes(i,:)))==0 )
       gdind = [gdind;i]; 
    end
end
U       = U(gdind,:);
xDes    = xDes(gdind,:);

%% Declare Matrices which will store the final outputs of the simulation
Hmat  = zeros(nmc,Nstim+1);       % Entropies
Emat  = zeros(nmc,Nstim+1);       % EMSE
Tmat  = zeros(nmc,Niter);         % Time

X1mat = zeros(nmc,Nstim);         % Stimuli - x1
X2mat = zeros(nmc,Nstim);         % Stimuli - x2
EVmat = zeros(nmc,4);             % Vector of final parameter estimates - Expectation Value
MLmat = zeros(nmc,4);             % Vector of final parameter estimates - Max Likelihood

%% Monte Carlo loop
for j=1:nmc
    %% Initialize Experiment
    
    % Initialize uniform prior over theta
    [NParticles,NThetaDim]  = size(thetaMat);
    PTheta                  = (1/NParticles)*ones(NParticles,1);
    nDes                    = length(xDes); 

    % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
    % (i.e. given the prior)
    thetaHatCur             = ((thetaMat')*PTheta)';
    HCur                    = CalcEnt(PTheta'); 
    ECur                    = norm(thetaHatCur-ThetaTrue);

    % Initialize plot variables
    HPlot                   = NaN*ones(Nstim + 1,1); 
    EPlot                   = NaN*ones(Nstim + 1,1);
    HPlot(1)                = HCur;
    EPlot(1)                = ECur;

    if(showPlots) 
        close all; 
        figure(1); 
        subplot(2,2,1); plot(0:Nstim,HPlot); axis square; ylabel('Entropy'); xlabel('Iteration');
        subplot(2,2,2); plot(0:Nstim,EPlot); axis square; ylabel('Error'); xlabel('Iteration'); 
    end

    tElapsedVec = zeros(Niter,1);
    Xall        = zeros(Nstim,4); % on each iteration you present 4 stimuli
    Rall        = zeros(Nstim,1); 

    %% Experiment loop
    for i=1:Niter
        % Find and present next design
        if(strcmp(sampMethod,'OED'))
            tic,
            % find best design
            uExp        = U*PTheta;
            [dum,mxind] = max(uExp);
            
            tElapsedVec(i)    = toc;
            
            xAug1       = [1,xDes(mxind,1:2),xDes(mxind,1)*xDes(mxind,2)];
            xAug2       = [1,xDes(mxind,3:4),xDes(mxind,3)*xDes(mxind,4)];
            xAug3       = [1,xDes(mxind,5:6),xDes(mxind,5)*xDes(mxind,6)];
            xAug4       = [1,xDes(mxind,7:8),xDes(mxind,7)*xDes(mxind,8)];
            
        else
            r1          = 5*rand;
            r2          = 5*rand;
            xAug1       = [1,r1,r2,r1*r2];
            
            r1          = 5*rand;
            r2          = 5*rand;
            xAug2       = [1,r1,r2,r1*r2];
            
            r1          = 5*rand;
            r2          = 5*rand;
            xAug3       = [1,r1,r2,r1*r2];
            
            r1          = 5*rand;
            r2          = 5*rand;
            xAug4       = [1,r1,r2,r1*r2];
        end

        % Stimulus 1
        
        % present stimuli + record responses
        if( rand <= gMat(xAug1,ThetaTrue))
            R1 = 1;   
        else
            R1 = 0;   
        end
        Rall(4*(i-1)+1,:) = R1;
        Xall(4*(i-1)+1,:) = xAug1;

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

        HPlot(4*(i-1) + 2)          = HCur;
        EPlot(4*(i-1) + 2)          = ECur;


        % Stimulus 2
        
        % present stimuli + record responses
        if( rand <= gMat(xAug2,ThetaTrue))
            R2 = 1;   
        else
            R2 = 0;   
        end
        Rall(4*(i-1)+2,:) = R2;
        Xall(4*(i-1)+2,:) = xAug2;

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

        HPlot(4*(i-1) + 3)          = HCur;
        EPlot(4*(i-1) + 3)          = ECur;


        % Stimulus 3
        
        % present stimuli + record responses
        if( rand <= gMat(xAug3,ThetaTrue))
            R3 = 1;   
        else
            R3 = 0;   
        end
        Rall(4*(i-1)+3,:) = R3;
        Xall(4*(i-1)+3,:) = xAug3;

        if(R3==1)
            PTheta = (gMat(xAug3,thetaMat)').*PTheta;
        else
            PTheta = ((1-gMat(xAug3,thetaMat))').*PTheta;
        end
        PTheta = (1/sum(PTheta))*PTheta;

        % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
        thetaHatCur = ((thetaMat')*PTheta)';
        HCur        = CalcEnt(PTheta'); 
        ECur        = norm(thetaHatCur-ThetaTrue);

        HPlot(4*(i-1) + 4)          = HCur;
        EPlot(4*(i-1) + 4)          = ECur;

        % Stimulus 4
        
        % present stimuli + record responses
        if( rand <= gMat(xAug4,ThetaTrue))
            R4 = 1;   
        else
            R4 = 0;   
        end
        Rall(4*(i-1)+4,:) = R4;
        Xall(4*(i-1)+4,:) = xAug4;
        
        if(R4==1)
            PTheta = (gMat(xAug4,thetaMat)').*PTheta;
        else
            PTheta = ((1-gMat(xAug4,thetaMat))').*PTheta;
        end
        PTheta = (1/sum(PTheta))*PTheta;

       % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
        thetaHatCur = ((thetaMat')*PTheta)';
        HCur        = CalcEnt(PTheta'); 
        ECur        = norm(thetaHatCur-ThetaTrue);

        HPlot(4*(i-1) + 5)          = HCur;
        EPlot(4*(i-1) + 5)          = ECur;

        % Plot
        if(showPlots)  
            figure(1); 
            subplot(2,2,1); plot(0:(4*Niter),HPlot); axis square;
            subplot(2,2,2); plot(0:(4*Niter),EPlot); axis square; ylabel('Error'); xlabel('Iteration'); 
            pause(0.1);
        end
               
 
    

        
    end     % end for i = 1:Niter

 
   
    
    
    %% Save experimental results
    Tmat(j,:)   = tElapsedVec'; 
    Hmat(j,:)    = HPlot';
    Emat(j,:)   = EPlot';
    X1mat(j,:)  = Xall(:,2)';
    X2mat(j,:)  = Xall(:,3)';
    
    % Optimize posterior to get final estimate
    options     = optimset('GradObj','on','Display','off');
    thetaHatML  = fminunc(@(theta) NegLogLikelihoodWithGrad(double(Xall),Rall,theta),double(thetaHatCur),options); 

    EVmat(j,:)  = thetaHatCur;
    MLmat(j,:)  = thetaHatML;

    disp(sprintf('Iteration: %d',j));



end % for j=1:nmc

%% Save outputs
savstr = strcat('PSI_Lookup_2D_Out_',sampMethod,'_',sysID,'_',num2str(NParticles),'.mat'); 
save(savstr,'Hmat','Emat','Tmat','X1mat','X2mat','EVmat','MLmat','nmc','Niter','sampMethod','sysID','ThetaTrue'); 

clear all; 

