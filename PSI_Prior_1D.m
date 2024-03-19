
% PSI_Prior_1D.m
% by
% Christopher DiMattina, PhD 
% Florida Gulf Coast University
%
% Brief Description: This program implements the PSI method for a 
%                    1-D psychophysical experiment with a set of
%                    particles drawn from a prior   
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
NParticles          = 1000;     % Number of supports for evolving posterior
nmc                 = 1; 
Niter               = 50;
sampMethod          = 'OED';    % set to 'IID' for method of constant stimuli
sysID               = computer; if(strcmp(sysID,'PCWIN')) sysID = 'PCWIN32'; end
showPlots           = 1;        % set to 0 to supress graphical output

%% Establish parameters of true underlying observer
LogBetaTrue         = 0; 
BetaTrue            = 10^LogBetaTrue; 
LambdaTrue          = 0; 
ThetaTrue           = [ -LambdaTrue*BetaTrue , BetaTrue ];  % convert to generic parameterization

if(exist('PSI_Prior_1D_Lookup.mat'))
    load('PSI_Prior_1D_Lookup.mat');
else
    
    LogBetaVec      = 0.5*randn(NParticles,1); 
    BetaVec         = 10.^LogBetaVec;
    LambdaVec       = randn(NParticles,1); 
    
    theta0vec       = -(LambdaVec.*BetaVec);
    theta1vec       = BetaVec;
    
    thetaMat        = [theta0vec, theta1vec];
    
    xMin            = - 7;
    xMax            = 7; 
    
    % Grid density
    NXGrid          = 51;
    xVec            = linspace(xMin,xMax,NXGrid)';

     % Compute the matrix P(R = 1 | x , theta)
    %
    % rows will be the different value of x
    % columns will be the values of theta
    %
    
    xAug    = [ ones(size(xVec)), xVec ];
    U       = xAug*(thetaMat');
    PR1XTh  = 1./(1 + exp(-U)); 
    PR0XTh  = 1-PR1XTh; 
    
    save('PSI_Prior_1D_Lookup.mat','PR1XTh','PR0XTh','xVec','xAug','ThetaTrue','thetaMat','theta0vec','theta1vec','NXGrid','NParticles');

end


%% Declare Matrices which will store the final outputs of the simulation
Hmat  = zeros(nmc,Niter+1);         % Entropies
Emat  = zeros(nmc,Niter+1);         % EMSE
Tmat  = zeros(nmc,Niter);           % Time for stimulus generation
Xmat  = zeros(nmc,Niter);           % Stimuli presented

Neffmat = zeros(nmc,Niter+1);       % Number of effective particles

EVmat = zeros(nmc,2);               % Vector of final parameter estimates - Expectation Value
MLmat = zeros(nmc,2);               % Vector of final parameter estimates - Max Likelihood

%% Monte Carlo loop

for j=1:nmc
    %% Initialize experiment
    % Initialize uniform prior over theta
    [NParticles,NThetaDim]  = size(thetaMat);
    PTheta                  = (1/NParticles)*ones(NParticles,1);
    nx                      = length(xVec); 

    % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
    % (i.e. given the prior)
    thetaHatCur             = ((thetaMat')*PTheta)';
    HCur                    = CalcEnt(PTheta'); 
    ECur                    = norm(thetaHatCur-ThetaTrue);
    
    % Initialize plot variables
    HPlot                   = NaN*ones(Niter+1,1); 
    EPlot                   = NaN*ones(Niter+1,1);
    NeffPlot                = NaN*ones(Niter+1,1);
    
    HPlot(1)                = HCur;
    EPlot(1)                = ECur;
    NeffPlot(1)             = 1./sum(PTheta.^2);
    
    if(showPlots)
        close all;
        figure(1);
        subplot(2,2,1); plot(0:Niter,HPlot); axis square; ylabel('Entropy'); xlabel('Iteration');
        subplot(2,2,2); plot(0:Niter,EPlot); axis square; ylabel('Error'); xlabel('Iteration');
        subplot(2,2,3); plot(ThetaTrue(1),ThetaTrue(2),'bo','MarkerSize',10,'LineWidth',3); hold on; axis square;
        plot(thetaHatCur(1),thetaHatCur(2),'k.');
        subplot(2,2,4); plot(0:Niter,NeffPlot); axis square; ylabel('N_e_f_f'); xlabel('Iteration');
    end

    
    tElapsedVec = zeros(Niter,1);
    Xall        = zeros(Niter,2);
    Rall        = zeros(Niter,1); 

    %% Experiment loop
    for i = 1:Niter

        % For each candidate stimulus, compute the expected value of the posterior entropy
        % P_n+1(Theta) if we present stimulus x using our current posterior P_n(Theta)
        
        if(strcmp(sampMethod,'OED'))
            tic,
            
            PThRep      = repmat(PTheta',nx,1);
            
            A1          = PR1XTh.*PThRep; 
            
            D1          = repmat(1./sum(A1,2),1,NParticles);
            Pn1Th_XR1   = D1.*A1;
            
            %D1 = diag(1./sum(A1,2));
            %Pn1Th_XR1   = D1*A1;
            
            A0          = PR0XTh.*PThRep; 
            
            D0          = repmat(1./sum(A0,2),1,NParticles);
            Pn1Th_XR0   = D0.*A0;
            
            %D0 = diag(1./sum(A0,2));
            %Pn1Th_XR0   = D0*A0;
            
            PR1X        = PR1XTh*PTheta;
            PR0X        = PR0XTh*PTheta;
            
            H1          = CalcEnt(Pn1Th_XR1);
            H0          = CalcEnt(Pn1Th_XR0);
            
            Hexp        = PR1X.*H1 + PR0X.*H0;
            
            % Present chosen stimulus + record response
            [dum,ind]   = min(Hexp);
            
            xNext       = xAug(ind,:);
            
            tElapsedVec(i)    = toc;
        elseif(strcmp(sampMethod,'IID') )
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
        Rall(i) = R;
        Xall(i,:) = xNext;
        
        % Update posterior density
        if(R==1)
            PTheta = (PR1XTh(ind,:)').*PTheta;
        else
            PTheta = (PR0XTh(ind,:)').*PTheta;
        end
        PTheta = (1/sum(PTheta))*PTheta;
        
        
        % Estimate Theta + Entropy (thetaHatCur, HCur) given current posterior
        thetaHatCur = ((thetaMat')*PTheta)';
        HCur        = CalcEnt(PTheta');
        ECur        = norm(thetaHatCur-ThetaTrue);
        
        % Plot
        HPlot(i+1)          = HCur;
        EPlot(i+1)          = ECur;
        NeffPlot(i+1)       = 1./sum(PTheta.^2);
        
        if(showPlots)
            figure(1);
            subplot(2,2,1); plot(0:Niter,HPlot); axis square;
            subplot(2,2,2); plot(0:Niter,EPlot); axis square; ylabel('Error'); xlabel('Iteration');
            subplot(2,2,3);  plot(ThetaTrue(1),ThetaTrue(2),'bo','MarkerSize',10,'LineWidth',3);  hold on; axis square;
            plot(thetaHatCur(1),thetaHatCur(2),'k.');
            subplot(2,2,4); plot(0:Niter,NeffPlot); axis square; ylabel('N_e_f_f'); xlabel('Iteration');
            pause(0.1);
        end
        
    end % for i = 1:Niter
    
    %% Save experimental results
    Tmat(j,:)   = tElapsedVec';
    Hmat(j,:)   = HPlot';
    Emat(j,:)   = EPlot';
    Xmat(j,:)   = Xall(:,2)';
    
    Neffmat(j,:)= NeffPlot'; 
    
    options     = optimset('GradObj','on','Display','off');
    thetaHatML  = fminunc(@(theta) NegLogLikelihoodWithGrad(Xall,Rall,theta),thetaHatCur,options); 

    EVmat(j,:)  = thetaHatCur;
    MLmat(j,:)  = thetaHatML;

if(showPlots)
    figure(1); subplot(2,2,3);  plot(thetaHatML(1),thetaHatML(2),'rd','MarkerSize',10,'LineWidth',3);
end
disp(sprintf('Iteration: %d',j));

end % for j=1:nmc

%% Save outputs
savstr = strcat('PSI_Prior_1D_Out_',sampMethod,'_',sysID,'.mat'); 
save(savstr,'Hmat','Emat','Tmat','Xmat','Neffmat','EVmat','MLmat','nmc','Niter','sampMethod','sysID','ThetaTrue'); 

clear all; 




