% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function [cSamples, cResults] = calcTMCMC(LogL, afLB, afUB, nSamplesTMCMC)
% Implementation based on Wolniak et al (2023), doi: 10.2139/ssrn.4648368 with minor changes
%
% This function calculates the transitional Markov chain Monte Carlo
% (TMCMC) method 
% 
% References: 
% Ching and Chen (2007), doi: 10.1061/(ASCE)0733-9399(2007)133:7(816)
% Lye et al. (2021), doi: 10.1016/j.ymssp.2021.107760

% INPUT: 
% - Problem function handle: LogL
% - Lower and upper bounds: afLB, afUB
% - Number of samples utilized: nSamplesTMCMC
% - Initial tempering parameter: pj_0

% OUTPUT: 
% - Generated samples: cSamples
% - Calculated objective function values: cResults

%%
% Generate Monte Carlo samples
nDim = size(afLB,2);
iIteration = 1;
pj(iIteration) = 0;
afSamples_0 = zeros(nSamplesTMCMC, nDim);
for iDim = 1:nDim
    afSamples_0(:,iDim) = unifrnd(afLB(iDim), afUB(iDim), nSamplesTMCMC, 1);
end

% Calculate first distribution in whole parameter space
parfor iSample = 1:nSamplesTMCMC 
    cCurrentResultsLogL = LogL(afSamples_0(iSample,:), iSample);
    fResults_0(iSample) = cCurrentResultsLogL;
end
cSamples{iIteration} = afSamples_0;
cResults{iIteration} = fResults_0;
fprintf('\t \t Iteration:\t%d\t\t pji:\t%f\n', 0, 0);

iCount = 0;
while pj(iIteration) < 1
    iCount = iCount + 1;

    % Evaluation stepsize
    wj = @(e) exp(abs(e) * cResults{iIteration});
    fmin = @(e) std(wj(e)) - 1 * mean(wj(e)) + realmin;
    e = abs(fzero(fmin, 0));
    pji = min(1, e + pj(iIteration));
    fprintf('\t \t Iteration:\t%d\t\t pji:\t%f\n', iIteration, pji);

    mu = zeros(1, nDim);
    a = (pji - pj(iIteration)) * (cResults{iIteration});
    wji = exp(a);
    wj_norm = wji./sum(wji);
    for iSample = 1:nSamplesTMCMC
        mu = mu + wj_norm(iSample) * cSamples{iIteration}(iSample,:);
    end

    % Calculate covariance matrix
    covGauss = zeros(nDim);
    for iSample = 1:nSamplesTMCMC
        tk_mu = cSamples{iIteration}(iSample,:) - mu;
        covGauss = covGauss + wj_norm(iSample)*(tk_mu'*tk_mu);
    end
    beta = 0.2; 
    covGauss = beta^2 * covGauss;

    % Sample generation
    samples_i_index = randsample(nSamplesTMCMC, nSamplesTMCMC, true, wj_norm);
    afSamples_i = cSamples{iIteration}(samples_i_index,:);
    fResults_i = cResults{iIteration}(samples_i_index);

    % MCMC with one step for new samples
    parfor iSample = 1:nSamplesTMCMC 
        thetaLead = afSamples_i(iSample,:);
        logLLead = fResults_i(iSample);

        % Candidate sample generation
        while true
            thetaCand = mvnrnd(thetaLead, covGauss);
            if all([(thetaCand > afLB),  (thetaCand < afUB)])
                break;
            end
        end

        % Sample calculation
        logLCand = LogL(thetaCand, iSample);

        % Acceptance / rejection step
        alpha = exp((logLCand - logLLead));
        if rand <= min(1, alpha)
            theta_j1(iSample, :) = thetaCand;
            logL_j1(iSample) = logLCand;
            thetaLead = thetaCand;
            logLLead = logLCand;
        else
            theta_j1(iSample, :) = thetaLead;
            logL_j1(iSample) = logLLead;
        end
    end
    iIteration = iIteration + 1;
    pj(iIteration) = pji;
    cSamples{iIteration} = theta_j1;
    cResults{iIteration} = logL_j1;
end

end