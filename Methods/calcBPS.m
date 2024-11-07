% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function  [mSamplesGPS,vResultsGPS,cStatesGPS,cSamplesTMCMC,cResultsTMCMC] = calcBPS (objective, vLB, vUB, nNDABC, nTDABC)
    % This function calculates the Bayesian pattern search (BPS) method
    % First step is based on code from the EngiO framework of Berger et al.
    % (2021), doi: 10.1016/j.advengsoft.2020.102959 with a new gridstrucure
    % Second step is based on code of Wolniak et al (2023), doi: 10.2139/ssrn.4648368
    
    % INPUT: 
    % - objective (function handle)
    %   objective function
    % - vLB and vUB (vector)
    %   Lower and upper bounds of design variables
    % - nNDABC (Number) 
    %   Muliplication factor of hall of fame, see N_{DABC} in paper
    % - nTDABC (Number)
    %   Number of tracked globally best samples, see T_{DABC} in paper
    %
    % Output:
    % - mSamplesGPS (matrix)
    %   Total samples of first step (GPS)
    % - vResultsGPS (vector)
    %   Objective function results of first step (GPS)
    % - cStatesGPS   (cell)
    %   Iteration details of GPS
    % - cSamplesTMCMC (cell)
    %   Total samples of second step (TMCMC)
    % - cResultsTMCMC (cell)
    %   Objective function results of first step (GPS)
    

    fprintf('Starting DABC \nStarting GPS \n');    
    % Start first step of DABC (GPS)
    nResolution=2^1023; 
        
    % Initialize sample and result history
    mSamplesGPS = [];
    vResultsGPS = [];
    % Initialize all states
    nDims = numel(vLB);

    states = struct('fBest', Inf, ...                  % Current best objective value
        'miSampleCache', [], ...           % Cache for sample coordinates
        'vfObjectiveValues', []);          % Cache for objective values

        % Calculate pattern
        vBases = eye(nDims);
        vLastPointDir = ones(nDims,1);

        fDistance = norm(vBases(:,1) - vBases(:,2));
        lambdafun = @(x) (abs( norm(vLastPointDir * x - vBases(:,1)) - fDistance ) );
        x = fminsearch(lambdafun, 0);

        vRawPattern = vBases;
        vRawPattern(:,nDims+1) = vLastPointDir * x;
        mfPattern = vRawPattern - mean(vRawPattern,2);

        states.mfTransform = mfPattern(1:nDims, 1:nDims);

        % Initialize resolution per dimension
        % Use given resolution for real dimensions
        states.viResolution = ones(1,nDims) * nResolution;

        % Initialize step sizes to half the search space
        states.viStep = floor((ones(1,nDims+1) * nResolution+1) / 2);
    
    % Allocate local variables
    numIters = 0;
    cStatesGPS ={};
    termFlag = 0;
    
    while ~termFlag
        % Iteration counter
        numIters = numIters+1;
        % Initialize and use all previous locations initially
            viBases = [zeros(1, numel(vLB)); states.miSampleCache];
            states.vfSortedMins = zeros(nTDABC,1);
            
            % If we have enough samples already, only use best ones
            if numel(states.miSampleCache) > 0
                if numel(states.miSampleCache(:, 1)) > nTDABC
                    % Sort by minima
                    [vfSortedMins, viSortedIndices] = sort(states.vfObjectiveValues);
                    viBases = states.miSampleCache(viSortedIndices(1:nTDABC),:);
                    states.vfSortedMins = vfSortedMins(1:nTDABC);
                    vfCoordinates = viBases ./ states.viResolution;
                    states.vfBases = (states.mfTransform * vfCoordinates')' .* (vUB-vLB) + (vUB+vLB)*.5;
                end
            end

        miSamples = []; % Stores sampling positions

            for iBase = 1:size(viBases,1)
                viBase = viBases(iBase, :);
                % Take care of initialization by adding the sample itself
                % This will be filtered out when already in cache
                miSamples = [miSamples ; viBase];
                
                % Move in axis directions
                for iDim = 1 : numel(vLB)
                    viSample = viBase;
                    viSample(iDim) = viSample(iDim) + states.viStep(iDim);
                    miSamples = [miSamples ; viSample];
                end
                
                % Move in last direction
                viSample = viBase - states.viStep(numel(vLB) + 1);
                miSamples = [miSamples ; viSample];
            end

        % Deduplicate samples using cache
        if numel(states.miSampleCache)
            miDedupSamples = setdiff(miSamples, states.miSampleCache, 'rows');
        else
            miDedupSamples = miSamples;
        end

        % Filter samples and output
        samples = [];
        for viSampleT = miDedupSamples'
            viSample = viSampleT';

            % Transform samples to floating point
            vfCoordinates = viSample ./ (states.viResolution);
            vfPosition = (states.mfTransform * vfCoordinates')' .* (vUB-vLB) + (vUB+vLB)*.5;

            if any(vfPosition > vUB) || any (vfPosition < vLB)
                continue
            end

            % We did not sample here before, put this sample to the cache
            states.miSampleCache = [states.miSampleCache; viSample];
            % Transform samples to floating point
            samples = [samples; vfPosition];
        end

        mSamplesGPS = [mSamplesGPS; samples]; 

        % Evaluate objectives
        % Initialize
        objectiveValues = zeros(size(samples, 1),1);

        parfor j = 1:size(samples, 1)
             objectiveValues(j,:) = -objective(samples(j,:), j); 
        end

        vResultsGPS = [vResultsGPS; objectiveValues]; 
        
        % Process results
        if termFlag == 0
            % Append results to objective value cache
            states.vfObjectiveValues = [states.vfObjectiveValues; objectiveValues];

            % Get minimum objective value
            fMin = min(objectiveValues);

            % Check if we have a new minimum
            if fMin < states.fBest
                states.fBest = fMin;
                terminate = false;
            else
                % Find out which dimension has the largest step
                [iMaxStep, iMaxDim] = max(states.viStep);
                terminate = false;
                % Reduce largest step by factor 2
                states.viStep(iMaxDim) = floor(states.viStep(iMaxDim)/ 2);
            end
            
            % Stopping cirterion inspired by the TMCMC
            if numIters>1 && ~sum(cStatesGPS{numIters-1}.vfSortedMins)==0
                ratio=std(exp(-states.vfSortedMins))/mean(exp(-states.vfSortedMins));
                if ratio < 1 && ratio > 0
                    termFlag = 3;
                end
                fprintf('\t \t Iteration:\t%d\t\t Samples:\t%d\t\t ratio:\t%d\t\t\n', numIters,size(mSamplesGPS,1),ratio );
            else
                fprintf('\t \t Iteration:\t%d\t\t Samples:\t%d\t\t ratio:\t--\t\t\n', numIters,size(mSamplesGPS,1));
            end

            if terminate
                termFlag = 1;
            end
        end

        % Store states to temporary struct
        cStatesGPS{numIters}= states;

        % Store sampling history
        cStatesGPS{numIters}.samples = samples;
        cStatesGPS{numIters}.objectiveValues = objectiveValues;
    end

    % Start second step od DABC (TMCMC)
    fprintf('Starting TMCMC \n');  
    afSamples_0=[];
    fResults_0=[];
    for i=1:nNDABC
        afSamples_0=[afSamples_0;cStatesGPS{end-1}.vfBases];
        fResults_0=[fResults_0;cStatesGPS{end-1}.vfSortedMins];
    end
    nSamplesTMCMC=size(fResults_0,1);
    nDim = numel(viBase);
    iIteration = 1;
    pj_TMCMC(iIteration) = 0;

    cSamplesTMCMC{iIteration} = afSamples_0;
    cResultsTMCMC{iIteration} = -fResults_0;
   

    while pj_TMCMC(iIteration) < 1
        % Evaluation stepsize
        wj = @(e) exp(abs(e) * cResultsTMCMC{iIteration});
        fmin = @(e) std(wj(e)) - 1 * mean(wj(e)) + realmin;
        e = abs(fzero(fmin, 0));
        pji_TMCMC = min(1, e + pj_TMCMC(iIteration));
  
        fprintf('\t \t Iteration:\t%d\t\t pji:\t%f\t Samples:\t%d\n', iIteration, pji_TMCMC,size(mSamplesGPS,1)+iIteration*nSamplesTMCMC);

        mu = zeros(1, nDim);
        a = (pji_TMCMC - pj_TMCMC(iIteration)) * (cResultsTMCMC{iIteration});
        wji = exp(a);
        wj_norm = wji./sum(wji);
        for iSample = 1:nSamplesTMCMC
            mu = mu + wj_norm(iSample) * cSamplesTMCMC{iIteration}(iSample,:);
        end

        % Calculate covariance matrix
        covGauss = zeros(nDim);
        for iSample = 1:nSamplesTMCMC
            tk_mu = cSamplesTMCMC{iIteration}(iSample,:) - mu;
            covGauss = covGauss + wj_norm(iSample)*(tk_mu'*tk_mu);
        end
        beta = 0.2; 
        covGauss = beta^2 * covGauss;

        % Sample generation
        samples_i_index = randsample(nSamplesTMCMC, nSamplesTMCMC, true, wj_norm);
        afSamples_i = cSamplesTMCMC{iIteration}(samples_i_index,:);
        fResults_i = cResultsTMCMC{iIteration}(samples_i_index);

        % MCMC with one step for new samples
        parfor iSample = 1:nSamplesTMCMC 
            thetaLead = afSamples_i(iSample,:);
            logLLead = fResults_i(iSample);

            % Candidate sample generation
            while true
                thetaCand = mvnrnd(thetaLead, covGauss);
                if all([(thetaCand > vLB),  (thetaCand < vUB)])
                    break;
                end
            end

            % Sample calculation
            cCurrentCand = -objective(thetaCand, iSample);
            % "Death Penalty"
            if isequal(cCurrentCand, Inf)
                cCurrentCand = -realmax; 
            end
            logLCand = -cCurrentCand;

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
        pj_TMCMC(iIteration) = pji_TMCMC;
        cSamplesTMCMC{iIteration} = theta_j1;
        cResultsTMCMC{iIteration} = logL_j1;
    end    
end