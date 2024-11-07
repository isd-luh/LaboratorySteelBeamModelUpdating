% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

close all 
clear all
clc

% runBeamProblem
% Model updating implementation for a laboratory steel cantilever beam with reversible damage mechanism by Wolniak et al. (doi: 10.1007/s13349-023-00701-9) 
% Including implementtion of: 
%   - TMCMC (transitional Markov chain Monte Carlo)
%   - BPS (Bayesian pattern search)
% The code is implemented using the Matlab programming syntax (MATLAB R2024a). 
%
% IMPORTANT
% To run the code, BayOMA identification results need to be
% downloaded from Wolniak et al. (2024),  https://doi.org/10.25835/r8pevw8m
%
% References:
% Dierksen et al. (2024), doi: *will be added after publication*
% Wolniak et al. (2024), doi: 10.25835/r8pevw8m
% Wolniak et al. (2023), doi: 10.1007/s13349-023-00701-9
% Wolniak et al (2023), preprint available at https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4991201
% Wolniak et al (2022), doi: 10.25835/123gy6gm
% Wernitz et al. (2022) doi: 10.1016/j.ymssp.2022.108808
% Hofmeister et al. (2019), doi: 10.1016/j.engstruct.2019.05.047
% Berger et al. (2021), doi: 10.1016/j.advengsoft.2020.102959
% Ching and Chen (2007), doi: 10.1061/(ASCE)0733-9399(2007)133:7(816)
% Lye et al. (2021), doi: 10.1016/j.ymssp.2021.107760

% Setting the individual methods 
nSettingTMCMC=2000;          % Number of samples per iteration, recommendation: 2000
vSettingBPS=[3,33];        % [Number of repetitions of the hall of fame ,Number of tracked globally best samples], recommendation: [3,33]

% Select damage scenarios
% Damage kinds: [Discrete, Gaussian distributed, uniformly distributed]
vDamageKind=1:3;
% Damage positions: [1:9], 1 close to fixed end, 9 close to the free end
vDamagePosition=1:9;

% No changes are required below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add problem folder to path
addpath('AnalyticalBeam');
% Add BayOMA Identification folder to path
addpath('AnalyticalBeam/BayOMAIdentification');
% Add methods folder to path
addpath('Methods');

rng(1)

% Running different methods for selected damage scenarios
for iDamageKind=vDamageKind       
    for iDamagePosition=vDamagePosition       
        problem=classBeam(iDamageKind,iDamagePosition);
        ObjFuncTMCMC = @(vTheta,iIndex) problem.ObjFuncTMCMC (vTheta,iIndex);
        % Run TMCMC
        if nSettingTMCMC > 0
            [cSamplesTMCMC{iDamageKind}{iDamagePosition},cObjValTMCMC{iDamageKind}{iDamagePosition}] = ...
                calcTMCMC(ObjFuncTMCMC, problem.vLB, problem.vUB, nSettingTMCMC);
        end
        % Run BPS
        if vSettingBPS(1) > 0
            [cSamplesBPS_GPS{iDamageKind}{iDamagePosition},~,~,cSamplesBPS_TMCMC{iDamageKind}{iDamagePosition},~] ...
                = calcBPS(ObjFuncTMCMC, problem.vLB, problem.vUB, vSettingBPS(1), vSettingBPS(2));
        end
        fprintf('\n\nDamage Kind %d, Damage Position %d\n\n',iDamageKind,iDamagePosition);  
    end
end

% Visualise results
vColors=parula(10);
cXLabel={'Damage position p in mm','Damage expansion e in mm','Damage intensity D'};
cYLabel={'F(p)','F(e)','F(D)'};
cTitle={'Gaussian distributed damage scenario','Uniformly distributed damage scenario','Discrete damage scenario'};
% Create three figures with three subplots
for iDamageKind=vDamageKind   
    figure()
    t=tiledlayout(2,2);
    for iPar=1:problem.nDim 
        if iPar==1
            ax1=nexttile(1,[1,2]);
        elseif iPar==2
            ax1=nexttile(3);
        elseif iPar==3
            ax1=nexttile(4);
            % Plot empty lines for legend
            iLegend=0;
            plot(ax1,[NaN NaN], [NaN NaN], 'Color', 'black','LineStyle', '--', 'LineWidth', 2);
            cLegend{iLegend+1}='TMCMC';
            hold on
            plot(ax1,[NaN NaN], [NaN NaN], 'Color', 'black','LineStyle', '-', 'LineWidth', 2);
            cLegend{iLegend+2}='BPS';
            hold on
            plot(ax1,[NaN NaN], [NaN NaN], 'Color', 'black','LineStyle', ':', 'LineWidth', 2);
            cLegend{iLegend+3}='Correct damage position';
            hold on
            for iPos=vDamagePosition
                plot(ax1,[NaN NaN], [NaN NaN], 'Color', vColors(iPos,:),'LineStyle', '-', 'LineWidth', 2);
                cLegend{iLegend+3+1}=['Damaged fishplate F',num2str(iPos)];
                iLegend=iLegend+1;
                hold on
            end
        end
        
        for iDamagePosition=vDamagePosition  
            % Plot TMCMC results
            if nSettingTMCMC > 0
                h2=cdfplot(cSamplesTMCMC{iDamageKind}{iDamagePosition}{end}(:,iPar));
                set(h2,'color',vColors(iDamagePosition,:),'LineStyle', '--');
                hold on
            end
            % Plot BPS results
            if vSettingBPS(1) > 0
                h3=cdfplot(cSamplesBPS_TMCMC{iDamageKind}{iDamagePosition}{end}(:,iPar));
                set(h3,'color',vColors(iDamagePosition,:),'LineStyle', '-');
                hold on
            end
            % Plot correct solutions
            if iPar==1
                xl=xline(problem.vTheaCor(iDamagePosition),'color',vColors(iDamagePosition,:),'LineStyle', ':');
                xl.LineWidth = 2;
                hold on
            end   
        end
        xlabel(cXLabel{iPar})
        ylabel(cYLabel{iPar})
        xlim([problem.vLB(iPar) problem.vUB(iPar)])
        title('Random title', 'Color','none')
    end
    title(t,cTitle{iDamageKind})
    lg  = legend(ax1,cLegend,'NumColumns',3,'Orientation','horizontal','FontSize',10); 
    lg.Layout.Tile = 'South';
end