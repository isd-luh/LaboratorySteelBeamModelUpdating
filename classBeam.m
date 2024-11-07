% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

classdef classBeam < handle
    % Class for performing Model updating implementation for a laboratory steel cantilever beam 
    %   with reversible damage mechanism
    % The implementation features a damage distribution function to reduce
    %   the dimension of the updating problem
    % References: 
    % Wolniak et al. (2023), doi: 10.1007/s13349-023-00701-9

    properties
        
        % Number of nodes in the FE model 
        nNodes=241; % nNodes=100 is recommended for high computational speed, nNodes=241 results in a discretisation of 5 mm and is recommended for more accurate results
        % Total length of the beam
        nBeamLength=1205; % (mm)
        % Select the modes taking into consideration (first four frequencies related to vertical bending modes)
        vModesMea=[1,2,3,4]; 
        vModesSim=[1,2,3,5];
        % Devide measurements in damaged and undamaged states
        vDamages=[1:2:18];
        vReferences=[2:2:18];
        
        % Lower und upper bound of the optimisation problem
        vLB=[];
        vUB=[];
        % Abbreviation M,S,D,R stand for measurement, simulation, damaged
        % and reference, respectively
        vMeanMD=[];
        vMeanMR=[];
        vStdMD=[];
        vStdMR=[];
        vSR=[];
        % Number of dimensions of the optimisation problem
        nDim=[];
        % Correct damage position
        vTheaCor=[];
    end
    
    methods
        % Constructor of class
        function this = classBeam(iDamageKind,iDamagePosition)
            % classBeam is the constructor of the class and selects the relevant OMA identification results and
            %   calculates the undamaged simulation state as a refernce
        
            % Load BayOMA identification results
            % Download data from https://doi.org/10.25835/r8pevw8m and save in
            % folder "BayOMAIdentification"

            strBayOMATitle{1}='GaussianDistributed';
            strBayOMATitle{2}='UniformlyDistributed';
            strBayOMATitle{3}='Discrete';
            for i=1:size(strBayOMATitle,2)
                tParquetFStd=parquetread([strBayOMATitle{i},'_','STD.parquet']);
                cFStdData{i}=cell2mat(table2array(tParquetFStd(:,2:end)));
            end
            for i=1:size(strBayOMATitle,2)
                tParquetF=parquetread([strBayOMATitle{i},'_','F.parquet']);
                cFData{i}=cell2mat(table2array(tParquetF(:,2:end)));
            end
            tParquetP=parquetread('DamagePositions.parquet');
            this.vTheaCor=table2array(tParquetP(:,2));
        
            % Set lower und upper bound of the optimisation problem (The second and third limits are based on experience and can theoretically be increased)
            this.vLB = [    0,                 0,                      0       ];
            this.vUB = [    this.nBeamLength,  0.5*this.nBeamLength,   0.5    ];
            % Set number of dimensions of the optimisation problem
            this.nDim=size(this.vLB,2);
            % Select corresponding measurement results
            this.vMeanMD=cFData{1, iDamageKind} (this.vDamages(iDamagePosition),this.vModesMea);
            this.vMeanMR=cFData{1, iDamageKind} (this.vReferences(iDamagePosition),this.vModesMea);
            this.vStdMD=cFStdData{1, iDamageKind} (this.vDamages(iDamagePosition),this.vModesMea);
            this.vStdMR=cFStdData{1, iDamageKind} (this.vReferences(iDamagePosition),this.vModesMea);
            % set stiffness variation for reference state
            vStiffVariationHealthy=ones(1,(this.nNodes)).';
            % calculate reference state
            [vFreqSR,~]= CalcFEBeam (vStiffVariationHealthy,this.nNodes);
            this.vSR=vFreqSR(this.vModesSim).';
        end
           % Objective function TMCMC
        function y=ObjFuncTMCMC (this,vTheta,iIndex)
            % ObjFuncTMCMC outputs the objective value for the TMCMC method depending on the
            % current sample vTheta
            %
            % INPUT:
            % - vTheta (vector of doubles)
            %   Parameters for damage distribution function
            % - iIndex (integer)
            %   Index for parallel computing
            %
            % Output:
            % - y (double)
            %   objective value

            nchains=size(vTheta,1);
            y=zeros(nchains,1);
            for n=1:nchains
                vFreqSD = this.CalcModel(vTheta(n,:));
                SD=vFreqSD(this.vModesSim);
                MD=this.vMeanMD.';
                MR=this.vMeanMR.';
                MDStd=this.vStdMD;
                MRStd=this.vStdMR;
                % Calculate the average normalised standard deviation
                nAllStd=mean(sqrt((MDStd./MD').^2+(MRStd./MR').^2));
                
                iObj=zeros(size(MR));
                % Calculate differece of simulation and observation of all
                % investigated modes
                for i=1:size(MR,1)
                    iObj(i)=-( (SD(i)-this.vSR(i))/this.vSR(i) - (MD(i)-MR(i))/MR(i) )^2/nAllStd^2;
                end
                y(n) = sum(iObj);
            end
            % Death penalty
            if ~isreal(vFreqSD) || y(n)==-Inf
                y=-realmax;
            end
        end
        
                % Objective function GPS
        function y=ObjFuncGPS (this,vTheta,iIndex)
            % ObjFuncGPS outputs the objective value for the GPS method depending on the
            % current sample vTheta
            %
            % INPUT:
            % - vTheta (vector of doubles)
            %   Parameters for damage distribution function
            % - iIndex (integer)
            %   Index for parallel computing
            %
            % Output:
            % - y (double)
            %   objective value

            if isempty(vTheta)
                y=1;
                return
            end
            % The GPS algorithm minimises the objective value (in contrast
            % to TMCMC), therefore the objective value of ObjFuncTMCMC
            % needs to be multiplied by (-1)
            y=ObjFuncTMCMC (this,vTheta,iIndex);safe
            y=y*(-1);
        end

             % calculate model
        function  vFreqSD = CalcModel(this,vDamageParameter)
            % CalcModel outputs the natural frequencies of the steel cantilever
            % beam based on the damage parameters for the damage
            % distribution function
            %
            % INPUT:
            % - vDamageParameter (vector of doubles)
            %   Parameters for damage distribution function
            %
            % Output:
            % - vFreqSD (vector of doubles)
            %   natural frequencies of the steel cantilever beam

            % Using damage distribution function
            vStiffVariation=CalcThetaFunc(vDamageParameter, this.nNodes,this.nBeamLength);
            % vStiffVariation2=CalcThetaFunc_auchalt([vDamageParameter(1)/this.nBeamLength*this.nNodes,vDamageParameter(2)/this.nBeamLength*this.nNodes,vDamageParameter(3)], this.nNodes,'Uniform');
            

            % Calculate FE model
            [vFreqSD,~] = CalcFEBeam (vStiffVariation,this.nNodes);
        end
    end
end