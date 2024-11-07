% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function [vFreq, vShape]= CalcFEBeam (vTheta,nResolution)
    % References: 
    % Wernitz et al. (2022) doi: 10.1016/j.ymssp.2022.108808
    %
    % CalcFEBeam outputs the natural frequencies and mode shapes of a
    % cantilever beam based on the damage parameters and the resolution.
    %
    % INPUT:
    % - vTheta (vector of doubles)
    %   scaling factors of the induvidual FEs
    % - nResolution (number)
    %   Number of FEs
    %
    % OUPUT:
    % - vFreq (vector of double)
    %   Natural frequencies
    % - vShape (matrix of double)
    %   Mass normalized eigenvectors
    
    % Calculate stiffness and mass distribution of cantilever beam (physical properties of the cantilever beam)
    [vRhoL, nBeamLength, mEI] = CantileverBeamMechanics (nResolution);
    
    % Calculate mass matrix
    M = CantileverBeamMass(vRhoL, nBeamLength);
    
    % Calculate stiffness matrix
    mEI_update = mEI .* vTheta';
    K = CantileverBeamFEM(mEI_update, nBeamLength);
     
    % Perform modal analysis
    vModeIndices = 1:2:size(M,1);
    [vFreq, vShape] = ModalAnalysis (M,K,vModeIndices);
    
end

