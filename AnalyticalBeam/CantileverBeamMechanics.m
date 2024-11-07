% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function [vRhoL, nBeamLength, vEI] = CantileverBeamMechanics(nResolution)
    % CantileverBeamMechanics outputs the density distribution, beam length
    % and stiffness distribution based in the resolution.
    %
    % INPUT:
    % - nResolution (number)
    %   Number of FEs
  
    % OUPUT:
    % - vRhoL (vector of doubles)
    %   specific masses along the beam
    % - nBeamLength (number)
    %   Total length of the beam 
    % - vEI (vector of doubles)
    %   sectional stiffness per element

    % Geometrical properties (m)
    nBeamLength =       1.205; % (m)
    nBeamWidth  =       .06; % (m)
    nBeamHeight =       .00515; % (m)
    nFishplateWidth =   .02; % (m)
    nFishplateHeight =  .00485; % (m)
    nSensorWeight =     .005; % (kg) 
    % Physical properties
    nRho = 7800; % (kg/m^3)
    nRhoFastener = 0.3; % (kg/m)   
    E = 2e+11 * 0.6352; % (N/m^2) 0.6352 was chosen to account for the deviation of the mechanical properties of the steel from literature values.

    nElementLength = (nBeamLength / nResolution);
    
    % Calculate area (m^2)
    nBeamA = nBeamWidth*nBeamHeight;
    nFishplateA = nFishplateWidth*nFishplateHeight;
    
    % Calculate moment of inertia (m^4)
    nBeamI = nBeamWidth*(nBeamHeight^3) / 12;
    nFishplateI = nFishplateWidth*(nFishplateHeight^3) / 12;
    
    % Calculate moment of inertia for single fishplate
    nCl = (nFishplateA * (nBeamHeight + nFishplateHeight)/2) / (nBeamA + nFishplateA);
    nSingleI = (nBeamI + nBeamA*nCl^2 + nFishplateI + nFishplateA*((nBeamHeight + nFishplateHeight)/2 - nCl)^2);
    nDoubleI = nBeamI + 2 * (nFishplateI + nFishplateA*((nBeamHeight + nFishplateHeight)/2)^2 );
    
    % Set geometrical properties
    afOverlapPos = .12 * (1:9);
    fOverlapLength = .01;
    afPropertyPos = [0, reshape([afOverlapPos; fOverlapLength+afOverlapPos],1,2*numel(afOverlapPos)), nBeamLength];
    
    bProperty = repmat([0,1],1,10);
    bProperty = bProperty(1:(end-1));
    
    % Specific mass distribution
    nRhoLSingle = nRho * (nBeamA + nFishplateA) + nRhoFastener; % (kg/m)
    nRhoLDouble = nRho * (nBeamA + 2*nFishplateA) + nRhoFastener; % (kg/m)
    
    vPosMass = diff(afPropertyPos).*bProperty * nRhoLDouble + diff(afPropertyPos).*(1-bProperty) * nRhoLSingle;
    vMassAccum = [0,cumsum(vPosMass)];
    
    % Stiffness distribution
    vStiffPos = diff(afPropertyPos).*bProperty * nDoubleI + diff(afPropertyPos).*(1-bProperty) * nSingleI;
    vStiffAccum = [0,cumsum(vStiffPos)];
    
    afElementPos = (0:nResolution) / nResolution * nBeamLength;
    vRhoL = diff(interp1(afPropertyPos, vMassAccum, afElementPos)) / nElementLength;
    vEI  = E * diff(interp1(afPropertyPos, vStiffAccum, afElementPos)) / nElementLength;
    
    % Add sensor weight
    vSensorIndices = CantileverBeamSensorIndices(nResolution,nBeamLength);
    vRhoL(vSensorIndices) = vRhoL(vSensorIndices) + nSensorWeight / nElementLength;
end