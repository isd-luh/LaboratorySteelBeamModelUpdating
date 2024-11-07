% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function vSensorIndices = CantileverBeamSensorIndices(nResolution,nBeamLength)
    % CantileverBeamSensorIndices outputs the indices of sensor positions
    % of cantilever beam experiment
    %
    % INPUT:
    % - nResolution (number)
    %   Number of FEs
    % - nBeamLength (number)
    %   Total length of the beam 
    %
    % OUPUT:
    % - vSensorIndices (matrix of doubles)
    %   Mass matrix
    
    vSensorStart=.01;
    vSensorDist=.075;
    nSensors=16;
    
    afSensorPos = vSensorStart + (0:nSensors-1) * vSensorDist;
    
    afElePos = (1:nResolution)/nResolution * nBeamLength;
    
    vSensorIndices = [];
    
    for iSensorPos=1:nSensors
        fSensor = afSensorPos(iSensorPos);
        [~,vSensorIndices(iSensorPos)] = min(abs(afElePos-fSensor));
    end
        
end
