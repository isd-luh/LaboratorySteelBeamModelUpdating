% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function vTheta=CalcThetaFunc(vDamageParameter,nResolution,nBeamLength)
    % CalcThetaFunc complies the damage distribution function and outputs
    %  the scaling factors of the induvidual FEs
    %
    % INPUT:
    % - vDamageParameter (vector of doubles)
    %   Parameters for damage distribution function
    % - nResolution (number)
    %   Number of FEs
    % - nBeamLength (number)
    %   Total beam length
    %
    % Output:
    % - vTheta (vector of doubles)
    %   scaling factors of the induvidual FEs
    
    vTheta=ones((nResolution),1);

    % In case the damage intensity and/or the damage width is equal to
    % zero
    if vDamageParameter(2) == 0 || vDamageParameter(3) == 0
        return
    end
    % Calculate first element that is scaled
    nNodeStart=(vDamageParameter(1)-vDamageParameter(2)/2)*((nResolution)/nBeamLength);
    if nNodeStart<0
        nNodeStart=0.5;
        vTheta(ceil(nNodeStart))=1-vDamageParameter(3);
    elseif nNodeStart==floor(nNodeStart)
        vTheta(nNodeStart+1)=1-vDamageParameter(3); 
    else
        vTheta(ceil(nNodeStart))=1-(vDamageParameter(3)*(ceil(nNodeStart)-nNodeStart));
    end

    % Calculate last element that is scaled
    nNodeEnd=(vDamageParameter(1)+vDamageParameter(2)/2)*((nResolution)/nBeamLength);
    if nNodeEnd>nResolution
        nNodeEnd=nResolution;
        vTheta(nResolution)=1-vDamageParameter(3);
    elseif nNodeEnd==floor(nNodeEnd)
        vTheta(nNodeEnd)=1-vDamageParameter(3); 
    else
        vTheta(ceil(nNodeEnd))=1-(vDamageParameter(3)*(nNodeEnd-floor(nNodeEnd)));
    end
    
    % If only one element is scaled
    if floor(nNodeStart) == floor(nNodeEnd)
        vTheta(ceil(nNodeStart))=1-(vDamageParameter(3)*(nNodeEnd-nNodeStart));
        return
    end
    % If multiple elements are scaled
    for n= ceil(nNodeStart)+1:floor(nNodeEnd)
        vTheta(n)=1-vDamageParameter(3);
    end
end


