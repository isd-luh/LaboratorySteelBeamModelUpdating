% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function M = CantileverBeamMass(mRhoL, L)
    % CantileverBeamMass outputs the mass matrix of a
    % cantilever beam based on its length and specific mass.
    %
    % INPUT:
    % - mRhoL (double)
    %   specific mass
    % - L (double)
    %   total length of the beam
    %
    % OUPUT:
    % - M (matrix of doubles)
    %   Mass matrix

    % calculate Mass matrix based on element properties

    % get number and length of elements
    nEl = numel(mRhoL);
    lenEl = L/nEl;

    % calculate beam element Mass
    Mel = (lenEl/420)*...                                                       
        [156 22*lenEl 54 -13*lenEl;...                                           
        22*lenEl 4*lenEl^2 13*lenEl -3*lenEl^2;...
        54 13*lenEl 156 -22*lenEl;...
        -13*lenEl -3*lenEl^2 -22*lenEl 4*lenEl^2];

    % assemble global Mass matrix
    M = zeros((nEl+1)*2);
    for iElem = 1:nEl
        M((iElem-1)*2+1:(iElem+1)*2,(iElem-1)*2+1:(iElem+1)*2) = ...
            M((iElem-1)*2+1:(iElem+1)*2,(iElem-1)*2+1:(iElem+1)*2) + Mel*mRhoL(iElem);        
    end

    % implement boundary conditions
    idxDiscard = [length(M)-1,length(M)];
    M(idxDiscard,:) = [];
    M(:,idxDiscard) = [];
end