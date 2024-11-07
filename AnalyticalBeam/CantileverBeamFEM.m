% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function K = CantileverBeamFEM(vEI, L)
    % CantileverBeamFEM outputs the stiffness matrix of a
    % cantilever beam based on its length and sectional bending stiffness.
    %
    % INPUT:
    % - vEI (vector of doubles)
    %   sectional stiffness per element
    % - L (double)
    %   total length of the beam
    %
    % OUPUT:
    % - K (matrix of doubles)
    %   stiffness matrix

    % calculate stiffness matrix based on element properties

    % get number and length of elements
    nEl = numel(vEI);
    lenEl = L/nEl;

    % calculate beam element stiffness
    Kel = (1/lenEl^3)*...                                                       
        [12 6*lenEl -12 6*lenEl;...                                           
        6*lenEl 4*lenEl^2 -6*lenEl 2*lenEl^2;...
        -12 -6*lenEl 12 -6*lenEl;...
        6*lenEl 2*lenEl^2 -6*lenEl 4*lenEl^2];

    % assemble global stiffness matrix
    K = zeros((nEl+1)*2);
    for iElem = 1:nEl
        K((iElem-1)*2+1:(iElem+1)*2,(iElem-1)*2+1:(iElem+1)*2) = ...
            K((iElem-1)*2+1:(iElem+1)*2,(iElem-1)*2+1:(iElem+1)*2) + Kel*vEI(iElem);        
    end

    % implement boundary conditions
    idxDiscard = [length(K)-1,length(K)];
    K(idxDiscard,:) = [];
    K(:,idxDiscard) = [];

end