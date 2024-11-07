% This file is part of LaboratorySteelBeamModelUpdating – A model updating implementation for damage localisation and quantification considering uncertainty from modal identification.
% Copyright (C) 2024 Leibniz Universität Hannover, Institut für Statik und Dynamik. Authors: Niklas Dierksen, Marlene Wolniak, Benedikt Hofmeister, Jasper Ragnitz, Clemens Hübler, Stefan Wernitz and Raimund Rolfes
% LaboratorySteelBeamModelUpdating is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% LaboratorySteelBeamModelUpdating is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Foobar. If not, see https://www.gnu.org/licenses/.

function [f0, Phi] = ModalAnalysis (M,K, indices)
    % ModelAnalysis outputs the natural frequencies and mode shapes of a cantilever beam based on its mass and stiffness matrix
    %
    % INPUT:
    % - M (matrix of doubles)
    %   Mass matrix
    % - K (matrix of doubles)
    %   stiffness matrix
    %
    % OUPUT:
    % - f0 (vector of double)
    %   Eigenfrequencies
    % - Phi (matrix of double)
    %   Mass normalized eigenvectors
    %

    if ~exist('indices')
        indices = [1:size(K,1)];
    end

    % extract eigenfrequencies
    [eigenVectors,eigenValues] = eig(K,M);
    f0 = sqrt(diag(eigenValues)) / (2*pi);
    
    [f0, ind] = sort(f0);
    eigenVectors = eigenVectors(:, ind);

    % extract modal matrix
    for iVector = 1:size(K, 1)
        % flip eigenvectors
        ev = eigenVectors(:,iVector);
        [~, maxIndex] = max(abs(ev ));
        ev = ev / eigenVectors(maxIndex,iVector);

        % normalize by mass
        Phi(:, iVector) = sqrt(1/(ev'*M*ev )) * ev;
    end

    %Reduce Mode Shape based on indices
    Phi = Phi(indices,:);
end


