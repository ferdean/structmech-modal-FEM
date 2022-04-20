function plotAxis3D(bbox)
%   Plots cartesian axis triad based on input bounding box
%
%   Copyright 2014-2019 The MathWorks, Inc.
% =========================================================================
% Modified by Ferran de Andrés (7.2021)
% =========================================================================

maxdim = max([abs(bbox(1,1)-bbox(1,2)), abs(bbox(2,1)-bbox(2,2)), abs(bbox(3,1)-bbox(3,2))]);
orgn = (ones(3,3).*(bbox(:,1) - 0.3*maxdim))';
vMag = .3*maxdim;
vec = vMag*eye(3);
quiver3(orgn(:,1),orgn(:,2),orgn(:,3),vec(:,1),vec(:,2),vec(:,3), 'Color', 'k','LineWidth',1.1);
labPts = (orgn + vec);
at = matlab.graphics.primitive.world.Text('VertexData', single(labPts)', 'String', {'x';'y';'z'}, 'HorizontalAlignment','center', 'Clipping','on','Interpreter','none', 'Parent', gca);
at.Font.Name = 'Helvetica';
at.Font.Size = 10;
end