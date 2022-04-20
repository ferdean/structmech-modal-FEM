function [visc_flag_nN,visc_flag_nE,...
    visc_lay_nN,visc_lay_nE] = materialdef(mesh,contour,hidemesh)
% Material definition of multi-material problems in 2D plates of material_1
% with an external layer of material_2 
% =========================================================================
% Created by:   Ferran de Andrés(3.2021) 
% =========================================================================
% INPUT: 
%   nodes      = Point matrix. The first and second rows contain x- and 
%                y-coordinates of the points in the mesh.
%   TOP        = Topology matrix
%              
%   model      = Struct containing the original p, e and t matrices of the
%                initial mesh. IMPORTANTE ACTUALIZAR. 
%   contour    = 4x1 array containing the limits of the internal body with
%                material_1.
%   hidemesh   = Flag for showing or not the mesh in the plot
%
% OUTPUT:
% visc_flag_nN = A flag is raised in the nodes in which material_2 is
%                applied. 1 x nN vector
% visc_flag_nE = A flag is raised in  the elements in which material_2 is
%                applied. 1 x nE vector
% visc_lay_nN  = Vector containing the nodes in which material_2 is
%                applied.
% visc_lay_nE  = Vector containing the elements in which material_2 is
%                applied.
% =========================================================================

%%% Input check (3 inputs)
narginchk(3,3);

%%% Program initialisation
visc_flag_nN    = zeros(1,size(mesh.nodes,1));
visc_flag_nE    = zeros(1,size(mesh.top  ,2));
visc_lay_nN     = [];
visc_lay_nE     = [];

%%% Nodes scrutineering
for ii = 1:size(mesh.nodes,1)
    if mesh.nodes(ii,1) < contour(1) || mesh.nodes(ii,1) > contour(2) || ...
       mesh.nodes(ii,2) > contour(3) || round(mesh.nodes(ii,2)*1e5)/1e5 < contour(4)
   
        visc_flag_nN(ii)   = 1;
        visc_lay_nN        = [visc_lay_nN ii];
        
    end
end

%%% Element scrutineering
for jj = 1:size(visc_lay_nN,2)
    for kk = 1:size(mesh.top,2)
        if mesh.top(1,kk) == visc_lay_nN(jj) || mesh.top(2,kk) == visc_lay_nN(jj) ||...
           mesh.top(3,kk) == visc_lay_nN(jj)
       
            visc_lay_nE = [visc_lay_nE kk];
            visc_flag_nE(kk)   = 1;
            
        end        
    end
end

%%% Solution postprocess
pdeplot(mesh.nodes',mesh.e,mesh.top,'xydata',visc_flag_nE,...
    'ColorMap','parula','XYStyle','flat');

t = mesh.top;
p = mesh.nodes';

hold on
if ~(hidemesh)
    for i = 1:size(t,2)
        plot([p(1,t(1,i)),p(1,t(2,i)),p(1,t(3,i)),p(1,t(1,i))],...
             [p(2,t(1,i)),p(2,t(2,i)),p(2,t(3,i)),p(2,t(1,i))],...
             'color',[0.5 0.5 0.5])
    end
end

set(gca,'fontname','times')
colorbar off
set(gcf,'units','normalized','outerposition',[0 0 1 1])
daspect([1 1 1]);
xlabel('Width [m]')
ylabel('Height [m]')
title({'Materials definition',...
    '\fontsize{12}\color{gray}Yellow: viscoelastic layer. Blue: steel core'},'fontsize',14)
end

