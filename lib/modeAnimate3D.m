function []=modeAnimate3D(mode,mesh,DOF,DOFfree,wnHz,phi,scaler,n_loops,original)
% Representation of the vibration modes.
% =========================================================================
% Created by:      Ferran de Andrés (2.2021)
% =========================================================================
% Versions:        1.0: (Ferran) First implementation of the function. 
%                                Steady representation of the mode.
%                  2.0: (Ferran) Adaptation for 3D models (tet4 elements)
% =========================================================================
% % INPUT:
%   mode      = Number (ordinary) of the represented mode
%               Scalar value
%   mesh      = Struct with all the meshing information listed below: 
%      nodes  = Point matrix. The first and second rows contain x- and 
%               y-coordinates of the points in the mesh.
%      top    = Triangle topology matrix. The first three rows contain 
%               indices to the corner points, given in counter clockwise 
%               order, and the fourth row contains the subdomain number.
%   DOF       = Matrix of dof's considering constained nodes
%               Size nn x Nn, with Nn nº dof's per node
%   DOFfree   = Matrix of dof's considering no constained nodes
%               Size nn x Nn, with Nn nº dof's per node
%   wnHz      = Natural frequencies [Hz]
%               Size tr x 1
%               Each row corresponds to mode in same column in phi
%   phi       = Modal amplitudes (vibration modes) scaled to unit mass
%               Size N x tr
%               Each column are modal amplitudes for one mode
%               Each row contains model amplitude of dof with global id (in
%               DOF) equal to that row
%   scaler    = Parameter controlling the plot scale
%               Scalar value
%   n_loops   = Number of loops
%               Scalar value
%   original  = Plot original undeformed mesh
%               Boolean value [0, 1]
% =========================================================================

%%% Input check (9 inputs)
narginchk(9,9);

%%% Phi* definition: Includes the constrained degrees of freedom in
%   the definition of the modal amplitudes (including them as zeros)

phi_star    = zeros(3*size(DOF,1),1); 

for dof   = 1:size(DOF,1)
   if DOF(dof,1) == 0 && DOF(dof,2) == 0 && DOF(dof,3) == 0
       
       phi_star(DOFfree(dof,1))     = 0;
       phi_star(DOFfree(dof,2))     = 0;
       phi_star(DOFfree(dof,3))     = 0;
   else
       
       phi_star(DOFfree(dof,1))     = phi(DOF(dof,1),mode);
       phi_star(DOFfree(dof,2))     = phi(DOF(dof,2),mode);  
       phi_star(DOFfree(dof,3))     = phi(DOF(dof,3),mode);  
     
   end  
end

%%% Nodal displacement computation

def         = zeros(size(mesh.nodes));
tot_disp    = def;

delta_t = [1:-0.075:-1, -1:0.075:1];

n_current_loops = 0;

h = figure;
filename = ['Mode_', num2str(mode),'.gif'];

while n_current_loops < n_loops
    
    for ii_t = 1:length(delta_t)
        
             
        for nn=1:size(mesh.nodes,1)

            def(nn,1)       = mesh.nodes(nn,1) + delta_t(ii_t) * phi_star(DOFfree(nn,1)) * scaler;
            def(nn,2)       = mesh.nodes(nn,2) + delta_t(ii_t) * phi_star(DOFfree(nn,2)) * scaler;
            def(nn,3)       = mesh.nodes(nn,3) + delta_t(ii_t) * phi_star(DOFfree(nn,3)) * scaler;
            tot_disp(nn)    = norm( [ delta_t(ii_t) * phi_star(DOFfree(nn,1)) ...
                                      delta_t(ii_t) * phi_star(DOFfree(nn,2)) ...
                                      delta_t(ii_t) * phi_star(DOFfree(nn,3))] );

        end

        %%% Intermediate auxiliar definitions
        modedef     = def.';
        nodeT       = mesh.nodes.';

        %%% Plot 
        plotFEM3D(modedef,mesh.top,'ColorMapData',tot_disp(:,1),'EdgeColor',[0.7 0.7 0.7],'Mesh','on')

        %%% Plot original mesh
        if original && (delta_t(ii_t) ~= 0)
            hold on
            plotFEM3D(nodeT,mesh.top(1:4,:),'EdgeColor',[0.9 0.9 0.9],'FaceAlpha',0.5)
        end 

        set(gca,'fontname','times')
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        daspect([1 1 1]);
        xlabel('Width [m]')
        ylabel('Height [m]')
        title({['Mode ', num2str(mode),'.    f = ',num2str(wnHz(mode),6),' Hz'],...
            '\fontsize{12}\color{gray}Colorbar: Modal displacement scaled to unit mass'},...
            'fontsize',14)
        caxis([0 0.12])
        xlim([min(mesh.nodes(:,1)) max(mesh.nodes(:,1))])
        ylim([min(mesh.nodes(:,2)) max(mesh.nodes(:,2))])
        zlim([min(mesh.nodes(:,3)) max(mesh.nodes(:,3))])
        hold off
        pause(.000001)
        
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if ii_t == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        end 
    end
    
    n_current_loops = n_current_loops + 1;
    
    
end
    
end

