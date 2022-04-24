function varargout = FRF(Mq,Kq,DOF,nDOF,phi,nN,node_ex,node_rep, verbose)
% Computes the frequency response function (FRF) given the modal
% transformed mass and stiffness matrices of the system in the selected
% observed/excited node.
% =========================================================================
% Created by:   Ferran de Andrés (7.2021) 
% =========================================================================
% INPUT: 
%   Mq        = Mass matrix after modal transformation (tr x tr)
%   Kq        = Stiffness matrix after modal transformation (tr x tr)
%   DOF       = Matrix of dof's considering constained nodes
%               Size nn x Nn, with Nn nº dof's per node
%   nDOF      = Number of active degrees of freedom
%               Scalar.
%   phi       = Modal amplitudes (vibration modes) scaled to unit mass
%               Size N x tr
%               Each column are modal amplitudes for one mode
%               Each row contains model amplitude of dof with global id (in
%               MGDL) equal to that row
%   nN        = 2D vector containing the ID of the excited node (position 1) 
%               and the observed node (position 2), such that
%               nN = [number_of_exited_nodem, number_of_observed_node]
%   node_i    = String specifying the direction (local) of exitation 
%                 * 'axial'
%                 * 'radial'
%                 * 'torsional'
% =========================================================================


nN_ex       = nN(1);
nN_obs      = nN(2);

i_ex_a      = DOF(nN_ex,1);     % Index of the axial coordinate at the contact node
i_ex_r      = DOF(nN_ex,2);     % Index of the radial coordinate at the contact node
i_ex_t      = DOF(nN_ex,3);     % Index of the torsional coordinate at the contact node

i_obs_a     = DOF(nN_obs,1);
i_obs_r     = DOF(nN_obs,2);
i_obs_t     = DOF(nN_obs,3);


% Frequency domain
wHz = 10:2:7000;

switch node_ex
    case 'axial'
        Fa      = 1 * ones(length(wHz),1);
        Fr      = 0 * ones(length(wHz),1);
        Ft      = 0 * ones(length(wHz),1);
        
    case 'radial'
        Fa      = 0 * ones(length(wHz),1);
        Fr      = 1 * ones(length(wHz),1);
        Ft      = 0 * ones(length(wHz),1);
    
    case 'torsional'
        Fa      = 0 * ones(length(wHz),1);
        Fr      = 0 * ones(length(wHz),1);
        Ft      = 1 * ones(length(wHz),1);       
end

switch node_rep
    case 'axial'
        i_rep = i_obs_a;
        
    case 'radial'
        i_rep = i_obs_r;
    
    case 'torsional'        
        i_rep = i_obs_t;     
end

uDOF    = zeros(nDOF,length(wHz));
% vDOF    = zeros(nDOF,length(wHz));


for i_w = 1:length(wHz)
    
        F_GDL = zeros(nDOF, 1);

        F_GDL(i_ex_a) = Fa(i_w);
        F_GDL(i_ex_r) = Fr(i_w);
        F_GDL(i_ex_t) = Ft(i_w);

        Qq = phi.' * F_GDL;
    
        w = wHz(i_w) * 2 * pi;

        qModos = (-w^2*Mq + Kq ) \ Qq; % nModos x 1
    
        uDOF(:,i_w) = phi * qModos; % N x 1
        
%         vDOF(:,i_w) = 1i * w * uDOF(:,i_w);
        
    
        if verbose && (i_w == 1 || mod(i_w, 100) == 0 || i_w == length(wHz))
            clc;    
            fprintf('Progress bar: %1.2f %%\n', i_w/length(wHz)*100)
        end
        
end

fprintf('Results are available')

if nargout == 2
    varargout{1} = uDOF;
    varargout{2} = wHz;
end


subplot(2,1,1)
semilogy(wHz,abs(uDOF(i_rep,:)),'k')
% legend({'Axial','Radial'},'Interpreter','latex')
haxis = gca; 
% Font size
FS    = 13; 
% Axis format
set(gca, 'FontName', 'TimesNewRoman');
hYLabel          = ylabel('$\mid  u(f) \mid$','Interpreter','latex');
hXLabel          = xlabel('Frecuencia (f) [Hz]','Interpreter','latex');
haxis.Box        = 'on';
haxis.TickDir    = 'out';
haxis.TickLength = [.02 .02];
haxis.XMinorTick = 'on';
haxis.YMinorTick = 'on';
haxis.YGrid      = 'on';
haxis.XMinorGrid = 'on';
haxis.YMinorGrid = 'on';
haxis.XGrid      = 'on';
haxis.XColor     = [.0 .0 .0];
haxis.YColor     = [.0 .0 .0];
haxis.LineWidth  = 0.75;
haxis.FontName   = 'Palatine';
haxis.FontSize   = 0.65*FS;
hXLabel.FontSize = 0.85*FS;
hYLabel.FontSize = 0.85*FS;
% ylim([min(abs(FRF_Y)) max(abs(FRF_Y))])
% xlim([10 max(w)])
% pbaspect([(1+sqrt(5))/2 1 1]) % Aurea proportion


subplot(2,1,2)
plot(wHz,angle(uDOF(i_rep,:)),'k')
%legend({'Sin capa','Modelo FEM con capa','Modelo analítico con capa'},'Interpreter','latex')
haxis = gca; 
% Font size
FS    = 13; 
% Axis format
set(gca, 'FontName', 'TimesNewRoman');
hYLabel          = ylabel('$\alpha(f)$','Interpreter','latex');
hXLabel          = xlabel('Frecuencia (f) [Hz]','Interpreter','latex');
haxis.Box        = 'on';
haxis.TickDir    = 'out';
haxis.TickLength = [.02 .02];
haxis.XMinorTick = 'on';
haxis.YMinorTick = 'on';
haxis.YGrid      = 'on';
haxis.XMinorGrid = 'on';
haxis.YMinorGrid = 'on';
haxis.XGrid      = 'on';
haxis.XColor     = [.0 .0 .0];
haxis.YColor     = [.0 .0 .0];
haxis.LineWidth  = 0.75;
haxis.FontName   = 'Palatine';
haxis.FontSize   = 0.65*FS;
hXLabel.FontSize = 0.85*FS;
hYLabel.FontSize = 0.85*FS;