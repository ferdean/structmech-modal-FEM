function [M, K] = freeMatrices3D(mesh,elementType,E,v,rho,eta,DOF)
%
% Mass and stifness matrix computation (sparse output)
% =========================================================================
% Created by:   Ferran de Andrés(7.2021) 
% Version: 1.0. First implementation for continuous media (same density in
%               all the geometry)
%          1.1  Addition of multi-domain studies for hex8 and tet4 problems
%          1.2  Implementation of quadratic element problems
% =========================================================================
% INPUT: 
%   mesh     = Struct with all the meshing information listed below: 
%      nodes = Point matrix.
%              nN x 3 matrix, where nN is the number of nodes 
%              X-Y-Z
%      top   = Topology matrix. The first nn_in_E rows contain indices to the
%              corner points (nn_in_E are the number of nodes per element)
%              given in counter clockwise order, and the (nn_in_E + 1) row 
%              contains the subdomain number.
%   elementtype    = Type of element
%                    Char value
%
%                    'tet4'   - Linear tetraedron (3D)
%                    'tet10'  - Quadratic tetraedron (3D)
%                    'hex8'   - Linear hexaedron (3D)
%                    'hex20'  - Quadratic hexaedron (3D)
%   E        = Young modulus [N/mm2]
%              Nº subdomains x 1 matrix 
%   v        = Poisson's ratio [-]
%              Nº subdomains x 1 matrix
%   rho      = Density of the material [kg/mm3]
%              Nº subdomains x 1 matrix
%   eta      = Damping loss factor [-]
%              Nº subdomains x 1 matrix
%   DOF      = Matrix of dof's considering no constained nodes
%              Size nn x Nn, with Nn nº dof's per node
%
% OUTPUT:
%   M        = Mass matrix (Dof x Dof)
%   K        = Stiffness matrix (Dof x Dof)
%   Kc       = Structural damping stiffness matrix (Dof x Dof)
% =========================================================================
% NOTE: the mass and stiffness matrices are defined for a problem in which 
% each node has 3 DOF's (displacement), in such a way that the displacement 
% vector has the following shape: 
% U        = [u1 v1 w1 u2 v2 w2 ...  un vn wn]'
% =========================================================================
 
%%% Input check (7 inputs)
narginchk(7,7);

%%% Topological information

nn          = size(mesh.nodes,1);           % Number of nodes  
nE          = size(mesh.top,2);             % Number of elements
nDOFn       = 3;                            % Number of DOF in each node 

switch elementType
    case "hex8"
        nn_in_E     = 8;                    % Number of nodes per element
        [N, dNdloc] = SF("hex8");           % Element shape functions
        
        nGauss      = 2;                    % Number of gauss points 
    
    case  "hex20"
        nn_in_E     = 20;                   
        [N, dNdloc] = SF("hex20");
        
        nGauss      = 3;
        
    case "tet4"
        nn_in_E     = 4; 
        [N, dNdloc] = SF("tet4");
        
        nGauss      = 2;
        
    case "tet10"
        nn_in_E     = 10; 
        [N, dNdloc] = SF("tet10");
        
        nGauss      = 3;        
end
        
nDOFe       = nDOFn*nn_in_E;                % Number of DOF in each element
nDOF_glob   = nDOFn*nn;                     % Number of global DOF

DOF_local   = reshape((1:nDOFe).',nDOFn,nn_in_E)';


%%% Gauss table
[xi_0, eta_0, tau_0, w]    = gaussTable(nGauss,elementType);
NPG                         = size(xi_0,1);

%%% Program initialisation
M = zeros(nDOF_glob);
K = zeros(nDOF_glob);


for e = 1:nE
    
    %%% Topological information of the element
    i_el            = mesh.top(1:nn_in_E,e);
    elementDOF      = DOF(i_el,:);
    node_coord      = mesh.nodes(i_el,:);
    
    if size(mesh.top,1) == (nn_in_E + 1) % If there's info about the material 
        E_eff   = E(mesh.top(nn_in_E + 1,e));
        v_eff   = v(mesh.top(nn_in_E + 1,e));
        rho_eff = rho(mesh.top(nn_in_E + 1,e));
        eta_eff = eta(mesh.top(nn_in_E + 1,e));
        
    else                               
        E_eff   = E(1);
        v_eff   = v(1);
        rho_eff = rho(1);
        eta_eff = eta(1);
        
    end
    
    %%% Element mass matrix
    me  = zeros(numel(elementDOF));
    ke  = zeros(numel(elementDOF));
    % kce = zeros(numel(elementDOF));
    
    %%% Material stiffness matrix
    k_D = E_eff/((1 + v_eff)*( 1 - 2* v_eff));
    
    D   = k_D * [1 - v_eff, v_eff, v_eff, 0, 0, 0;...
                 v_eff, 1 - v_eff, v_eff  0, 0, 0;...
                 v_eff, v_eff, 1 - v_eff, 0, 0, 0;...
                 0, 0, 0, 1/2*(1 - 2*v_eff), 0, 0;...
                 0, 0, 0, 0, 1/2*(1 - 2*v_eff), 0;...
                 0, 0, 0, 0, 0, 1/2*(1 - 2*v_eff)];
             
    for it = 1:NPG      % Gauss quadrature
    
        %%% Shape functions and its local derivatives
        N_0         = N(xi_0(it), eta_0(it), tau_0(it));
        dNdloc_0    = dNdloc(xi_0(it), eta_0(it), tau_0(it));
        
        %%% Local to global coordinates
        [Jaco,~]    = jacobian(node_coord,dNdloc_0);
        J           = det(Jaco);
        
        dNdglob     = Jaco\dNdloc_0;  % Global derivatives
        
        %%% Element deformation matrix (B)
        B1 = [dNdglob(1,1), 0, 0;...
              0, dNdglob(2,1), 0;...
              0, 0, dNdglob(3,1);...
              dNdglob(2,1), dNdglob(1,1), 0;...
              dNdglob(3,1), 0, dNdglob(1,1);...
              0, dNdglob(3,1), dNdglob(2,1)];
          
        B2 = [dNdglob(1,2), 0, 0;...
              0, dNdglob(2,2), 0;...
              0, 0, dNdglob(3,2);...
              dNdglob(2,2), dNdglob(1,2), 0;...
              dNdglob(3,2), 0, dNdglob(1,2);...
              0, dNdglob(3,2), dNdglob(2,2)];
          
        B3 = [dNdglob(1,3), 0, 0;...
              0, dNdglob(2,3), 0;...
              0, 0, dNdglob(3,3);...
              dNdglob(2,3), dNdglob(1,3), 0;...
              dNdglob(3,3), 0, dNdglob(1,3);...
              0, dNdglob(3,3), dNdglob(2,3)];
          
        B4 = [dNdglob(1,4), 0, 0;...
              0, dNdglob(2,4), 0;...
              0, 0, dNdglob(3,4);...
              dNdglob(2,4), dNdglob(1,4), 0;...
              dNdglob(3,4), 0, dNdglob(1,4);...
              0, dNdglob(3,4), dNdglob(2,4)];
        
        B = [B1 B2 B3 B4];
        
        if elementType == "hex8"   
        
            B5 = [dNdglob(1,5), 0, 0;...
                  0, dNdglob(2,5), 0;...
                  0, 0, dNdglob(3,5);...
                  dNdglob(2,5), dNdglob(1,5), 0;...
                  dNdglob(3,5), 0, dNdglob(1,5);...
                  0, dNdglob(3,5), dNdglob(2,5)];

            B6 = [dNdglob(1,6), 0, 0;...
                  0, dNdglob(2,6), 0;...
                  0, 0, dNdglob(3,6);...
                  dNdglob(2,6), dNdglob(1,6), 0;...
                  dNdglob(3,6), 0, dNdglob(1,6);...
                  0, dNdglob(3,6), dNdglob(2,6)];

            B7 = [dNdglob(1,7), 0, 0;...
                  0, dNdglob(2,7), 0;...
                  0, 0, dNdglob(3,7);...
                  dNdglob(2,7), dNdglob(1,7), 0;...
                  dNdglob(3,7), 0, dNdglob(1,7);...
                  0, dNdglob(3,7), dNdglob(2,7)];


            B8 = [dNdglob(1,8), 0, 0;...
                  0, dNdglob(2,8), 0;...
                  0, 0, dNdglob(3,8);...
                  dNdglob(2,8), dNdglob(1,8), 0;...
                  dNdglob(3,8), 0, dNdglob(1,8);...
                  0, dNdglob(3,8), dNdglob(2,8)];
              
            B = [B B5 B6 B7 B8];
            
        elseif elementType == "tet10"
            
            B5 = [dNdglob(1,5), 0, 0;...
                  0, dNdglob(2,5), 0;...
                  0, 0, dNdglob(3,5);...
                  dNdglob(2,5), dNdglob(1,5), 0;...
                  dNdglob(3,5), 0, dNdglob(1,5);...
                  0, dNdglob(3,5), dNdglob(2,5)];

            B6 = [dNdglob(1,6), 0, 0;...
                  0, dNdglob(2,6), 0;...
                  0, 0, dNdglob(3,6);...
                  dNdglob(2,6), dNdglob(1,6), 0;...
                  dNdglob(3,6), 0, dNdglob(1,6);...
                  0, dNdglob(3,6), dNdglob(2,6)];

            B7 = [dNdglob(1,7), 0, 0;...
                  0, dNdglob(2,7), 0;...
                  0, 0, dNdglob(3,7);...
                  dNdglob(2,7), dNdglob(1,7), 0;...
                  dNdglob(3,7), 0, dNdglob(1,7);...
                  0, dNdglob(3,7), dNdglob(2,7)];


            B8 = [dNdglob(1,8), 0, 0;...
                  0, dNdglob(2,8), 0;...
                  0, 0, dNdglob(3,8);...
                  dNdglob(2,8), dNdglob(1,8), 0;...
                  dNdglob(3,8), 0, dNdglob(1,8);...
                  0, dNdglob(3,8), dNdglob(2,8)];
              
            B9 = [dNdglob(1,9), 0, 0;...
                  0, dNdglob(2,9), 0;...
                  0, 0, dNdglob(3,9);...
                  dNdglob(2,9), dNdglob(1,9), 0;...
                  dNdglob(3,9), 0, dNdglob(1,9);...
                  0, dNdglob(3,9), dNdglob(2,9)];
            
            B10= [dNdglob(1,10), 0, 0;...
                  0, dNdglob(2,10), 0;...
                  0, 0, dNdglob(3,10);...
                  dNdglob(2,10), dNdglob(1,10), 0;...
                  dNdglob(3,10), 0, dNdglob(1,10);...
                  0, dNdglob(3,10), dNdglob(2,10)];
              
            B = [B B5 B6 B7 B8 B9 B10];
            
            
        end   
          
        me = me +  w(it) * ( N_0.'*N_0) * J * rho_eff;
        ke = ke +  w(it) * ((B.' * D * B ) * J);  
          
    end
    
       %%% Matrix assembly
   
      for i_1 = 1 : nn_in_E
        
        dof1_local       = DOF_local(i_1, :);        
        dof1_global      = elementDOF(i_1, :);       
        
        for i_2 = 1 : nn_in_E
            dof2_local       = DOF_local(i_2, :);    
            dof2_global      = elementDOF(i_2, :);   
            
            %%% Mass matrix
            M(dof1_global, dof2_global) = ...
                M(dof1_global, dof2_global) + ...
                me(dof1_local, dof2_local);
          
            %%% Stiffness matrix
            K(dof1_global, dof2_global) = ...
                K(dof1_global, dof2_global) + ...
                ke(dof1_local, dof2_local);

%             Kc(dof1_global, dof2_global) = ...
%                 Kc(dof1_global, dof2_global) + ...
%                 kce(dof1_local, dof2_local);
            
        end
      end     
    clc;
    fprintf('Progress bar: %1.3f %%\n', e/nE*100); 
    
end

M = sparse(M);
K = sparse(K);

end

