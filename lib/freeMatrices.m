function [M,K,T1,T2,T3]=freeMatrices(E,v,rho,eta,t,N,DSF,mesh,dof)
% Mass and stifness matrix computation
% =========================================================================
% Created by:   Ferran de Andr�s(1.2021) 
% Version: 1.0. First implementation for continuous media (same density in
%               all the geometry)
% =========================================================================
% INPUT: 
%   E        = Young modulus [N/mm2]
%              Scalar value
%   v        = Poisson's ratio [-]
%              Scalar value
%   rho      = Density of the material [kg/mm3]
%              Scalar value
%   eta      = Damping loss factor [-]
%              Scalar value
%   t        = Thickness (2D problems) [mm]
%              Scalar value
%   SF       = Shape Functions
%              Function handle
%   DSF      = Natural derivatives of shape functions (wrt local
%              coordinates)
%   mesh     = Struct with all the meshing information listed below: 
%      nodes = Point matrix. The first and second rows contain x- and 
%              y-coordinates of the points in the mesh.
%      top   = Triangle topology matrix. The first three rows contain 
%              indices to the corner points, given in counter clockwise 
%              order, and the fourth row contains the subdomain number.
%      e     = Mesh edges, returned as a 7-by-Ne matrix, where Ne is the 
%              number of boundary edges in the mesh. 
%              An edge is a pair of points in p containing a boundary 
%              between subdomains, or containing an outer boundary
%   dof      = Matrix of dof's considering no constained nodes
%              Size nn x Nn, with Nn n� dof's per node
% OUTPUT:
%   M        = Mass matrix (Dof x Dof)
%   K        = Stiffness matrix (Dof x Dof)
%   H        = Structural damping matrix (Dof x Dof)
%   Ti       = Computation times (1 x nE)
% =========================================================================
% NOTE: the mass and stiffness matrices are defined for a problem in which 
% each node has 2 Dof (horizontal and vertical displacement), in such a way 
% that the displacement vector has the following shape: 
% U        = [u1 v1 u2 v2 ...  un vn]'
% =========================================================================

% syms eta xi
%%% Input check (9 inputs)
narginchk(9,9);

%%% Mesh characteristics
nn          = size(mesh.nodes,1);
nE          = size(mesh.top,2);
nDof_glob   = 1*nn;
nn_in_E     = 3;            % N� of nodes in an element
NPG         = 2;            % Number of integration points

nDofn       = 1;
nodes_el    = 3;
nDofEl      = nDofn*nodes_el;

dof_local   = reshape((1:nDofEl).',nDofn,nodes_el)';
    
%%% Gauss table
[xi_0, eta_0, w] = gaussTable_Tri(NPG);

%%% Program initialisation
M = zeros(nDof_glob);
K = zeros(nDof_glob);

if size(mesh.top,1) == 4
    multimaterial = 1;
else
    multimaterial = 0;
end



T1 = zeros(nE,1);
T2 = T1;
T3 = T2;



if ~ multimaterial % If there's info about the material 
        E_eff   = E(1);
        v_eff   = v(1);
        rho_eff = rho(1);
        eta_eff = eta(1);

        D = E_eff/(1-v_eff^2)  * [1 v_eff 0 ;  v_eff 1 0; 0 0 (1- v_eff)/2];
        
end


%%% Mass matrix
for e=1:nE
    tic
    
    i           = mesh.top(1:3,e);
    elementDof  = dof(i,:);
    node_coord  = mesh.nodes(i,:);

     if multimaterial 
        E_eff   = E(mesh.top(nn_in_E + 1,e));
        v_eff   = v(mesh.top(nn_in_E + 1,e));
        rho_eff = rho(mesh.top(nn_in_E + 1,e));
        eta_eff = eta(mesh.top(nn_in_E + 1,e));
        
        D = E_eff/(1-v_eff^2)  * [1 v_eff 0 ;  v_eff 1 0; 0 0 (1- v_eff)/2];
     end
     
     
    %%% Element mass matrix
    me  = zeros(numel(elementDof));
    ke  = zeros(numel(elementDof));
    
    for it = 1:NPG      % Gauss quadrature
        
        if isa(DSF,'function_handle')
            dNdloc_0=DSF(eta_0(it),xi_0(it));
        else
            dNdloc_0=DSF;
        end
        
        N_0         = N(eta_0(it),xi_0(it));
        
        
        [Jaco,~]    = jacobian(node_coord,dNdloc_0);
        J           = det(Jaco);
        
        dNdglob     = Jaco\dNdloc_0; % Global derivatives
        
%         B           = zeros(3,(2*nn_in_E));
%         
%         B(1,1:2:(2*nn_in_E)) = dNdglob(1,:);
%         B(2,2:2:(2*nn_in_E)) = dNdglob(2,:);
%         B(3,1:2:(2*nn_in_E)) = dNdglob(2,:);
%         B(3,2:2:(2*nn_in_E)) = dNdglob(1,:);
        
        % [dN1/dx    0     dN2/dx     0    dN3/dx    0     ;
        %    0     dN1/dy    0     dN2/dy    0     dN3/dy  ;
        %  dN1/dy  dN1/dx  dN2/dy  dN2/dx  dN3/dy  dN3/dy ];

        T1(e)=toc;
        
        me = me +  w(it) * ( N_0.'*N_0) * J * rho_eff * t;
%         ke = ke +  w(it) * ((B.' * D * B ) * J * t);
    end
    
%         ke = ke .* (1 +  1i*eta_eff);
    
        T2(e)=toc;

%    %%% Matrix assembly
%    for i_1 = 1 : nodes_el
%         
%         dof1_local       = dof_local(i_1, :);        
%         dof1_global      = elementDof(i_1, :);       
%         
%         for i_2 = 1 : nodes_el
%             dof2_local       = dof_local(i_2, :);    
%             dof2_global      = elementDof(i_2, :);   
%             
%             %%% Mass matrix
%             M(dof1_global, dof2_global) = ...
%                 M(dof1_global, dof2_global) + ...
%                 me(dof1_local, dof2_local);
%           
% %             %%% Stiffness matrix
% %             K(dof1_global, dof2_global) = ...
% %                 K(dof1_global, dof2_global) + ...
% %                 ke(dof1_local, dof2_local);
% 
%         end
%     end
    
    %%% Matrix assembly [ref]

    T   = meshgrid(1:nDofEl).';
    Tt  = T.';
    ii  = T(:);
    jj  = Tt(:);
    kk  = (1:nDofEl ^2).';
    kkk = kk + nDofEl ^2*(e-1);

    DOF_glob   = reshape(elementDof.', nDofEl, 1);
    Ig(kkk)    = DOF_glob(ii);
    Jg(kkk)    = DOF_glob(jj);
    Mg(kkk)    = me(:);
%     Kg(kkk)    = ke(:);


    T3(e)=toc;
end
MFree = sparse(Ig, Jg, Mg, nn, nn);
% KFree = sparse(Ig, Jg, Kg, nDOF_glob, nDOF_glob);
%%% Mass matrix check: it should be symmetric positive definite. 
%   
%   The most efficient method to check whether a matrix is symmetric 
%   positive definite is to simply attempt to use chol on the matrix. 
%   If the factorization fails, then the matrix is not symmetric positive 
%   definite. This method does not require the matrix to be symmetric for a 
%   successful test (if the matrix is not symmetric, then the factorization 
%   fails).

% try chol(M);
% catch
%     error('The computed mass  matrix is not symmetric positive definite')
% end
% try chol(K);
% catch
%     error('The computed stiffness matrix is not symmetric positive definite')
% end
end