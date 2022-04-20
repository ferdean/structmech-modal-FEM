function [M, K, DOF, N, DOFfree] = fastMatrices(mesh,E,v,rho,eta,constrNodes, verbose)
%
% Fast mass and stifness matrix computation for tetrahedral elements
%
% [ref]  Cuvelier, François, Caroline Japhet, and Gilles Scarella. 
%        "An efficient way to perform the assembly of finite element 
%        matrices in Matlab and Octave." 
%        arXiv preprint arXiv:1305.3122 (2013) 
% =========================================================================
% Child functions: computeB
% =========================================================================
% Created by:   Ferran de Andrés(8.2021) 
% Version: 1.0. First implementation for tetrahedral elements
% =========================================================================
% INPUT: 
%   mesh        = Struct with all the meshing information listed below: 
%      nodes    = Point matrix.
%                 nN x 3 matrix, where nN is the number of nodes 
%                 X-Y-Z
%      top      = Topology matrix. The first nn_in_E rows contain indices 
%                 to the corner points (nn_in_E are the number of nodes per 
%                 element) given in counter clockwise order, and the 
%                 (nn_in_E + 1) row contains the subdomain number.
%   E           = Young modulus [N/mm2]
%                 Nº subdomains x 1 matrix 
%   v           = Poisson's ratio [-]
%                 Nº subdomains x 1 matrix
%   rho         = Density of the material [kg/mm3]
%                 Nº subdomains x 1 matrix
%   eta         = Damping loss factor [-]
%                 Nº subdomains x 1 matrix
%   constrNodes = Id of constrained nodes
%                 nConstr x 1
%
% OUTPUT:
%   M           = Mass matrix (Dof x Dof)
%   K           = Stiffness matrix (Dof x Dof)
%   DOF         = Matrix of dof's considering no constained nodes
%                 Size nn x Nn, with Nn nº dof's per node
%   N           = Number of dof's in solid
%                 Scalar value
%   DOFfree     = Matrix of dof's considering no constained nodes
%                 Size nn x Nn, with Nn nº dof's per node
% =========================================================================
% NOTE: the mass and stiffness matrices are defined for a problem in which 
% each node has 3 DOF's (displacement), in such a way that the displacement 
% vector has the following shape: 
% U        = [u1 v1 w1 u2 v2 w2 ...  un vn wn]'
% =========================================================================
 
%%% Input check (7 inputs)
narginchk(7, 7);

%%% Topological information
nn          = size(mesh.nodes,1);           % Number of nodes  
nE          = size(mesh.top,2);             % Number of elements
nDOFn       = 3;                            % Number of DOF in each node 

%%% Element type (updated ONLY for tetrahedrons)
if size(mesh.top,1) == 4
        elementType     = "tet4";           % Linear or quadratic element
        
        nn_in_E         = 4;                % Nodes per element
        nGauss          = 2;                % Integration points per dimension
        multimaterial   = 0;                

elseif size(mesh.top,1) == 5
        elementType     = "tet4";
        
        nn_in_E         = 4;
        nGauss          = 2;
        multimaterial   = 1;
        
elseif size(mesh.top,1) == 10
        elementType     = "tet10";
        
        nn_in_E         = 10;
        nGauss          = 3;
        multimaterial   = 0;
        
elseif size(mesh.top,1) == 11
        elementType     = "tet10";
        
        nn_in_E         = 10;
        nGauss          = 3;
        multimaterial   = 1;
        
elseif size(mesh.top,1) == 8
        elementType     = "hex8";
        
        nn_in_E         = 8;
        nGauss          = 2;
        multimaterial   = 0;
        
elseif size(mesh.top,1) == 9
        elementType     = "hex8";
        
        nn_in_E         = 8;
        nGauss          = 2;
        multimaterial   = 1;
else
    error('The mesh is not formed by tetrahedral (not linear nor quadratic) elements')
end

%%% Element shape functions
[N, dNdloc] = SF(elementType);

%%% Element topological information and DOF matrices
nDOFe       = nDOFn*nn_in_E;                % Number of DOF in each element
nDOF_glob   = nDOFn*nn;                     % Number of global DOF

DOFfree     = reshape((1:nDOF_glob).',nDOFn,nn).';

%%% Gauss table
[xi_0, eta_0, tau_0, w]    = gaussTable(nGauss,elementType);
NPG                         = size(xi_0,1);

%%% Program initialisation (Fast assembly)

Ig  = zeros(nDOFe^2+nDOFe^2*(nE-1), 1);
Jg  = Ig;

Mg  = zeros(nDOFe^2+nDOFe^2*(nE-1), 1);
Kg  = Mg;

%%% Material data
if ~ multimaterial % If there's info about the material 
        E_eff   = E(1);
        v_eff   = v(1);
        rho_eff = rho(1);
        eta_eff = eta(1);
        
        %%% Material stiffness matrix
        k_D = E_eff/((1 + v_eff)*( 1 - 2* v_eff));

        D   = k_D * [1 - v_eff, v_eff, v_eff, 0, 0, 0;...
                     v_eff, 1 - v_eff, v_eff  0, 0, 0;...
                     v_eff, v_eff, 1 - v_eff, 0, 0, 0;...
                     0, 0, 0, 1/2*(1 - 2*v_eff), 0, 0;...
                     0, 0, 0, 0, 1/2*(1 - 2*v_eff), 0;...
                     0, 0, 0, 0, 0, 1/2*(1 - 2*v_eff)];
end

%%% Main mesh loop
for e = 1:nE
    
    %%% Topological information of the element
    i_el            = mesh.top(1:nn_in_E,e);
    elementDOF      = DOFfree(i_el,:);
    node_coord      = mesh.nodes(i_el,:);
    
    if multimaterial 
        E_eff   = E(mesh.top(nn_in_E + 1,e));
        v_eff   = v(mesh.top(nn_in_E + 1,e));
        rho_eff = rho(mesh.top(nn_in_E + 1,e));
        eta_eff = eta(mesh.top(nn_in_E + 1,e));
        
        %%% Material stiffness matrix
        k_D = E_eff/((1 + v_eff)*( 1 - 2* v_eff));

        D   = k_D * [1 - v_eff, v_eff, v_eff, 0, 0, 0;...
                     v_eff, 1 - v_eff, v_eff  0, 0, 0;...
                     v_eff, v_eff, 1 - v_eff, 0, 0, 0;...
                     0, 0, 0, 1/2*(1 - 2*v_eff), 0, 0;...
                     0, 0, 0, 0, 1/2*(1 - 2*v_eff), 0;...
                     0, 0, 0, 0, 0, 1/2*(1 - 2*v_eff)];
    end

    
    %%% Element mass matrix
    me  = zeros(numel(elementDOF));
    ke  = zeros(numel(elementDOF));
    % kce = zeros(numel(elementDOF));
    

    %%% Numerical integration loop  (Gauss quadrature)       
    for it = 1:NPG      
    
        %%% Shape functions and its local derivatives
        N_0         = N(xi_0(it), eta_0(it), tau_0(it));
        dNdloc_0    = dNdloc(xi_0(it), eta_0(it), tau_0(it));
        
        %%% Local to global coordinates
        [Jaco,~]    = jacobian(node_coord,dNdloc_0);
        J           = det(Jaco);
        
        dNdglob     = Jaco\dNdloc_0;  % Global derivatives
        
        %%% Element deformation matrix (B)
        B           = computeB(dNdglob,elementType);   
        
        %%% Element matrices (global)
        me = me +  w(it) * ( N_0.'*N_0) * J * rho_eff;
        ke = ke +  w(it) * ((B.' * D * B ) * J) * (1 + 1i * eta_eff);  
        
    end
    
    %%% Matrix assembly [ref]

    T   = meshgrid(1:nDOFe).';
    Tt  = T.';
    ii  = T(:);
    jj  = Tt(:);
    kk  = (1:nDOFe^2).';
    kkk = kk + nDOFe^2*(e-1);

    DOF_glob   = reshape(elementDOF.', nDOFe, 1);
    Ig(kkk)    = DOF_glob(ii);
    Jg(kkk)    = DOF_glob(jj);
    Mg(kkk)    = me(:);
    Kg(kkk)    = ke(:);

    %%% Progress bar
    if verbose && (e == 1 || mod(e, 250) == 0 || e == nE)
        clc;    
        fprintf('Assembly progress: %1.2f %%\n', e/nE*100);
    end
end

MFree = sparse(Ig, Jg, Mg, nDOF_glob, nDOF_glob);
KFree = sparse(Ig, Jg, Kg, nDOF_glob, nDOF_glob);

%%% Boundary conditions
%   DOF real
DOF                 = zeros(nn, nDOFn);
freeNodes           = setdiff( (1:nn).', constrNodes);
nFreeNodes          = length(freeNodes);
N                   = nDOFn*nFreeNodes;
DOF(freeNodes, :)   = reshape( (1:N).', nDOFn, nFreeNodes).';
                            
%   K and M real
freeDofsPrev        = sort(reshape(DOFfree(freeNodes, :).', [], 1));
K                   = KFree(freeDofsPrev, freeDofsPrev);
M                   = MFree(freeDofsPrev, freeDofsPrev);
end

function B = computeB(dNdglob,elementType)

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
end

