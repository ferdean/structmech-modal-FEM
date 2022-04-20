function [M, K, DOF, N, DOFfree] = computeAxisymMatrices(mesh,E,v,rho,eta,constrNodes,nodalDiam)
% Computes axisymmetric 2D matrices
% =========================================================================
% Created by:   Ferran de Andrés(5.2021) 
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
%              Size nn x Nn, with Nn nº dof's per node
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

%%% Input check (7 inputs)
narginchk(7,7);

%%% Topological information
nn          = size(mesh.nodes,1);           % Number of nodes  
nE          = size(mesh.top,2);             % Number of elements
nDOFn       = 3;                            % Number of DOF in each node 

%%% Element type (updated ONLY for tetrahedrons)
if size(mesh.top,1) == 3
        elementType     = "tri3";           % Linear or quadratic element
        
        nn_in_E         = 3;                % Nodes per element
        nGauss          = 3;                % Integration points per dimension
        multimaterial   = 0;                

elseif size(mesh.top,1) == 4
        elementType     = "tri3";
        
        nn_in_E         = 3;
        nGauss          = 3;
        multimaterial   = 1;
        
elseif size(mesh.top,1) == 6
        elementType     = "tri6";
        
        nn_in_E         = 6;
        nGauss          = 3;
        multimaterial   = 0;
        
elseif size(mesh.top,1) == 7
        elementType     = "tri6";
        
        nn_in_E         = 6;
        nGauss          = 3;
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
[xi_0, eta_0, w]            = gaussTable_Tri_Victor(nGauss);
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
        k_D_1 = E_eff/((1 + v_eff)*( 1 - 2* v_eff));
        k_D_2 = E_eff/(2*(1 + v_eff));

        D_1   = k_D_1 * [1 - v_eff, v_eff, v_eff, 0;...
                       v_eff, 1 - v_eff, v_eff  0;...
                       v_eff, v_eff, 1 - v_eff, 0;...
                       0, 0, 0, 1/2*(1 - 2*v_eff)];
                   
        D_2   = k_D_2  * eye(2);

end

% Note: 
% xi    = z = y
% eta   = r = x

for e = 1:nE
    
    %%% Topological information of the element
    i_el            = mesh.top(1:nn_in_E,e);
    elementDOF      = DOFfree(i_el,:);
    node_coord      = mesh.nodes(i_el,:);
    node_area       = 1/2 * det([node_coord(1,1) node_coord(1,2) 1;...
                                 node_coord(2,1) node_coord(2,2) 1;...
                                 node_coord(3,1) node_coord(3,2) 1]);
    
    if multimaterial 
        E_eff   = E(mesh.top(nn_in_E + 1,e));
        v_eff   = v(mesh.top(nn_in_E + 1,e));
        rho_eff = rho(mesh.top(nn_in_E + 1,e));
        eta_eff = eta(mesh.top(nn_in_E + 1,e));
        
        %%% Material stiffness matrix
        k_D_1 = E_eff/((1 + v_eff)*( 1 - 2* v_eff));
        k_D_2 = E_eff/(2*(1 + v_eff));

        D_1   = k_D_1 * [1 - v_eff, v_eff, v_eff, 0;...
                       v_eff, 1 - v_eff, v_eff  0;...
                       v_eff, v_eff, 1 - v_eff, 0;...
                       0, 0, 0, 1/2*(1 - 2*v_eff)];
        D_2 =  k_D_2  * eye(2);
    end
    
    %%% Element mass matrix
    me  = zeros(numel(elementDOF));
    ke  = zeros(numel(elementDOF));
    
    for it = 1:NPG      % Gauss quadrature
                
        %%% Shape functions and its local derivatives
        N_0         = N(xi_0(it), eta_0(it));
        
        if isa(dNdloc,'function_handle')
            dNdloc_0 = dNdloc(eta_0(it),xi_0(it));
        else
            dNdloc_0 = dNdloc;
        end
        
        
        %%% Local to global coordinates
        [Jaco,~]    = jacobian(node_coord,dNdloc_0);
        J           = det(Jaco);
        
        dNdglob     = Jaco\dNdloc_0;  % Global derivatives
        
        %%% Element deformation matrix (B)
        [B1, B2] = computeB_2D(dNdglob,N_0,xi_0(it),elementType,nodalDiam);   
              
        %%% Element matrices (global)
        me = me +  pi * node_area * xi_0(it) * w(it) * ( N_0.'*N_0) * J * rho_eff;
        ke = ke +  pi * node_area * xi_0(it) * w(it) * ((B1.' * D_1 * B1 + B2.' * D_2 * B2) * J) * (1 + 1i * eta_eff) ;  
        
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

end

MFree = sparse(Ig, Jg, Mg, nDOF_glob, nDOF_glob);
KFree = sparse(Ig, Jg, Kg, nDOF_glob, nDOF_glob);

%%% Boundary conditions
%   DOF real
DOF                 = zeros(nn, nDOFn);
freeNodes           = setdiff( (1:nn).', constrNodes);
nFreeNodes          = length(freeNodes);
nDOF                = nDOFn*nFreeNodes;
DOF(freeNodes, :)   = reshape( (1:nDOF).', nDOFn, nFreeNodes).';
                            
%   K and M real
freeDofsPrev        = sort(reshape(DOFfree(freeNodes, :).', [], 1));
K                   = KFree(freeDofsPrev, freeDofsPrev);
M                   = MFree(freeDofsPrev, freeDofsPrev);
end

function [B1, B2] = computeB_2D(dNdglob,N,r,elementType,nodalDiam)

           elementType;

           B11 = [dNdglob(1,1), 0, 0;...
                  N(1,1)/r, nodalDiam*N(1,1)/r, 0;...
                  0, 0, dNdglob(2,1);...
                  dNdglob(2,1), 0, dNdglob(1,1)]; 
           
           B12 = [dNdglob(1,2), 0, 0;...
                  N(2,2)/r, nodalDiam*N(2,2)/r, 0;...
                  0, 0, dNdglob(2,2);...
                  dNdglob(2,2), 0, dNdglob(1,2)]; 
                  
           B13 = [dNdglob(1,3), 0, 0;...
                  N(3,3)/r, nodalDiam*N(3,3)/r, 0;...
                  0, 0, dNdglob(2,3);...
                  dNdglob(2,3), 0, dNdglob(1,3)];
              
             
           B21 = [ -nodalDiam*N(1,1)/r,  (dNdglob(1,1)-N(1,1)/r),  0;...
                     0, dNdglob(2,1), -nodalDiam*N(1,1)/r];
           
           B22 = [ -nodalDiam*N(2,2)/r,  (dNdglob(1,2)-N(2,2)/r),  0;...
                 0, dNdglob(2,2), -nodalDiam*N(2,2)/r];
           
           B23 = [ -nodalDiam*N(3,3)/r,  (dNdglob(1,3)-N(3,3)/r),  0;...
                 0, dNdglob(2,3), -nodalDiam*N(3,3)/r];

           
                  
           B1 = [B11 B12 B13];
           B2 = [B21 B22 B23];

end






