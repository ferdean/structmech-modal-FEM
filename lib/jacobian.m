function[Jacobian,InverseJacobian]=jacobian(node_coord,DSF)
% Jacobian matrix computation
% =========================================================================
% Created by:   Ferran de Andrés(2.2021) 
% =========================================================================
% INPUT: 
%   node_coord      = Coordinates of the nodes of the element
%                     3x2 matrix (if triangular elements)
%   DSF             = Natural derivatives of shape functions (wrt local
%                     coordinates)
% OUTPUT:
%   Jacobian        = Jacobian matrix 
%                     nxn matrix, being n the dimensions of the problem
%   InverseJacobian = Inverse of the Jacobian matrix
% =========================================================================
%%% Input check (2 inputs)
narginchk(2,2);
%%% Ensure inputs are well defined
% node_coord must be a 1x6 matrix with real elements
    if ~(size(node_coord,1)==1||size(node_coord,2)==6||isreal(node_coord))
       error('Nodes coordinates are not well defined') 
    end

%%% Check if symbolic and Jacobian calculation
if isa(DSF,'function_handle')
    syms eta xi
    Jacobian=DSF(eta,xi)*node_coord;
else
    Jacobian=DSF*node_coord;
end

%%% Inverse matrix
InverseJacobian=inv(Jacobian); % A more efficient routine will be 
                               % implemented in further versions

end