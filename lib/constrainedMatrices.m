function [K, M, MGDL, N] = constrainedMatrices(...
    KFree,MFree, MGDLFree, constrNodes)
% Includes constrained dofs in matrices
% =========================================================================
% Child functions: 
% Created by:      Víctor (11.2020) ; Modified by:      Ferran (2.2021)
% =========================================================================
% Versions:        1.0: (Víctor) Function implementation. 
%                  1.1: (Ferran) Generalisation on the scope of the 
%                                functions, comments and minor changes.
% =========================================================================
% INPUT:
%   KFree       = Stiffnes matrix
%                 Size NFree x NFree, with NFree number of dof's
%   MFree       = Mass matrix
%                 Size NFree x NFree
%   MGDLFree    = Matrix of dof's considering no constained nodes
%                 Size nn x Nn, with Nn nº dof's per node
%   constrNodes = Id of constrained nodes
%                 Size nConstr x 1
% OUTPUT:
%   K           = Stiffnes matrix
%                 Size N x N, with N number of dof's
%   M           = Mass matrix
%                 Size N x N
%   MGDL        = Matrix of dof's
%                 Size nn x Nn, with Nn nº dof's per node
%                 0 indicates with no dof

% =========================================================================

%%% Prev definition about dof's
[nn, NperNode]   = size(MGDLFree);
NinFree          = nn * NperNode;

%%% DOF's constrained:
% Assuming that all dof's of constrNodes are constrained
dofConstrained   = MGDLFree(constrNodes, :);
dofConstrained   = sort(dofConstrained(:));
% Not repeated dofConstrained must be
if any(unique(dofConstrained) ~= dofConstrained)
    error('Repeated dofs constrained');
end

% Free dofs
freeDofsPrev        = setdiff((1:NinFree).', dofConstrained);
N                   = length(freeDofsPrev);

%%% MGDL real
MGDL                = zeros(NinFree, 1);
% Fill non zero
MGDL(freeDofsPrev)  = (1:N).';
% Reshape
MGDL                = reshape(MGDL.', NperNode, nn).';

%%% K and M 
K                   = KFree(freeDofsPrev, freeDofsPrev);  
% Kc                  = KcFree(freeDofsPrev, freeDofsPrev); 
M                   = MFree(freeDofsPrev, freeDofsPrev);
                        
end


