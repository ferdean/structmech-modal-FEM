function [wnHz, phi, eta_r] = modal(K, M, tr)
% Modal solution 
% =========================================================================
% Child functions:
% Created by:      Víctor (11.2020) ; Modified by: Ferran (2.2021)
% =========================================================================
% Versions:        1.0: (Víctor) Modal solution and first implementation of
%                       the function. Main objetives of the functions are
%                       covered. 
%                  1.1: (Ferran) Numerical near-zero imaginary results
%                       chopping.
%                  1.2: (Ferran) Adaptation for structural damping studies.
%                 
% =========================================================================
% INPUT:
%   K         = Stiffness matrix [Pa]
%               Size N x N, with N nº dof's
%   M         = Mass matrix [kg]
%               Size N x N
%   tr        = Nº vibration modes extracted [#]
%               Scalar value
% OUTPUT:
%   phi       = Modal amplitudes (vibration modes) scaled to unit mass
%               Size N x tr
%               Each column are modal amplitudes for one mode
%               Each row contains model amplitude of dof with global id (in
%               MGDL) equal to that row
%   wnHz      = Natural frequencies [Hz]
%               Size tr x 1
%               Each row corresponds to mode in same column in phi
% =========================================================================

%%% Input check (3 inputs)
narginchk(3,3);

%%% Ensure inputs are well defined
% K and M square matrices with same order 
if ~ismatrix(K) || ~ismatrix(M) || size(K,1)~=size(K,2) || ...
        size(M,1)~=size(M,2) || size(K,1)~=size(M,1)
    error('Stiffness/mass matrix bad defined');
end
% tr must be scalar and real and natural and non-zero and 
% lesser or equal than Nº dofs
if ~(isscalar(tr) && imag(tr)==0 && mod(tr,1)==0 && tr>0) || tr>size(K, 1)
    error('Truncation tr bad defined');
end

%%% Previous calc
% Nº dofs
N      = size(K, 1);

%%% Eigen calc
[V, D] = eigs(sparse(K), sparse(M), tr, 'sm');

%%% Vibration modes scaled to unit mass
aux    = diag( V.'*(M*V) ).';
phi    = V./repmat(sqrt(aux), [N 1]);

%%% Natural frequencies
wn      = sqrt(real(diag(D)));    % [rad/s]
eta_r   = imag(diag(D))./wn.^2;   % [-]
wnHz    = wn/2/pi;                % [Hz]

%%% Natural frequencies negatives/imaginary are not possible
%%% Numerically it can be negative or imaginary as long as very small value
wnHz(wnHz<0 & wnHz>-1e-6 & wnHz<1e-6) = 0;
wnHz(imag(wnHz)~=0 & abs(imag(wnHz))<1e-3)=...
    real(wnHz(imag(wnHz)~=0 & abs(imag(wnHz))<1e-3)); 


%%% Sort natural frequencies
[wnHz, ord] = sort(wnHz, 'ascend');
phi         = phi(:, ord);

%%% Must be accomplished:
% phi.' * M * phi -> unit matrix
% phi.' * K * phi -> wn^2 (D)


end

