%%% Material data (CASE 1 - flat plate)

Lx                  = 460e-3;       % [m]
Ly                  = 100e-3;       % [m]
thick               = 3e-3;         % [m]

E_core              = 2.1e11;       % [N/m2]
v_core              = 0.285;        % [-] Poisson's ratio
rho_core            = 7850;         % [kg/m3] Average steel density
eta_core            = 5e-4;         % [-] Damping loss coefficient

% % EVA (Ethylene-vinyl acetate)
% E_viscoelastic      = 0.025e9;      %[N/m2]
% rho_viscoelastic    = 950;          %[kg/m3]
% v_viscoelastic      = 0.48;         %[-]
% eta_viscoelastic    = 0.83;         % [-] Damping loss coefficient

% Parametric study material
E_viscoelastic      = 0.002e9;        %[N/m2]
rho_viscoelastic    = 950;          %[kg/m3]
v_viscoelastic      = 0.4;          %[-]
eta_viscoelastic    = 0.8;          %[-] Damping loss coefficient


% % IIR (Butyl rubber)
% E_viscoelastic      = 0.0011e9;     %[N/m2]
% rho_viscoelastic    = 930;          %[kg/m3]
% v_viscoelastic      = 0.5;          %[-]
% eta_viscoelastic    = 0.9;          % [-] Damping loss coefficient


E_alum              = 72e9;         %[N/m2]
rho_alum            = 2700;         %[kg/m3]
v_alum              = 0.335;        %[-]
eta_alum            = 0.0011;       % [-] Damping loss coefficient




E       = [E_core,      E_viscoelastic,     E_alum];
v       = [v_core,      v_viscoelastic,     v_alum];
rho     = [rho_core,    rho_viscoelastic,   rho_alum];
eta     = [eta_core,    eta_viscoelastic,   eta_alum];



%%% Steel core      = 1. High carbon steel (0.7 - 1.7 % carbon)
%   Viscous layer   = 2. Two options: 
%                           - EVA. Ethylene Vinyl Acetate
%                           - IIR: Butyl rubber
%   Al restrictor   = 3. 
