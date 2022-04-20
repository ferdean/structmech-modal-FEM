function [xi, eta, tau, w] = gaussTable(n,elementType)
% Generates gaussian quadrature points and weights for hexahedrons in
% normalized coordinates
% Reference: Petyt, M. Finite element vibration analysis (Second edition,
% page 157, Table 5.1)
% =================================================================
% Child functions: 
% Created by:      Ferran de Andrés (7.2021)
% =================================================================
% INPUT:
%   n         = Nº of integration points in each direction
%               Must be 2, 3 
% OUTPUT:
%   xi        = Normalized x coordinates
%               Size n^3 x 1
%   eta       = Normalized y coordinates
%               Size n^3 x 1
%   tau       = Normalized z coordinates
%               Size n^3 x 1
%   w         = Integration weights
%               Size n^3 x 1
% =================================================================

if n~=2 && n~=3
    error('Number of gauss points in each direction must be 2 or 3.')
end

% if elementType == 'tet4 '
%     elementType = 'tet4';
% elseif elementType == 'hex8 '
%     elementType = 'hex8';
% end

switch elementType
    case "hex8"
        %%% Table
        if n==2
            k = 1/sqrt(3);

            xi  = [-k; -k; -k; -k;  k;  k;  k;  k];
            eta = [-k; -k;  k;  k; -k; -k;  k;  k];
            tau = [-k;  k; -k;  k; -k;  k; -k;  k];

            w   = ones(n^3,1);

        elseif n==3
            %%% 
            fprintf('\n To be defined')
        end
    case "tet4" 
        if n==1
            xi  = 1/4;
            eta = 1/4;
            tau = 1/4; 
            
            w   = 1/6;
            
        elseif n == 2
            k_1 = 1/20 * (5 - sqrt(5));
            k_2 = 1/20 * (5 + 3*sqrt(5));
            
            xi  = [k_1; k_1; k_1; k_2];
            eta = [k_1; k_1; k_2; k_1];
            tau = [k_1; k_2; k_1; k_1];
            
            w   = [1/24; 1/24; 1/24; 1/24];
            
        elseif n==3
            %%% 
            fprintf('\n To be defined')
        end
    case "tet10"
        if n==1
            xi  = 1/4;
            eta = 1/4;
            tau = 1/4; 
            
            w   = 1/6;
            
        elseif n==2
            k_1 = 1/20 * (5 - sqrt(5));
            k_2 = 1/20 * (5 + 3*sqrt(5));
            
            xi  = [k_1; k_1; k_1; k_2];
            eta = [k_1; k_1; k_2; k_1];
            tau = [k_1; k_2; k_1; k_1];
            
            w   = [1/24; 1/24; 1/24; 1/24];
            
        elseif n==3
            k_1 = 0.25;
            k_2 = 1/6;
            k_3 = 0.5; 
            
            xi  = [k_1; k_2; k_2; k_2; k_3];
            eta = [k_1; k_2; k_2; k_3; k_2];
            tau = [k_1; k_2; k_3; k_2; k_2];
            
            w   = [-2/15; 3/40; 3/40; 3/40; 3/40];
        end
end