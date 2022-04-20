function [xi, eta, w] = gaussTable_Tri(n)
% Generates gaussian quadrature points and weights for triangles in
% normalized axisymmetric coordinates
%
% Reference: Petyt, M. Finite element vibration analysis (Second edition,
% page 157, Table 5.1)
%
% =================================================================
% Child functions: 
% Created by:      Víctor (10.2020) ; Modified by: Ferran (1.2021)
% =================================================================
% INPUT:
%   n         = Formula number 
%               Must be [1, 3]
% OUTPUT:
%   xi        = Normalized x coordinates
%               Size n x 1
%   eta       = Normalized y coordinates
%               Size n x 1
%   zeta      = Normalized z coordinates
%               Size n x 1
%   w         = Integration weights
%               Size n x 1
% =================================================================

if n~=1 && n~=2 && n~=3
    error('Number of gauss points must be 1, 2 or 3. If more integration points are needed, please see the references at the documentation of the function')
end

%%% Table
if n == 1
    xi      = 1/3;
    eta     = 1/3;
%     zeta    = 1/3;
    w       = 1;

elseif n == 2
    xi      = [1/6; 2/3; 1/6];
    eta     = [1/6; 1/6; 2/3];
%     zeta    = [1/6; 1/6; 2/3];
    w       = [1/6; 1/6; 1/6];

elseif n == 3
    xi      = [1/3; 3/5; 1/5; 1/5];
    eta     = [1/3; 1/5; 3/5; 1/5];
%     zeta    = [1/3; 1/5; 1/5; 3/5];
    w       = [-27/48; 25/48; 25/48; 25/48];    
    
end
    


end