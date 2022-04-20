function [xi, eta, w] = gaussTable_Tri_Victor(n)
% Generates gaussian quadrature points and weights for triangles in
% normalized coordinates
% =================================================================
% Child functions: 
% Created by:      Víctor (10.2020) ; Modified by: Ferran (1.2021)
% =================================================================
% INPUT:
%   n         = Nº of integration points
%               Must be -3, 1, 3, 4
% OUTPUT:
%   xi        = Normalized x coordinates
%               Size n x 1
%   eta       = Normalized y coordinates
%               Size n x 1
%   w         = Integration weights
%               Size n x 1
% =================================================================

if n~=-3 && n~=1 && n~=3 && n~=4
    error('Number of gauss points must be -3, 1, 3 or 4.')
end

%%% Table
if n==1
    xi  = 1/3;
    eta = 1/3;
    w   = 1/2;

elseif n==3
    xi  = [1/6; 2/3; 1/6];
    eta = [1/6; 1/6; 2/3];
    w   = [1/6; 1/6; 1/6];

elseif n==-3
    xi  = [1/2; 1/2; 0];
    eta = [0;   1/2; 1/2];
    w   = [1/6; 1/6; 1/6];

elseif n==4
    xi  = [1/3;    1/5;   3/5;   1/5];
    eta = [1/3;    1/5;   1/5;   3/5];
    w   = [-27/96; 25/96; 25/96; 25/96];
    
    
end
    


end