function [N, dNdloc] = SF(elementtype)
% Definition of shape functions for several element types used in FEM
% analysis
% =========================================================================
% Created by:   Ferran de Andrés(7.2021) 
% Version: 1.0. First implementation and definition of linear triangular
%               elements shape functions
%          1.1  Introduction of linear shape function derivatives
%          2.0  Introduction of quadratic elements
%          3.0  Implementation of 3D elements (linear)
%          3.1  Implementation of 3D elements (quadratic)
% =========================================================================
% INPUT: 
%   elementtype    = Type of element
%                    Char value
%
%                    'tri3'   - Triangular linear element (2D)
%                    'tri6'   - Triangular quadratic element (2D)
%                    'tet4'   - Linear tetraedron (3D)
%                    'tet10'  - Quadratic tetraedron (3D)
%                    'hex8'   - Linear hexaedron (3D)
%                    'hex20'  - Quadratic hexaedron (3D)
% OUTPUT:
%   N              = Shape Functions 
%  dNdloc          = Differentiation of shape functions in local coordinates
%                       
% =========================================================================

%%% Input check
narginchk(1,1);

if ~(isstring(elementtype))
       error('Element type is not well defined') 
end

%%% Shape functions and derivatives

switch elementtype

    case "test"
        N1= @(xi,eta)  1-eta-xi;
        N2= @(xi,eta)  xi;
        N3= @(xi,eta)  eta;
        N = @(xi,eta)  [N1(xi,eta) N2(xi,eta) N3(xi,eta);
                        N1(xi,eta) N2(xi,eta) N3(xi,eta)];



        %%% Derivatives
        dN1dxi = -1;
        dN2dxi =  1;
        dN3dxi =  0;
        dN1deta  = -1;
        dN2deta  =  0;
        dN3deta  =  1;

        dNdloc = [dN1dxi dN2dxi dN3dxi;...
                  dN1deta  dN2deta  dN3deta]; 

    case "tri3"
        
        N1= @(xi,eta)  1-eta-xi;
        N2= @(xi,eta)  eta;
        N3= @(xi,eta)  xi;
        N = @(xi,eta)  [N1(xi,eta) 0 N2(xi,eta) 0 N3(xi,eta) 0;
                        0 N1(xi,eta) 0 N2(xi,eta) 0 N3(xi,eta)];
                        % 0 0 N1(eta,xi) 0 0 N2(xi,eta) 0 0 N3(xi,eta)];


        %%% Derivatives
        dN1dxi = -1;
        dN2dxi =  0;
        dN3dxi =  1;
        dN1deta  = -1;
        dN2deta  =  1;
        dN3deta  =  0;

        dNdloc = [dN1dxi dN2dxi dN3dxi;...
                  dN1deta  dN2deta  dN3deta]; 
   
    case "tri6"
        
        %%% Shape functions
        N1 = @(xi,eta) (1-xi-eta)*(1-2*xi-2*eta);
        N2 = @(xi,eta) -xi*(1-2*xi);
        N3 = @(xi,eta) -eta*(1-2*eta);
        N4 = @(xi,eta) 4*xi*(1-xi-eta);
        N5 = @(xi,eta) 4*xi*eta;
        N6 = @(xi,eta) 4*eta*(1-xi-eta);

        %%% Matrix
        I  = eye(2);
        N  = @(xi,eta) [N1(xi,eta)*I,  N2(xi,eta)*I, ...
                        N3(xi,eta)*I,  N4(xi,eta)*I, ...
                        N5(xi,eta)*I,  N6(xi,eta)*I];

        %%% Derivative
        dN1dxi  = @(xi,eta) 4*eta + 4*xi - 3;
        dN1deta = @(xi,eta) 4*eta + 4*xi - 3;

        dN2dxi  = @(xi,eta) 4*xi - 1;
        dN2deta = @(xi,eta) 0;

        dN3dxi  = @(xi,eta) 0;
        dN3deta = @(xi,eta) 4*eta - 1;

        dN4dxi  = @(xi,eta) 4 - 8*xi - 4*eta;
        dN4deta = @(xi,eta) -4*xi;

        dN5dxi  = @(xi,eta) 4*eta;
        dN5deta = @(xi,eta) 4*xi;

        dN6dxi  = @(xi,eta) -4*eta;
        dN6deta = @(xi,eta) 4 - 4*xi - 8*eta;



        dNdloc = @(xi,eta) ...
            [dN1dxi(xi,eta),    dN2dxi(xi,eta),    dN3dxi(xi,eta), ...
             dN4dxi(xi,eta),    dN5dxi(xi,eta),    dN6dxi(xi,eta); ...
             dN1deta(xi,eta)    dN2deta(xi,eta),   dN3deta(xi,eta), ...
             dN4deta(xi,eta),   dN5deta(xi,eta),   dN6deta(xi,eta)];    

    case "hex8"
        
        %%% Shape functions
        N1 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 - eta) * (1 - tau);
        N2 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 - eta) * (1 - tau);
        N3 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 + eta) * (1 - tau);
        N4 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 + eta) * (1 - tau);
        N5 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 - eta) * (1 + tau);
        N6 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 - eta) * (1 + tau);
        N7 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 + eta) * (1 + tau);
        N8 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 + eta) * (1 + tau);
        
        %%% Matrix
        I  = eye(3);
        
        N = @(xi,eta,tau) [N1(xi,eta,tau)*I, N2(xi,eta,tau)*I,...
                           N3(xi,eta,tau)*I, N4(xi,eta,tau)*I,...
                           N5(xi,eta,tau)*I, N6(xi,eta,tau)*I,...
                           N7(xi,eta,tau)*I, N8(xi,eta,tau)*I];
        
        %%% Derivative
        % N1 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 - eta) * (1 - tau);
        dN1dxi  = @(xi,eta,tau) - 1/8 * (1 - eta) * (1 - tau);
        dN1deta = @(xi,eta,tau) - 1/8 * (1 - xi)  * (1 - tau);
        dN1dtau = @(xi,eta,tau) - 1/8 * (1 - xi)  * (1 - eta);
        
        % N2 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 - eta) * (1 - tau);
        dN2dxi  = @(xi,eta,tau)   1/8 * (1 - eta) * (1 - tau);
        dN2deta = @(xi,eta,tau) - 1/8 * (1 + xi)  * (1 - tau);
        dN2dtau = @(xi,eta,tau) - 1/8 * (1 + xi)  * (1 - eta);
        
        % N3 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 + eta) * (1 - tau);
        dN3dxi  = @(xi,eta,tau)   1/8 * (1 + eta) * (1 - tau);
        dN3deta = @(xi,eta,tau)   1/8 * (1 + xi)  * (1 - tau);
        dN3dtau = @(xi,eta,tau) - 1/8 * (1 + xi)  * (1 + eta);
        
        % N4 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 + eta) * (1 - tau);
        dN4dxi  = @(xi,eta,tau) - 1/8 * (1 + eta) * (1 - tau);
        dN4deta = @(xi,eta,tau)   1/8 * (1 - xi)  * (1 - tau);
        dN4dtau = @(xi,eta,tau) - 1/8 * (1 - xi)  * (1 + eta);
        
        % N5 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 - eta) * (1 + tau);
        dN5dxi  = @(xi,eta,tau) - 1/8 * (1 - eta) * (1 + tau);
        dN5deta = @(xi,eta,tau) - 1/8 * (1 - xi)  * (1 + tau);
        dN5dtau = @(xi,eta,tau)   1/8 * (1 - xi)  * (1 - eta);

        % N6 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 - eta) * (1 + tau);
        dN6dxi  = @(xi,eta,tau)   1/8 * (1 - eta) * (1 + tau);
        dN6deta = @(xi,eta,tau) - 1/8 * (1 + xi)  * (1 + tau);
        dN6dtau = @(xi,eta,tau)   1/8 * (1 + xi)  * (1 - eta);
        
        % N7 = @(xi,eta,tau) 1/8 * (1 + xi) * (1 + eta) * (1 + tau);
        dN7dxi  = @(xi,eta,tau)   1/8 * (1 + eta) * (1 + tau);
        dN7deta = @(xi,eta,tau)   1/8 * (1 + xi)  * (1 + tau);
        dN7dtau = @(xi,eta,tau)   1/8 * (1 + xi)  * (1 + eta);
        
        % N8 = @(xi,eta,tau) 1/8 * (1 - xi) * (1 + eta) * (1 + tau);
        dN8dxi  = @(xi,eta,tau) - 1/8 * (1 + eta) * (1 + tau);
        dN8deta = @(xi,eta,tau)   1/8 * (1 - xi)  * (1 + tau);
        dN8dtau = @(xi,eta,tau)   1/8 * (1 - xi)  * (1 + eta);
        
        
        dNdxi  = @(xi,eta,tau) [dN1dxi(xi,eta,tau), dN2dxi(xi,eta,tau),...
                                dN3dxi(xi,eta,tau), dN4dxi(xi,eta,tau),...
                                dN5dxi(xi,eta,tau), dN6dxi(xi,eta,tau),...
                                dN7dxi(xi,eta,tau), dN8dxi(xi,eta,tau)];
                            
                            
        dNdeta  = @(xi,eta,tau) [dN1deta(xi,eta,tau), dN2deta(xi,eta,tau),...
                                 dN3deta(xi,eta,tau), dN4deta(xi,eta,tau),...
                                 dN5deta(xi,eta,tau), dN6deta(xi,eta,tau),...
                                 dN7deta(xi,eta,tau), dN8deta(xi,eta,tau)];
                            
                            
        dNdtau  = @(xi,eta,tau) [dN1dtau(xi,eta,tau), dN2dtau(xi,eta,tau),...
                                 dN3dtau(xi,eta,tau), dN4dtau(xi,eta,tau),...
                                 dN5dtau(xi,eta,tau), dN6dtau(xi,eta,tau),...
                                 dN7dtau(xi,eta,tau), dN8dtau(xi,eta,tau)];
        
        dNdloc  = @(xi,eta,tau) [dNdxi(xi,eta,tau);...
                                 dNdeta(xi,eta,tau);...
                                 dNdtau(xi,eta,tau)];
    case "tet4"
        
        %%% Shape functions
        N1 = @(xi,eta,tau) eta;
        N2 = @(xi,eta,tau) tau;
        N3 = @(xi,eta,tau) 1 - xi - eta - tau;
        N4 = @(xi,eta,tau) xi;
        
        %%% Matrix
        I  = eye(3);
        
        N = @(xi,eta,tau) [N1(xi,eta,tau)*I, N2(xi,eta,tau)*I,...
                           N3(xi,eta,tau)*I, N4(xi,eta,tau)*I];
               
        %%% Derivative
        
        dN1dxi  = @(xi,eta,tau)   0;
        dN1deta = @(xi,eta,tau)   1;
        dN1dtau = @(xi,eta,tau)   0;
        
        dN2dxi  = @(xi,eta,tau)   0;
        dN2deta = @(xi,eta,tau)   0;
        dN2dtau = @(xi,eta,tau)   1;
        
        dN3dxi  = @(xi,eta,tau) - 1;
        dN3deta = @(xi,eta,tau) - 1;
        dN3dtau = @(xi,eta,tau) - 1;
        
        dN4dxi  = @(xi,eta,tau)   1;
        dN4deta = @(xi,eta,tau)   0;
        dN4dtau = @(xi,eta,tau)   0;
        
        
        dNdxi  = @(xi,eta,tau) [dN1dxi(xi,eta,tau), dN2dxi(xi,eta,tau),...
                                dN3dxi(xi,eta,tau), dN4dxi(xi,eta,tau)];
                            
                            
        dNdeta  = @(xi,eta,tau) [dN1deta(xi,eta,tau), dN2deta(xi,eta,tau),...
                                 dN3deta(xi,eta,tau), dN4deta(xi,eta,tau)];
                            
                            
        dNdtau  = @(xi,eta,tau) [dN1dtau(xi,eta,tau), dN2dtau(xi,eta,tau),...
                                 dN3dtau(xi,eta,tau), dN4dtau(xi,eta,tau)];
        
        dNdloc  = @(xi,eta,tau) [dNdxi(xi,eta,tau);...
                                 dNdeta(xi,eta,tau);...
                                 dNdtau(xi,eta,tau)];
                             
    case "tet10"
        
        %%% Shape functions
        N1  = @(xi,eta,tau) eta * (2*eta - 1);
        N2  = @(xi,eta,tau) tau * (2*tau - 1);
        N3  = @(xi,eta,tau) (1 - xi - eta - tau) * (1 - 2*xi - 2*eta - 2*tau);
        N4  = @(xi,eta,tau) xi  * (2*xi  - 1);
        N5  = @(xi,eta,tau) 4 * eta * tau;
        N6  = @(xi,eta,tau) 4 * tau * (1 - xi - eta - tau);
        N7  = @(xi,eta,tau) 4 * eta * (1 - xi - eta - tau);
        N8  = @(xi,eta,tau) 4 * xi  * eta;
        N9  = @(xi,eta,tau) 4 * xi  * tau; 
        N10 = @(xi,eta,tau) 4 * xi  * (1 - xi - eta - tau);
        
        %%% Matrix
        I  = eye(3);
        
        N = @(xi,eta,tau) [N1(xi,eta,tau)*I, N2(xi,eta,tau)*I,...
                           N3(xi,eta,tau)*I, N4(xi,eta,tau)*I,...
                           N5(xi,eta,tau)*I, N6(xi,eta,tau)*I,...
                           N7(xi,eta,tau)*I, N8(xi,eta,tau)*I,...
                           N9(xi,eta,tau)*I, N10(xi,eta,tau)*I];
                       
        %%% Derivative
        dN1dxi  = @(xi,eta,tau)   0;
        dN1deta = @(xi,eta,tau) - 1 + 4 * eta;
        dN1dtau = @(xi,eta,tau)   0;
        
        dN2dxi  = @(xi,eta,tau)   0;
        dN2deta = @(xi,eta,tau)   0;
        dN2dtau = @(xi,eta,tau) - 1 + 4 * tau;
        
        dN3dxi  = @(xi,eta,tau) - 3 + 4 * xi + 4 * eta + 4 * tau;
        dN3deta = @(xi,eta,tau) - 3 + 4 * xi + 4 * eta + 4 * tau;
        dN3dtau = @(xi,eta,tau) - 3 + 4 * xi + 4 * eta + 4 * tau;
        
        dN4dxi  = @(xi,eta,tau) - 1 + 4 * xi;
        dN4deta = @(xi,eta,tau)   0;
        dN4dtau = @(xi,eta,tau)   0;
        
        dN5dxi  = @(xi,eta,tau)   0;
        dN5deta = @(xi,eta,tau)   4 * tau;
        dN5dtau = @(xi,eta,tau)   4 * eta;
        
        dN6dxi  = @(xi,eta,tau) - 4 * tau;
        dN6deta = @(xi,eta,tau) - 4 * tau;
        dN6dtau = @(xi,eta,tau) - 4 * (- 1 + xi + eta + 2 * tau);
        
        dN7dxi  = @(xi,eta,tau) - 4 * eta;
        dN7deta = @(xi,eta,tau) - 4 * (- 1 + xi + 2 * eta + tau);
        dN7dtau = @(xi,eta,tau) - 4 * eta;
        
        dN8dxi  = @(xi,eta,tau)   4 * eta;
        dN8deta = @(xi,eta,tau)   4 * xi;
        dN8dtau = @(xi,eta,tau)   0;
        
        dN9dxi  = @(xi,eta,tau)   4 * tau;
        dN9deta = @(xi,eta,tau)   0;
        dN9dtau = @(xi,eta,tau)   4 * xi;
        
        dN10dxi  = @(xi,eta,tau) - 4 * (- 1 + 2 * xi + eta + tau);
        dN10deta = @(xi,eta,tau) - 4 * xi;
        dN10dtau = @(xi,eta,tau) - 4 * xi;
        
        dNdxi  = @(xi,eta,tau) [dN1dxi(xi,eta,tau), dN2dxi(xi,eta,tau),...
                                dN3dxi(xi,eta,tau), dN4dxi(xi,eta,tau),...
                                dN5dxi(xi,eta,tau), dN6dxi(xi,eta,tau),...
                                dN7dxi(xi,eta,tau), dN8dxi(xi,eta,tau),...
                                dN9dxi(xi,eta,tau), dN10dxi(xi,eta,tau)];
                            
                            
        dNdeta  = @(xi,eta,tau) [dN1deta(xi,eta,tau), dN2deta(xi,eta,tau),...
                                 dN3deta(xi,eta,tau), dN4deta(xi,eta,tau),...
                                 dN5deta(xi,eta,tau), dN6deta(xi,eta,tau),...
                                 dN7deta(xi,eta,tau), dN8deta(xi,eta,tau),...
                                 dN9deta(xi,eta,tau), dN10deta(xi,eta,tau)];
                            
                            
        dNdtau  = @(xi,eta,tau) [dN1dtau(xi,eta,tau), dN2dtau(xi,eta,tau),...
                                 dN3dtau(xi,eta,tau), dN4dtau(xi,eta,tau),...
                                 dN5dtau(xi,eta,tau), dN6dtau(xi,eta,tau),...
                                 dN7dtau(xi,eta,tau), dN8dtau(xi,eta,tau),...
                                 dN9dtau(xi,eta,tau), dN10dtau(xi,eta,tau)];
        
        dNdloc  = @(xi,eta,tau) [dNdxi(xi,eta,tau);...
                                 dNdeta(xi,eta,tau);...
                                 dNdtau(xi,eta,tau)];
        
        
end
        
        
end

