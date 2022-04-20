function [I,nit] = simpson2(func,a,b,c,d,maxsteps,tol)
% Numerical integration: Simpson's double integral in a region defined by
% functions
% =========================================================================
% Created by:   Ferran de Andrés(2.2021)
% =========================================================================
%  INPUT: 
%   func     = Function to be integrated
%              Handle function
%   [a,b]    = Integration limit of variable 1
%              Scalar values
%   [c,d]    = Integration limit of variable 2
%              Handle functions
%   maxsteps = Maximum iterations [-]
%              Scalar value
%   tol      = Numerical approach tolerance [-]
%              Scalar value
%
% OUTPUT:
%   I        = Defined integral result [-]
%              Scalar value
%   nit      = Number of iterations to meet the tolerance
%              Scalar value
% =========================================================================

steps=0;
In1=tol+1;
In=0;
while (abs(max(max((In1-In))))>tol)&&steps<maxsteps
    steps=steps+1;
    h=(b-a)/(2*steps);
    In=In1;

    p=zeros(1,(2*steps+1));
    for i=1:length(p)
        if i==1||i==length(p)
            p(i)=1;
        elseif mod(i,2)==0
            p(i)=4;
        else
            p(i)=2;
        end
    end
    clear i

    x=a:h:b;
    cx=feval(c,x);
    dx=feval(d,x);
    kx=(dx-cx)/(2*steps);
    y=ones(2*steps+1,1)*cx+(1:2*steps+1)'*kx;

    I1=0;I2=0;

    for i=1:(2*steps+1)
        for j=1:(2*steps+1)
            I1=I1+p(j)*feval(func,x(i),y(i,j));
        end
        I2=I2+p(i)*kx(i)/3*I1;
        I1=0;
    end

    In1=I2*h/3;   
end
I=In1;
nit=steps;

% if nit==maxsteps
%     warning('Maximum number of steps have been used')
% end
end

