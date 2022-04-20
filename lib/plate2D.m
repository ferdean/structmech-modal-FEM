function []=plate2D(Lx,Ly,name,plot)
% Geometry definition of a 2D plate
% =========================================================================
% Created by:   Ferran de Andrés(1.2021) 
% =========================================================================
% INPUT: 
%   Lx      = Width of the plate [mm]
%             Scalar value
%   Ly      = Height of the plate [mm]
%             Scalar value
%   name    = Name of the auxiliar .txt document [-]
%             Str value
%   plot    = Defines if the output includes a plot [-]
%             Binary value (0,1)
% OUTPUT:
%   name.txt archive with the coordinates of the vertices of the 2D plate
% =========================================================================

%%% Input check (4 inputs)
narginchk(4,4);

%%% Ensure inputs are well defined
% Lengths must be scalar and real and greater than 0
if ~(isscalar(Lx) && imag(Lx)==0 && Lx>0)
    error('Length Lx error');
end
if ~(isscalar(Ly) && imag(Ly)==0 && Ly>0)
    error('Length Ly error');   
end
% name must be a character array
if ~(isstring(name))
    error('Name of the geometry document error');
end
% Plot decision should be binary
if (plot~=1)&&(plot~=0)
    error('Plot variable error');
end

%%% Vertices definition
vertices=[0  0; Lx  0; Lx Ly; 0 Ly; 0 0];

%%% Geometry archive 
geom=fopen(join([name,".txt"],""),'w');
dlmwrite(join([name,".txt"],""),vertices);
fclose(geom);

%%% Geometry plot
if plot
    x=vertices(:,1);
    y=vertices(:,2);
    patch(x,y,'k')
    set(gca,'fontname','times')
    daspect([1 1 1]);
    xlim([min(x)*0.5 max(x)*1.5])
    xlabel('Height [mm]')
    ylim([min(y)*0.5 max(y)*1.5])
    ylabel('Width [mm]')
    title('Geometry plot','fontsize',14)    
end

end


