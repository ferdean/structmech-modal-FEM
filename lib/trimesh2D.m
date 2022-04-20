function[mesh,dof,q]=trimesh2D(geom,acc,plotting)
% Mesh definition for 2D plate
% =========================================================================
% Created by:   Ferran de Andrés(1.2021) 
% =========================================================================
% INPUT: 
%   geom     = Name of the document with the vertices of the geometry
%              Str value
%   acc      = Accuracy of the meshing
%              Scalar value (0,inf)
%   plotting = Defines if the output includes a plot [-]
%              Binary value (0,1)
% OUTPUT:
%   mesh     = Struct with all the meshing information listed below: 
%      nodes = Point matrix. The first and second rows contain x- and 
%              y-coordinates of the points in the mesh.
%      top   = Triangle topology matrix. The first three rows contain 
%              indices to the corner points, given in counter clockwise 
%              order, and the fourth row contains the subdomain number.
%      e     = Mesh edges, returned as a 7-by-Ne matrix, where Ne is the 
%              number of boundary edges in the mesh. 
%              An edge is a pair of points in p containing a boundary 
%              between subdomains, or containing an outer boundary
%
%   dof      = Matrix with degrees of freedom corresponding to each node
%   q        = Quality of the elements in the mesh. The quality
%              criterion measures the likeness of the elment to the
%              reference. Max: q=1. Unacceptable: q<0.6 
% =========================================================================

%%% Input check (3 inputs)
narginchk(3,3);

%%% Ensure inputs are well defined
% geom must be an string array
if ~(isstring(geom))
    error('Name of the geometry document error');
end
% acc must be an scalar value
if ~(isscalar(acc))&&acc>0
     error('Accuracy definition error');
end
% Plotting decision should be binary
if (plotting~=1)&&(plotting~=0)
    error('Plot variable error');
end

if      size(geom,1) == 1

    %%% Geometry input
    % [File,Path] = uigetfile('.txt');
    geometry    = importdata(join([geom, ".txt"],""));
    geo         = [2                   
                   size(geometry,1)-1      % Number of segments
                   geometry(1:end-1,1)     % x-coordinates
                   geometry(1:end-1,2)];   % y-coordinates
    [g,~]       = decsg(geo);

elseif  size(geom,1)  == 2
    geometry    = importdata(join([geom(1), ".txt"],""));
    geo         = [2                   
                   size(geometry,1)-1      % Number of segments
                   geometry(1:end-1,1)     % x-coordinates
                   geometry(1:end-1,2)];   % y-coordinates
    geometry2    = importdata(join([geom(2), ".txt"],""));
    geo2         = [2                   
                   size(geometry2,1)-1      % Number of segments
                   geometry2(1:end-1,1)     % x-coordinates
                   geometry2(1:end-1,2)];   % y-coordinates

    ns=(char('in','out'))';
    sf='out + in ';

    geo3=[geo2,geo];

    g       = decsg(geo3,sf,ns);
    
elseif size(geom,1) == 3
    geometry    = importdata(join([geom(1), ".txt"],""));
    geo         = [2                   
                   size(geometry,1)-1      % Number of segments
                   geometry(1:end-1,1)     % x-coordinates
                   geometry(1:end-1,2)];   % y-coordinates
    geometry2    = importdata(join([geom(2), ".txt"],""));
    geo2         = [2                   
                   size(geometry2,1)-1      % Number of segments
                   geometry2(1:end-1,1)     % x-coordinates
                   geometry2(1:end-1,2)];   % y-coordinates
    geometry3    = importdata(join([geom(3), ".txt"],""));
    geo3         = [2                   
                   size(geometry3,1)-1      % Number of segments
                   geometry3(1:end-1,1)     % x-coordinates
                   geometry3(1:end-1,2)];   % y-coordinates           
    
    ns=(char('undamped','damped','restrictor'))';
    sf='undamped + damped + restrictor';
    
    geo4 = [geo3, geo2, geo];
    
    g = decsg(geo4,sf,ns);
    % pdegplot(g,'FaceLabels','on')
end
    
%%% Meshing
r = min(range(geometry(:,1)),range(geometry(:,2)));
[p,e,t] = initmesh(g,'Hgrad',1.2,'Hmax',r/acc,'Init','off');

mesh.nodes = p'; mesh.e = e; mesh.top = t; 

% IMPORTANT: initmesh creates only 'linear' elements

dof=( zeros( size(p) ) )';


for n=1:size(p,2)
   dof(n,:)=[2*n-1 2*n]; 
end


%%% Check mesh quality
x1 = p(1,t(1,:)); x2 = p(1,t(2,:)); x3 = p(1,t(3,:));
y1 = p(2,t(1,:)); y2 = p(2,t(2,:)); y3 = p(2,t(3,:));

areas   = (abs(x1.*y2+x2.*y3+x3.*y1-x1.*y3-x2.*y1-x3.*y2)/2)';

l(:,1)  = sqrt((x1-x3).^2+(y1-y3).^2);
l(:,2)  = sqrt((x1-x2).^2+(y1-y2).^2);
l(:,3)  = sqrt((x2-x3).^2+(y2-y3).^2);

q       = zeros(size(areas,1),1);

for i = 1:size(areas,1)
  q(i) = 4*sqrt(3)*areas(i)/(l(i,1)^2+l(i,2)^2+l(i,3)^2);
end

if min(q)<0.6
   warning('Unacceptable mesh quality') 
end

%%% Plotting
if plotting
    % Mesh quality
    pdeplot(p,e,t,'xydata',q);
    
    hold on
    
    % Mesh
    for i=1:size(t,2)
        plot([p(1,t(1,i)),p(1,t(2,i)),p(1,t(3,i)),p(1,t(1,i))],...
             [p(2,t(1,i)),p(2,t(2,i)),p(2,t(3,i)),p(2,t(1,i))],...
             'color',[0.5 0.5 0.5])
    end
    
    % Geometry
    if      size(geom,1) == 1
        plot(geometry(:,2),geometry(:,1),'k')
    else
        plot(geometry(:,1),geometry(:,2),'k')
        plot(geometry2(:,1),geometry2(:,2),'k')   
    end

% Plot configuration
set(gca,'fontname','times')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
daspect([1 1 1]);
xlabel('Width [m]')
xlim([0 1.5*max(x2)])
ylabel('Height [m]')
title({['Mesh plot. Number of nodes: ',num2str(size(p,2))],...
    '\fontsize{12}\color{gray}Colorbar: Modal displacement scaled to unit mass'},'fontsize',14)
hold off
end

