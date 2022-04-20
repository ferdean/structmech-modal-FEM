function [q] = elementquality(mesh)
% Element quality plot
% =========================================================================
% Created by:   Ferran de Andrés(3.2021) 
% =========================================================================
% INPUT: 
%   mesh     = Struct with all the meshing information listed below: 
%      nodes = Point matrix. The first and second rows contain x- and 
%              y-coordinates of the points in the mesh.
%      top   = Triangle topology matrix. The first three rows contain 
%              indices to the corner points, given in counter clockwise 
%              order, and the fourth row contains the subdomain number.
% OUTPUT:
%   quality_histogram_plot 
%   q        = Quality of the elements in the mesh. The quality
%              criterion measures the likeness of the elment to the
%              reference. Max: q=1. Unacceptable: q<0.6
%              1 x nE matrix 
% =========================================================================

%%% Input check (1 inputs)
narginchk(1,1);

if ~(isstruct(mesh))
    error('Mesh definition error');
end

%%% Compute mesh quality
p = mesh.nodes;
t = mesh.top;

x1 = p(t(1,:),1); x2 = p(t(2,:),1); x3 = p(t(3,:),1);
y1 = p(t(1,:),2); y2 = p(t(2,:),2); y3 = p(t(3,:),2);
y1 = p(t(1,:),2); y2 = p(t(2,:),2); y3 = p(t(3,:),2);

areas   = (abs(x1.*y2+x2.*y3+x3.*y1-x1.*y3-x2.*y1-x3.*y2)/2);

l(:,1)  = sqrt((x1-x3).^2+(y1-y3).^2);
l(:,2)  = sqrt((x1-x2).^2+(y1-y2).^2);
l(:,3)  = sqrt((x2-x3).^2+(y2-y3).^2);

q = 4*sqrt(3).*areas .* ((l(:,1).^2+l(:,2).^2+l(:,3).^2)).^(-1);

if min(q)<0.6
   warning('Unacceptable mesh quality') 
end

%%% Solution representation
zero = num2cell(zeros(1,9));
[q_09,q_08,q_07,q_06,q_05,q_04,q_03,q_02,q_01] = deal(zero{:});

for i   = 1:length(q)
    if q(i) > 0.9
        q_09    = q_09+1;
    elseif q(i) < 0.9 &&  q(i) > 0.8 
        q_08    = q_08+1;
    elseif q(i) < 0.8 &&  q(i) > 0.7
        q_07    = q_07+1;
    elseif q(i) < 0.7 &&  q(i) > 0.6 
        q_06    = q_06+1;
    elseif q(i) < 0.6 &&  q(i) > 0.5 
        q_05    = q_05+1;
    elseif q(i) < 0.5 &&  q(i) > 0.4
        q_04    = q_04+1;    
    elseif q(i) < 0.4 &&  q(i) > 0.3
        q_03    = q_03+1;   
    elseif q(i) < 0.3 &&  q(i) > 0.2
        q_02    = q_02+1;
    else
        q_01    = q_01+1;
    end
end
q_qual_y  = [q_01,q_02,q_03,q_04,q_05,q_06,q_07,q_08,q_09];
q_qual_x = 0.1:0.1:1;



h   = histogram('BinEdges',q_qual_x,'BinCounts',q_qual_y);

h.FaceColor = [0, 0, 0];
% set(gca,'fontname',)

haxis.YAxis.MinorTick = 'on';
axis tight
grid minor
% minHeight = min(PdB(:)) - 0.1*(max(PdB(:)) - min(PdB(:)));
% Font size
FS    = 14; 
% Axis format
set(gca, 'FontName', 'Palatino');
hYLabel          = ylabel('Number of elements [-]','Interpreter','latex');
hXLabel          = xlabel('Mesh quality [-]','Interpreter','latex');
haxis.Box        = 'on';
haxis.TickDir    = 'out';
haxis.TickLength = [.02 .02];
haxis.XMinorTick = 'off';
haxis.YMinorTick = 'on';
haxis.YGrid      = 'on';
haxis.XMinorGrid = 'on';
haxis.YMinorGrid = 'on';
haxis.XGrid      = 'on';
haxis.XColor     = [.0 .0 .0];
haxis.YColor     = [.0 .0 .0];
haxis.LineWidth  = 0.75;
haxis.FontName   = 'Palatino';
haxis.FontSize   = 0.65*FS;
hXLabel.FontSize = 0.85*FS;
hYLabel.FontSize = 0.85*FS;
% ylim([180 250])
% xlim([10 7000])
pbaspect([(1+sqrt(5))/2 1 1]) % Aurea proportion