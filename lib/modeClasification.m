function [] = modeClasification(nodes, fixednodes, phi, wnHz)
% Classification of he modes in terms of their components (axial, radial
% and torsional). This is specific for the study of wheels/discs.
% =========================================================================
% Created by:      Ferran de Andrés (2.2021)
% =========================================================================

nodesConstrained   = nodes;
nodesConstrained(fixednodes,:) = [];

tr = size(phi,2);

angleGlobLoc   = atan(nodesConstrained(:,3)./nodesConstrained(:,2));
phi_modal_loc  = zeros(size(nodesConstrained,1),3);

axial_perc      = zeros(tr,1);
radial_perc     = zeros(tr,1); 
torsional_perc  = zeros(tr,1); 

for i_mode = 1:tr
    
%   phi_nodal_glob = reshape(phi(:,mode), [size(phi,1)/3 3]);
    phi_modal_glob = [phi(1:3:end,i_mode) phi(2:3:end,i_mode) phi(3:3:end,i_mode)];
    
  for iR = 1:size(nodesConstrained,1)
    
    R_glob2loc = [1   0                      0;
                  0   cos(angleGlobLoc(iR))  sin(angleGlobLoc(iR));
                  0  -sin(angleGlobLoc(iR))  cos(angleGlobLoc(iR))]; 
    
    phi_modal_loc(iR,:) = ( R_glob2loc * (phi_modal_glob(iR,:)).').';
    
  end
  
  axialPart     = abs(mean(phi_modal_loc(:,1)));
  radialPart    = abs(mean(phi_modal_loc(:,2)));
  torsionalPart = abs(mean(phi_modal_loc(:,3)));
  
  tot     = axialPart + radialPart + torsionalPart;

  axial_perc(i_mode)        = axialPart/tot*100;
  radial_perc(i_mode)       = radialPart/tot*100;
  torsional_perc(i_mode)    = torsionalPart/tot*100;   
end

figure(1)
yyaxis right
modeplot    = scatter(1:tr,real(wnHz)/1000,'k','x');
modeplot.LineWidth = 1.2;
ylabel('Natural frequency [kHz]','Interpreter','latex')

haxis = gca; 
haxis.Box        = 'on';
haxis.TickDir    = 'out';
haxis.TickLength = [.02 .02];
haxis.XMinorTick = 'on';
haxis.YMinorTick = 'on';
haxis.YGrid      = 'on';
haxis.XMinorGrid = 'on';
haxis.YMinorGrid = 'on';
haxis.XGrid      = 'on';
haxis.XColor     = [.0 .0 .0];
haxis.YColor     = [.0 .0 .0];


yyaxis left
barplot     = bar([axial_perc,radial_perc,torsional_perc],'stacked');
ylim([0 100])
ytickformat('percentage')
barplot(1).FaceColor = [73 76 80]/255;
barplot(2).FaceColor = [205 205 205]/255;
barplot(3).FaceColor = [255 255 255]/255;
ylabel('Modal participation [-]','Interpreter','latex')

xlim([0 tr])
xlabel('Mode [-]','Interpreter','latex')
pbaspect([(1+sqrt(5))/2 1 1]) % Aurea proportion

haxis = gca; 
haxis.Box        = 'on';
haxis.TickDir    = 'out';
haxis.TickLength = [.02 .02];
haxis.XMinorTick = 'on';
haxis.YMinorTick = 'on';
haxis.YGrid      = 'on';
haxis.XMinorGrid = 'on';
haxis.YMinorGrid = 'on';
haxis.XGrid      = 'on';
haxis.XColor     = [.0 .0 .0];
haxis.YColor     = [.0 .0 .0];