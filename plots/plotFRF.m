subplot(2,1,1)
semilogy(wHz,abs(HaaW),'k','LineWidth',1.2)

haxis = gca; 
% Font size
FS    = 13; 
% Axis format
set(gca, 'FontName', 'TimesNewRoman');
hYLabel          = ylabel('$\mid  u(f) \mid$','Interpreter','latex');
hXLabel          = xlabel('Frecuencia (f) [Hz]','Interpreter','latex');
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
haxis.LineWidth  = 0.75;
haxis.FontName   = 'Palatine';
haxis.FontSize   = 0.65*FS;
hXLabel.FontSize = 0.85*FS;
hYLabel.FontSize = 0.85*FS;
% ylim([min(abs(FRF_Y)) max(abs(FRF_Y))])
% xlim([10 max(w)])
% pbaspect([(1+sqrt(5))/2 1 1]) % Aurea proportion


subplot(2,1,2)
plot(wHz,angle(HaaW),'k','LineWidth',1.2)

%legend({'Sin capa','Modelo FEM con capa','Modelo analítico con capa'},'Interpreter','latex')
haxis = gca; 
% Font size
FS    = 13; 
% Axis format
set(gca, 'FontName', 'TimesNewRoman');
hYLabel          = ylabel('$\alpha(f)$','Interpreter','latex');
hXLabel          = xlabel('Frecuencia (f) [Hz]','Interpreter','latex');
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
haxis.LineWidth  = 0.75;
haxis.FontName   = 'Palatine';
haxis.FontSize   = 0.65*FS;
hXLabel.FontSize = 0.85*FS;
hYLabel.FontSize = 0.85*FS;