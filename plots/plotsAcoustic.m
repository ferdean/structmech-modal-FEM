stairs(cf(2:(end-2)),pAdBradial(2:(end-2),1),'Color','k','LineWidth',1.2)
hold on
% stairs(cfA(2:(end-2)),pAdBvisc_anal(2:(end-2),1),'Color',[180, 182, 177]/255,'LineWidth',1.2,'LineStyle','--')
% hold on
stairs(cfA(2:(end-2)),pAdBradialvisc(2:(end-2),1),'Color',[180, 182, 177]/255,'LineWidth',1.2)

% 'Color',[180, 182, 177]/255

cfLabels = [27.34, 54.6875, 109.375, 218.75, 437.5, 875,1750,3500,7000];
haxis = gca;
haxis.XScale = 'log';
haxis.XTick = cfLabels;
haxis.YAxis.MinorTick = 'on';
axis tight
grid minor
% minHeight = min(PdB(:)) - 0.1*(max(PdB(:)) - min(PdB(:)));
% Font size
FS    = 13; 
% Axis format
set(gca, 'FontName', 'TimesNewRoman');
hYLabel          = ylabel('Potencia acustica emitida [dB re 1pW]','Interpreter','latex');
hXLabel          = xlabel('Frecuencia (f) [Hz]','Interpreter','latex');
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
haxis.FontName   = 'Palatine';
haxis.FontSize   = 0.65*FS;
hXLabel.FontSize = 0.85*FS;
hYLabel.FontSize = 0.85*FS;
ylim([180 250])
xlim([100 5000])
pbaspect([(1+sqrt(5))/2 1 1]) % Aurea proportion