clear all, 
close all,

% Load data
load sample_ppg.mat

figure('Position', [10, 10, 1000, 800]);
ms=mean(ppg);

% Get PPG, PPG', PPG", and PPG'" signals
ppg=ppg+0.5*ms;

dt=1/Fs;
vpg = savitzky_golay(ppg, 1, 9)./dt;
vpg = (vpg/max((vpg)-min(vpg))*(max(ppg)-min(ppg)))/2+0.5*ms;

apg = savitzky_golay(vpg, 1, 9)./dt;
apg = (apg/max((apg)-min(apg))*(max(ppg)-min(ppg)))/2;

jpg = savitzky_golay(apg, 1, 9)./dt;
jpg = (jpg/max((jpg)-min(jpg))*(max(ppg)-min(ppg)))/2-0.5*ms;

% Plot signals
plot(ppg,'k','LineWidth',2),hold on;
plot(vpg,'k--','LineWidth',1.5),hold on;
plot(apg,'k:','LineWidth',2),hold on;
plot(jpg,'k-.','LineWidth',2),hold on;

ylim([-0.2 1.2])
xlim([0 1.3*Fs])
ylim([min(jpg) max(ppg)])
xticks(0:Fs/2:Fs)
xticklabels({'0','500','1000'})

% Set legend texts
text(1, ppg(1)+0.03,'PPG','FontSize', 15);
text(1, vpg(1)+0.03,'PPG''','FontSize', 15);
text(1, apg(1)+0.03,'PPG''''','FontSize', 15);
text(1, jpg(1)+0.03,'PPG''''''','FontSize', 15);

% Set labels
xlabel('Time [ms]','FontSize', 20)
set(gca,'FontSize',20);
ylabel('Amplitude [nu]','FontSize', 20)
set(gca,'ytick',[])

% Get the current date and time
currentDate = datestr(now, 'yyyy-mm-dd');

% Specify the filename
filename = ['Sample_Signals_', currentDate, '.svg'];

% Save the figure as SVG
saveas(gcf, filename, 'svg');