function [Dedrift] = Dedrift_v3b(Result,folder,dname)
%DEDRIFT_V3B Summary of this function goes here
%   Detailed explanation goes here
% This script will call the dedrifting

% File = 'D:\_Gamma35_Results\2018_01_30_Gamma35\Result_Gamma35_Phi2_dt5\_complete_tracked';
% File = 'D:\_DatenJohn\Result_x159_1_dt3_v3\_complete_tracked';
% File = 'D:\_Gamma35_Results\2017_12_07_Gamma35\Result_Gamma35_Phi1_dt1\_complete_tracked';
% 
% Data = dlmread([File '.txt'],'\t',1,0);

Data = Result;
Data(:,5) = Data(:,5)+1;
%%
timecol= 5; IDcol = 6;
winsize = 30;
[ driftx, drifty, driftxsmooth, driftysmooth] = FindDrift_v2( Data, timecol,IDcol, winsize);
[ Data_driftkorr ] = CorrectDrift_v2( Data, timecol, driftxsmooth, driftysmooth);
[ driftxneu, driftyneu,driftxneusmooth,driftyneusmooth] = FindDrift_v2( Data_driftkorr, timecol,IDcol,winsize);
%%
t = 1:length(driftx);
fig = figure('visible','off','Color',[1 1 1]);
subplot(2,1,1)
plot(t, driftx, 'b-', 'LineWidth', 1)
hold on
plot(t, drifty, 'r-', 'LineWidth', 1)
hold on 
plot(t, driftxsmooth, 'b-', 'LineWidth', 2)
hold on
plot(t, driftysmooth, 'r-', 'LineWidth', 2)
ylabel('Drift [µm]', 'Fontsize',12, 'FontWeight','bold','FontName','Helvetica' )
xlabel('timestep [{\tau} = 1]', 'Fontsize',12, 'FontWeight','bold', 'FontName','Helvetica' )
set(gca, 'Fontsize',12)
set(gca, 'FontWeight','bold')
set(gca,'LineWidth',1)
legend('X direction','Y direction')

subplot(2,1,2)
plot(t, driftxneu, 'b-', 'LineWidth', 1)
hold on
plot(t, driftyneu, 'r-', 'LineWidth', 1)
hold on
plot(t, driftxneusmooth, 'b-', 'LineWidth', 2)
hold on
plot(t, driftyneusmooth, 'r-', 'LineWidth', 2)
ylabel('Drift [µm]', 'Fontsize',12, 'FontWeight','bold','FontName','Helvetica' )
xlabel('timestep [{\tau} = 1]', 'Fontsize',12, 'FontWeight','bold', 'FontName','Helvetica' )
set(gca, 'Fontsize',12)
set(gca, 'FontWeight','bold')
set(gca,'LineWidth',1)
legend('X direction','Y direction')
output = strcat (folder,dname);
export_fig([output '_Dedrift'],fig,'-pdf','-tiff','-r300','-painters')
%% write down dedrifted data
Dedrift = Data_driftkorr;
dlmwrite([folder dname '_dedrift_v3.txt'], 'X Y Ux Uy T ID', 'delimiter', ' ');
dlmwrite([folder dname '_dedrift_v3.txt'], Data_driftkorr, '-append', 'delimiter', '\t','precision', 9);

a = [t' driftx drifty ];
dlmwrite([folder dname '_drift_raw_v3.txt'],a,'\t');
% 
% set(fig, 'Position', get(0, 'Screensize'));
% savefig(fig,[folder outfile(1:(end-5)) '_drift_v3.fig'])
