function [Result] = StartTracking(Alltracks,folder,dname)
%STARTTRACKING Summary of this function goes here
%   Detailed explanation goes here

% Directory = 'D:\_DatenJohn\Result_x1_3_dt3_v4\';
% Directory = 'D:\_Gamma35Directory = 'D:\_Gamma35_Results\2018_01_30_Gamma35\Result_Gamma35_Phi2_dt5\';
% Results\2018_03_14_Gamma35\Result_Gamma35_Phi3_Overview2\';
% File = '_complete';

% Data = dlmread('C:\Users\Du\Desktop\Promotion\LM-10\Messungen\2018_10_30\Singleparticle_1_500_NaCl\M1 2018-10-30 11-33-21_AllTracks.csv',',',1,0);
Data = Alltracks;
N=length(Data(:,1));
% Old Script by Jörg
% Data = [Data(:,5) Data(:,6) zeros(N,1) zeros(N,3) Data(:,4)];
% NEW Script by Christoph:
Data = [Data(:,5) Data(:,6) zeros(N,2)  Data(:,4)];

% Data = dlmread([Directory File '.txt'],'\t',1,0);
Data = sortrows(Data,5);
%% for dt1 and dt2 rmax = 2, thetamax = 30°, for dt5 rmax = 3, thetamax = 60
maxdisp = 10; %maximum displacement (not squared)
maxrot = 90; %maximum rotation in degree (90° means obviously inactive)
    
tic
[Result, links_active] = trackJR(Data, maxdisp, maxrot);
toc

% return
%%
Titles = 'X Y Ux Uy t ID';
output = strcat(folder,dname);
dlmwrite([output '_tracked.txt'] ,Result,'-append', 'delimiter', '\t', 'precision', 9)
% dlmwrite([Directory File '_tracked.txt'], Titles, '\t')
% dlmwrite([Directory File '_tracked.txt'],Result, '-append', 'delimiter', '\t', 'precision', 9)
% disp('file is written...')

%% calculate Displacements and save them!
disp('Calculate Displacement for tau=1...');
[trans] = Displ_v2(Result,1,10000);
disp('...done')
disp('Calculate Tracklifetime')
life = Tracklifetime_v2(Result,10000);
disp('...done')
%%
fontsize1 = 16;
fontsize2 = 14;
Fig = figure('visible','off','Color',[1 1 1]);
% figure(1) 
subplot(1,2,1)
hist(trans, 1000)
set(gco, 'LineWidth',3)
title('Translation','Fontsize',fontsize1, 'FontWeight','bold','FontName','Helvetica' )
axis([0 400 -0.1 1])
xlabel('\Deltar [µm]', 'Fontsize',fontsize1, 'FontWeight','bold','FontName','Helvetica' )
ylabel('P(\Deltar)', 'Fontsize',fontsize1, 'FontWeight','bold', 'FontName','Helvetica' )
zlabel('Z-coord. mod 1', 'Fontsize',8, 'FontWeight','normal', 'FontName','Helvetica' )
set(gca, 'Fontsize',fontsize2)
set(gca, 'FontWeight','bold')
set(gca,'LineWidth',3)
%daspect([1 1 1])

% subplot(1,3,2)
% hist(rot*57,1000)
% set(gco, 'LineWidth',3)
% title('Rotation','Fontsize',fontsize1, 'FontWeight','bold','FontName','Helvetica' )
% axis([0 400 -0.1 1])
% xlabel('\Delta\theta [deg]', 'Fontsize',fontsize1, 'FontWeight','bold','FontName','Helvetica' )
% ylabel('P(\Delta\theta)', 'Fontsize',fontsize1, 'FontWeight','bold', 'FontName','Helvetica' )
% zlabel('Z-coord. mod 1', 'Fontsize',8, 'FontWeight','normal', 'FontName','Helvetica' )
% set(gca, 'Fontsize',fontsize2)
% set(gca, 'FontWeight','bold')
% set(gca,'LineWidth',3)
% % daspect([1 1 1])

subplot(1,2,2)
len = max(Result(:,5))-min(Result(:,5))+1;
hist(life(:,1),len)
set(gco, 'LineWidth',3)
title('Tracklength','Fontsize',fontsize1, 'FontWeight','bold','FontName','Helvetica' )
axis([0 400 -0.1 1])
xlabel('Number of steps', 'Fontsize',fontsize1, 'FontWeight','bold','FontName','Helvetica' )
ylabel('P(T)', 'Fontsize',fontsize1, 'FontWeight','bold', 'FontName','Helvetica' )
zlabel('Z-coord. mod 1', 'Fontsize',8, 'FontWeight','normal', 'FontName','Helvetica' )
set(gca, 'Fontsize',fontsize2)
set(gca, 'FontWeight','bold')
set(gca,'LineWidth',3)
%daspect([1 1 1])



%daspect([1 1 1])
%print(figure(1),[Directory File '_Displacements'],'-dpng','-r300')
%savefig(Fig,[output '_Displacements.fig'])
export_fig([output '_Displacements'],Fig,'-pdf','-tiff','-r300','-painters')

dlmwrite([output '_Displ_trans.txt'],trans,'\t');
% dlmwrite([output '_Displ_rot.txt'],rot,'\t');
dlmwrite([output '_Lifetime.txt'],life,'\t');
