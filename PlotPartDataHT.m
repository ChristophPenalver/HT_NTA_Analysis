function [] = PlotPartDataHT(ParticleData,dname,folder)
%DATAPLOTTER Summary of this function goes here
%   Detailed explanation goes here

%% TRACKLENGTH VS SIZE Gernerating Data for Plotting 
SortTL = sort(ParticleData(:,6));
s_comb = sortrows ([ParticleData(:,2),ParticleData(:,6)],2);
TLSize = [];
TLNumbers = [];
for tl = 1:1:100
    idxt = s_comb(:,2) == tl;
    s_idxt = sum (idxt);
    if s_idxt ~= 0
        MeanV = mean(s_comb(find(idxt),1));
        MaxV = max (s_comb(find(idxt),1));
        MinV = min (s_comb(find(idxt),1)); 
        TLs = [tl,MeanV,MaxV,MinV];
        TLSize = [TLSize;TLs];
%         f_nr = SortTL(idxt);
        [sf_nr,~] = size(SortTL(idxt));
        TLn = [tl,sf_nr];
        TLNumbers = [TLNumbers;TLn];
    end
end 
[leng,~] = size(TLSize);
Mean = ones(leng,1) .* mean(TLSize(:,2));
Mode = ones(leng,1) .* mode(TLSize(:,2));
Residue = abs(TLSize(:,2) - Mean);
%-------------------------------------------------------------------------%
%% TRACKLENGTH VS SIZE Creating PLOT 
% Create figure
h(1) = figure('visible','off','Color',[1 1 1]);
axis ([0 101 0 (leng+10)])
plot (TLSize(:,1),Mean,'-k','LineWidth',2)
hold on 
plot (TLSize(:,1),Mode,'-r','LineWidth',2)
hold on 
plot (TLSize(:,1),TLSize(:,2),'bx','Marker','x')
hold on
plot(TLSize(:,1),TLSize(:,3),'rx','Marker','^')
hold on
plot(TLSize(:,1),TLSize(:,4),'gx','Marker','v')
xlabel ('Tracklength [steps]','Fontsize',12, 'FontName','Helvetica')
ylabel ('Diameter [nm]','Fontsize',12 ,'FontName','Helvetica')
set (gca ,'LineWidth' ,1) 
legend ('Mean Value','Mode Value','Average Size','Max Value','Min Value',...
'Most frequent','Box','off','Location','northeast');
legend boxoff
%% TRACKLENGTH VS SIZE Save figure
figdes = '_SZ&TL';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,h,'-pdf','-tiff','-r300','-painters')
close(h)
TLSZ = array2table([TLSize,Mean,Mode]);
headers = {'Tracklength' 'AverageSize' 'Max' 'Min' 'Mean' 'Mode'};
TLSZ.Properties.VariableNames = headers;
writetable(TLSZ,foldname);
%% RESIDUAL Creating PLOT
i(1) = figure('visible','off','Color',[1 1 1]);
bar (TLSize(:,1),Residue)
xlabel ('Tracklength [steps]','Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel ('Absolute residue [nm]','Fontsize'  ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1) 
%% RESIDUAL Save figure
figdes = '_RESIDUE';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,i,'-pdf','-tiff','-r300','-painters')
close(i)
ResDue = array2table([TLSize(:,1),Residue]);
headers = {'Tracklength' 'AbsoluteResidue'};
ResDue.Properties.VariableNames = headers;
writetable(ResDue,foldname);
%% TRACKLENGTH VS NUMBEROFTRACKS Creating PLOT
g(1) = figure('visible','off','Color',[1 1 1]);
plot (TLNumbers(:,1),TLNumbers(:,2),'-k')
xlabel ('Tracklength [steps]','Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel ('Number of Tracks [counts]','Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1) 
hold off
%% TRACKLENGTH VS NUMBEROFTRACKS Save the figures
figdes = '_TL&NOT';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,g,'-pdf','-tiff','-r300','-painters')
close(g)
TLNOP = array2table(TLNumbers);
headers = {'Tracklength' 'NOT'};
TLNOP.Properties.VariableNames = headers;
writetable(TLNOP,foldname);
end

