function [x_max,x_maxDiff] = PlotAllTrHT(data, dname_alltr, folder_alltr)
%PLOTALLTRHT Summary of this function goes here
%   Detailed explanation goes here
%% SIZE DISTRIBUTION -- DATA COLLECTION
r_s = ceil(sort(data(find(diff([data(:,2);max(data(:,2))+1])),2)));
min_vals = min(r_s);
max_vals = max(r_s);
bins = (max_vals - min_vals);
[N,edges] = histcounts(r_s,bins);
edges = edges(1:(end-1));
%% SIZE DISTRIBUTION -- Calculation of MEAN VALUE
options = fitoptions('gauss2', 'Lower', [-Inf -Inf 0 -Inf -Inf 0]);
fitobject = fit(edges',N','gauss2',options);
%% SIZE DISTRIBUTION -- Creating Bar FIGURE
h(1) = figure ('visible','off','Color',[1 1 1]);
plot (edges,N,':k','LineWidth',0.5);
axis ([0 max(edges) 0 (max(N)+10)])
hold on
f = plot(fitobject);
f.LineWidth = 2.5;
[Peak, PeakIdx] = findpeaks(f.YData); 
x_max = f.XData(PeakIdx); 
% PeakVal = x_max;
text(f.XData(PeakIdx), Peak, sprintf('Max = %6.3f', round(x_max,2)),...
'VerticalAlignment','bottom')
hold on 

xlabel('Diameter [nm]', 'Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth',1)
legend (dname_alltr,'Gaussian fit')
legend boxoff
%% SIZE DISTRIBUTION -- SAVE the FIGURE
figdes = '_sizeDistr';
foldname = strcat(folder_alltr,dname_alltr,figdes);
export_fig(foldname,h,'-pdf','-tiff','-r300','-painters')
close (h)
SizeDistribution = array2table([edges',N']);
headers = {'Size' 'ParticleCounts'};
SizeDistribution.Properties.VariableNames = headers;
writetable(SizeDistribution,foldname);
%% DIFFUSION COEFFICIENT -- DATA COLLECTION
d_single = round(data(find(diff([data(:,3);max(data(:,3))+1])),3),-4);
min_valD = min(d_single);
max_valD = max(d_single);
bins = (max_valD - min_valD)/10000;
[N,edges] = histcounts(d_single,bins);
edges = edges(1:(end-1));

options = fitoptions('gauss2', 'Lower', [-0.1 -Inf 0 -0.1 -Inf 0]);
fitobject = fit(edges',N','gauss2',options);
%% DIFFUSION COEFFICIENT -- Creating Bar FIGURE
H(1) = figure ('visible','off','Color',[1 1 1]);
plot (edges,N,':k','LineWidth',0.5);
axis ([0 max(edges) 0 (max(N)+10)])
hold on
%% DIFFUSION COEFFICIENT -- FITTING GAUSS
f_D = plot(fitobject);
f_D.LineWidth = 2.5;
[Peak_D, PeakIdx_D] = findpeaks(f_D.YData); 
x_maxDiff = f_D.XData(PeakIdx_D);
x_maxD = (f_D.XData(PeakIdx_D)/1000000); 
text(f_D.XData(PeakIdx_D), Peak_D, sprintf('Max = %6.3f', x_maxD),...
'VerticalAlignment','bottom')
hold on 
%% DIFFUSION COEFFICIENT -- Generating FIGURE
xlabel('D [nm^2/s]', 'Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1) 
legend (dname_alltr,'Gaussian fit')
legend boxoff
%% DIFFUSION COEFFICIENT -- SAVE the FIGURE
figdes = '_DifcoDistr';
foldname_Dc = strcat(folder_alltr,dname_alltr,figdes);
export_fig(foldname_Dc,H,'-pdf','-tiff','-r300','-painters')
close(H)
DiffDistribution = array2table([edges',N']);
headers = {'DiffusionCoefficient_nm2s' 'Intensity'};
DiffDistribution.Properties.VariableNames = headers;
writetable(DiffDistribution,foldname_Dc);
end

