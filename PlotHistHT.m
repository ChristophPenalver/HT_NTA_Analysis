function [Allspecies,OvPer] = PlotHistHT(data_ev,appear_part,maximalval,maximalD,dname,folder,NumberOfParticles)
%PLOTHISTHT Summary of this function goes here
%   Detailed explanation goes here
%% !!!SIZE!!!
%% Compare the DiffCoeff or Size of the Clusterelements
% Look at Cluster and define the species of a oligomer
[len_cluster,~] = size (appear_part);
AddRow = [zeros(len_cluster,1),appear_part];
if len_cluster > 0
    [~,Code] = SpeciesDetermine(data_ev,AddRow,maximalval);
    Check4Cell = iscell(Code);
    Spec = zeros(length(Code),1);
    for xx = 1:length(Code)
        Test4Row = isnumeric(cell2mat(Code(xx,1)));
        if Test4Row == 0
            Convert = cell2mat(Code(xx,1));
            LengthConvert = length(Convert);
            if LengthConvert == 6
                Spec(xx,1) = min([str2num(Convert(1,1)),str2num(Convert(1,6))]);
            else 
                Spec(xx,1) = min([str2num(Convert(1:3)),str2num(Convert(7:end))]);
            end
        else
            Spec(xx,1) = cell2mat(Code(xx,1));
        end
    end
end 
Numbercodes = (1.00:0.50:8.00)';
Allspecies = {'1-Mer','2-Mer','3-Mer','4-Mer','Unknown'};
Specnames = {'Monomer','Pot. Monomer','Dimer','Pot. Dimer',...
    'Trimer Linear','Pot. Trimer Linear','Trimer Triangle',...
    'Pot. Trimer Triangle','Tetramer Linear','Pot. Tetramer Linear',...
    'Tetramer Square','Pot. Tetramer Square','Tetramer Tetrahedron',...
    'Pot. Tetramer Tetrahedron','Unknown'};
Percentages = zeros(length(Numbercodes),1);
OverallPercentages = zeros(5,1);
Extrema = zeros(15,2);
g = figure('visible','off','Color',[1 1 1]);
f = 1.5; % factor for dividing bins
runspec = {}; 
CM = jet(length(Extrema)); 
for xx = 1:15
%     IsSpec = ismember(Spec,Numbercodes(xx));
    SumSpec = sum(ismember(Spec,Numbercodes(xx)));
    RadiiSpec = ceil(AddRow(ismember(Spec,Numbercodes(xx)),3));
    if sum(RadiiSpec) >=1
        Extrema(xx,1) = min(RadiiSpec);
        Extrema(xx,2) = max(RadiiSpec);
    else
        Extrema(xx,1) = 0;
        Extrema(xx,2) = 0;
    end
    bins = ceil((Extrema(xx,2)/f - Extrema(xx,1)/f));
    if bins ~= 0 && xx ~= 15
        [Values,Edges] = histcounts(RadiiSpec,bins);
        Edges = Edges(1:(end-1));
        bar (Edges, Values,'FaceColor',CM(xx,:));
        runspec{1,xx} = cell2mat(Specnames(xx));
        hold on
    else 
        runspec{1,xx} = cell2mat(Specnames(xx));
    end
    Percentages(xx,1) = SumSpec/length(Spec);
    % Calculationg the overall percentage of each species.
    if xx<=2
        OverallPercentages(1,1) = OverallPercentages(1,1)+Percentages(xx,1);
    elseif xx<=4 && xx>=3
        OverallPercentages(2,1) = OverallPercentages(2,1)+Percentages(xx,1);
    elseif xx<=8 && xx>=5
        OverallPercentages(3,1) = OverallPercentages(3,1)+Percentages(xx,1);
    elseif xx<=14 && xx>=9
        OverallPercentages(4,1) = OverallPercentages(4,1)+Percentages(xx,1);
    else
        OverallPercentages(5,1) = OverallPercentages(5,1)+Percentages(xx,1);
    end
end
OverallPercentages = OverallPercentages.*100;
OvPer = OverallPercentages';
% Plotting the colored Histogram in a 
xlabel('Diameter [nm]', 'Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica')
% set (gca ,'LineWidth' ,1) 
legend (runspec,'Location','northeast','FontSize',9);
legend boxoff
hold off 
figdes = '_SpecHist';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,g,'-pdf','-tiff','-r300','-painters')
close(g)
%% Unknown Bar Plot:
h = figure('visible','off','Color',[1 1 1]);
bins = floor((Extrema(xx,2)/f - Extrema(xx,1)/f));
if bins ~= 0
    [Values,Edges] = histcounts(RadiiSpec,bins);
    Edges = Edges(1:(end-1));
    bar (Edges, Values,'FaceColor',[1 1 0]);
    runspec{1,xx} = cell2mat(Specnames(xx));
    hold on
end
xlabel('Diameter [nm]', 'Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1) 

figdes = '_UknSpecHist';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,h,'-pdf','-tiff','-r300','-painters')
close(h)
UknBar = array2table([Edges',Values']);
headers = {'Size' 'NumberOfParticles'};
UknBar.Properties.VariableNames = headers;
writetable(UknBar,foldname);
%%
% OverallPercentages = OverallPercentages.*100;
% % IsSpec = ismember(Spec,Unkwn);
% SumSpec = sum(ismember(Spec,Numbercodes(15,1)));
% Percentages = 100 .* SumSpec/length(Spec);
% OverallPercentages = [OverallPercentages;Percentages];
% OutputStrg = cell(5,1);
% for xx = 1:5
%     OutputStrg(xx,1) = {[cell2mat(Allspecies(1,xx)) ': ' num2str(OverallPercentages(xx,1)) ' %']};
% end
foldname = strcat(folder,dname,'_Results');
fileID = fopen([foldname '.txt'],'w');
fprintf(fileID,'%6s \r\n','Percentages of different oligomers:');
fprintf(fileID,'%12s \r\n',num2str(cell2mat(Allspecies)));
fprintf(fileID,'%6s \r\n',num2str(OvPer));
fclose(fileID);
% for ii = 1:1:5
%     fprintf(fileID,'%8s \r\n',cell2mat(OutputStrg(ii,:))');
% end
% fclose(fileID);

J = figure('visible','off','Color',[1 1 1]);
setfigure = categorical(Spec,Numbercodes, Specnames);
specplot = histogram(setfigure,'BarWidth',0.5,'Normalization','probability',...
    'FaceColor', [1 1 1 ]);
NoP = num2str(NumberOfParticles);
pt = [dname,NoP,' particles in total'];
title (pt, 'Fontsize', 16, 'FontName', 'Helvetica')
xlabel('Species', 'Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel('Occurance', 'Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1)

figdes = '_Hist';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,J,'-pdf','-tiff','-r300','-painters')
ExpData = array2table([Spec,Numbercodes,Specnames]);
headers = {'Species' 'NumberCode' 'SpeciesNames'};
writetable(ExpData,foldname);
close(J)
end

