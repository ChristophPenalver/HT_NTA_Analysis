function [] = Int_check(data,dname,folder)
%INT_CHECK Summary of this function goes here
%   Detailed explanation goes here

%% SIZE VS INTENSITY  -- DATA COLLECTION
% Preparing a matrix with the maximal intensity of the particle
IntensityData = [];
for idx = 1:max(data(:,1))
    int = ismember (data(:,1),idx);
    s_int = sum (int);
    if s_int ~= 0
        TempData = [max(data(int,1)),max(data(int,2)),max(data(int,7))];
        IntensityData = [IntensityData;TempData];
    end
end
%% SIZE VS INTENSITY  -- CREATING FIGURE
% plot intensity vs size of all particles
H = figure ('visible', 'off','Color',[1 1 1]);
plot (IntensityData(:,2),IntensityData(:,3),'bx');
xlabel ('Diameter [nm]','Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel ('Adjusted Intensity [a.u.]','Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1) 
%% AVERAGE INTENSITY VS SIZE -- DATA COLLECTION
AvIntensitySize =[];
new_sz = 0;
s_rwsz= sortrows(IntensityData,2);
for chk_sz = min(s_rwsz(:,2)):1:max(s_rwsz(:,2))
    new_sz = s_rwsz(:,2) <= (chk_sz) & s_rwsz(:,2) >= (chk_sz-1);
    AvIntensitySize = [AvIntensitySize;[chk_sz,mean(s_rwsz(new_sz,3))]];
end
%% AVERAGE INTENSITY VS SIZE -- CREATING FIGURE
% plot averrage intensity vs size
G = figure ('visible', 'off','Color',[1 1 1]);
plot (AvIntensitySize(:,1),AvIntensitySize(:,2),':kx');
xlabel ('Diameter [nm]','Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel ('Average Adjusted Intensity [a.u.]','Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1) 
%% AVERAGE SIZE VS INTENSITY -- DATA COLLECTION
IntensityAvSize = [];
new_row = 0;
s_row = sortrows(IntensityData,3);
for chk = 4:0.05:11
    s_int = s_row(:,3);
        new_row = s_int <= chk & s_int >= (chk-0.05);
        IntensityAvSize = [IntensityAvSize;[chk,mean(s_row(new_row,2))]];
end
[leng,~] = size(IntensityAvSize);
Mean = ones(leng,1).* nanmean(IntensityAvSize(:,2));
res = IntensityAvSize(:,2) - Mean;
Mode = ones(leng,1).* nanmean(IntensityAvSize(:,2));

% Elimination of NaN values
% ck = isnan(combdat(:,1));
% ck_not = ~ck;
% new_combdat = combdat(ck_not,:);
% x_pol = new_combdat(:,2);
% y_pol = new_combdat(:,1);
%% AVERAGE SIZE VS INTENSITY -- CREATING FIGURE
% Plot averrage intensity vs size 
F = figure('visible', 'off','Color',[1 1 1]);
plot (IntensityAvSize(:,1),IntensityAvSize(:,2),'kx');
hold on
%Plot mean
plot (IntensityAvSize(:,1),Mean,'-r')
hold on
%% AVERAGE SIZE VS INTENSITY -- Definition of parameters 
% axesLimits1 = xlim(axes1);
% fig = polyfit (combdat(:,2),combdat(:,1),3);
% x_pol_2 = linspace(min(combdat(:,2)), max(combdat(:,2)));
% y_pol_2 = polyval(fig,x_pol_2);
% plot(x_pol_2,y_pol_2,'r')
% hold on
%% AVERAGE SIZE VS INTENSITY -- PLOTTING FIGURE
xlabel ('Adjusted Intensity [a.u.]','Fontsize' ,12 , 'FontName', 'Helvetica')
ylabel ('Average Diameter [nm]','Fontsize' ,12 , 'FontName', 'Helvetica')
set (gca ,'LineWidth' ,1) 
legend (dname,'Average Size','Location','best');
legend boxoff
%% SIZE VS INTENSITY  -- SAVE FIGURE
figdes = '_IN-SZ';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,H,'-png','-tiff','-pdf','-q101','-eps','-r600','-painters')
IntensityData = array2table(IntensityData(:,2:3));
headers = {'Size' 'Intensity'};
IntensityData.Properties.VariableNames = headers;
writetable(IntensityData,foldname);
%% AVERAGE INTENSITY VS SIZE -- SAVE FIGURE
figdes = '_AvIN-SZ';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,G,'-png','-tiff','-pdf','-q101','-eps','-r600','-painters')
AvIntensitySize = array2table(AvIntensitySize);
headers = {'Size' 'AverageIntensity'};
AvIntensitySize.Properties.VariableNames = headers;
writetable(AvIntensitySize,foldname);
%% AVERAGE SIZE VS INTENSITY -- SAVE FIGURE
figdes = '_IN-AvSZ';
foldname = strcat(folder,dname,figdes);
export_fig(foldname,F,'-png','-tiff','-pdf','-q101','-eps','-r600','-painters')
IntensityAvSize = array2table(IntensityAvSize);
headers = {'Intensity' 'AverageSize'};
IntensityAvSize.Properties.VariableNames = headers;
writetable(IntensityAvSize,foldname);

end

