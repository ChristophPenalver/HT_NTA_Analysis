function [AllTracks,AllTracksTable] = Trackcalculator_v2(AlltrackData,dname,folder,Temperature)
%TRACKCALCULATOR_V2 Summary of this function goes here
%   Detailed explanation goes here
%% Script for Track Analysis
% Calculation of particle movement with data given by the PTM.
% Alltracks Import
% [infile,inpath] = uigetfile ({'*.csv', 'CSV Files (*.csv)';...
%         '*.*', 'ALL Files(*.*)'},'Import an Alltracks file','MultiSelect', 'on');
% [outfile,folder] = uiputfile ({'*.xlsx', 'EXCEL Files (*.xlsx)';...
%         '*.*', 'ALL Files(*.*)'},'Select an export file for Data Information and Size Distributions');
% dataname = outfile_alltr(1:(end-5));
output = strcat (folder,dname);
% Ex = 5;
% AllTracks = csv_Parser_alltr (infile_alltr,inpath_alltr,Ex); 
% tof = Alltracks(:,8);
AllTracks = AlltrackData((AlltrackData(:,8)==1),1:7);
% Retracking of particles with tracker of Jörg Roller
[Result] = StartTracking(AllTracks,folder,dname);
% Driftcorrection of data
[Dedrift] = Dedrift_v3b(Result,folder,dname);
NTASize = AllTracks((diff([AllTracks(:,2);max(AllTracks(:,2)+1)]) ~=0),2);
NTADt = AllTracks((diff([AllTracks(:,2);max(AllTracks(:,2)+1)]) ~=0),3);
clearvars  Result Alltracks;
% convert to usable data in Script
Dedrift = sortrows(Dedrift,6);
%% Definition of particle information. Look at every cluster
% How many particles in data set?
% Check for multiple calls of ID, take first row of occurance
DedriftValues = [1,2,5,6];
Allpart = Dedrift((diff([Dedrift(:,6);max(Dedrift(:,6)+1)]) ~=0),DedriftValues);
[Num_part,~] = size(Dedrift((diff([Dedrift(:,6);max(Dedrift(:,6)+1)]) ~=0),DedriftValues));
%% Constants for calculations:
% Boltzmann constant in J/K = N*m/K
kb = 1.38064852*10^-23;
% % Temperature input by user
% prompt = {'Enter a value of temperature (in °C)'};
% dlgtitle = 'Temperature in °C';
% definput = {'25'};
% opts.Interpreter = 'tex';
% answer = inputdlg(prompt,dlgtitle,[1 40],definput,opts);
% Temperature = str2num(cell2mat(answer));
% Temperature in Kelvin
T = Temperature + 273.15;
% Viscosity of water at 25 °C in Pa*s = N*s*m^-2 = 0.000891;
% Viscosity of water at 21.8 °C in Pa*s = N*s*m^-2 = 0.000929;
[viscosity] = VisCalc(Temperature);
viscosity = viscosity/1000;

%% Evaluating data and creating MSD Data
MSD1 = figure('visible','off','Color',[1 1 1]);
w=waitbar(0, 'Evaluating...'); 
% particle_info = zeros(Num_part,9);
Sizes = zeros(Num_part,3);
% PdeltaX =zeros(Num_part,3);
yy = 1;
frame = 0.0333;
% tvalue = tic;
AllID = [];
AllSizes = [];
AllDt = [];
AllFrame = [];
AllX = [];
AllY = [];
for particle_idx = 1:max(Dedrift(:,6))
    part_member = ismember(particle_idx,Dedrift(:,6));
%     trav_nm = 0;
    msd2d = zeros(length(Dedrift((Dedrift(:,6)==particle_idx))),1); %Prealloc msd-vector
    stat = zeros(length(Dedrift((Dedrift(:,6)==particle_idx))),1); %Prealloc statistic values
    MSD2D = zeros(length(Dedrift((Dedrift(:,6)==particle_idx))),1);
    dist = zeros(length(Dedrift((Dedrift(:,6)==particle_idx))),1);
    if part_member >= 1 && length(dist)>1
        % Use coordinates of the particle 
        place = Dedrift(Dedrift(:,6)==particle_idx,1:2);
        [steps,~] = size (place);
        Tracklength = steps;
        if Tracklength >=5 % Just particles above a certain step length are evaluated
            for ni=1:1:steps-1
                for mi=ni+1:steps
                    pos1=place(ni,1:2);
                    pos0=place(mi,1:2);
                    if isempty(pos1)==false && isempty(pos0)==false
                       msd2d(mi-ni,1)= msd2d(mi-ni,1)+ norm(pos1-pos0)^2;
                       dist(mi-ni,1)=dist(mi-ni,1)+ norm(pos1-pos0);%sum for msd
                       stat(mi-ni,1)=stat(mi-ni,1)+1; %Counting statistic, therefore indiv. N
                    end
                end
            end

            Distance=sum(dist);
            MSD2D = msd2d./stat; % in pixel 
            MSD2D = 166^2.*MSD2D(1:end-1); % convert to nm 
%             MSD2D = MSD2D(~isnan(MSD2D))'; % delete NaN values
            stat = stat(1:end-1);
            Step = 1:steps-1;
            % plot MSD2D in plot vs tau 
            plot((Step*frame),MSD2D);
            hold on

            % Calculations with travelled distance
            trav_nm = sum(Distance).*166^2;
            % Dt in nm^2/s
%             D_nm = (trav_nm^2)/(4*(Tracklength/30));
            % Dt in m^2/s
            D_m = (trav_nm*10^-9)^2/(4*(Tracklength/30));
            % Radius in m
            r = ((kb.*T)./(D_m.*6.*pi.*viscosity));
            % Diameter in nm
            Size = 2.*r.*10^9;

            % Calculations with MSD2D
            D_MSD = MSD2D./(4.*1*frame);
            % Calculation of D in m²/s
            D_MSD = D_MSD.*10^-18;
            % Calculation of Particle diameter by stokes Einstein Equation
            Diameter_MSD_Einstein = 2.*10^9.*((kb.*T)./(D_MSD.*6.*pi.*viscosity));
            
            % Calculation of Size with equation (1) from John G Walker:
            % Improved nano-particle tracking analysis (2012) Meas. Sci.
            % Technol. 23(2012) 
            Diameter_MSD = 2.*10^9.*((2.*kb.*T.*frame)/(3.*pi.*viscosity.*(MSD2D(1,1)*10^-18)));
            
            % Generation of Allrtacks-Data:
            Entries = find(Dedrift(:,6) == particle_idx);
            NumEntr = length(Entries);
            ParticleID = zeros(NumEntr,1);
            ParticleSize = zeros(NumEntr,1);
            ParticleDt = zeros(NumEntr,1);
            for xx = 1:NumEntr
                ParticleID(xx,1) = particle_idx;
                ParticleSize(xx,1) = Diameter_MSD_Einstein(1,1);
                ParticleDt(xx,1) = D_MSD(1,1)*1E18;
            end
            AllID = [AllID;ParticleID];
            AllSizes = [AllSizes;ParticleSize];
            AllDt = [AllDt;ParticleDt];
            AllFrame = [AllFrame;Dedrift(Dedrift(:,6)==particle_idx,5)];
            AllX = [AllX;Dedrift(Dedrift(:,6)==particle_idx,1)];
            AllY = [AllY;Dedrift(Dedrift(:,6)==particle_idx,2)];
            
            Sizes(yy,:) = [Size,Diameter_MSD(1,1),Diameter_MSD_Einstein(1,1)];        
            yy = yy+1;
            % Display the progress of the script in the command window
            % (time consuming)
%             disp([datestr(now) ' - ' ' Evaluated particle: ' num2str(particle_idx) ' of ' num2str(max(Dedrift(:,8)))]);
        end
    end
%     tElapsed = toc(tvalue);
    % update waitbar
    if rem(particle_idx, 50) == 0 
        waitbar(particle_idx/max(Dedrift(:,6)));
    end
end
xlabel('\tau [s]', 'Fontsize',16 ,'FontName', 'Helvetica','FontWeight','bold')
ylabel('MSD(\tau) [nm²]', 'Fontsize',16 ,'FontName','Helvetica','FontWeight','bold')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
set (gca ,'LineWidth' ,1.5); 
set(gcf,'color','w');
close(w)
export_fig([output '_MSD'],MSD1,'-pdf','-tiff','-r300','-painters')

% Size Distribution of data
Diameters = array2table(Sizes);
headers = {'Calc_Size' 'Size_MSD' 'Size_MSD_Einstein'};
Diameters.Properties.VariableNames = headers;

%% Generation of NTA Size Distribution
H1 = figure('visible','off','Color',[1 1 1]);
if max(NTASize)<= 1000
    edges = 0:0.5:max(NTASize);
else
    edges = 0:0.5:1000;
end
[h1,~] = hist (NTASize,edges);
hist (NTASize,edges);
h1smoothed = smoothdata(h1,'sgolay',20);
xlabel('Diameter [nm]', 'Fontsize' ,12 , 'FontName', 'Helvetica');
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica');
% Sets font size and line strength of axis (gca)
axis([0 max(edges) 0 (max(h1)+15)]);
set (gca ,'LineWidth' ,1.5); 
set(gcf,'color','w');
output1 = [output '_NTASize'];
Name = 'NTA';
[Res1,GOF1,Loc1] = Gauss2FitSize(edges,h1smoothed,output1,Name);
%% Generation of Einstein MSD Size Distribution
H2 = figure('visible','off','Color',[1 1 1]);
if max(Diameters.Size_MSD_Einstein)<= 1000
    edges = 0:0.5:max(Diameters.Size_MSD_Einstein);
else
    edges = 0:0.5:1000;
end
[h2,~] = hist (Diameters.Size_MSD_Einstein,edges);
hist (Diameters.Size_MSD_Einstein,edges);
h2smoothed = smoothdata(h2,'sgolay',20);
xlabel('Diameter [nm]', 'Fontsize' ,12 , 'FontName', 'Helvetica');
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica');
% Sets font size and line strength of axis (gca)
axis([0 max(edges) 0 (max(h2)+15)]);
set (gca ,'LineWidth' ,1.5); 
set(gcf,'color','w');
output2 = [output '_EinsteinSize'];
Name = 'Einstein';
[Res2,GOF2,Loc2] = Gauss2FitSize(edges,h2smoothed,output2,Name);
%% Generation of NTA DiffusionCoefficient Distribution
D1 = figure('visible','off','Color',[1 1 1]);
if max(NTADt)<= 2E7
    edges = 1E5:1E4:max(NTADt);
else
    edges = 0:0.5:1000;
end
[d1,~] = hist(NTADt,edges);
hist (NTADt,edges);
d1smoothed = smoothdata(d1,'sgolay',20);
xlabel('Dt [nm²/s]', 'Fontsize' ,12 , 'FontName', 'Helvetica');
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica');
% Sets font size and line strength of axis (gca)
axis([0 max(edges) 0 (max(d1)+15)]);
set (gca ,'LineWidth' ,1.5); 
set(gcf,'color','w');
output3 = [output '_NTADt'];
Name = 'NTA-D_t';
[Res3,GOF3,Loc3] = Gauss1FitDt(edges,d1smoothed,output3,Name);
%% Generation of Einstein DiffusionCoefficient Distribution
D2 = figure('visible','off','Color',[1 1 1]);
DtEinst = AllDt((diff([AllDt(:,1);max(AllDt(:,1)+1)]) ~=0),1);
if max(DtEinst)<= 2E7
    edges = 1E5:1E4:max(DtEinst);
else
    edges = 0:0.5:1000;
end
[d2,~] = hist(DtEinst,edges);
hist (DtEinst,edges);
d2smoothed = smoothdata(d2,'sgolay',20);
xlabel('Dt [nm²/s]', 'Fontsize' ,12 , 'FontName', 'Helvetica');
ylabel('Particle Counts', 'Fontsize' ,12 , 'FontName', 'Helvetica');
% Sets font size and line strength of axis (gca)
axis([0 max(edges) 0 (max(d2)+15)]);
set (gca ,'LineWidth' ,1.5); 
set(gcf,'color','w');
output4 = [output '_EinsteinDt'];
Name = 'Einstein-D_t';
[Res4,GOF4,Loc4] = Gauss1FitDt(edges,d2smoothed,output4,Name);
%% Output of Files for further analysis
% 1: ParticleID
% 2: Size in nm
% 3: Diffusion coefficient
% 4: Frame 
% 5: X in pixels
% 6: Y in pixels
% 7: Intensity
% 8: included in distribution

Intensity = zeros(length(AllY),1);
Included = ones(length(AllY),1);
AllTracks = [AllID,AllSizes,AllDt,AllFrame,AllX,AllY,Intensity,Included];
AllTracksTable = array2table(AllTracks);
headers = {'ParticleID' 'Diameter' 'D_t' 'Frame' 'X' 'Y' 'Intensity' 'Included'};
AllTracksTable.Properties.VariableNames = headers;

end

