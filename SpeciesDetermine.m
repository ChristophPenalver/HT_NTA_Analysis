function [Clusterspecies,Code,SpeciesCode] = SpeciesDetermine(largedata,Cluster,maximalvalue)
%Species Determine: Determines species of particle signals 
%   Checking the Diffusion coefficient of the particle signal to check if
%   the particle is a Monomer or a smaller cluster species. The script
%   compares Dt of the particle with the calculated particle data.

[Clusterlength,~] = size (Cluster);

Percentage = 0.02;
PotPercentage = 0.04;

% To calculate the used diffusion coefficient one has to use:            
siz_sp = floor(largedata (:,2)); 
% Calculation of lower Monomer limit
ll = ceil(maximalvalue - Percentage *maximalvalue);
mem = ismember(siz_sp,ll);
% If the searched size is not present in the data set, search for the
% next bigger size and take the MINIMUM of the Diffusion coefficient
% and evaluate further
if sum(mem) >= 1
    Dll = max(largedata(mem,3));
else
    while sum(mem) <=1
        ll = ll+1;
        mem = ismember(siz_sp,ll);
    end
    Dll = min(largedata(mem,3));
end
% Calculation of potential lower Monomer limit
potll = ceil(maximalvalue - PotPercentage*maximalvalue);
mem = ismember(siz_sp,potll);
% If the searched size is not present in the data set, search for the
% next bigger size and take the MINIMUM of the Diffusion coefficient
% and evaluate further
if sum(mem) >= 1
    potDll = max(largedata(mem,3));
else
    while sum(mem) <=1
        potll = potll+1;
        mem = ismember(siz_sp,potll);
    end
    potDll = min(largedata(mem,3));
end
%     Dll = a/ll;
%     potDll = a/potll;
% Calculation of higher Monomer limit
hl = floor(maximalvalue + Percentage*maximalvalue);
mem = ismember(siz_sp,hl);
% If the searched size is not present in the data set, search for the
% next smaller size and take the MAXIMUM of the Diffusion coefficient
% and evaluate further
if sum(mem) >= 1
    Dhl = min(largedata(mem,3));
else
    while sum(mem) <=1
        hl = hl-1;
        mem = ismember(siz_sp,hl);
    end
    Dhl = max(largedata(mem,3));
end
% Calculation of potential higher Monomer limit
pothl = floor(maximalvalue + PotPercentage*maximalvalue);
mem = ismember(siz_sp,pothl);
% If the searched size is not present in the data set, search for the
% next smaller size and take the MAXIMUM of the Diffusion coefficient
% and evaluate further
if sum(mem) >= 1
    potDhl = min(largedata(mem,3));
else
    while sum(mem) <=1
        pothl = pothl-1;
        mem = ismember(siz_sp,pothl);
    end
    potDhl = max(largedata(mem,3));
end
%% Species Evaluation
% One characterisitc number for each species (1 = Monomer, 1.5 =Pot
% Monomer, 2 = Dimer, 3 = Linear Trimer, 4 = Trigonal Trimer, 5 = Linear
% Tetramer, 6 = Square Tetramer, 7 = Tetragonal Tetramer, 8 = Unknown;
Numbercodes = (1.00:0.50:7.50)';
SpecCode = [1,1,2,2,3,3,3,3,4,4,4,4,4,4]';
Unkwn = 8.00;
UnkwnSpecies = {'Unknown'};
Omega = [1;1;0.735;0.735;0.589;0.589;0.651;0.651;0.5000;0.5000;0.578;0.578;0.605;0.605];
% Omega = [1;0.735;0.589;0.651;0.5000;0.578;0.605];

MaximaSpecies = zeros(1,length(Numbercodes));
MinimaSpecies = zeros(1,length(Numbercodes));
% PotMaximaSpecies = zeros(1,length(Numbercodes));
% PotMinimaSpecies = zeros(1,length(Numbercodes));
for spec = 1:1:length(Numbercodes)
    if rem(spec,2) == 1
        MaximaSpecies(1,spec) = Dhl * Omega(spec,1);
        MinimaSpecies(1,spec) = Dll * Omega(spec,1);
    else
        MaximaSpecies(1,spec) = potDhl * Omega(spec,1);
        MinimaSpecies(1,spec) = potDll * Omega(spec,1);
    end 
end
SpeciesNames = {'1-Mer';'Pot. 1-Mer';'2-Mer';'Pot. 2-Mer';'3_L-Mer';...
    'Pot. 3_L-Mer';'3_T-Mer';'Pot. 3_T-Mer';'4_L-Mer';'Pot. 4_L-Mer';...
    '4_S-Mer';'Pot. 4_S-Mer';'4_T-Mer';'Pot. 4_T-Mer';};
Clusterspecies = [];
Code = [];
SpeciesCode =[];
% Fitting = zeros(length(Numbercodes),1);
% DiffC = Cluster(Clusterlength,4);
for variable = 1:1:Clusterlength       
%     DiffC = Cluster(variable,3);
    DiffC = Cluster(variable,4);
    FitSpec = DiffC <= MinimaSpecies & DiffC >= MaximaSpecies;
    Rows = find(FitSpec);
    Even = rem(Rows,2)==0;
    HoldSpecies = {};
    NumberCode = 0;
    if sum(rem(Rows,2))>=2
        Multispec = Rows(find(rem(Rows,2)));
        for xx = 1:sum((rem(Rows,2)))
            HoldSpecies(xx,1) = SpeciesNames(Multispec(1,xx));
            HoldCode(xx,1) = Numbercodes(Multispec(1,xx));
            HoldSpecCode(xx,1) = SpecCode(Multispec(1,xx));
        end
        Species = {[cell2mat(HoldSpecies(1,1)) ' or ' cell2mat(HoldSpecies(2,1))]};
        NumberCode = {[num2str(HoldCode(1,1)) ' or ' num2str(HoldCode(2,1))]};
        SpCd = {[num2str(HoldSpecCode(1,1)) ' or ' num2str(HoldSpecCode(2,1))]};
    elseif sum(Even)>=2
        Multispec = Rows(find(Even));
        for xx = 1:sum(Even)
            HoldSpecies(xx,1) = SpeciesNames(Multispec(1,xx));
            HoldCode(xx,1) = Numbercodes(Multispec(1,xx));
            HoldSpecCode(xx,1) = SpecCode(Multispec(1,xx));
        end
        Species = {[cell2mat(HoldSpecies(1,1)) ' or ' cell2mat(HoldSpecies(2,1))]};
        NumberCode = {[num2str(HoldCode(1,1)) ' or ' num2str(HoldCode(2,1))]};
        SpCd = {[num2str(HoldSpecCode(1,1)) ' or ' num2str(HoldSpecCode(2,1))]};
    elseif sum(rem(Rows,2)) == 1
        HoldSpecies = SpeciesNames{Rows(find(rem(Rows,2))),1};
        NumberCode = Numbercodes(Rows(find(rem(Rows,2))),1);
        SpCd = SpecCode(Rows(find(rem(Rows,2))),1);
        NumberCode = {NumberCode};
        Species = {HoldSpecies};
        SpCd = {SpCd};
    elseif sum(Even) == 1 
        HoldSpecies = SpeciesNames{Rows(find(Even)),1};
        NumberCode = Numbercodes(Rows(find(Even)),1);
        SpCd = SpecCode(Rows(find(Even)),1);
        NumberCode = {NumberCode};
        Species = {HoldSpecies};
        SpCd = {SpCd};
    else   
        Species = {UnkwnSpecies};
        NumberCode = 8;
        NumberCode = {NumberCode};
        SpCd = 8;
        SpCd = {SpCd};
    end
    Clusterspecies = [Clusterspecies;Species];
    Code = [Code;NumberCode];
    SpeciesCode = [SpeciesCode;SpCd];
end
% Code = Code';
end


