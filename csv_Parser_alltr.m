function [ data_alltr] = csv_Parser_alltr( inpath_alltr,infile_alltr,Ex )
%CSV_PARSER Parse xlsx file to array
%   Loads the input Excel File(important convert to xlsx beforehand) 

%% Initialize variables.
filename_alltr = [infile_alltr,inpath_alltr];
delimiter = ',';
if Ex == 5 
    startRow = 2;
else
    startRow = 1;
end
%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: categorical (%C)
% For more information, see the TEXTSCAN documentation.
if Ex == 5
    formatSpec = '%f%f%f%f%f%f%f%C%[^\n\r]';
else
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
end

%% Open the text file.
fileID = fopen(filename_alltr,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
if Ex ==5 
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1,...
    'ReturnOnError', false, 'EndOfLine', '\r\n');
else
    dataArray = textscan(fileID, formatSpec, 'Delimiter',...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN,...
    'ReturnOnError', false);
end
%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
if Ex == 5
    res = table(dataArray{1:end-1}, 'VariableNames', {'ParticleID','Sizenm',...
        'Diffusioncoefficientnm2s1','Frame','Xpixels','Ypixels',...
        'LnAdjustedintensityAU','Includedindistribution'});
    vars = {'ParticleID','Sizenm',...
        'Diffusioncoefficientnm2s1','Frame','Xpixels','Ypixels',...
        'LnAdjustedintensityAU'};
    data_ind = table2array(res(1:end,vars));
    true = res.Includedindistribution;
    exc = true == 'True';
    data_alltr = [data_ind,exc];
else
    res = table(dataArray{1:end-1}, 'VariableNames', {'ParticleID','Sizenm',...
        'Diffusioncoefficientnm2s1','Frame','Xpixels','Ypixels',...
        'LnAdjustedintensityAU'});
    vars = {'ParticleID','Sizenm',...
        'Diffusioncoefficientnm2s1','Frame','Xpixels','Ypixels',...
        'LnAdjustedintensityAU'};
    data_ind = table2array(res(1:end,vars));
    data_alltr = data_ind;
end
%% Clear temporary variables
clearvars filename_alltr delimiter startRow formatSpec fileID dataArray ans;

end