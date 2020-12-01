% add proteomics subfolders to path 
folder = strsplit(fileparts(which(mfilename)), filesep);
folder = folder(1:end-2);
addpath(genpath(fullfile(folder{:})));


%% wild-type QC data
clear, clc
filename = 'QC_wild-type.txt';

delimiter = '\t';

% read columns of data as strings:
formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';

% open the text file.
fileid = fopen(filename,'r');

% read columns of data according to format string.
dataArray = textscan(fileid, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

% close the text file.
fclose(fileid);

% convert the contents of columns containing numeric strings to numbers.
% replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[3,4,5,6,7]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% split data into numeric and cell columns.
rawNumericColumns = raw(:, [3,4,5,6,7]);
rawCellColumns = raw(:, [1,2]);

% replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

% create output variable
qc_wild_type = raw;

% clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

% construct wild-type data structure
for i=2:size(qc_wild_type,1)
    c1 = strsplit(qc_wild_type{i,1},'_');
    qc_wild_type_struct(i-1).pro = c1{1};
    qc_wild_type_struct(i-1).pep = c1{2};
    qc_wild_type_struct(i-1).rep = c1{4};
    qc_wild_type_struct(i-1).conf = qc_wild_type{i,2};
    qc_wild_type_struct(i-1).qc1 = qc_wild_type{i,3};
    qc_wild_type_struct(i-1).qc2 = qc_wild_type{i,4};
    qc_wild_type_struct(i-1).qc3 = qc_wild_type{i,5};
    qc_wild_type_struct(i-1).qc4 = qc_wild_type{i,6};
    qc_wild_type_struct(i-1).qc5 = qc_wild_type{i,7};
    
    tmp = [qc_wild_type{i,3},qc_wild_type{i,4},qc_wild_type{i,5},qc_wild_type{i,6},qc_wild_type{i,7}];
    
    qc_wild_type_struct(i-1).mean = mean(tmp);
    qc_wild_type_struct(i-1).std = std(tmp);
    qc_wild_type_struct(i-1).cov = qc_wild_type_struct(i-1).std/qc_wild_type_struct(i-1).mean;
end

dir = 'C:\Users\Francis Motta\Documents\GraduateSchool\Research\Biochronicity\proteomics\data\';
writetable(struct2table(qc_wild_type_struct), [dir, 'QC_wild-type_cov.txt'],'Delimiter','\t','WriteRowNames',true)

% =======================================================================

%% clb1-6 QC data
clear
clc
filename = 'QC_clb1-6.txt';

delimiter = '\t';

% read columns of data as strings:
formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';

% open the text file.
fileid = fopen(filename,'r');

% read columns of data according to format string.
dataArray = textscan(fileid, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

% close the text file.
fclose(fileid);

% convert the contents of columns containing numeric strings to numbers.
% replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[3,4,5,6,7]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% split data into numeric and cell columns.
rawNumericColumns = raw(:, [3,4,5,6,7]);
rawCellColumns = raw(:, [1,2]);

% replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

% create output variable
qc_clb16 = raw;

% clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

% construct clb1-6 data structure
for i=2:size(qc_clb16,1)
    c1 = strsplit(qc_clb16{i,1},'_');
    qc_clb16_struct(i-1).pro = c1{1};
    qc_clb16_struct(i-1).pep = c1{2};
    qc_clb16_struct(i-1).rep = c1{4};
    qc_clb16_struct(i-1).conf = qc_clb16{i,2};
    qc_clb16_struct(i-1).qc1 = qc_clb16{i,3};
    qc_clb16_struct(i-1).qc2 = qc_clb16{i,4};
    qc_clb16_struct(i-1).qc3 = qc_clb16{i,5};
    qc_clb16_struct(i-1).qc4 = qc_clb16{i,6};
    qc_clb16_struct(i-1).qc5 = qc_clb16{i,7};
    
    tmp = [qc_clb16{i,3},qc_clb16{i,4},qc_clb16{i,5},qc_clb16{i,6},qc_clb16{i,7}];
    
    qc_clb16_struct(i-1).mean = mean(tmp);
    qc_clb16_struct(i-1).std = std(tmp);
    qc_clb16_struct(i-1).cov = qc_clb16_struct(i-1).std/qc_clb16_struct(i-1).mean;
end

dir = 'C:\Users\Francis Motta\Documents\GraduateSchool\Research\Biochronicity\proteomics\data\';
writetable(struct2table(qc_clb16_struct), [dir, 'QC_clb1-6_cov.txt'],'Delimiter','\t','WriteRowNames',true)
