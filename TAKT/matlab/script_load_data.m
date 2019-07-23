%% Load proteomics data from 01/04/17 email exchange with Tina
% Data from the 5 quality control (QC) samples (cell samples from different
% time points mixed together and run through 5 technical replicates) from 
% wild-type yeast cells and from clb1-6 delete cells. 
% Replicate experiments are designated by "r1" and "r2". The data points 
% are ratios between the "heavy" spike-in peptides and the "light" 
% endogenous peptides.

% add proteomics subfolders to path 
clear
clc
folder = strsplit(fileparts(which(mfilename)), filesep);
folder = folder(1:end-2);
addpath(genpath(fullfile(folder{:})));

% -----------------------------------------------------------------------

%% wild-type QC data
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

% plot historgram of coefficient of variation
figure('color','w')
hist([qc_wild_type_struct.cov],30)
xlabel('CoV')
ylabel('counts')
[h,p] = lillietest([qc_wild_type_struct.cov],'Distr','exp')


% -----------------------------------------------------------------------

%% clb1-6 QC data
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

% plot historgram of coefficient of variation
figure('color','w')
hist([qc_clb16_struct.cov],30)
xlabel('CoV')
ylabel('counts')
[h,p] = lillietest([qc_clb16_struct.cov],'Distr','exp')

% -----------------------------------------------------------------------

%% wild-type time-series data
filename = 'wild-type_linearSpline_time-series.txt';

delimiter = '\t';

% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% open the text file.
fileID = fopen(filename,'r');

% read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

% close the text file.
fclose(fileID);

% convert the contents of columns containing numeric strings to numbers.
% replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
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
            % Convert numeric strings to numbers.
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
rawNumericColumns = raw(:, [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]);
rawCellColumns = raw(:, 1);

% create output variable
wildtype_timeseries = raw;
% clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;

% construct wild-type time series data structure
wildtype_timeseries_struct = struct('pro',[],'type',[],'pep',[],'rep',[],'ts',[]);
for i=2:size(wildtype_timeseries,1)
    c1 = strsplit(wildtype_timeseries{i,1},'_');
    wildtype_timeseries_struct(i-1).pro = c1{1};
    if length(c1) == 1
        wildtype_timeseries_struct(i-1).type = 'rna';
    else
       if strcmp(c1{2},'wb')
           wildtype_timeseries_struct(i-1).type = 'wb';
           wildtype_timeseries_struct(i-1).rep = c1{4};
       else
           wildtype_timeseries_struct(i-1).type = 'ms';
           wildtype_timeseries_struct(i-1).pep = c1{2};
           wildtype_timeseries_struct(i-1).rep = c1{4};
       end
    end
    wildtype_timeseries_struct(i-1).ts = [wildtype_timeseries{i,2:end}];
end

% find unique proteins
proteins = unique({wildtype_timeseries_struct.pro});

% for each protein compare replicates 1 and 2 
% mass spec data only
pro_pep_analysis_struct = struct('pro',[],'pep',[],'data',[]);
k=1;
for i=1:length(proteins)
    tmp1 = wildtype_timeseries_struct( strcmp({wildtype_timeseries_struct.pro},proteins{i}) & strcmp({wildtype_timeseries_struct.type},'ms'));
    peps = unique({tmp1.pep});
    for j=1:length(peps)
        tmp2 = tmp1( strcmp({tmp1.pep},peps{j}) );
        if length(tmp2) == 2
            pro_pep_analysis_struct(k).pro = proteins{i};
            pro_pep_analysis_struct(k).pep = peps{j};
            pro_pep_analysis_struct(k).data = [tmp2(1).ts;tmp2(2).ts];
            k=k+1;
        end
    end
end


proteins = unique({pro_pep_analysis_struct.pro});
pro_analysis_struct = struct('pro',[],'me',[]);
for i=1:length(proteins)
    tmp1 = pro_pep_analysis_struct( strcmp({pro_pep_analysis_struct.pro},proteins{i}) );
    peps = unique({tmp1.pep});
    tmp2 = [];
    for j=1:length(peps)
        tmp2 = [tmp2,tmp1(j).data];
    end
    
    for j=1:size(tmp2,2)
        if tmp2(1,j)+tmp2(2,j) == 0
            tmp3(j) = 0;
        else
            tmp3(j) = sqrt(2)*abs(tmp2(1,j)-tmp2(2,j))/((tmp2(1,j)+tmp2(2,j))/2);
        end
    end
    
    pro_analysis_struct(i).pro = proteins{i};
    pro_analysis_struct(i).me = mean(tmp3);
end

[~,idx] = sort([pro_analysis_struct.me]);
pro_analysis_struct = pro_analysis_struct(idx);

% histogram of mean normalized absolute error
figure('color','w')
hist([pro_analysis_struct.me],10)
mean([pro_analysis_struct.me])

% CDC28
figure('color','w')
maxlim = 1.05*max([pro_pep_analysis_struct(6).data(1,:),pro_pep_analysis_struct(6).data(2,:),pro_pep_analysis_struct(7).data(1,:),pro_pep_analysis_struct(7).data(2,:)]);
scatter(pro_pep_analysis_struct(6).data(1,:),pro_pep_analysis_struct(6).data(2,:),'filled')
hold
scatter(pro_pep_analysis_struct(7).data(1,:),pro_pep_analysis_struct(7).data(2,:),'filled')
[h,~] = legend({'pep 1', 'pep 2'});
axis square
xlim([0,maxlim])
xlim([0,maxlim])
line([0,maxlim],[0,maxlim])
xlabel('rep 1')
ylabel('rep 2')
title('CDC28')

% CLN2
figure('color','w')
maxlim = 1.05*max([pro_pep_analysis_struct(15).data(1,:),pro_pep_analysis_struct(15).data(2,:),pro_pep_analysis_struct(16).data(1,:),pro_pep_analysis_struct(16).data(2,:)]);
scatter(pro_pep_analysis_struct(15).data(1,:),pro_pep_analysis_struct(15).data(2,:),'filled')
hold
scatter(pro_pep_analysis_struct(16).data(1,:),pro_pep_analysis_struct(16).data(2,:),'filled')
[h,~] = legend({'pep 1', 'pep 2'});
axis square
xlim([0,maxlim])
xlim([0,maxlim])
line([0,maxlim],[0,maxlim])
xlabel('rep 1')
ylabel('rep 2')
title('CLN2')

% ELM1
figure('color','w')
maxlim = 1.05*max([pro_pep_analysis_struct(18).data(1,:),pro_pep_analysis_struct(18).data(2,:),pro_pep_analysis_struct(19).data(1,:),pro_pep_analysis_struct(19).data(2,:)]);
scatter(pro_pep_analysis_struct(18).data(1,:),pro_pep_analysis_struct(18).data(2,:),'filled')
hold
scatter(pro_pep_analysis_struct(19).data(1,:),pro_pep_analysis_struct(19).data(2,:),'filled')
[h,~] = legend({'pep 1', 'pep 2'});
axis square
xlim([0,maxlim])
xlim([0,maxlim])
line([0,maxlim],[0,maxlim])
xlabel('rep 1')
ylabel('rep 2')
title('ELM1')


k=1;
for i=1:size(pro_pep_analysis_struct,2)
    tmpdata = pro_pep_analysis_struct(i).data;
    for j=1:size(tmpdata,2)
        if tmpdata(1,j)+tmpdata(2,j) == 0
            alldata(k,1) = 0;
            alldata(k,2) = 0;
        else
            alldata(k,1) = abs(tmpdata(1,j)-tmpdata(2,j))/((tmpdata(1,j)+tmpdata(2,j))/2);
            alldata(k,2) = (tmpdata(1,j)+tmpdata(2,j))/2;
        end
        k = k+1;
    end
end
    
figure('color','w')
scatter(alldata(:,1),alldata(:,2),'filled')
axis square
xlabel('mean value of replicates')
ylabel('mean normalized absolute replicate difference')