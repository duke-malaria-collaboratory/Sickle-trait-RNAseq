%% Compute mean normalized replicate error of time series data
% actually computing sqrt(2)*CoV(X) since this is equal to |x-y|/((x+y)/2)
% when the sample X = [x,y].

% load data

% add proteomics subfolders to path 
folder = strsplit(fileparts(which(mfilename)), filesep);
folder = folder(1:end-2);
addpath(genpath(fullfile(folder{:})));

% -----------------------------------------------------------------------
clear, clc
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
clearvars delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;

% -----------------------------------------------------------------------

% construct wild-type cell array
wildtype_timeseries_struct = struct('pro',[],'type',[],'pep',[],'rep',[],'ts',[]);
for i=2:size(wildtype_timeseries,1)
    c1 = strsplit(wildtype_timeseries{i,1},'_');
    wildtype_timeseries_struct(i-1).pro = c1{1};
    if length(c1) == 1
        wildtype_timeseries_struct(i-1).type = 'rna';
    else
       if strcmp(c1{2},'wb')
           wildtype_timeseries_struct(i-1).type = 'wb';
           wildtype_timeseries_struct(i-1).pep = '0';
           wildtype_timeseries_struct(i-1).rep = c1{4};
       else
           wildtype_timeseries_struct(i-1).type = 'ms';
           wildtype_timeseries_struct(i-1).pep = c1{2};
           wildtype_timeseries_struct(i-1).rep = c1{4};
       end
    end
    wildtype_timeseries_struct(i-1).ts = [wildtype_timeseries{i,2:end}];
end

% remove rna data
wildtype_timeseries_struct(strcmp({wildtype_timeseries_struct.type},'rna')) = [];

% find unique proteins
proteins = unique({wildtype_timeseries_struct.pro});

% for each protein compare replicates 1 and 2 
% mass spec data only
pro_pep_analysis_struct = struct('pro',[],'type',[],'pep',[],'data',[]);
k=2;
pro_pep_ts_cov{1,1} = 'pro';
pro_pep_ts_cov{1,2} = 'tec';
pro_pep_ts_cov{1,3} = 'pep';
pro_pep_ts_cov = {pro_pep_ts_cov{:},wildtype_timeseries{1,2:31}};
times = [wildtype_timeseries{1,2:31}]-min([wildtype_timeseries{1,2:31}]);
T = max(times);
pro_pep_ts_cov{1,34} = 'avg_cov';
pro_pep_ts_cov{1,35} = 'simp_snr';
pro_pep_ts_cov{1,36} = 'comp_snr';
for i=1:length(proteins)
    tmp1 = wildtype_timeseries_struct( strcmp({wildtype_timeseries_struct.pro},proteins{i}) );
    peps = unique({tmp1.pep});
    for j=1:length(peps)
        tmp2 = tmp1( strcmp({tmp1.pep},peps{j}) );
        
        tstmp = [];
        for t=1:length(tmp2(1).ts)
            datatmp = [];
            for l=1:length(tmp2)
                datatmp = [datatmp;tmp2(l).ts(t)];
            end
            tstmp = [tstmp, datatmp];
            
            if length(tmp2) == 1
                pro_pep_ts_cov{k,t+3} = NaN;
            elseif mean(datatmp) == 0
                pro_pep_ts_cov{k,t+3} = 0;
            else
                pro_pep_ts_cov{k,t+3} = std(datatmp)/mean(datatmp);
            end
        end
        
        Pn = 0;
        signal = mean(tstmp);
        noise = std(tstmp);
        for l=1:length(tmp2)
            Pn = Pn + 1/T*trapz(times,(tstmp(l,:)-signal).^2);
        end
        Pn = Pn/(length(tmp2)-1);
        Ps = 1/T*trapz(times,mean(tstmp.^2,1)) - Pn/(length(tmp2));
        
        pro_pep_ts_cov{k,1} = proteins{i};
        pro_pep_ts_cov{k,2} = tmp2(1).type;
        pro_pep_ts_cov{k,3} = tmp2(1).pep;
        pro_pep_ts_cov{k,34} = mean([pro_pep_ts_cov{k,4:31}]);
        pro_pep_ts_cov{k,35} = mean(signal)/std(noise);
        pro_pep_ts_cov{k,36} = sqrt(Ps/Pn);
        k=k+1;
    end
end

dir = 'C:\Users\Francis Motta\Documents\GraduateSchool\Research\Biochronicity\proteomics\data\';
fileID = fopen([dir , filename(1:end-4), '_rep_cov.txt'], 'w');
formatSpec = '%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\n';
fprintf(fileID,formatSpec,pro_pep_ts_cov{1,:});
formatSpec = '%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n';

[nrows,ncols] = size(pro_pep_ts_cov);
for row = 2:nrows
    fprintf(fileID,formatSpec,pro_pep_ts_cov{row,:});
end

fclose(fileID);


% -----------------------------------------------------------------------
clear, clc
filename = 'clb1-6_linearSpline_time-series.txt';

delimiter = '\t';

% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

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

for col=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
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
rawNumericColumns = raw(:, [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]);
rawCellColumns = raw(:, 1);

% create output variable
clb16_timeseries = raw;
% clear temporary variables
clearvars delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;

% -----------------------------------------------------------------------
% construct clb1-6 cell array
clb16_timeseries_struct = struct('pro',[],'type',[],'pep',[],'rep',[],'ts',[]);
for i=2:size(clb16_timeseries,1)
    c1 = strsplit(clb16_timeseries{i,1},'_');
    clb16_timeseries_struct(i-1).pro = c1{1};
    if length(c1) == 1
        clb16_timeseries_struct(i-1).type = 'rna';
    else
       if strcmp(c1{2},'wb')
           clb16_timeseries_struct(i-1).type = 'wb';
           clb16_timeseries_struct(i-1).pep = '0';
           clb16_timeseries_struct(i-1).rep = c1{4};
       else
           clb16_timeseries_struct(i-1).type = 'ms';
           clb16_timeseries_struct(i-1).pep = c1{2};
           clb16_timeseries_struct(i-1).rep = c1{4};
       end
    end
    clb16_timeseries_struct(i-1).ts = [clb16_timeseries{i,2:end}];
end

% remove rna data
clb16_timeseries_struct(strcmp({clb16_timeseries_struct.type},'rna')) = [];

% find unique proteins
proteins = unique({clb16_timeseries_struct.pro});

% for each protein compare replicates 1 and 2 
% mass spec data only
pro_pep_analysis_struct = struct('pro',[],'type',[],'pep',[],'data',[]);
k=2;
pro_pep_ts_cov{1,1} = 'pro';
pro_pep_ts_cov{1,2} = 'tec';
pro_pep_ts_cov{1,3} = 'pep';
pro_pep_ts_cov = {pro_pep_ts_cov{:},clb16_timeseries{1,2:29}};
times = [clb16_timeseries{1,2:29}]-min([clb16_timeseries{1,2:29}]);
T = max(times);
pro_pep_ts_cov{1,32} = 'avg_cov';
pro_pep_ts_cov{1,33} = 'simp_snr';
pro_pep_ts_cov{1,34} = 'comp_snr';
for i=1:length(proteins)
    tmp1 = clb16_timeseries_struct( strcmp({clb16_timeseries_struct.pro},proteins{i}) );
    peps = unique({tmp1.pep});
    for j=1:length(peps)
        tmp2 = tmp1( strcmp({tmp1.pep},peps{j}) );
        tstmp = [];
        for t=1:length(tmp2(1).ts)
            datatmp = [];
            for l=1:length(tmp2)
                datatmp = [datatmp;tmp2(l).ts(t)];
            end
            tstmp = [tstmp, datatmp];
            
            if length(tmp2) == 1
                pro_pep_ts_cov{k,t+3} = NaN;
            elseif mean(datatmp) == 0
                pro_pep_ts_cov{k,t+3} = 0;
            else
                pro_pep_ts_cov{k,t+3} = std(datatmp)/mean(datatmp);
            end
        end
        
        Pn = 0;
        signal = mean(tstmp);
        noise = std(tstmp);
        for l=1:length(tmp2)
            Pn = Pn + 1/T*trapz(times,(tstmp(l,:)-signal).^2);
        end
        Pn = Pn/(length(tmp2)-1);
        Ps = 1/T*trapz(times,mean(tstmp.^2,1)) - Pn/(length(tmp2));
        
        pro_pep_ts_cov{k,1} = proteins{i};
        pro_pep_ts_cov{k,2} = tmp2(1).type;
        pro_pep_ts_cov{k,3} = tmp2(1).pep;
        pro_pep_ts_cov{k,32} = mean([pro_pep_ts_cov{k,4:31}]);
        pro_pep_ts_cov{k,33} = mean(signal)/std(noise);
        pro_pep_ts_cov{k,34} = sqrt(Ps/Pn);
        k=k+1;
    end
end

dir = 'C:\Users\Francis Motta\Documents\GraduateSchool\Research\Biochronicity\proteomics\data\';
fileID = fopen([dir , filename(1:end-4), '_rep_cov.txt'], 'w');
formatSpec = '%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\n';
fprintf(fileID,formatSpec,pro_pep_ts_cov{1,:});
formatSpec = '%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n';

[nrows,ncols] = size(pro_pep_ts_cov);
for row = 2:nrows
    fprintf(fileID,formatSpec,pro_pep_ts_cov{row,:});
end

fclose(fileID);
