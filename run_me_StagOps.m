% Explanation:
% This script organizes the data of 3 arrays into 1 and run the
% identification of the same neuron across days
% One annoying thing is that they don't consider neurons that disappear and
% reappears so I had to do dates separately
% addpath('/Users/jiaxintu/Documents/MATLAB/General-Helper-Functions')
clear
DataDir = '/Users/jiaxintu/Documents/Data/StagOpts_wrapper/neuralData/';
% DataDir = '/Volumes/My Passport/sorted/TrackingSameNeurons/New folder';
n_channels_per_array = 32;
n_arrays = 2; %
areaIDs = ['F','G'];

fileNames = dir([DataDir,'/*StagOps*']);
dates = regexp([fileNames(:).name],'(\d{6,8})','tokens');
dates = unique(vertcat(dates{:}));
n_dates = length(dates);
allcomb = [1,2;2,3;1,3];
% allcomb = [1,2;1,3;1,4;1,5;2,3;2,4;2,5;3,4;3,5;4,5]; % covers all possibility, I don't know how to be more intelligent on this
survived_neurons = cell(length(allcomb),3);
for irun = 1:length(allcomb)
    whichdates = allcomb(irun,:); % change me for different dates (either 2 any day, or 3 consecutive days)
    [channel,unit,spiketimes,wmean] = deal(cell(length(whichdates),n_arrays));
    
    for arrayID = areaIDs
        j = double(arrayID)-areaIDs(1)+1;
        array_offset = n_channels_per_array*(j-1);
        fileNames = dir([DataDir,'*',arrayID,'_StagOps*']);
        %     n_dates = length(fileNames);
        for iDate = whichdates
            neural_data = load(fullfile(DataDir,fileNames(iDate).name));
            fs = fieldnames(neural_data);
            expr = ['.*','(?<channel>\d{3})','(?<unit>[a-d]{1})','$'];
            [unitstr,units] = regexp(fs,expr,'match','tokens');
            units = [units{:}]';
            unitstr = [unitstr{:}]';
            units = vertcat(units{:});
            
            % remove <500spikes
            totalspk = arrayfun(@(i)length(neural_data.(unitstr{i})),1:length(unitstr));
            units = units(totalspk>=500,:);
            unitstr = unitstr(totalspk>=500);
            
            channel{iDate,j} = str2num(vertcat(units{:,1}))+ array_offset;
            unit{iDate,j} = double(vertcat(units{:,2}))-'a'+1;
            
            waveforms = find(contains(fs,'wf'));
            waveforms = waveforms(totalspk>=500);
            for i = 1:length(waveforms)
                %        neural_data.(fs{waveforms(i)}) = mean(neural_data.(fs{waveforms(i)}),2);
                wmean{iDate,j}{i} =  mean(neural_data.(fs{waveforms(i)}),2);
                spiketimes{iDate,j}{i} = neural_data.(unitstr{i});
            end
        end
    end
    %% Format the data into required input form
    for iDate = whichdates
        unit{iDate,1} = vertcat(unit{iDate,:});
        wmean{iDate,1} = horzcat(wmean{iDate,:});
        channel{iDate,1} = vertcat(channel{iDate,:});
        spiketimes{iDate,1} = horzcat(spiketimes{iDate,:});
    end
    
    unit = unit(:,1);
    wmean = wmean(:,1);
    channel = channel(:,1);
    spiketimes = spiketimes(:,1);
    
    empty_cells = cellfun(@(i)isempty(i),spiketimes);
    unit = unit(~empty_cells);
    wmean = wmean(~empty_cells);
    channel = channel(~empty_cells);
    spiketimes = spiketimes(~empty_cells);
    %%
    [survival, score, corrscore, wavescore, autoscore, basescore, correlations] = ...
        unitIdentification(channel,unit, spiketimes, wmean);
    %% What are the cell indices?
    % survived_neurons = cell(length(survival),2);
    % N.B.: survived_neurons is a cell containing matrix of 2 cols
    % first col is cell_idx the day before, and the second col is the cell_idx of the current day
    % for i = 1:length(survival)
    %     temp = find(survival{i});
    %     [survivedi,survivedj] = ind2sub(size(survival{i}),temp);
    %     survived_neurons{i,1} = [survivedi,survivedj] ;
    %     survived_neurons{i,2} = sprintf('%s,%s',dates{whichdates(i)}(5:8),dates{whichdates(i+1)}(5:8));
    % end
    
    temp = find(survival{1});
    [survivedi,survivedj] = ind2sub(size(survival{1}),temp);
    survived_neurons{irun,1} = [survivedi,survivedj] ;
    survived_neurons{irun,2} = sprintf('%s,%s',dates{whichdates(1)}(5:8),dates{whichdates(2)}(5:8));
    survived_neurons{irun,3} = cellfun(@length,unit)';
end

save(['Survived_neurons3.mat'],'survived_neurons')
%%
repeatedID = survived_neurons{1};
%% How many cells survived for n sessions?
% % cells that survived for all two sessions
% n_both = numel(intersect(survived_neurons{1}(:,2),survived_neurons{2}(:,1)));
% % cells that survived for either session
% n_either = arrayfun(@(i)size(survived_neurons{i},1),1:length(survived_neurons))-n_both;
%
% n_total = arrayfun(@(i)numel(unit{i}),1:length(unit)-1);
%
