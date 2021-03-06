
%% FSAReward, analysis ssVEP

close all; clear all; clc;

% experiment settings
rng(7); % seed for pseudo-random number generator ('A scuppetta!)
addpath(genpath('E:\Experiments\Grahek_Ivan\FSAReward\repo\EEG_preprocessing\')); % add path to EEGLAB + plugins & various functions
expname = 'FSAReward'; % experiment name
pathexp = ['E:\Experiments\Grahek_Ivan\' expname '\']; % main directory
pathdata = [pathexp 'data\EEG\']; % where to get the .bdf files
pathanalysis = [pathexp 'analysis\EEG\']; % where to store the logfile of artifact rejection
Avg.channs = [pathexp 'repo\EEG_preprocessing\toolboxes_functions\BioSemi68.locs']; % channel locations
begin_epoch = .5; % begin epoch (in seconds)
end_epoch = 3.25; % end epoch (in seconds)

% TRIGGERS:
% 1: baseline, red attended
% 2: baseline, blue attended
% 3: acquisition, red attended
% 4: acquisition, blue attended
% 5: extinction, red attended
% 6: extinction, blue attended
% 11: baseline, red attended, movement
% 12: baseline, blue attended, movement
% 13: acquisition, red attended, movement
% 14: acquisition, blue attended, movement
% 15: extinction, red attended, movement
% 16: extinction, blue attended, movement
Avg.trig = [1 2 3 4 5 6 11 12 13 14 15 16]; % select triggers

% final sample (N = 44):
Avg.subjects = [1 2 3 4 6 7 8 9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48];

Avg.pathssj = [pathanalysis 'preproc\']; % where to retrieve the preprocessed individual subject data (for grand average)
Avg.pathGM = [pathanalysis 'mo\mean\']; % where to save the grand average

%% grand averages: create & save separate grand averages for each condition

[Trials] = MeanData_FSAReward(Avg);

%% calculate and plot topography and frequency amplitude spectrum

% set variable for topography
topo = struct('numchans', 64, ... % number of scalp electrodes
    'freq', [10 12], ... % frequencies of interest
    'samprate', 512, ... % sampling rate
    'time', [begin_epoch end_epoch], ... % epoch length (in seconds)
    'timeplot', 1, ... % random time point (within the epoch length) used as reference to extract topography (in seconds)
    'maptype', 'hot', ... % color map (see help topoplot)
    'elec_coords', Avg.channs, ... % electrode coordinates
    'labels', 'numbers', ... % labels to use in plot (numbers or labels)
    'subjects', 1:numel(Avg.subjects), ... % index of all participants
    'subjnum', Avg.subjects, ... % all participants
    'conds', 1:numel(Avg.trig), ... % all conditions
    'condnum', Avg.trig, ... % all conditions
    'pathin', Avg.pathGM,... % where to take the grand averages
    'pathout', [pathexp 'repo\EEG_preprocessing\movement\'], ... % where to save the output
    'amp_topomaplim', [0 1], ... % min/max amplitude value in topography (common across frequencies)
    'amp_spectramaplim', [0 1.6]); % min/max amplitude value in spectrum (common across frequencies)

eventlabels = {'BslnRedAttended' 'BslnBlueAttended' 'AcqRedAttended' 'AcqBlueAttended' 'ExtRedAttended' 'ExtBlueAttended' 'BslnRedAttendedMov' 'BslnBlueAttendedMov' 'AcqRedAttendedMov' 'AcqBlueAttendedMov' 'ExtRedAttendedMov' 'ExtBlueAttendedMov'}; % condition labels

% plot topography
[elec_rank] = plot_topo_FSAReward(topo);

% select electrodes for spectral analysis
elec_choice = 1;
switch elec_choice % select method for electrode selection
    case 1 % select electrodes with max amplitude for each participant
        elec = findmaxelec(elec_rank(:, :, 1), 4, 0); % select 4 electrodes with highest amplitude
        elec = repmat(elec, numel(topo.subjects), 1); % select electrodes
    case 2 % manually select electrodes
        elec = repmat([28], numel(topo.subjects), 1);
end

% fields to add to the topo structure (for graph of amplitude spectra)
topo.channels = elec; % selected electrodes
topo.bound = [0 16]; % range of frequencies to plot

% plot frequency amplitude spectrum, extract & save amplitude values
[amplitudes] = spec_plot_FSAReward(topo, eventlabels);

%% extract single-subject, single-trial ssVEPs

filenames = dir([pathdata '*.bdf']); % read file names in folder path (puts it in a structure called filenames)

st_data = []; temp_st_data = [];
for isub = 1:numel(filenames) % loop through participants
    dispname = filenames(isub).name(1:end-4);
    EEG = pop_loadset([dispname '_elist_be_artrej.set'], Avg.pathssj);
    EEG = eeg_detrend_widmann(EEG); % detrend
    EEG = pop_select(EEG, 'channel', topo.channels(isub, :)); % select 4 electrodes with highest amplitude
    sprange = eegF_Time2Sp(EEG, topo.time(1)):eegF_Time2Sp(EEG, topo.time(2)); % time range of interest (i.e., excluding the initial ERP response)
    freqs = ((0:(numel(sprange) * 8) - 1) / (numel(sprange) * 8)) * topo.samprate; % frequencies
    [~, freq1] = min(abs(freqs - topo.freq(1))); % sampling point corresponding to first frequency of interest
    [~, freq2] = min(abs(freqs - topo.freq(2))); % sampling point corresponding to second frequency of interest
    % FFT zeropadded to 8 times the length of the data, in order to:
    % - increase frequency resolution and be able to visualize exactly 10 Hz and 12 Hz
    % - visually magnify possible differences between conditions
    temp_st_data = abs(fft(EEG.data(:, sprange, :), numel(sprange) * 8, 2)) * 2 / size(sprange, 2); % normalized
    temp_st_data = squeeze(mean(temp_st_data, 1)); % average across electrodes
    temp_st_data = temp_st_data([freq1 freq2], :)';  % extract frequencies of interest
    temp_st_data = num2cell(temp_st_data(:)); % put column 2 below column 1 and convert to cell
    % create matrix ready to be saved as .csv
    temp_ssj = repelem(cellstr(dispname), numel(temp_st_data))'; % repeat participant number
    temp_freq = num2cell(repelem(topo.freq, numel(temp_st_data) / 2))'; % repeat frequency
    temp_conds = num2cell(repmat([EEG.event(:).type], 1, 2))'; % copy condition
    temp_st_data = [temp_ssj temp_freq temp_conds temp_st_data]; % concatenate
    st_data = [st_data; temp_st_data]; % add to all data
end

st_data_header = {'participant' 'frequency' 'condition' 'amplitude'};
st_data = [st_data_header; st_data];
cell2csv([topo.pathout 'singleTrial_amplitudes.csv'], st_data);

%%
