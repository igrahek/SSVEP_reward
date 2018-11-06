
%% FSAReward, analysis ssVEP

% add path to EEGLAB + plugins & various functions
addpath(genpath('E:\Experiments\Grahek_Ivan\FSAReward\scripts\EEG\toolboxes_functions\'));

close all; clear all; clc;

% experiment settings
rng(7); % seed for pseudo-random number generator ('A scuppetta!)
expname = 'FSAReward'; % experiment name
pathexp = ['E:\Experiments\Grahek_Ivan\' expname]; % main directory
pathanalysis = 'analysis\EEG\Exp1'; % where to store the analysis
pathpreproc = 'preproc'; % where to store the preprocessed single-subject data (for grand average)
pathGA = 'mean'; % where to store the grand-averages
pathblocks = 'control_analyses\blocks'; % where to store the preprocessed single-subject data separated into blocks
Avg.channs64 = [pathexp '\scripts\EEG\Antonio_BioSemi64.locs']; % channel locations without external channels
Avg.prefix = 'VP'; % prefix of data files
begin_epoch = .5; % begin epoch (in seconds)
end_epoch = 3.25; % end epoch (in seconds)

% Participants:
%  * technical problems (no EEG recording from the beginning of the experiment, problems with battery):
%     - VP05, VP10, VP13, VP27
%  * percentage of rejected epochs (ssVEP preprocessing pipeline) in one or more conditions is >= 33%:
%     - VP08 [17.50, 12.50, 15.00, 7.50, 27.50, 22.50, 5.00, 6.67, 18.33, 18.33, 25.00, 35.00]
%     - VP35 [20.00, 15.00, 27.50, 32.50, 42.50, 35.00, 8.33, 18.33, 41.67, 23.33, 38.33, 30.00]
%     - VP39 [10.00, 7.50, 12.50, 10.00, 45.00, 40.00, 15.00, 6.67, 5.00, 3.33, 36.67, 25.00]
%     - VP44 [2.50, 7.50, 30.00, 37.50, 15.00, 22.50, 5.00, 5.00, 23.33, 33.33, 10.00, 18.33]
%     - VP46 [27.50, 27.50, 5.00, 15.00, 0.00, 5.00, 35.00, 36.67, 20.00, 21.67, 10.00, 5.00]

% final sample (N = 39):
Avg.subjects = [1 2 3 4 6 7 9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 36 37 38 40 41 42 43 45 47 48];
Avg.pathin = [pathexp '\' pathanalysis '\' pathpreproc '\']; % where to retrieve the preprocessed individual subject data (for grand average)
Avg.pathout = [pathexp '\' pathanalysis '\' pathGA '\']; % where to save the grand average

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
    'elec_coords', Avg.channs64, ... % electrode coordinates
    'labels', 'numbers', ... % labels to use in plot (numbers or labels)
    'subjects', Avg.subjects, ... % number of datasets
    'pathin', Avg.pathout, ... % where to take the data
    'amp_topomaplim', [0 1], ... % min/max amplitude value in topography (common across frequencies)
    'amp_spectramaplim', [0 1.6]); % min/max amplitude value in spectrum (common across frequencies)

% select participant's group for analysis
topo.subjnum = Avg.subjects; % all participants
topo.condnum = Avg.trig; % all conditions
topo.subjects = 1:numel(Avg.subjects); % index of all participants
topo.conds = 1:numel(Avg.trig); % index of all conditions
eventlabels = {'BslnRedAttended' 'BslnBlueAttended' 'AcqRedAttended' 'AcqBlueAttended' 'ExtRedAttended' 'ExtBlueAttended' 'BslnRedAttendedMovement' 'BslnBlueAttendedMovement' 'AcqRedAttendedMovement' 'AcqBlueAttendedMovement' 'ExtRedAttendedMovement' 'ExtBlueAttendedMovement'}; % condition labels
csvname = 'rewardBoth'; % name of csv file containing amplitude values at specified frequencies

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

% plot frequency amplitude spectrum & extract amplitude values
[amplitudes] = spec_plot_FSAReward(topo, eventlabels);
cell2csv([pathexp '\' pathanalysis '\amplitudes_' csvname '.csv'], amplitudes); % save amplitude values in a .csv file

%%
