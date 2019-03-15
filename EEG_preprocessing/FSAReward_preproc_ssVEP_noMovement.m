
%% FSAReward, preprocessing for analysis ssVEP

close all; clear all; clc;

% experiment settings
rng(7); % seed for pseudo-random number generator ('A scuppetta!)
addpath(genpath('E:\Experiments\Grahek_Ivan\FSAReward\scripts\EEG\toolboxes_functions\')); % add path to EEGLAB + plugins & various functions
expname = 'FSAReward'; % experiment name
pathexp = ['E:\Experiments\Grahek_Ivan\' expname '\']; % main directory
pathdata = [pathexp 'data\EEG\']; % where to get the .bdf files
pathanalysis = [pathexp 'analysis\EEG\']; % where to store the logfile of artifact rejection
Avg.channs = [pathexp '\scripts\EEG\BioSemi68.locs']; % channel locations
% Avg.channs64 = 'E:\Experiments\EmoSSR_faces\toolboxes_functions\leipzig64_Christopher.locs'; % channel locations for plot_topo.m
begin_epoch = 0; % begin epoch (in seconds)
end_epoch = 3.25; % end epoch (in seconds)

% reference = 'mastoids'; % change to 'average' if you want to reference to the average of all scalp electrodes

% TRIGGERS (only conditions in which no dot movement occurred):
% 1: baseline, red attended; 2: baseline, blue attended;
% 3: acquisition, red attended; 4: acquisition, blue attended;
% 5: extinction, red attended; 6: extinction, blue attended.
Avg.trig=[1 2 3 4 5 6]; % select triggers

filenames = dir([pathdata '*.bdf']); % read file names in folder path (puts it in a structure called filenames)
% discarded participants due to  technical problems (no EEG recording from the beginning of the experiment, problems with battery):
% - VP05, VP10, VP13, VP27
% 
% % participants from which we connected 64 electrodes but mistakenly recorded 128 channels
% % (the signal from the second electrode bundle is flat)
% VPN128 = ['VP37' 'VP38' 'VP39' 'VP40' 'VP41' 'VP42' 'VP43' 'VP44' 'VP45' 'VP46' 'VP47' 'VP48'];

logfile = []; elapsed = []; % preallocate logfile and matrix with elapsed time (in seconds) of single-subject preprocessing

%% PREPROCESSING

% start preprocessing
for isub = 1:numel(filenames) % loop through files
    
    tic % start timer
    
    % file name
    dispname = filenames(isub).name(1:end-4);
    
    % display participant number
    disp('***********************')
    disp(['Processing ' dispname '...'])
    disp('***********************')
    
    % due to problems recording mastoids in several participants, we opted for average reference
    EEG = pop_biosig([pathdata dispname '.bdf'], 'channels', 1:68, 'ref', 1:64, 'refoptions', {'keepref' 'on'});
    
    %     pop_eegplot(EEG,1,1,1); % check data
    
    EEG = pop_chanedit(EEG, 'load', {Avg.channs, 'filetype', 'autodetect'}); % assign channel locations
    %     topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); % see channel locations
    EEG = pop_rmbase(EEG, 'Warning', 'off'); % remove DC offset
    
    %     pop_eegplot(EEG,1,1,1); % check data
    
    EEG = pop_epoch(EEG, num2cell(Avg.trig), [begin_epoch end_epoch]); % epoching
    
    %     pop_eegplot(EEG, 1, 1, 1); % check data

    % for log file
    for irej = 1:numel(EEG.epoch)
        trigtotepochs{irej, :} = EEG.epoch(1, irej).eventtype; % create vector with the sum of total epochs for each condition
    end
    trigtotepochs = cell2mat(trigtotepochs)'; % convert this vector from cell to number (otherwise the "sum" function won't work) and transpose
    for k = 1:numel(Avg.trig)
        totepochs(:, k) = sum(trigtotepochs == Avg.trig(k)); % find all triggers belonging to each condition and sum them
    end
    
    %%%%%%%%%%%%%%%%%%
    %%% ARTIFACT REJECTION %%%
    %%%%%%%%%%%%%%%%%%
    
    % channel interpolation & artifact rejection using FASTER
    cfg = []; % create Fieldtrip-like structure
    cfg.datachan = 1:64; % select scalp electrodes
    cfg.thresh = [3 3 0 3 3 12]; % see help eegF_FASTER for a description of each number (lower numbers are more conservative)
    [gen_bad_chans, EEG, trials2removeFASTER] = eegF_FASTER(cfg, EEG); % find & interpolate noisy channels; find trials contaminated by artifacts (does not eliminate them!)
    EEG = pop_reref(EEG, 1:64,'keepref','on'); % re-reference to average (the eegF_FASTER function re-references to Cz, hence the need for re-referencing here)
    
    % detect residual blinks and horizontal eye movements with SCADS
    EEG=eegF_Bipolarize(EEG);  % create bipolar vEOG and hEOG channels
    [EEG,trials2removeBlinks]=SCADS_RemoveBlinks(EEG,65); % find trials contaminated by blinks (does not eliminate them!)
    [EEG,trials2removeEyeMov]=SCADS_RemoveEyeMovements(EEG,66); % find trials contaminated by horizontal eye movements (does not eliminate them!)
    EEG=pop_select(EEG,'nochannel',[65 66]); % remove bipolar eye channels
    
    TotTrials2Remove=unique([trials2removeFASTER,trials2removeBlinks,trials2removeEyeMov]); % total trials to remove according to all artifact rejection methods
    
    EEG=pop_select(EEG,'notrial',TotTrials2Remove); % remove epochs contaminated by artifacts
    
    %     pop_eegplot(EEG,1,1,1); % check data
    
    pop_saveset(EEG,'filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d')],'filepath',[pathexp '\' pathanalysis '\' pathpreproc '\']); % save (for grand averages)
    % EEG=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '.set'],'filepath',[pathexp '\' pathanalysis '\' pathpreproc '\']); % load file
    
    % for log file
    for irej2=1:numel(EEG.epoch)
        trigepochsleft{irej2,:}=EEG.epoch(1,irej2).eventtype(1,1); % create vector with the sum of epochs left for each condition (after artifact rejection)
    end
    trigepochsleft=cell2mat(trigepochsleft)'; % convert this vector from cell to number, otherwise the "sum" function won't work. Transpose.
    for j=1:numel(Avg.trig)
        epochsleft(:,j)=sum(trigepochsleft==Avg.trig(j)); % find all triggers belonging to each condition and sum them
        rejepxcond(:,j)=100-(epochsleft(:,j)*100)/totepochs(:,j); % calculate percentage of rejected epochs separately for each condition
    end
    
    % save log file in a structure variable
    logfile(:,isub)=struct('subject',num2str(Avg.subjects(1,isub),'%.2d'),... % participant number
        'interp_chans',gen_bad_chans',... % interpolated channels
        'perc_rejepochs',rejepxcond); % percentage of rejected epochs per condition
    
    % clear variables before starting the preprocessing of a new participant
    clear EEG EEG_block trigtotepochs totepochs trigepochsleft epochsleft rejepxcond
    
end

% save log file
cd([pathexp '\' pathanalysis '\']) % set current directory
save log_preproc.mat logfile % save log file

%%
