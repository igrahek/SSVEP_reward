
%% FSAReward, preprocessing for analysis ssVEP

% add path to EEGLAB + plugins & various functions
addpath(genpath('E:\Experiments\Grahek_Ivan\FSAReward\scripts\EEG\toolboxes_functions\'));

close all; clear all; clc;

% experiment settings
rng(7); % seed for pseudo-random number generator ('A scuppetta!)
expname='FSAReward'; % experiment name
pathexp=['E:\Experiments\Grahek_Ivan\' expname]; % main directory
pathdata='data\EEG'; % where to get the .bdf files
pathanalysis='analysis\EEG'; % where to store the analysis
pathpreproc='preproc'; % where to store the preprocessed single-subject data (for grand average)
Avg.channs=[pathexp '\scripts\EEG\Antonio_BioSemi70.locs']; % channel locations
Avg.prefix='VP'; % prefix of data files
begin_epoch=0; % begin epoch (in seconds)
end_epoch=3.25; % end epoch (in seconds)
nblocks=12; % number of trial blocks

% discarded participants due to  technical problems (no EEG recording from the beginning of the experiment, problems with battery):
% - VP05, VP10, VP13, VP27
% final sample
Avg.subjects=[1 2 3 4 6 7 8 9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48];

% TRIGGERS (only conditions in which no dot movement occurred):
% 1: baseline, red attended; 2: baseline, blue attended;
% 3: acquisition, red attended; 4: acquisition, blue attended;
% 5: extinction, red attended; 6: extinction, blue attended.
Avg.trig=[1 2 3 4 5 6]; % select triggers

%% PREPROCESSING

for isub=1:numel(Avg.subjects) % index of participant number
    
    % display participant number
    disp(['Processing ' Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '...'])
    
    % import .bdf files
    if isub>=33 % in the latest batch of participants (from VP37 to VP48), we connected 64 electrodes but mistakenly recorded 128 channels (the signal from the second electrode bundle is flat). Here we fix the mistake by selecting only the right electrodes.
        EEG=pop_biosig([pathexp '\' pathdata '\' Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '.bdf'],'ref',[133 134] ,'refoptions',{'keepref' 'on'}); % referenced to average mastoids
        if isub>=37 % due to a mistake, VP41 to VP48 have 7 exceeding external electrodes
            EEG=pop_select(EEG,'nochannel',[65:128 135:143]); % remove channels with no signal
        else  EEG=pop_select(EEG,'nochannel',[65:128 135:136]); % remove channels with no signal
        end
    else
        EEG=pop_biosig([pathexp '\' pathdata '\' Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '.bdf'],'channels',[1:70],'ref',[69 70],'refoptions',{'keepref' 'on'}); % referenced to average mastoids
    end
    
    %     pop_eegplot(EEG,1,1,1); % check data
    
    EEG=pop_chanedit(EEG,'load',{Avg.channs,'filetype','autodetect'}); % assign channel locations
    %     topoplot([],EEG.chanlocs,'style','blank', 'electrodes','labelpoint','chaninfo',EEG.chaninfo); % see channel locations
    EEG=pop_rmbase(EEG,'Warning','off'); % remove DC
    
    %     pop_eegplot(EEG,1,1,1); % check data
    
    EEG=pop_epoch(EEG,num2cell(Avg.trig),[begin_epoch end_epoch]); % epoching
    
    % for log file
    for irej=1:numel(EEG.epoch)
        trigtotepochs{irej,:}=EEG.epoch(1,irej).eventtype(1,1); % create vector with the sum of total epochs for each condition
    end
    trigtotepochs=cell2mat(trigtotepochs)'; % convert this vector from cell to number, otherwise the "sum" function won't work. Transpose.
    for k=1:numel(Avg.trig)
        totepochs(:,k)=sum(trigtotepochs==Avg.trig(k)); % find all triggers belonging to each condition and sum them
    end
    
    % channel interpolation & artifact rejection using FASTER
    cfg=[]; % create Fieldtrip-like structure
    cfg.datachan=1:64; % select scalp electrodes
    cfg.thresh=[3 3 0 3 3 12]; % see help eegF_FASTER for a description of each number (lower numbers are more conservative)
    [gen_bad_chans,EEG,trials2removeFASTER]=eegF_FASTER(cfg,EEG); % find & interpolate noisy channels; find trials contaminated by artifacts (does not eliminate them!)
    EEG=pop_reref(EEG,[69 70],'keepref','off'); % re-reference to mastoids (the eegF_FASTER function re-references to Cz, hence the need for re-referencing here) and discard them afterwards
    
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
