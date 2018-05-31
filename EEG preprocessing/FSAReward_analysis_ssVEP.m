
%% FSAReward, analysis ssVEP

% add path to EEGLAB + plugins & various functions
addpath(genpath('E:\Experiments\Grahek_Ivan\FSAReward\scripts\EEG\toolboxes_functions\'));

close all; clear all; clc;

% experiment settings
rng(7); % seed for pseudo-random number generator ('A scuppetta!)
expname='FSAReward'; % experiment name
pathexp=['E:\Experiments\Grahek_Ivan\' expname]; % main directory
pathanalysis='analysis\EEG\Exp1'; % where to store the analysis
pathpreproc='preproc'; % where to store the preprocessed single-subject data (for grand average)
pathGA='mean'; % where to store the grand-averages
pathblocks='control_analyses\blocks'; % where to store the preprocessed single-subject data separated into blocks
Avg.channs64=[pathexp '\scripts\EEG\Antonio_BioSemi64.locs']; % channel locations without external channels
Avg.prefix='VP'; % prefix of data files
begin_epoch=.5; % begin epoch (in seconds)
end_epoch=3.25; % end epoch (in seconds)

% Participants:
%  * technical problems (no EEG recording from the beginning of the experiment, problems with battery):
%     - VP05, VP10, VP13, VP27
%  * percentage of rejected epochs (ssVEP preprocessing pipeline) in one or more conditions exceeds 35%:
%     - VP02 [27.5,25,15,17.5,27.5,40]
%     - VP35 [20,15,27.5,32.5,42.5,35]
%     - VP39 [10,7.5,12.5,10,45,40]
%     - VP44 [2.5,7.5,30,37.5,15,22.5]
%  * poor behavioral performance (cut-off criterion: 60% hit rate across all conditions):
%     - VP04, VP06, VP14, VP15, VP17, VP20, VP24, VP25, VP26, VP31, VP34, VP38, VP41, VP45
%  * good behavioral performance, good EEG signal (N = 26):
%     - VP01, VP03, VP07, VP08, VP09, VP11, VP12, VP16, VP18, VP19, VP21, VP22, VP23, VP28, VP29, VP30, VP32, VP33, VP36, VP37, VP40, VP42, VP43, VP46, VP47, VP48

% participants with poor behavioral performance (N=14):
% Avg.subjects=[4 6 14 15 17 20 24 25 26 31 34 38 41 45];
% Avg.pathin=[pathexp '\' pathanalysis '\' pathpreproc '\poor_performance\']; % where to retrieve the preprocessed individual subject data (for grand average)
% Avg.pathout=[pathexp '\' pathanalysis '\' pathGA '\poor_performance\']; % where to save the grand average

% participants with good behavioral performance  (N = 26):
% Avg.subjects=[1 3 7 8 9 11 12 16 18 19 21 22 23 28 29 30 32 33 36 37 40 42 43 46 47 48];
% Avg.pathin=[pathexp '\' pathanalysis '\' pathpreproc '\good_performance\']; % where to retrieve the preprocessed individual subject data (for grand average)
% Avg.pathout=[pathexp '\' pathanalysis '\' pathGA '\good_performance\']; % where to save the grand average

% all participants (N=40):
Avg.subjects=[1 3 4 6 7 8 9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 36 37 38 40 41 42 43 45 46 47 48];
Avg.pathin=[pathexp '\' pathanalysis '\' pathpreproc '\']; % where to retrieve the preprocessed individual subject data (for grand average)
Avg.pathout=[pathexp '\' pathanalysis '\' pathGA '\']; % where to save the grand average

% TRIGGERS (only conditions in which no dot movement occurred):
% 1: baseline, red attended; 2: baseline, blue attended;
% 3: acquisition, red attended; 4: acquisition, blue attended;
% 5: extinction, red attended; 6: extinction, blue attended.
Avg.trig=[1 2 3 4 5 6]; % select triggers

%% grand averages: create & save separate grand averages for each condition

[Trials]=MeanData_FSAReward(Avg);

%% calculate and plot topography and frequency amplitude spectrum

% set variable for topography
topo=struct('numchans',64,... % number of scalp electrodes
    'freq',[10 12],... % frequencies of interest
    'samprate',512,... % sampling rate
    'time',[begin_epoch end_epoch],... % epoch length (in seconds)
    'timeplot',1,... % random time point (within the epoch length) used as reference to extract topography (in seconds)
    'maptype','hot',... % color map (see help topoplot)
    'elec_coords',Avg.channs64,... % electrode coordinates
    'labels','numbers',... % labels to use in plot (numbers or labels)
    'subjects',Avg.subjects,... % number of datasets
    'pathin',Avg.pathout,... % where to take the data
    'amp_topomaplim',[0 1],... % min/max amplitude value in topography (common across frequencies)
    'amp_spectramaplim',[0 1.6]); % min/max amplitude value in spectrum (common across frequencies)

% select participant's group for analysis
plotgroup=3;
switch plotgroup
    case 1 % odd numbers --> red was rewarded
        topo.subjnum=Avg.subjects(mod(Avg.subjects,2)==1); % participants with odd numbers
        topo.condnum=Avg.trig(mod(Avg.trig,2)==1); % conditions with odd numbers
        topo.subjects=find(mod(Avg.subjects,2)==1); % index of selected participants
        topo.conds=find(mod(Avg.trig,2)==1); % index of selected conditions
        eventlabels={'BslnRedAttended' 'AcqRedAttended' 'ExtRedAttended'}; % condition labels
        csvname='rewardRed'; % name of csv file containing amplitude values at specified frequencies
    case 2 % even numbers --> blue was rewarded
        topo.subjnum=Avg.subjects(mod(Avg.subjects,2)==0); % participants with even numbers
        topo.condnum=Avg.trig(mod(Avg.trig,2)==0); % conditions with even numbers
        topo.subjects=find(mod(Avg.subjects,2)==0); % index of selected participants
        topo.conds=find(mod(Avg.trig,2)==0); % index of selected conditions
        eventlabels={'BslnBlueAttended' 'AcqBlueAttended' 'ExtBlueAttended'}; % condition labels
        csvname='rewardBlue'; % name of csv file containing amplitude values at specified frequencies 
    case 3 % all numbers --> all reward and attention conditions
        topo.subjnum=Avg.subjects; % all participants
        topo.condnum=Avg.trig; % all conditions
        topo.subjects=1:numel(Avg.subjects); % index of all participants
        topo.conds=1:numel(Avg.trig); % index of all conditions
        eventlabels={'BslnRedAttended' 'BslnBlueAttended' 'AcqRedAttended' 'AcqBlueAttended' 'ExtRedAttended' 'ExtBlueAttended'}; % condition labels
        csvname='rewardBoth'; % name of csv file containing amplitude values at specified frequencies
    case 4 % subset of participants, selected based on behavioral performance
        topo.subjnum=[1 3 4 6 7 8 9 12 16 19 21 22 23 24 25]; % participants selected based on behavioral performance
        topo.condnum=Avg.trig; % all conditions
        topo.subjects=find(ismember(Avg.subjects,topo.subjnum)); % index of selected participants
        topo.conds=1:numel(Avg.trig); % index of all conditions
        eventlabels={'BslnRedAttended' 'BslnBlueAttended' 'AcqRedAttended' 'AcqBlueAttended' 'ExtRedAttended' 'ExtBlueAttended'}; % condition labels
        csvname='rewardBoth_goodBehavior'; % name of csv file containing amplitude values at specified frequencies
end

% plot topography
[elec_rank]=plot_topo_FSAReward(topo,plotgroup);

% select electrodes for spectral analysis
elec_choice=1;
switch elec_choice % select method for electrode selection
    case 1 % select electrodes with max amplitude for each participant
        elec=findmaxelec(elec_rank(:,:,1),4,0); % select 4 electrodes with highest amplitude
        elec=repmat(elec,numel(topo.subjects),1); % select electrodes
    case 2 % manually select electrodes
        elec=repmat([28],numel(topo.subjects),1);
end

% fields to add to the topo structure (for graph of amplitude spectra)
topo.channels=elec; % selected electrodes
topo.bound=[0 16]; % range of frequencies to plot

% plot frequency amplitude spectrum & extract amplitude values
[amplitudes]=spec_plot_FSAReward(topo,eventlabels,plotgroup);
cell2csv([pathexp '\' pathanalysis '\amplitudes_' csvname '.csv'],amplitudes); % save amplitude values in a .csv file

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FSAReward, blockwise extraction ssVEP amplitudes %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fft_data=[]; fft_amps=[]; % initialize matrices
for isub=1:numel(Avg.subjects) % loop through participants
    
    % load files
    EEG_block01=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block1.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 1
    EEG_block02=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block2.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 2
    EEG_block03=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block3.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 3
    EEG_block04=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block4.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 4
    EEG_block05=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block5.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 5
    EEG_block06=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block6.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 6
    EEG_block07=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block7.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 7
    EEG_block08=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block8.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 8
    EEG_block09=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block9.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 9
    EEG_block10=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block10.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 10
    EEG_block11=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block11.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 11
    EEG_block12=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block12.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % block 12
 
    % detrending
    EEG_block01=eeg_detrend_widmann(EEG_block01); % block 1
    EEG_block02=eeg_detrend_widmann(EEG_block02); % block 2
    EEG_block03=eeg_detrend_widmann(EEG_block03); % block 3
    EEG_block04=eeg_detrend_widmann(EEG_block04); % block 4
    EEG_block05=eeg_detrend_widmann(EEG_block05); % block 5
    EEG_block06=eeg_detrend_widmann(EEG_block06); % block 6
    EEG_block07=eeg_detrend_widmann(EEG_block07); % block 7
    EEG_block08=eeg_detrend_widmann(EEG_block08); % block 8
    EEG_block09=eeg_detrend_widmann(EEG_block09); % block 9
    EEG_block10=eeg_detrend_widmann(EEG_block10); % block 10
    EEG_block11=eeg_detrend_widmann(EEG_block11); % block 11
    EEG_block12=eeg_detrend_widmann(EEG_block12); % block 12
    
    sprange=eegF_Time2Sp(EEG_block01,topo.time(1)):eegF_Time2Sp(EEG_block01,topo.time(2)); % time range of interest (i.e., excluding the initial ERP response)
    totfreqs=((0:(numel(sprange)*2)-1)/(numel(sprange)*2))*topo.samprate; % frequencies from which to extract amplitudes
    
    % merge blocks in pairs, so that EEG data will have the same structure as behavioral data (i.e., 2 parts for each experimental phase)
    EEG_bsln_1stPart=pop_mergeset(EEG_block01,EEG_block02); % first part of baseline
    EEG_bsln_2ndPart=pop_mergeset(EEG_block03,EEG_block04); % second part of baseline
    
    EEG_acq_1stPart=pop_mergeset(EEG_block05,EEG_block06); % first part of acquisition
    EEG_acq_2ndPart=pop_mergeset(EEG_block07,EEG_block08); % second part of acquisition
    
    EEG_ext_1stPart=pop_mergeset(EEG_block09,EEG_block10); % first part of extinction
    EEG_ext_2ndPart=pop_mergeset(EEG_block11,EEG_block12); % second part of extinction
    
    % baseline
    for icond=1:2 % loop through baseline triggers
        % first part
        EEG_bsln_1stPart_icond=pop_select(EEG_bsln_1stPart,'trial',find([EEG_bsln_1stPart.event(:).type]==icond)); % select trials of current condition
        % normalized FFT of EEG signal (averaged across trials and electrodes with maximal amplitude)
        % (zeropadded to double the length of the data, in order to increase frequency resolution and be able to visualize exactly 10 Hz and 12 Hz)
        % matrix dimensions: participants x blocks x conditions x frequencies
        fft_data(isub,1,icond,:)=abs(fft(mean(mean(EEG_bsln_1stPart_icond.data(topo.channels(isub,:),sprange,:),3),1),numel(sprange)*2,2))*2/size(sprange,2);
        for ifreq=1:numel(topo.freq) % loop through frequencies of interest
            [~,spoi]=min(abs(totfreqs-topo.freq(ifreq))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
            % extract amplitude (use 2 dimensions for easy export to .csv)
            temp_fft_amps=fft_data(isub,1,icond,spoi); % extract amplitude
            fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
        end
        % second part
        EEG_bsln_2ndPart_icond=pop_select(EEG_bsln_2ndPart,'trial',find([EEG_bsln_2ndPart.event(:).type]==icond)); % select trials of current condition
        fft_data(isub,2,icond,:)=abs(fft(mean(mean(EEG_bsln_2ndPart_icond.data(topo.channels(isub,:),sprange,:),3),1),numel(sprange)*2,2))*2/size(sprange,2); % normalized FFT of EEG signal (averaged across trials and electrodes with maximal amplitude)
        for ifreq=1:numel(topo.freq) % loop through frequencies of interest
            [~,spoi]=min(abs(totfreqs-topo.freq(ifreq))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
            temp_fft_amps=fft_data(isub,2,icond,spoi); % extract amplitude
            fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
        end
    end

    % acquisition
    for icond=3:4 % loop through acquisition triggers
        % first part
        EEG_acq_1stPart_icond=pop_select(EEG_acq_1stPart,'trial',find([EEG_acq_1stPart.event(:).type]==icond)); % select trials of current condition
        fft_data(isub,3,icond,:)=abs(fft(mean(mean(EEG_acq_1stPart_icond.data(topo.channels(isub,:),sprange,:),3),1),numel(sprange)*2,2))*2/size(sprange,2); % normalized FFT of EEG signal (averaged across trials and electrodes with maximal amplitude)
        for ifreq=1:numel(topo.freq) % loop through frequencies of interest
            [~,spoi]=min(abs(totfreqs-topo.freq(ifreq))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
            temp_fft_amps=fft_data(isub,3,icond,spoi); % extract amplitude
            fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
        end
        % second part
        EEG_acq_2ndPart_icond=pop_select(EEG_acq_2ndPart,'trial',find([EEG_acq_2ndPart.event(:).type]==icond)); % select trials of current condition
        fft_data(isub,4,icond,:)=abs(fft(mean(mean(EEG_acq_2ndPart_icond.data(topo.channels(isub,:),sprange,:),3),1),numel(sprange)*2,2))*2/size(sprange,2); % normalized FFT of EEG signal (averaged across trials and electrodes with maximal amplitude)
        for ifreq=1:numel(topo.freq) % loop through frequencies of interest
            [~,spoi]=min(abs(totfreqs-topo.freq(ifreq))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
            temp_fft_amps=fft_data(isub,4,icond,spoi); % extract amplitude
            fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
        end
    end
    
    % extinction
    for icond=5:6 % loop through extinction triggers
        % first part
        EEG_ext_1stPart_icond=pop_select(EEG_ext_1stPart,'trial',find([EEG_ext_1stPart.event(:).type]==icond)); % select trials of current condition
        fft_data(isub,5,icond,:)=abs(fft(mean(mean(EEG_ext_1stPart_icond.data(topo.channels(isub,:),sprange,:),3),1),numel(sprange)*2,2))*2/size(sprange,2); % normalized FFT of EEG signal (averaged across trials and electrodes with maximal amplitude)
        for ifreq=1:numel(topo.freq) % loop through frequencies of interest
            [~,spoi]=min(abs(totfreqs-topo.freq(ifreq))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
            temp_fft_amps=fft_data(isub,5,icond,spoi); % extract amplitude
            fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
        end
        % second part
        EEG_ext_2ndPart_icond=pop_select(EEG_ext_2ndPart,'trial',find([EEG_ext_2ndPart.event(:).type]==icond)); % select trials of current condition
        fft_data(isub,6,icond,:)=abs(fft(mean(mean(EEG_ext_2ndPart_icond.data(topo.channels(isub,:),sprange,:),3),1),numel(sprange)*2,2))*2/size(sprange,2); % normalized FFT of EEG signal (averaged across trials and electrodes with maximal amplitude)
        for ifreq=1:numel(topo.freq) % loop through frequencies of interest
            [~,spoi]=min(abs(totfreqs-topo.freq(ifreq))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
            temp_fft_amps=fft_data(isub,6,icond,spoi); % extract amplitude
            fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
        end
    end
end

% save data to .csv
tot_fft_amps=[num2cell(repelem(Avg.subjects,numel(eventlabels)*numel(topo.freq)*2))',... % participants (repeat for number of conditions x number of frequencies x number of parts per experimental phase)
    repmat(repelem(eventlabels,4),1,numel(Avg.subjects))',... % conditions
    repmat([repelem({'1st'},numel(topo.freq)) repelem({'2nd'},numel(topo.freq))],1,numel(fft_amps)/4)',... % part of experiment phase
    num2cell(repmat(topo.freq,1,numel(fft_amps)/numel(topo.freq))'),... % frequencies
    num2cell(fft_amps)]; % amplitudes
tot_fft_amps_header={'participant','condition','part_expPhase','frequency','amplitude'}; % header of big matrix
tot_fft_amps=[tot_fft_amps_header;tot_fft_amps]; % concatenate header and data
cell2csv([pathexp '\' pathanalysis '\' pathblocks '\amplitudes_part_expPhase.csv'],tot_fft_amps); % save amplitude values in a .csv file

%%

% nblocks=12; % number of trial blocks
% 
% fft_data=[]; fft_amps=[]; % initialize matrices
% for isub=1:numel(Avg.subjects) % loop through participants
%     for iblockcount=1:nblocks % loop through blocks
%         EEG_block=pop_loadset('filename',[Avg.prefix num2str(Avg.subjects(1,isub),'%.2d') '_block' num2str(iblockcount) '.set'],'filepath',[pathexp '\' pathanalysis '\' pathblocks '\']); % load file
%         EEG_block=eeg_detrend_widmann(EEG_block); % detrend
%         sprange=eegF_Time2Sp(EEG_block,topo.time(1)):eegF_Time2Sp(EEG_block,topo.time(2)); % time range of interest (i.e., excluding the initial ERP response)
%         totfreqs=((0:(numel(sprange)*2)-1)/(numel(sprange)*2))*topo.samprate; % frequencies from which to extract amplitudes
%         if iblockcount<5
%             i_cond=1:2;
%         elseif all(iblockcount > 4 & iblockcount < 9)
%             i_cond=3:4;
%         elseif iblockcount>8
%             i_cond=5:6;
%         end
%         for icond=1:numel(i_cond) % loop through conditions
%             trials_icond=find([EEG_block.event(:).type]==i_cond(icond)); % select trials of current condition
%             EEG_block2=pop_select(EEG_block,'trial',trials_icond); % ... select trials
%             if numel(size(EEG_block2.data))==3 % if the data contains at least one trial
%                 % normalized FFT of EEG signal (averaged across trials and electrodes with maximal amplitude)
%                 % (zeropadded to double the length of the data, in order to increase frequency resolution and be able to visualize exactly 10 Hz and 12 Hz)
%                 % matrix dimensions: participants x blocks x conditions x frequencies
%                 fft_data(isub,iblockcount,i_cond(icond),:)=abs(fft(mean(mean(EEG_block2.data(topo.channels(isub,:),sprange,:),3),1),numel(sprange)*2,2))*2/size(sprange,2);
%                 for ifreq=1:numel(topo.freq) % loop through frequencies of interest
%                     [~,spoi]=min(abs(totfreqs-topo.freq(ifreq))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
%                     % fft_amps(isub,iblockcount,icond,ifreq)=fft_data(isub,iblockcount,icond,spoi); % extract amplitude
%                     % instead of saving results as matrix, only use 2 dimensions (to save as .csv)
%                     temp_fft_amps=fft_data(isub,iblockcount,i_cond(icond),spoi); % extract amplitude
%                     fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
%                 end
%             else
%                 for ifreq=1:numel(topo.freq) % loop through frequencies of interest
%                     temp_fft_amps=NaN;
%                     fft_amps=[fft_amps;temp_fft_amps]; % concatenate in big matrix
%                 end
%             end
%         end
%     end
% end
% 
% % save data to .csv
% tot_fft_amps=[num2cell(repelem(Avg.subjects,1,nblocks*(numel(i_cond)*numel(topo.freq)))'),... % participants
%     num2cell(repmat(repelem(1:nblocks,1,numel(i_cond)*numel(topo.freq))',numel(Avg.subjects),1)),... % blocks
%     repmat([repmat(repelem(eventlabels(1:2),1,2),1,4) repmat(repelem(eventlabels(3:4),1,2),1,4) repmat(repelem(eventlabels(5:6),1,2),1,4)]',numel(Avg.subjects),1),... % conditions
%     num2cell(repmat(topo.freq,1,numel(fft_amps)/2)'),... % frequencies
%     num2cell(fft_amps)]; % amplitudes
% tot_fft_amps_header={'participant','block','condition','frequency','amplitude'}; % header of big matrix
% tot_fft_amps=[tot_fft_amps_header;tot_fft_amps]; % concatenate header and data
% cell2csv([pathexp '\' pathanalysis '\' pathblocks '\amplitudes_blocks.csv'],tot_fft_amps); % save amplitude values in a .csv file

%%
