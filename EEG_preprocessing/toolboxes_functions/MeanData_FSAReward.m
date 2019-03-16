
%% Generate one mean file for each condition containing all subjects

function [Trials] = MeanData_FSAReward(Avg)

disp('---------------------------');
disp('Averaging...');
disp('---------------------------');

% Calculate Mean
tempfiles = dir([Avg.pathssj '*.set']);

for iSub = 1:numel(tempfiles)
    FileName = tempfiles(iSub).name;
    disp(FileName);
    
    % Load Data
    EEG = pop_loadset(FileName, Avg.pathssj);
%     EEG = eeg_detrend_widmann(EEG); % detrending
    
    if iSub == 1
        NbChannels = EEG.nbchan;
        NbSamples = EEG.pnts;
        AllSub = zeros(numel(Avg.trig), NbChannels, NbSamples, numel(Avg.subjects));
    end
    
    for iCond = 1:numel(Avg.trig)
        EEG2 = pop_selectevent(EEG, 'type', Avg.trig(iCond), 'deleteevents', 'off', 'deleteepochs', 'on');
        AllSub(iCond, :, :, iSub) = mean(EEG2.data, 3);
        TrialsAveraged(iCond, iSub) = size(EEG2.data, 3);
    end
    clear EEG2
end

clear EEG

% Store Mean in Different Files
EEG = pop_loadset(FileName, Avg.pathssj);

EEG.trials = numel(Avg.subjects);
EEG.event = EEG.event(1:numel(Avg.subjects));
EEG.epoch = EEG.epoch(1:numel(Avg.subjects));

for iCond = 1:numel(Avg.trig);
    EEG.data = squeeze(AllSub(iCond, :, :, :));
    EEG.setname = ['Mean-C' int2str(iCond)];
    pop_saveset(EEG, [EEG.setname '.set'], Avg.pathGM); % save
end

% display number of averaged epochs for each condition
Trials = TrialsAveraged';

end
