
%% FSAReward: preprocessing behavioral data
%       Ivan Grahek & Antonio Schettino

clear all; clc;
addpath(genpath('E:\Experiments\Grahek_Ivan\FSAReward\scripts\EEG\toolboxes_functions\functions_AS\')); % add path to various functions
%   pathdata='C:\Users\igrahek\Google Drive\Work computer\FSAReward\data\behavior\Exp1\'; % directory with behavioral data
%   pathdest='C:\Users\igrahek\Google Drive\Work computer\FSAReward\analysis\behavior\Exp1\'; % where to save individual csv files  
pathdata='E:\Experiments\Grahek_Ivan\FSAReward\data\behavior\Exp1\'; % directory with behavioral data
pathdest='E:\Experiments\Grahek_Ivan\FSAReward\analysis\Behavior\Exp1\'; % where to save individual csv files  
files=dir([pathdata '*.mat']); % list of all data files in the directory
cd(pathdata) % go to data directory
datatot=[]; % preallocate matrix containing datasets of all participants

for i=1:numel(files) % cycle through participants
    load(files(i,1).name) % load participant's data
    data=dataset({zeros(720,7) 'ParticipantNo', 'Trial', 'AttendedColor', 'Reward','EventType', 'Response', 'RT'}); % initialize data matrix (number of events: 720)
    NumberOfTrials=600; % total number of trials
    TempNoEvents=0; % event counter for current event
    for TrialBeingAnalyzed=1:NumberOfTrials % cycle through trials
        if numel(VPR.Trial(TrialBeingAnalyzed).EvType)~=0 % check if there is an event, otherwise go to the next trial
            for EventNr=1:numel(VPR.Trial(TrialBeingAnalyzed).EvType) % cycle through every event in the trial
                data.Trial(TempNoEvents+EventNr)=TrialBeingAnalyzed; % record trial number
                data.AttendedColor(TempNoEvents+EventNr)=VPR.TrgEvPerCond{VPR.Trial(TrialBeingAnalyzed).Condition}; % record attended color
                data.Reward(TempNoEvents+EventNr)=VPR.Trial(TrialBeingAnalyzed).EvRewd(EventNr); % record reward info
                data.EventType(TempNoEvents+EventNr)=VPR.Trial(TrialBeingAnalyzed).EvType(EventNr); % record event type info (which color moved?)
                data.ParticipantNo(TempNoEvents+EventNr)=S.SubjNr; % record participant number
                if numel(VPR.Trial(TrialBeingAnalyzed).ReType)~=0 % check if there are responses in this trial
                    for ResponseNr=1:size(VPR.Trial(TrialBeingAnalyzed).ReType,2) % cycle through every response in the trial
                        if  VPR.Trial(TrialBeingAnalyzed).EvTime(EventNr)+VPR.MinDelay<=VPR.Trial(TrialBeingAnalyzed).ReTime(ResponseNr) && ... % reaction time is above minimal delay - This makes sure that we're analyzing a response to a particular event
                                VPR.Trial(TrialBeingAnalyzed).EvTime(EventNr)+VPR.MaxDelay>=VPR.Trial(TrialBeingAnalyzed).ReTime(ResponseNr) && ... % reaction time is below maximal delay - This makes sure that we're analyzing a response to a particular event
                                data.Response(TempNoEvents+EventNr)==0 % there was no previous response on the same event
                            if VPR.TrgEvPerCond{VPR.Trial(TrialBeingAnalyzed).Condition}==VPR.Trial(TrialBeingAnalyzed).EvType(EventNr) % if the to-be-attended color of dots equals the color of dots which moved in the event...
                                data.RT(TempNoEvents+EventNr)=VPR.Trial(TrialBeingAnalyzed).ReTime(ResponseNr)-VPR.Trial(TrialBeingAnalyzed).EvTime(EventNr); % ... calculate the RT and write it down in the appropriate row
                                data.Response(TempNoEvents+EventNr)=1; % record response as hit (1): the attended color moved and participants responded within the time window
                            else  % if the to-be-attended color of dots does NOT equal the color of dots which moved in the event...
                                data.RT(TempNoEvents+EventNr)=VPR.Trial(TrialBeingAnalyzed).ReTime(ResponseNr)-VPR.Trial(TrialBeingAnalyzed).EvTime(EventNr); % calculate RT (originally NaN)
                                data.Response(TempNoEvents+EventNr)=2; % record response as false alarm (2): the attended color did not move but participants responded
                            end
                        else
                            data.RT(TempNoEvents+EventNr)=NaN; % mark RT as NaN
                            data.Response(TempNoEvents+EventNr)=0; % record response as miss (0): the attended color moved, participants responded but outside the time window
                        end
                    end
                elseif VPR.TrgEvPerCond{VPR.Trial(TrialBeingAnalyzed).Condition}==VPR.Trial(TrialBeingAnalyzed).EvType(EventNr) % if the to-be-attended color of dots equals the color of dots which moved in the event...
                    data.RT(TempNoEvents+EventNr)=NaN; % mark RT as NaN
                    data.Response(TempNoEvents+EventNr)=0; % record response as miss (0): the attended color moved but participants did not respond
                else
                    data.RT(TempNoEvents+EventNr)=NaN; % mark RT as NaN
                    data.Response(TempNoEvents+EventNr)=3; % record response as correct rejection (3): the attended color did not move and participants correctly withheld the response
                end
            end
        end
        TempNoEvents=TempNoEvents+numel(VPR.Trial(TrialBeingAnalyzed).EvType);
    end
    cell2csv([pathdest 'ssj' num2str(S.SubjNr) '.csv'],dataset2cell(data)); % save individual data in a .csv file
    datatot=[datatot;data]; % save all participants in one big dataset
end

cell2csv([pathdest 'allssj.csv'],dataset2cell(datatot)); % save all data in a .csv file

%%
