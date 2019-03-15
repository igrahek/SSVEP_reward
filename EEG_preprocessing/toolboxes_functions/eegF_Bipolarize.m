% EEG=eegF_Bipolarize(EEG)
%
% Assumes that:
%   HEOG = EXG1 - EXG2
%   VEOG = EXG3 - EXG4
%
% (c) 2005,2008 - S.Andersen
function[EEG]=eegF_Bipolarize(EEG)
disp('Replacing: VEOG = EXG1-EXG2, HEOG = EXG3-EXG4');
ChannelNames=['EXG1';'EXG2';'EXG3';'EXG4'];
ChannelPos=[];

for i=1:size(EEG.chanlocs,2)
    for j=1:size(ChannelNames,1)
        if strcmp(EEG.chanlocs(i).labels,ChannelNames(j,:))
            ChannelPos(j)=i;
        end
    end
end

if size(ChannelPos,2)==size(ChannelNames,1)
    EEG.chanlocs(ChannelPos(1)).labels='VEOG';
    EEG.data(ChannelPos(1),:,:)=EEG.data(ChannelPos(1),:,:)-EEG.data(ChannelPos(2),:,:);
    EEG.chanlocs(ChannelPos(3)).labels='HEOG';
    EEG.data(ChannelPos(3),:,:)=EEG.data(ChannelPos(3),:,:)-EEG.data(ChannelPos(4),:,:);
    EEG=pop_select(EEG,'nochannel',[ChannelPos(2) ChannelPos(4)]);
else
    error('Failed to find correct channels for subtraction !')
end