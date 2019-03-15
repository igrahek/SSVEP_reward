% [EEG,Trials2Remove]=SCADS_RemoveEyeMovements(EEG,EOGChannelNrs,TimeRange,Threshold)
%
% Finds (and optionally removes) trials with Eye Movements exceeding Threshold microvolot within
% the specified timerange. TimeRange is a vector with two values: the
% start-time and the end-time of the to-be considered timerange. If
% TimeRange is omitted the entire epoch is used for rejection.
% - Use eegF_Bipolarize to generate a bipolar VEOG & HEOG-Channels first
%
% (c) 2007 - S.Andersen
function[EEG,Trials2Remove]=SCADS_RemoveEyeMovements(EEG,EOGChannelNrs,TimeRange,Threshold);
if nargin<2 || nargin>4 || (nargin>2 && ~ismember(numel(TimeRange),[0 2]))
   help SCADS_RemoveEyeMovements;
   return;
end
if nargin<3 || numel(TimeRange)==0
    TimeRange(1)=EEG.xmin;
    TimeRange(2)=EEG.xmax;
end
if nargin<4
    Threshold=25;
end
SPRange=[eegF_Time2Sp(EEG,TimeRange(1)):eegF_Time2Sp(EEG,TimeRange(2))];
MaxMoveDur = 20;    %Time for an eye movement in millisenconds
FiltWidth  =100;    %Width of Median Filter in milliseconds

MaxMoveDur=round(MaxMoveDur*EEG.srate/1000);   %Time for an eye movement in sampling points
FiltWidth=round(0.5*FiltWidth*EEG.srate/1000); %Width of Median Filter in sampling points
% Median filter eye channels
EOG=EEG.data(EOGChannelNrs,:,:);
N=size(EOG,2);
for n=1:N
    EOG_filt(:,n,:)=median(EOG(:,max(1,n-FiltWidth):min(N,n+FiltWidth),:),2);
end
EOG_filt=EOG_filt(:,SPRange,:);
N=size(EOG_filt,2);
% Find and delete trials with eye movements
y=squeeze(max(max(abs(EOG_filt(:,1:N-MaxMoveDur,:)-EOG_filt(:,1+MaxMoveDur:N,:)),[],2),[],1));
if median(y)+10>Threshold
    Threshold=median(y)+10;
    disp(['EOG is very noisy - threshold raised to ' num2str(Threshold) ' microvolt.'])
    disp('Please recheck eye movements manually for this participant!!!');
end
Trials2Remove=find(y>Threshold)';
disp(['Trials with eye movements: ' num2str(Trials2Remove,'%4d')]);
% EEG=pop_selectevent(EEG,'omitepoch',Trials2Remove,'deleteevents','off','deleteepochs','on');