%[EEG,Trials2Remove]=SCADS_RemoveBlinks(EEG,VEOGChannelNr,TimeRange)
%
% Finds (and optionally removes) trials with blinks within the specified timerange. Narrow spikes
% from electrical artifacts are ignored (hopefully). TimeRange is a vector
% with two values: the start-time and the end-time of the to-be considered
% time-range. If TimeRange is omitted the entire epoch is used for rejection.

% - Use eegF_Bipolarize to generate a bipolar VEOG-Channel first
% - Results do not depend on whether eegF_Detrend was used prior to
%   calling this function
%
% (c) 2004 - Gruber, Hassler & 2005,2006,2007,2008 - S.Andersen
function[EEG,Trials2Remove]=SCADS_RemoveBlinks(EEG,VEOGChannelNr,TimeRange)
if nargin<2 || nargin>3 || (nargin>2 && ~ismember(numel(TimeRange),[0 2]))
   help(mfilename)
   return
end
if nargin<3 || numel(TimeRange)==0, TimeRange=[EEG.xmin EEG.xmax]; end
MinWidth=20;                % Minimum Width of a blink
Lambda=20;                  % smaller Lambda = more strict control
Gamma=round(MinWidth*EEG.srate/1000);
SPRange=eegF_Time2Sp(EEG,TimeRange(1)):eegF_Time2Sp(EEG,TimeRange(2));
%define threshold
VEOGChannelAmpl=abs(detrend(squeeze(EEG.data(VEOGChannelNr,SPRange,:))));
VEOGChannelAmpl=sort(VEOGChannelAmpl,1,'descend');
Max_VEOGChannelAmpl=VEOGChannelAmpl(Gamma,:);
Threshold=median(Max_VEOGChannelAmpl) + Lambda*median(abs(Max_VEOGChannelAmpl-median(Max_VEOGChannelAmpl)));
%remove trials
Trials2Remove=find(Max_VEOGChannelAmpl>Threshold);
disp(['Trials with blinks: ' num2str(Trials2Remove,'%4d')]);
% EEG=pop_selectevent(EEG,'omitepoch',Trials2Remove,'deleteevents','off','deleteepochs','on');