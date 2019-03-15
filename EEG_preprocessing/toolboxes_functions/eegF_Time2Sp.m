% [SamplingPoint] = eegF_Time2Sp(EEG,Time);
%
% Returns the sampling point(s) corresponding to the timepoint(s) given in
% Time (in seconds)
%
% (c) 2004,2006,2007,2008 - S.Andersen
function[SamplingPoint]=eegF_Time2Sp(EEG,Time)
if nargin ~=2, help(mfilename), return, end
Time(Time<EEG.xmin)=EEG.xmin;
Time(Time>EEG.xmax)=EEG.xmax;
SamplingPoint=floor((Time-EEG.xmin)/(EEG.xmax-EEG.xmin)*(EEG.pnts-1))+1;