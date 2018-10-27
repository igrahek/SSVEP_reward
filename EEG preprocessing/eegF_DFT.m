% [EEG]=eegF_DFT(EEG,Frequency,TimeRange,Output)
%
% Calculates a digital fourier tronsformation at the specified frequency
% for a certain timerange. TimeRange is a vector with two values: the
% start-time and the end-time of the to-be considered time-range.
%   Output: 'amplitude' (default) or 'complex'
%
% (c) 2005,2007,2008 - S.Andersen
function[EEG]=eegF_DFT(EEG,Frequency,TimeRange,Output)
if nargin<2 || nargin>4 || (nargin>2 && ~ismember(numel(TimeRange),[0 2]))
    help(mfilename)
    return
end
if nargin<3 || numel(TimeRange)==0, TimeRange=[EEG.xmin EEG.xmax]; end
if nargin<4, Output='amplitude'; end
disp('Calculating Fourier Transformation...');
SPRange=eegF_Time2Sp(EEG,TimeRange(1)):eegF_Time2Sp(EEG,TimeRange(2));
x=exp(-i*(0:(size(SPRange,2)-1))*2*pi*Frequency/EEG.srate);%Complexe Schwingung
x=x*2/size(SPRange,2); %Normalisierung
x=repmat(x,[size(EEG.data,1) 1 1]);
for Trial=1:size(EEG.data,3)
    xM=x.*EEG.data(:,SPRange,Trial); %e-Funktion mit Daten multiplizieren
    xM=repmat(sum(xM,2),[1 numel(SPRange) 1]);
    switch lower(Output)
        case 'amplitude'
            EEG.data(:,SPRange,Trial)=abs(xM); %Betrag in EEG.Data zurückschreiben
        case 'complex'
            EEG.data(:,SPRange,Trial)=xM;
        otherwise
            error([mfilename ': Unknown value ''' Output ''' for ''output''.'])
    end
end