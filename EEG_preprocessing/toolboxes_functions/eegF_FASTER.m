% function [gen_bad_chans,EEG,trials2remove,comps2remove]=eegF_FASTER(cfg,EEG)
%
% calls FASTER routines
%
% example call:
% cfg=[];
% cfg.datachan=1:64;
% [EEG,trials2remove,comps2remove]=eegF_FASTER(cfg,EEG)
%
% REQUIREMENTS:
%       1. EEG.data has to be a 3D-matrix with dimorder: channels x points x epochs
%          -> data has to be epoched
%       2. Epochs should be zero-mean (either by mean subtraction or highpass filter)
%       3. Channel locations must be included in EEG (see EEG.chanlocs)
%
% NOTES:
%       This function does not actually remove epochs from the data set
%       (although verbose might suggest otherwise). It yields a vector of
%       suspicious epochs that can be used for subsequent removal,
%       e.g. executing: 
%           EEG=pop_select(EEG,'notrial',trials2remove);
%       after FASTER has been run.
%       Same applies to independent components (step 3):
%           EEG=pop_subcomp(EEG,comps2remove,0);
% 
% INPUT:
%       cfg.(...)   - fieldtrip-like configuration structure containing
%                     the following fields:
%       mandatory fields:
%       .datachan   - vector of scalp data channel indices,
%                     for 64-electrode recordings: cfg.datachan=1:64;
%       optional:
%       .timerange  - select time range within epochs,
%                     provide as two-element vector in seconds (e.g. [-1 2]),
%                     if omitted FASTER considers entire epochs
%       .thresh     - determine FASTER rejection criteria,
%                     provide as six-element vector (e.g. [3 3 3 3 3 12]).
%                     First five numbers are z-scores, one for each step
%                     of the faster algorithm. Higher scores are more liberal,
%                     lower scores more conservative. The last number determines
%                     the maximal number of electrodes that may be interpolated
%                     in a given epoch before it is marked for rejection.
%                     If omitted, default values of
%                     cfg.thresh=[3 3 3 3 3 round(n_electrodes/6)] are used.
%                     Artifacts are detected and corrected in five aspects
%                     of the EEG data (see paper for details):
%                     1) CHANNELS: channels with variance, mean correlation,
%                        and Hurst exponent larger than the specified Z-value
%                        (e.g., 3) are interpolated;
%                     2) EPOCHS: mean across channels is computed for each epoch.
%                        If amplitude range, variance, and channel deviation 
%                        exceed the specified Z-value (e.g., 3), the epoch is removed;
%                     3) INDEPENDENT COMPONENTS (ICs): this step requires
%                        that an ICA has been run on the dataset and is
%                        automatically skipped if no ICA weights are found
%                        in the input EEG structure. If you have run an ICA
%                        but nonetheless want FASTER to skip this step,
%                        enter a 0 at the third position in cfg.thresh (e.g., [3 3 0 3 3 12]).
%                        If correlation with EOG chans, spatial kurtosis,
%                        slope in filter band, Hurst exponent, and median gradients
%                        are larger than the specified Z-value (e.g., 3),
%                        the component is subtracted;
%                     4) SINGLE-CHANNEL, SINGLE-EPOCH ARTIFACTS: within each epoch,
%                        channels with variance, median gradient, amplitude range,
%                        and channel deviation larger than the specified Z-value
%                        (e.g., 3) are interpolated; 
%                     5) GRAND-AVERAGE: grand-averages with amplitude range,
%                        variance, channel deviation, and maximum EOG value
%                        larger than the specified Z-value (e.g., 3) are removed.
%       .eyechan    - you may specify indices of eye channels to facilitate
%                     the detection of blink components in step 3, e.g.,
%                     cfg.eyechan=[65 66 67 68].
%                     This requires ICA weights (see above).          
%       .badchan    - provide list of known bad channels,
%                     this option can be used to make sure that FASTER
%                     interpolates generally contaminated channels or
%                     flatliners noted during data recordings.
%                     Provide as an n-element vector (e.g. [3 17 38]).
%
%       EEG         - common EEGLAB data struct containing epoched EEG.data
%     
% OUTPUT:
%       EEG         - EEGLAB data struct with the same number of epochs
%                     as the input dataset yet interpolated where necessary
%                     and containing marks for the to-be-rejected trials in
%                     EEG.reject.rejmanual
%       trials2remove - list of epoch indices that have been marked for rejection
%       comps2remove  - list of components that have been marked for subtraction
%                       (is empty if step 3 is omitted)
%
% based on (cite when using):
%       Nolan, Whelan & Reilly (2010) FASTER: Fully Automated Statistical
%       Thresholding for EEG artifact Rejection. J Neurosci Methods
%
% download the FASTER toolbox from:
%       http://www.mee.tcd.ie/neuraleng/Research/Faster
% and put it in the plugin folder of your EEGLAB.
%
% (c) 2012 - Craddock & Keitel
% 2013 - Updated by Craddock for Schettino
% 2015 - Schettino updates help

function [gen_bad_chans,EEG,trials2remove,comps2remove]=eegF_FASTER(cfg,EEG)

% % do some input checking
if nargin~=2, help(mfilename); return; end
if any([~isstruct(cfg) ~isstruct(EEG)]), help(mfilename); return; end
if ~isfield(cfg,'thresh'), Thresh=[3 3 3 3 3 floor(EEG.nbchan./6)]; else Thresh=cfg.thresh; end
if ~isfield(cfg,'timerange'), TimeRange=[]; else TimeRange=cfg.timerange; end
if ~isfield(cfg,'badchan'), BadChans=[]; else BadChans=cfg.badchan; end
if ~isfield(cfg,'eyechan'), EyeChans=[]; else EyeChans=cfg.eyechan; end
if ~isfield(cfg,'datachan')
    error('Vector of data channels must be provided as cfg.datachan=1:x!');
else
    datachan=cfg.datachan;
end
if ~isfield(EEG.chanlocs,'X'), error('Channel locations must be read before running FASTER!'); end

if ~ismember(numel(TimeRange),[0 2]), help(mfilename); return; end
if numel(Thresh)~=6, help(mfilename); return; end

fprintf('### running FASTER\n\n')

% select time-range data if requested
if ~isempty(TimeRange)
    TMP=pop_select(EEG,'time',TimeRange);
else
    TMP=EEG;
end

fprintf('STEP 1: determine and interpolate generally contaminated channels...\n')
% find vertex electrode...
refc=1;while ~strcmpi(EEG.chanlocs(refc).labels,'cz');refc=refc+1;end
% and use it as a standard reference
list_props=channel_properties(TMP,datachan,refc);
gen_bad_chans=find(min_z(list_props,prep_rej_opt(list_props,Thresh(1))));
TMP=h_eeg_interp_spl(TMP,unique([gen_bad_chans' BadChans]));

% determine contaminated epochs and
% store indices of these for subsequent removal
fprintf('STEP 2: determine and remove contaminated epochs...\n')
list_props=epoch_properties(TMP,datachan);
trials2remove=find(min_z(list_props,prep_rej_opt(list_props,Thresh(2))));

% trial bookkeeping
[TMP,trials]=tbookkeep(TMP,true(1,TMP.trials),trials2remove);

fprintf('STEP 3: determine and subtract artifact components...\n')
if ~isempty(TMP.icaweights) && Thresh(3)
    list_props=component_properties(TMP,EyeChans);
    comps2remove=find(min_z(list_props,prep_rej_opt(list_props,Thresh(3))));
    fprintf('Found %d suspicious component(s).\n',numel(comps2remove));
else
    comps2remove=[];
    if ~isfield(TMP,'icaweights')
        fprintf('No ICA performed on this dataset. Proceeding...\n');
    else
        fprintf('Skipping step as requested. Proceeding...\n');
    end
end

% determine contaminated channels per epoch
fprintf('STEP 4: determine and interpolate contaminated channels per epoch...\n')
for iepoch=1:TMP.trials
    list_props=single_epoch_channel_properties(TMP,iepoch,datachan);
    loc_bad_chans{iepoch}=find(min_z(list_props,prep_rej_opt(list_props,Thresh(4))));
    chan_thresh(iepoch,1)=numel(loc_bad_chans{iepoch})+numel(gen_bad_chans);
end
% interpolate
TMP=h_epoch_interp_spl(TMP,loc_bad_chans);
% extra step: reject epochs that exceed the limit of contaminated channels
trials2remove=find(chan_thresh>Thresh(6));

% trial bookkeeping
[TMP,trials]=tbookkeep(TMP,trials,trials2remove);

% determine epochs that contaminate the GA and
% store indices of these for subsequent removal
fprintf('STEP 5: determine and remove epochs based on GA properties...\n')
list_props=GA_properties(TMP,datachan,[]);
trials2remove=find(min_z(list_props,prep_rej_opt(list_props,Thresh(5))));

% trial bookkeeping
[TMP,trials]=tbookkeep(TMP,trials,trials2remove);

% reintegrate interpolated data in original EEG dataset
datapoints=time2sp(EEG,TMP.xmin):time2sp(EEG,TMP.xmax);
if numel(datapoints)>size(TMP.data,2)
    datapoints=datapoints(1:end-1);
end
EEG.data(:,datapoints,trials)=TMP.data;

% get final vector of trials 2 remove
trials2remove=find(~trials);
% ... and store it in the output EEG struct as well
EEG.reject.rejmanual=~trials;

% indicate that FASTER has been run on this dataset
EEG.reject.artifact_rej_method='FASTER';

%% subfunctions
function sp=time2sp(EEG,toi)
% convert times to sampling points
toi(toi<EEG.xmin)=EEG.xmin; toi(toi>EEG.xmax)=EEG.xmax;
sp=floor((toi-EEG.xmin)/(EEG.xmax-EEG.xmin)*(EEG.pnts-1))+1;

function [TMP,trials]=tbookkeep(TMP,trials,trials2remove)
% trial bookkeeping
subtrials=true(1,TMP.trials);
TMP=pop_select(TMP,'notrial',trials2remove);
subtrials(trials2remove)=false;
trials(trials)=subtrials;

function rej_opt=prep_rej_opt(list_props,thresh)
% prepare input to min_z function
% !assumes that evaluations within one step use the same threshold!
rej_opt.measure=ones(1,size(list_props,2));
rej_opt.z=thresh*ones(1,size(list_props,2));
