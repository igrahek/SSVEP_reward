% function [gen_bad_chans,EEG_ERP]=eegF_FASTER_chaninterp(cfg,EEG_ERP)
%
% calls FASTER routines (for details, see help file of eegF_FASTER.m)
% ONLY FOR CHANNEL INTERPOLATION
% ONLY FOR FSAReward EXPERIMENT (not tested with any other experiment)
% 
% 2017 - adapted by Schettino (from eegF_FASTER.m function, authored by Craddock & Keitel)
% 
function [gen_bad_chans,EEG_ERP]=eegF_FASTER_chaninterp(cfg,EEG_ERP)

% % do some input checking
if ~isfield(cfg,'thresh'), Thresh=3; else Thresh=cfg.thresh; end
if ~isfield(cfg,'badchan'), BadChans=[]; else BadChans=cfg.badchan; end
if ~isfield(cfg,'datachan'), error('Vector of data channels must be provided as cfg.datachan=1:x!'); else datachan=cfg.datachan; end
if ~isfield(EEG_ERP.chanlocs,'X'), error('Channel locations must be read before running FASTER!'); end

TMP=EEG_ERP;

fprintf('Determine and interpolate generally contaminated channels...\n');
% find vertex electrode...
refc=1;while ~strcmpi(EEG_ERP.chanlocs(refc).labels,'cz');refc=refc+1;end
% and use it as a standard reference
list_props=channel_properties(TMP,datachan,refc);
gen_bad_chans=find(min_z(list_props));
TMP=h_eeg_interp_spl(TMP,unique([gen_bad_chans' BadChans])); % interpolate

EEG_ERP=TMP;
%%
