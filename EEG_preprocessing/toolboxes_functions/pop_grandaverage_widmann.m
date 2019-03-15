% pop_grandaverage - Grand average epoched EEGLAB EEG datasets
%
% Usage:
%   >> [EEG, com] = pop_grandaverage(); % pop-up window mode
%   >> [EEG, com] = pop_grandaverage({'dataset1.set' 'dataset2.set' ...
%                       'datasetn.set'}, 'key1', value1, 'key2', ...
%                       value2, 'keyn', valuen);
%   >> [EEG, com] = pop_grandaverage(ALLEEG); % pop-up window mode
%   >> [EEG, com] = pop_grandaverage(ALLEEG, 'key1', value1, 'key2', ...
%                       value2, 'keyn', valuen);
%
% Inputs:
%   ALLEEG        - vector of EEGLAB EEG structures OR cell array of
%                   strings with dataset filenames
%   'datasets'    - vector datasets to average
%
% Optional inputs:
%   'pathname'    - char array path name {default '.'}
%   'eventtype'   - new event type {default: first time locking event of
%                   each dataset}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%   com           - history string
%
% Note:
%   Use with separate datasets for each subject and event type. All trials
%   are averaged per dataset and the averages are collected in the third
%   (trial) dimension of the output dataset. The event type of each average
%   is defined by the first time locking event. In file mode all datasets
%   should reside in a single directory. Electrode locations are NOT
%   averaged. Electrode locations are preserved if all datasets contain the
%   same or only a single dataset contains electrode location information.
%
% Author: Andreas Widmann, University of Leipzig, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id$

function [EEG, com] = pop_grandaverage(ALLEEG, varargin)

EEG = [];
com = '';

% Convert args to struct
args = struct(varargin{:});

if nargin < 1
    % Pop up file selection ui
    [args.filenames args.pathname] = uigetfile('*.set', 'Select datasets -- pop_grandaverage', 'Multiselect', 'on');
    if isnumeric(args.filenames) && args.filenames == 0, return, end
    if ischar(args.filenames)
        args.filenames = {args.filenames};
    end
elseif iscell(ALLEEG)
    args.filenames = ALLEEG;
    clear ALLEEG;
else
    % Pop up dataset selection ui
    if ~isfield(args, 'datasets')
        drawnow;
        uigeom = {1 1};
        uilist = {{'style' 'text' 'string' 'Select datasets:'} ...
                  {'style' 'listbox' 'string' {ALLEEG(:).setname} 'max' length(ALLEEG)}};
        result = inputgui(uigeom, uilist, 'pophelp(''pop_grandaverage'')', 'Grand average loaded datasets -- pop_grandaverage()', [], 'normal', [1 length(ALLEEG)]);
        if isempty(result), return, end

        args.datasets = result{1};
    end
    ALLEEG = ALLEEG(args.datasets);
end

EEG = eeg_emptyset;

% Load datasets
if isfield(args, 'filenames')
    if ~isfield(args, 'pathname')
        args.pathname = '.';
    end
    for file = 1:length(args.filenames)
        ALLEEG(file) = pop_loadset('filename', args.filenames{file}, 'filepath', args.pathname);
    end
end

% Grand average data
for dataset = 1:length(ALLEEG)
    EEG.data = cat(3, EEG.data, mean(ALLEEG(dataset).data, 3));

    % Get event type from first time locking event
    if isfield(args, 'eventtype')
        EEG.event(dataset).type = args.eventtype;
    elseif isfield(ALLEEG(dataset).event, 'type') && isfield(ALLEEG(dataset).event, 'latency')
        tleLatency = (0 - ALLEEG(dataset).xmin) * ALLEEG(dataset).srate + 1;
        [foo, tleIndex] = min([ALLEEG(dataset).event.latency] - tleLatency);
        EEG.event(dataset).type = ALLEEG(dataset).event(tleIndex).type;
    else
        EEG.event(dataset).type = 'TLE'; % default
    end
end

% pnts, nbchan, srate, xmin, xmax
for iField = {'pnts' 'nbchan' 'srate' 'xmin' 'xmax'}
    EEG.(iField{:}) = unique([ALLEEG.(iField{:})]);
    if length(EEG.(iField{:})) ~= 1
        error(['Field EEG.' iField{:} ' not uniform.']);
    end
end

% Reference
try
    EEG.ref = unique(cat(1, ALLEEG.ref), 'rows');
    if size(EEG.ref, 1) > 1
        error('Field EEG.ref not uniform.');
    end
catch
    error('Field EEG.ref not uniform.');
end

% Trials
EEG.trials = size(EEG.data, 3);

% Channel labels and locations
Chanlocs = cat(1, ALLEEG.chanlocs);
fieldArray = fieldnames(Chanlocs);

for iField = 1:length(fieldArray)
    for iChan = 1:EEG.nbchan

        EEG.chanlocs(iChan).(fieldArray{iField}) = unique(cat(1, Chanlocs(:, iChan).(fieldArray{iField})), 'rows');

        if size(EEG.chanlocs(iChan).(fieldArray{iField}), 1) > 1 && strcmpi(fieldArray{iField}, 'labels')
            error('Structure EEG.chanlocs.labels not uniform.');
        elseif size(EEG.chanlocs(iChan).(fieldArray{iField}), 1) > 1
            warning(['Structure EEG.chanlocs not uniform. Grand-averaging of electrode locations not yet implemented. Skipping field ' fieldArray{iField} '.']);
            EEG.chanlocs = rmfield(EEG.chanlocs, fieldArray{iField});
            break
        end

    end
end

% Events
tmp = num2cell(round((0 - EEG.xmin) * EEG.srate + 1) + (0:length(EEG.event) - 1) * EEG.pnts);
[EEG.event(:).latency] = deal(tmp{:});
tmp = num2cell(1:length(EEG.event));
[EEG.event(:).epoch] = deal(tmp{:});
[EEG.event(:).trials] = deal(ALLEEG.trials);
[EEG.event(:).setname] = deal(ALLEEG.setname);

% Check dataset
EEG = eeg_checkset(EEG, 'eventconsistency');

% History string
if isfield(args, 'filenames')
    com = sprintf('EEG = %s({', mfilename);
    for file = 1:length(args.filenames) - 1
        com = [com sprintf('''%s'' ', args.filenames{file})];
    end
    com = [com sprintf('''%s''}', args.filenames{end})];
    args = rmfield(args, 'filenames');
else
    com = sprintf('EEG = %s(ALLEEG', mfilename);
end
for c = fieldnames(args)'
    if ischar(args.(c{:}))
        com = [com sprintf(', ''%s'', ''%s''', c{:}, args.(c{:}))];
    else
        com = [com sprintf(', ''%s'', %s', c{:}, mat2str(args.(c{:})))];
    end
end
com = [com ');'];
