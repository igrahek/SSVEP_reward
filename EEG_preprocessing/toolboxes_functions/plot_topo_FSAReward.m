
function [elec_rank] = plot_topo_FSAReward(topo)

data = zeros(topo.numchans, round(topo.time(2) * topo.samprate), numel(topo.subjects), numel(topo.conds), size(topo.freq, 2)); % preallocate big matrix: electrodes x timepoints x participants x conditions x frequencies of interest
amp_order = zeros(topo.numchans, numel(topo.subjects), size(topo.freq, 2)); % preallocate matrix containing sorted amplitude values for each electrode
elec_rank = zeros(topo.numchans, numel(topo.subjects), size(topo.freq, 2)); %  preallocate matrix containing sorted electrode indices
for ifreq = 1:numel(topo.freq) % loop through frequencies
    for icond = 1:numel(topo.conds) % loop through conditions
        filename = ['Mean-C' num2str(topo.conds(icond)) '.set']; % which grand average to load
        EEG = pop_loadset(filename, topo.pathin); % load grand average
        EEG = eeg_detrend_widmann(EEG); % detrend
        EEG = eegF_DFT(EEG, topo.freq(ifreq), topo.time); % discrete Fourier transform
        for  isub = 1:numel(topo.subjects) % loop through participants
            data(:, :, isub, icond, ifreq) = EEG.data(:, :, topo.subjects(isub)); % save transformed data in big matrix
        end
    end
    elecdata = squeeze(mean(data(:, eegF_Time2Sp(EEG, topo.timeplot), :, :, ifreq), 4)); % average across conditions and display the selected time point (topo.timeplot)
    [temp_amp_order, temp_elec_rank] = sort(elecdata, 1, 'descend'); % sort electrodes (highest amplitude on top)
    amp_order(:, :, ifreq) = temp_amp_order; % sorted amplitude values (electrodes x number of participants x frequencies of interest)
    elec_rank(:, :, ifreq) = temp_elec_rank; % sorted electrode indices (electrodes x number of participants x frequencies of interest)
end
elec_rank = permute(elec_rank, [2 1 3]); % rearrange matrix dimensions for compatibility with findmaxelec.m (participants x electrodes x frequencies)

%% PLOT TOPO
EEG = pop_chanedit(EEG, 'load', {topo.elec_coords, 'filetype', 'autodetect'}); % assign channel locations
topodata = squeeze(mean(mean(data(:, eegF_Time2Sp(EEG, topo.timeplot), :, :, :), 4), 3)); % average across conditions and participants, and display the selected time point (topo.timeplot)
disp('****************************')
disp('Saving data...')
disp('****************************')
topodata_R = [[{'10 Hz'} {'12 Hz'}]; num2cell(topodata)];
cell2csv([topo.pathout 'grandAverage_topos.csv'], topodata_R); % save values in a .csv file (for topographies in R)

figure
for ifreq = 1:numel(topo.freq) % loop through frequencies    
    subplot(1, numel(topo.freq), ifreq); % number of subplots according to frequencies
    topoplot(topodata(:, ifreq), EEG.chanlocs, 'maplimits', 'maxmin', 'style', 'fill', 'conv', 'on', 'electrodes', topo.labels); % plot topography
    title(['Avg ' num2str(topo.freq(ifreq)) ' Hz'], 'fontsize', 15); % title subplot
    colormap(topo.maptype); % color map
    color_bar = colorbar; xlabel(color_bar, '\muV', 'fontsize', 15); set(gca, 'CLim', topo.amp_topomaplim, 'fontsize', 15); % add title and limits to colorbar
    hold on
end
suptitle('both colors') % reward is ignored
hold off

end
