
function [amplitudes] = spec_plot_FSAReward(topo, eventlabels)

fft_data = [];
for ifreq = 1:numel(topo.freq) % loop through frequencies
    for icond = 1:numel(topo.conds) % loop through conditions
        filename = ['Mean-C' num2str(topo.conds(icond)) '.set']; % which grand average to load
        EEG = pop_loadset(filename, topo.pathin); % load grand average
        EEG = pop_eegfiltnew(EEG, [], 1, 1690, true, [], 0); % Hamming windowed sinc FIR filter, passband edge 1 Hz, filter order 1691 (estimated), transition bandwidth 1 Hz, cutoff frequency (-6 dB) 0.5 Hz
        sprange = eegF_Time2Sp(EEG, topo.time(1)):eegF_Time2Sp(EEG, topo.time(2)); % time range of interest (i.e., excluding the initial ERP response)
        EEG = eeg_detrend_widmann(EEG); % detrend
        for isub = 1:numel(topo.subjects) % loop through participants
            fft_data(:, :, isub, icond, ifreq) = abs(fft(EEG.data(topo.channels(isub, :), sprange, topo.subjects(isub)), numel(sprange) * 8, 2)) * 2 / size(sprange, 2); % fast Fourier transform (normalized)
            % zeropadded to 8 times the length of the data, to increase frequency resolution and be able to visualize exactly 10 Hz and 12 Hz (as well as visually magnify possible differences between conditions)
        end
    end
end

% save data to plot in R
freqs = ((0:(numel(sprange) * 8) - 1) / (numel(sprange) * 8)) * topo.samprate; % frequencies
% extract amplitude until 30 Hz (more than enough for plotting)
[~, stopfreq] = min(abs(freqs - 30)); % look for the sampling point (on x-axis) corresponding to 30 Hz
data = squeeze(mean(fft_data(:, 1:stopfreq, :, :, :), 1)); % average across electrodes
data_plot = [];
counter = 1;
for i_freq = 1:size(data, 4) % loop through frequencies
    for i_cond = 1:size(data, 3) % loop through conditions
        for i_ssj = 1:size(data, 2) % loop through participants
            temp = data(:, i_ssj, i_cond, i_freq)';
            data_plot(counter, :) = [topo.subjnum(i_ssj), topo.freq(i_freq), topo.conds(i_cond), temp];
            counter = counter + 1;
        end
    end
end
data_header = ['participant' 'frequency' 'condition' num2cell(freqs(1:stopfreq))];
data_plot = [data_header; num2cell(data_plot)];
disp('****************************')
disp('Saving data...')
disp('Please be patient.')
disp('****************************')
cell2csv([topo.pathout 'grandAverage_spectra.csv'], data_plot); % slow

%% PLOT SPECTRA

allplots = [];
figure
freq_axis = ((0:(numel(sprange) * 8) - 1) / (numel(sprange) * 8)) * topo.samprate; % frequencies from which to extract amplitudes
for ifreq = 1:numel(topo.freq) % loop through frequencies
    for icond = 1:numel(topo.conds) % loop through conditions
        spectradata = squeeze(mean(mean(fft_data(:, :, :, :, ifreq), 3), 1)); % average across participants and selected electrodes
        allplots(icond) = plot(freq_axis, spectradata(:, icond), 'LineWidth', 2.5); % plot spectra
        axis([topo.bound topo.amp_spectramaplim]); % set boundaries of x- axis and y-axis
        xlabel('Frequency (Hz)', 'fontsize', 15); % label x-axis
        ylabel ('Amplitude(\muV)', 'fontsize', 15); % label y-axis
        hold all
    end
end
suptitle(['Grand AVG ' num2str(topo.freq) ' Hz']) % see 10 Hz and 12 Hz
legend(allplots, eventlabels); % condition names

% look for point corresponding to frequency of interest
for freqnum = 1:numel(topo.freq)
    [~,spoi] = min(abs(freq_axis - topo.freq(freqnum))); % look for the sampling point (on x-axis) corresponding to frequencies of interest
    amps(:, :, freqnum) = squeeze(mean(fft_data(:, spoi, :, :, freqnum), 1)); % extract amplitude value of specified frequency averaged across electrodes
end

% permute amps so that they can be saved as .csv
amplitudes = [];
for freqnum = 1:numel(topo.freq) % loop through frequencies
    temp = amps(:, :, freqnum); % extract amplitude value of specified frequency averaged across electrodes
    amplitudes = [amplitudes; temp]; % append to big matrix
    clear temp
end
subject = [topo.subjnum, topo.subjnum]'; % repeat participants number twice (for 10 and 12 Hz, respectively)
% subject = kron(topo.subjnum, ones(1, numel(topo.freq)))'; % repeat participants number as many times as there are frequencies
frequency = [];
for freqnum2 = 1:numel(topo.freq) % loop through frequencies
    tempfreq = [repmat(topo.freq(freqnum2), 1, numel(topo.subjnum))]'; % repeat frequency once for each participant
    frequency = [frequency; tempfreq];
end
amplitudes = cat(1, [{'Subject'}, {'Frequency'}, eventlabels], num2cell(cat(2, subject, frequency, amplitudes))); % concatenate

cell2csv([topo.pathout '\grandAverage_amplitudes.csv'], amplitudes); % save amplitude values in a .csv file

end
