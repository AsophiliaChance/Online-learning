%% Detect vocal onset/offset from ACCx/ACCy/ACCz and flag EEG alpha trials
% Purpose:
%   - Use one of ACCx/ACCy/ACCz as laryngeal vibration channel, without knowing in advance which axis is best.
%   - Detect vocal onset/offset after each S111 trigger.
%   - Add VON_*/VOFF_* events to the EEGLAB .set file.
%   - Export trial-level QC table indicating whether a fixed listening-window alpha segment is speech-free.
%
% Recommended use:
%   1) Run this script on *_ica_clean.set files that still contain ACCx/ACCy/ACCz.
%   2) In the later alpha-power script, keep only trials with AlphaUsable == true.
%   3) Extract a fixed window, e.g., S111 + 0.30 to 1.50 s, rather than variable-length windows.
%
% Notes:
%   - ACC channels are NOT used as EEG channels.
%   - This script does NOT delete EEG data.
%   - This script only adds vocal events and saves QC tables.
%   - It handles labels such as ACCx/ACCy/ACCz and EACCx/EACCy/EACCz.

clear; clc;

%% ===================== User settings =====================
datadir = 'D:\online_learning\training_psychophysiology';
filt    = '*_ica_clean.set';

% Candidate ACC channels. The script will search both ACCx and EACCx variants.
acc_candidate_labels = {'ACCx','ACCy','ACCz'};

% Training trigger
trig_train = 'S111';   % accepts S111 or S 111 in data because spaces are removed internally

% Expected vocalization search window relative to S111
search_start_sec = 0.50;
search_end_sec   = 3.50;

% Fixed EEG alpha window relative to S111
% Trials are usable only if vocal onset occurs after alpha_end_sec + speech_buffer_sec.
alpha_start_sec = 0.30;
alpha_end_sec   = 1.50;
speech_buffer_sec = 0.20;

% Optional block indexing, if each training block contains 63 S111 trials.
trials_per_block = 63;

% ACC preprocessing for laryngeal vibration
acc_bp_low  = 70;      % Hz
acc_bp_high = 400;     % Hz; automatically limited by Nyquist
filter_order = 4;
env_smooth_sec = 0.030;  % 30 ms moving average envelope smoothing

% Thresholding
thr_mad_k = 5;              % threshold = baseline median + k * MAD
min_vocal_score = 5;        % candidate peak z-like score must be >= this value
min_on_dur_sec  = 0.060;    % sustained above threshold for at least 60 ms
min_off_dur_sec = 0.150;    % below threshold for at least 150 ms to define vocal offset

% Baseline window for ACC noise threshold, relative to S111.
% This should be before expected vocalization.
base_start_sec = 0.05;
base_end_sec   = 0.45;

% Output
out_suffix = '_with_vocal_events';
qc_csv_name = 'vocal_onset_qc_from_ACCxyz.csv';
thr_csv_name = 'vocal_ACC_thresholds.csv';

% Save per-file QC table too?
save_per_file_qc = true;

%% ===================== Start EEGLAB =====================
eeglab; close(gcf);

files = dir(fullfile(datadir, filt));
if isempty(files)
    error('No files found: %s', fullfile(datadir, filt));
end

ALL = table();
THR = table();

%% ===================== File loop =====================
for curfile = 1:numel(files)

    filename = files(curfile).name;
    [~, nam, ~] = fileparts(filename);

    fprintf('\n==============================\n');
    fprintf('Processing %s\n', filename);
    fprintf('==============================\n');

    EEG = pop_loadset('filename', filename, 'filepath', datadir);
    EEG = eeg_checkset(EEG, 'eventconsistency');

    fs = EEG.srate;
    fprintf('Sampling rate: %.3f Hz\n', fs);

    %% ----------------------------------------------------------
    % 1) Find ACCx/ACCy/ACCz channels
    %% ----------------------------------------------------------
    raw_labels = {EEG.chanlocs.labels};
    canon_labels = cellfun(@canonical_label, raw_labels, 'UniformOutput', false);

    acc_idx = [];
    acc_labels_found = {};
    acc_axis_found = {};

    for a = 1:numel(acc_candidate_labels)
        axis_label = acc_candidate_labels{a};
        canon_target1 = canonical_label(axis_label);          % e.g., ACCX
        canon_target2 = canonical_label(['E' axis_label]);    % e.g., EACCX

        idx = find(strcmpi(canon_labels, canon_target1) | strcmpi(canon_labels, canon_target2), 1);
        if ~isempty(idx)
            acc_idx(end+1) = idx; %#ok<SAGROW>
            acc_labels_found{end+1} = EEG.chanlocs(idx).labels; %#ok<SAGROW>
            acc_axis_found{end+1} = axis_label; %#ok<SAGROW>
        end
    end

    if isempty(acc_idx)
        warning('No ACCx/ACCy/ACCz channels found in %s. Skip this file.', filename);
        continue;
    end

    fprintf('Found ACC channels: %s\n', strjoin(acc_labels_found, ', '));

    %% ----------------------------------------------------------
    % 2) Preprocess each ACC axis and compute envelope
    %% ----------------------------------------------------------
    bp_high = min(acc_bp_high, fs/2 - 1);
    bp_low  = acc_bp_low;

    if bp_low >= bp_high
        error('Invalid ACC band-pass range [%.1f %.1f] for fs=%.1f. Lower acc_bp_low or check data.', ...
            bp_low, bp_high, fs);
    end

    [b, a] = butter(filter_order, [bp_low bp_high] / (fs/2), 'bandpass');
    smooth_n = max(1, round(env_smooth_sec * fs));

    nACC = numel(acc_idx);
    acc_env_all = nan(nACC, EEG.pnts);

    for chA = 1:nACC
        x = double(EEG.data(acc_idx(chA), :));
        x = x - mean(x, 'omitnan');
        x = detrend(x);

        x_filt = filtfilt(b, a, x);
        env = abs(hilbert(x_filt));
        env = movmean(env, smooth_n);

        acc_env_all(chA, :) = env;
    end

    %% ----------------------------------------------------------
    % 3) Find S111 events
    %% ----------------------------------------------------------
    ev_type = cell(numel(EEG.event), 1);
    for e = 1:numel(EEG.event)
        ev_type{e} = canonical_event_type(EEG.event(e).type);
    end

    trig_clean = canonical_event_type(trig_train);
    idx_train = find(strcmp(ev_type, trig_clean));

    if isempty(idx_train)
        warning('No %s events found in %s. Skip.', trig_train, filename);
        continue;
    end

s111_lat = [EEG.event(idx_train).latency];
s111_lat = round(s111_lat(:));
    s111_lat = round(s111_lat(:));
    fprintf('Found %d %s events.\n', numel(s111_lat), trig_train);

    %% ----------------------------------------------------------
    % 4) Compute axis-specific robust thresholds from early non-speech windows
    %% ----------------------------------------------------------
    acc_thr = nan(nACC, 1);
    acc_noise_med = nan(nACC, 1);
    acc_noise_mad = nan(nACC, 1);
    acc_base_n = nan(nACC, 1);

    for chA = 1:nACC
        base_values = [];

        for tr = 1:numel(s111_lat)
            lat = s111_lat(tr);
            b1 = lat + round(base_start_sec * fs);
            b2 = lat + round(base_end_sec   * fs);
            b1 = max(1, b1);
            b2 = min(EEG.pnts, b2);

            if b2 > b1
                base_values = [base_values; acc_env_all(chA, b1:b2)']; %#ok<AGROW>
            end
        end

        base_values = base_values(isfinite(base_values));
        acc_base_n(chA) = numel(base_values);

        if isempty(base_values)
            warning('No baseline values for %s in %s. This ACC axis will be ignored.', ...
                acc_labels_found{chA}, filename);
            continue;
        end

        acc_noise_med(chA) = median(base_values, 'omitnan');
        acc_noise_mad(chA) = mad(base_values, 1);

        if acc_noise_mad(chA) <= 0 || ~isfinite(acc_noise_mad(chA))
            warning('MAD is zero/invalid for %s in %s. This ACC axis may be unreliable.', ...
                acc_labels_found{chA}, filename);
            acc_noise_mad(chA) = eps;
        end

        acc_thr(chA) = acc_noise_med(chA) + thr_mad_k * acc_noise_mad(chA);

        fprintf('%s threshold = %.6g | median=%.6g | MAD=%.6g | Nbase=%d\n', ...
            acc_labels_found{chA}, acc_thr(chA), acc_noise_med(chA), acc_noise_mad(chA), acc_base_n(chA));

        tmp_thr = table();
        tmp_thr.File = string(filename);
        tmp_thr.Subject = string(nam);
        tmp_thr.ACC_Axis = string(acc_axis_found{chA});
        tmp_thr.ACC_Label = string(acc_labels_found{chA});
        tmp_thr.BandpassLow_Hz = bp_low;
        tmp_thr.BandpassHigh_Hz = bp_high;
        tmp_thr.EnvelopeSmooth_sec = env_smooth_sec;
        tmp_thr.ThresholdMADk = thr_mad_k;
        tmp_thr.NoiseMedian = acc_noise_med(chA);
        tmp_thr.NoiseMAD = acc_noise_mad(chA);
        tmp_thr.Threshold = acc_thr(chA);
        tmp_thr.NBaselineSamples = acc_base_n(chA);

        THR = [THR; tmp_thr]; %#ok<AGROW>
    end

    %% ----------------------------------------------------------
    % 5) Detect vocal onset/offset trial by trial
    %% ----------------------------------------------------------
    min_on_samp  = max(1, round(min_on_dur_sec  * fs));
    min_off_samp = max(1, round(min_off_dur_sec * fs));

    rows = table();
    new_events = EEG.event;

    for tr = 1:numel(s111_lat)

        lat = s111_lat(tr);
        block_id = ceil(tr / trials_per_block);

        search_start = lat + round(search_start_sec * fs);
        search_end   = lat + round(search_end_sec   * fs);
        search_start = max(1, search_start);
        search_end   = min(EEG.pnts, search_end);

        onset_candidates = [];
        offset_candidates = [];
        onset_scores = [];
        onset_channels = {};
        onset_axes = {};

        if search_end > search_start

            for chA = 1:nACC

                if ~isfinite(acc_thr(chA))
                    continue;
                end

                env = acc_env_all(chA, :);
                thr = acc_thr(chA);

                seg_env = env(search_start:search_end);
                above = seg_env > thr;

                [burst_start, burst_end] = find_boolean_bursts(above);
                if isempty(burst_start)
                    continue;
                end

                burst_dur = burst_end - burst_start + 1;
                valid = find(burst_dur >= min_on_samp);
                if isempty(valid)
                    continue;
                end

                % Choose first sustained burst in this ACC axis.
                first_burst = valid(1);
                onset_local = burst_start(first_burst);
                onset_samp  = search_start + onset_local - 1;

                % Estimate offset as the first sustained below-threshold period after onset.
                above_after_onset = above(onset_local:end);
                first_below_run = find_first_sustained_false(above_after_onset, min_off_samp);

                if isempty(first_below_run)
                    offset_samp = search_end;
                else
                    % first_below_run is local to above_after_onset; offset is just before sustained below period.
                    offset_local = onset_local + first_below_run - 2;
                    offset_samp = search_start + max(offset_local, onset_local) - 1;
                end

                % Reliability score: peak envelope in the first burst relative to this axis baseline noise.
                peak_val = max(seg_env(burst_start(first_burst):burst_end(first_burst)));
                score = (peak_val - acc_noise_med(chA)) / max(acc_noise_mad(chA), eps);

                if score >= min_vocal_score
                    onset_candidates(end+1) = onset_samp; %#ok<SAGROW>
                    offset_candidates(end+1) = offset_samp; %#ok<SAGROW>
                    onset_scores(end+1) = score; %#ok<SAGROW>
                    onset_channels{end+1} = acc_labels_found{chA}; %#ok<SAGROW>
                    onset_axes{end+1} = acc_axis_found{chA}; %#ok<SAGROW>
                end
            end
        end

        % Conservative choice: earliest reliable vocal onset among ACCx/ACCy/ACCz.
        if isempty(onset_candidates)
            voice_onset_samp  = NaN;
            voice_offset_samp = NaN;
            voice_onset_sec   = NaN;
            voice_offset_sec  = NaN;
            vocal_channel_used = "none";
            vocal_axis_used = "none";
            vocal_score = NaN;
            vocal_detected = false;
        else
            [voice_onset_samp, best_i] = min(onset_candidates);
            voice_offset_samp = offset_candidates(best_i);

            voice_onset_sec  = (voice_onset_samp  - lat) / fs;
            voice_offset_sec = (voice_offset_samp - lat) / fs;

            vocal_channel_used = string(onset_channels{best_i});
            vocal_axis_used = string(onset_axes{best_i});
            vocal_score = onset_scores(best_i);
            vocal_detected = true;

            % Add EEG events for detected vocalization.
new_events(end+1).type = 'VON';
new_events(end).latency = double(voice_onset_samp);
new_events(end).duration = [];

new_events(end+1).type = 'VOFF';
new_events(end).latency = double(voice_offset_samp);
new_events(end).duration = [];
        end

        %% Decide whether fixed alpha window is speech-free
        alpha_start_samp = lat + round(alpha_start_sec * fs);
        alpha_end_samp   = lat + round(alpha_end_sec   * fs);
        speech_safe_until = alpha_end_samp + round(speech_buffer_sec * fs);

        if alpha_start_samp < 1 || alpha_end_samp > EEG.pnts || alpha_end_samp <= alpha_start_samp
            alpha_usable = false;
            reject_reason = "alpha_window_out_of_range";
        elseif ~vocal_detected
            alpha_usable = true;
            reject_reason = "no_vocal_onset_detected";
        elseif voice_onset_samp >= speech_safe_until
            alpha_usable = true;
            reject_reason = "speech_after_alpha_window";
        else
            alpha_usable = false;
            reject_reason = "early_vocal_onset";
        end

        % Save trial row.
        tmp = table();
        tmp.File = string(filename);
        tmp.Subject = string(nam);
        tmp.Trial = tr;
        tmp.Block = block_id;
        tmp.S111_latency_sample = lat;
        tmp.Fs = fs;
        tmp.VocalDetected = vocal_detected;
        tmp.VocalAxisUsed = vocal_axis_used;
        tmp.VocalChannelUsed = vocal_channel_used;
        tmp.VocalScore = vocal_score;
        tmp.VocalOnset_sample = voice_onset_samp;
        tmp.VocalOffset_sample = voice_offset_samp;
        tmp.VocalOnset_sec_from_S111 = voice_onset_sec;
        tmp.VocalOffset_sec_from_S111 = voice_offset_sec;
        tmp.AlphaWindowStart_sec = alpha_start_sec;
        tmp.AlphaWindowEnd_sec = alpha_end_sec;
        tmp.SpeechBuffer_sec = speech_buffer_sec;
        tmp.AlphaUsable = alpha_usable;
        tmp.RejectReason = reject_reason;

        rows = [rows; tmp]; %#ok<AGROW>
    end

    %% ----------------------------------------------------------
    % 6) Save EEG with vocal events
    %% ----------------------------------------------------------
    EEG.event = new_events;
    EEG = eeg_checkset(EEG, 'eventconsistency');

    out_set = [nam out_suffix];
    pop_saveset(EEG, 'filename', out_set, 'filepath', datadir);

    %% Save per-file QC
    if save_per_file_qc
        per_file_csv = fullfile(datadir, [nam '_vocal_onset_qc_from_ACCxyz.csv']);
        writetable(rows, per_file_csv);
    end

    %% Append all rows
    ALL = [ALL; rows]; %#ok<AGROW>

    %% Print summary
    fprintf('Saved EEG with vocal events: %s\n', out_set);
    fprintf('Trials total: %d\n', height(rows));
    fprintf('Vocal detected: %d / %d\n', sum(rows.VocalDetected), height(rows));
    fprintf('Alpha usable: %d / %d\n', sum(rows.AlphaUsable), height(rows));
    fprintf('Early vocal onset rejected: %d\n', sum(rows.RejectReason == "early_vocal_onset"));

    % Channel usage summary
    if any(rows.VocalDetected)
        fprintf('Vocal axis usage among detected trials:\n');
        disp(groupsummary(rows(rows.VocalDetected,:), 'VocalAxisUsed'));
    end
end

%% ===================== Save group-level QC =====================
out_qc = fullfile(datadir, qc_csv_name);
writetable(ALL, out_qc);

out_thr = fullfile(datadir, thr_csv_name);
writetable(THR, out_thr);

fprintf('\nAll done.\n');
fprintf('Group-level QC table saved to:\n%s\n', out_qc);
fprintf('ACC threshold table saved to:\n%s\n', out_thr);

%% ===================== Local helper functions =====================
function s = canonical_label(x)
    % Normalize channel labels: remove spaces/underscores/hyphens and upper-case.
    if isstring(x), x = char(x); end
    if ~ischar(x), x = char(string(x)); end
    s = upper(regexprep(x, '[\s_\-]', ''));
end

function s = canonical_event_type(x)
    % Normalize event labels, e.g., 'S 111' -> 'S111'.
    if ischar(x) || isstring(x)
        s = char(x);
    elseif isnumeric(x)
        s = num2str(x);
    else
        s = '';
    end
    s = upper(regexprep(strtrim(s), '\s+', ''));
end

function [burst_start, burst_end] = find_boolean_bursts(x)
    % Return start/end indices of true-runs in a logical vector.
    x = logical(x(:));
    d = diff([false; x; false]);
    burst_start = find(d == 1);
    burst_end   = find(d == -1) - 1;
end

function idx = find_first_sustained_false(x, min_len)
    % In a logical vector x, find first index where x becomes false and stays false
    % for at least min_len samples. Return [] if not found.
    x = logical(x(:));
    below = ~x;

    if isempty(below) || numel(below) < min_len
        idx = [];
        return;
    end

    % Find false-runs in x, i.e., true-runs in below.
    [run_start, run_end] = find_boolean_bursts(below);
    run_dur = run_end - run_start + 1;
    valid = find(run_dur >= min_len, 1, 'first');

    if isempty(valid)
        idx = [];
    else
        idx = run_start(valid);
    end
end