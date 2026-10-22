%% Compute resting-state posterior alpha baseline from 2-s EEG epochs
% Purpose:
%   - Compute resting-state alpha power baseline for the VOFF alpha LMM.
%   - Output alpha_rest_baseline_2s.csv with required columns:
%       Subject, AlphaPower_baseline
%
% Processing:
%   1) Load resting-state EEGLAB .set files.
%   2) Resample temporarily to 250 Hz.
%   3) Apply average reference using EEG channels only.
%      Auxiliary channels ECG/EDA/ACCx/ACCy/ACCz are retained but excluded
%      from rereferencing and spectral estimation.
%   4) Extract resting segments from S111 + 10 s to S111 + 100 s.
%      By default, the first two S111 markers are used.
%   5) Divide each resting segment into non-overlapping 2-s epochs.
%   6) Compute posterior ROI alpha power, 8-12 Hz, using a Hamming-tapered PSD estimate.
%   7) Save epoch-level, segment-level, and subject-level baseline CSV files.
%
% Notes:
%   - This script does NOT overwrite original .set files.
%   - The output alpha_rest_baseline_2s.csv is automatically used by
%     Alpha_VOFF_LMM_TypeIII.m if placed in the same base_dir.

clear; clc;

%% ===================== User settings =====================
datadir = 'D:\online_learning\resting_state';

% Candidate file patterns for resting-state cleaned files.
% The script will combine all matching files and remove duplicates.
file_patterns = {'*_ica_clean.set'};

% Rest trigger and segment definition
rest_trigger = 'S111';
use_first_n_rest_triggers = 2;    % use first two S111 triggers
rest_start_sec = 10;              % segment starts 10 s after S111
rest_end_sec   = 100;             % segment ends 100 s after S111, duration = 90 s

% Resampling
target_srate = 250;

% Average reference
% Use EEG channels only. Auxiliary channels are retained but excluded from reference computation.
apply_average_reference = true;

% Alpha settings
alpha_band = [8 12];
use_log10_power = true;

% Posterior alpha ROI; keep identical to training VOFF alpha analysis.
alpha_roi = {'E70','E71','E72','E75','E76','E77','E67','E83'};

% Auxiliary channels to retain but exclude from EEG analysis/rereference
aux_labels = {'ECG','EDA','ACCx','ACCy','ACCz','EACCx','EACCy','EACCz'};

% Resting PSD epoch settings
epoch_len_sec = 2.0;
epoch_overlap_sec = 0.0;

% PSD settings for each 2-s epoch
% Use one Hamming-tapered window spanning the full epoch.
welch_window_sec = [];   % [] = full epoch length
welch_overlap_frac = 0;
nfft_min = 512;

% Optional output of resampled/rereferenced .set copy
save_preprocessed_set = false;
preprocessed_suffix = '_rest_alpha250_avgref.set';

% Output files
out_epoch_csv   = fullfile(datadir, 'alpha_rest_epoch_level_2s.csv');
out_segment_csv = fullfile(datadir, 'alpha_rest_segment_level_2s.csv');
out_baseline_csv = fullfile(datadir, 'alpha_rest_baseline_2s.csv');
out_channel_csv = fullfile(datadir, 'alpha_rest_channel_level_2s.csv');

%% ===================== Start EEGLAB =====================
eeglab; close(gcf);

files = collect_files(datadir, file_patterns);
if isempty(files)
    error('No resting-state .set files found in %s using the specified patterns.', datadir);
end

EPOCH_ALL = table();
CHAN_ALL = table();

%% ===================== File loop =====================
for curfile = 1:numel(files)

    filename = files(curfile).name;
    [~, nam, ~] = fileparts(filename);

    fprintf('\n==============================\n');
    fprintf('Processing %s\n', filename);
    fprintf('==============================\n');

    EEG = pop_loadset('filename', filename, 'filepath', datadir);
    EEG = eeg_checkset(EEG, 'eventconsistency');

    %% ----------------------------------------------------------
    % 1) Resample to 250 Hz
    %% ----------------------------------------------------------
    if abs(EEG.srate - target_srate) > 1e-6
        fprintf('Resampling from %.3f Hz to %.3f Hz...\n', EEG.srate, target_srate);
        EEG = pop_resample(EEG, target_srate);
        EEG = eeg_checkset(EEG, 'eventconsistency');
    else
        fprintf('Sampling rate already %.3f Hz. Skip resampling.\n', EEG.srate);
    end

    fs = EEG.srate;

    %% ----------------------------------------------------------
    % 2) Identify EEG channels and exclude auxiliary channels
    %% ----------------------------------------------------------
    raw_labels = {EEG.chanlocs.labels};
    canon_labels = cellfun(@canonical_label, raw_labels, 'UniformOutput', false);
    canon_aux = cellfun(@canonical_label, aux_labels, 'UniformOutput', false);

    is_aux = ismember(canon_labels, canon_aux);

    % EEG channels are defined as standard EGI electrodes labeled E<number>.
    % This prevents ECG/EDA/ACC or other non-EEG channels from entering the
    % average reference or spectral analysis even if they are not listed above.
    is_eeg_label = false(size(canon_labels));
    for ch = 1:numel(canon_labels)
        is_eeg_label(ch) = ~isempty(regexp(canon_labels{ch}, '^E\d+$', 'once'));
    end

    eeg_idx = find(is_eeg_label & ~is_aux);

    % Fallback for datasets whose EEG labels are not E<number>.
    if isempty(eeg_idx)
        warning('No E<number> EEG labels found in %s. Falling back to all non-auxiliary channels.', filename);
        eeg_idx = find(~is_aux);
    end

    if isempty(eeg_idx)
        warning('No EEG channels left after excluding auxiliary channels in %s. Skip.', filename);
        continue;
    end

    fprintf('EEG channels used for average reference and alpha analysis: %d\n', numel(eeg_idx));
    fprintf('Auxiliary/non-EEG channels retained but excluded: ');
    disp(raw_labels(~ismember(1:numel(raw_labels), eeg_idx)));

    %% ----------------------------------------------------------
    % 3) Average reference using EEG channels only
    %% ----------------------------------------------------------
    if apply_average_reference
        eeg_data = double(EEG.data(eeg_idx, :));
        avg_ref = mean(eeg_data, 1, 'omitnan');
        EEG.data(eeg_idx, :) = eeg_data - avg_ref;
        EEG = eeg_checkset(EEG);
        fprintf('Applied average reference using EEG channels only.\n');
    else
        fprintf('Average reference skipped because apply_average_reference=false.\n');
    end

    %% ----------------------------------------------------------
    % 4) Identify posterior alpha ROI channels
    %% ----------------------------------------------------------
    canon_roi = cellfun(@canonical_label, alpha_roi, 'UniformOutput', false);
    roi_idx = find(ismember(canon_labels, canon_roi) & ~is_aux);

    if isempty(roi_idx)
        warning('No posterior alpha ROI channels found in %s. Skip.', filename);
        continue;
    end

    fprintf('Posterior alpha ROI channels found (%d): %s\n', ...
        numel(roi_idx), strjoin(raw_labels(roi_idx), ', '));

    %% Optional save of resampled/rereferenced copy
    if save_preprocessed_set
        pop_saveset(EEG, 'filename', [nam preprocessed_suffix], 'filepath', datadir);
    end

    %% ----------------------------------------------------------
    % 5) Find resting S111 events
    %% ----------------------------------------------------------
    ev_type = cell(numel(EEG.event), 1);
    for e = 1:numel(EEG.event)
        ev_type{e} = canonical_event_type(EEG.event(e).type);
    end

    rest_clean = canonical_event_type(rest_trigger);
    idx_rest = find(strcmp(ev_type, rest_clean));

    if isempty(idx_rest)
        warning('No %s events found in %s. Skip.', rest_trigger, filename);
        continue;
    end

    rest_lat = [EEG.event(idx_rest).latency];
    rest_lat = round(rest_lat(:));
    rest_lat = sort(rest_lat);

    if numel(rest_lat) < use_first_n_rest_triggers
        warning('Only %d %s events found in %s; expected at least %d. Using all available.', ...
            numel(rest_lat), rest_trigger, filename, use_first_n_rest_triggers);
        rest_lat = rest_lat(:);
    else
        rest_lat = rest_lat(1:use_first_n_rest_triggers);
    end

    fprintf('Using %d resting %s segment(s).\n', numel(rest_lat), rest_trigger);

    %% ----------------------------------------------------------
    % 6) Extract rest segments, split into 2-s epochs, compute alpha power
    %% ----------------------------------------------------------
    epoch_len_samp = round(epoch_len_sec * fs);
    step_samp = round((epoch_len_sec - epoch_overlap_sec) * fs);
    if step_samp <= 0
        error('Invalid epoch overlap. Step size must be > 0.');
    end

    subject_epoch_rows = table();
    subject_chan_rows = table();

    subj_num = extract_subject_number(string(nam));

    for seg_id = 1:numel(rest_lat)

        lat = rest_lat(seg_id);
        seg_start = lat + round(rest_start_sec * fs);
        seg_end   = lat + round(rest_end_sec   * fs) - 1;

        seg_start = max(1, seg_start);
        seg_end   = min(EEG.pnts, seg_end);

        if seg_end <= seg_start || (seg_end - seg_start + 1) < epoch_len_samp
            warning('Rest segment %d in %s is too short or out of range. Skip.', seg_id, filename);
            continue;
        end

        ep_starts = seg_start:step_samp:(seg_end - epoch_len_samp + 1);

        for ep = 1:numel(ep_starts)

            ep_start = ep_starts(ep);
            ep_end = ep_start + epoch_len_samp - 1;

            seg = double(EEG.data(roi_idx, ep_start:ep_end));

            [alpha_by_chan_raw, alpha_by_chan_out] = compute_alpha_for_channels( ...
                seg, fs, alpha_band, use_log10_power, welch_window_sec, ...
                welch_overlap_frac, nfft_min);

            roi_alpha_raw = mean(alpha_by_chan_raw, 'omitnan');
            roi_alpha = mean(alpha_by_chan_out, 'omitnan');
            n_roi_chan_used = sum(isfinite(alpha_by_chan_out));

            for c = 1:numel(roi_idx)
                tmpc = table();
                tmpc.File = string(filename);
                tmpc.Subject = string(nam);
                tmpc.SubjectNum = subj_num;
                tmpc.RestSegment = seg_id;
                tmpc.Epoch = ep;
                tmpc.Channel = string(raw_labels{roi_idx(c)});
                tmpc.EpochStart_sample = ep_start;
                tmpc.EpochEnd_sample = ep_end;
                tmpc.EpochLength_sec = epoch_len_sec;
                tmpc.AlphaPower_raw = alpha_by_chan_raw(c);
                tmpc.AlphaPower = alpha_by_chan_out(c);
                tmpc.AlphaMetric = string(metric_name(use_log10_power));
                subject_chan_rows = [subject_chan_rows; tmpc]; %#ok<AGROW>
            end

            tmpe = table();
            tmpe.File = string(filename);
            tmpe.Subject = string(nam);
            tmpe.SubjectNum = subj_num;
            tmpe.RestSegment = seg_id;
            tmpe.Epoch = ep;
            tmpe.RestTrigger = string(rest_trigger);
            tmpe.RestTrigger_latency_sample = lat;
            tmpe.RestStart_sec_from_trigger = rest_start_sec;
            tmpe.RestEnd_sec_from_trigger = rest_end_sec;
            tmpe.EpochStart_sample = ep_start;
            tmpe.EpochEnd_sample = ep_end;
            tmpe.EpochLength_sec = epoch_len_sec;
            tmpe.Fs = fs;
            tmpe.N_ROI_Channels = n_roi_chan_used;
            tmpe.AlphaPower_raw = roi_alpha_raw;
            tmpe.AlphaPower = roi_alpha;
            tmpe.AlphaMetric = string(metric_name(use_log10_power));
            tmpe.Usable = isfinite(roi_alpha) && n_roi_chan_used > 0;

            subject_epoch_rows = [subject_epoch_rows; tmpe]; %#ok<AGROW>
        end
    end

    fprintf('Usable resting alpha epochs: %d / %d\n', ...
        sum(subject_epoch_rows.Usable), height(subject_epoch_rows));

    EPOCH_ALL = [EPOCH_ALL; subject_epoch_rows]; %#ok<AGROW>
    CHAN_ALL = [CHAN_ALL; subject_chan_rows]; %#ok<AGROW>
end

%% ===================== Save epoch-level CSV =====================
writetable(EPOCH_ALL, out_epoch_csv);
fprintf('\nEpoch-level output saved to:\n%s\n', out_epoch_csv);

%% ===================== Save channel-level CSV =====================
writetable(CHAN_ALL, out_channel_csv);
fprintf('Channel-level output saved to:\n%s\n', out_channel_csv);

%% ===================== Segment-level summary =====================
if isempty(EPOCH_ALL)
    warning('No resting epoch rows produced. No baseline file created.');
    return;
end

Euse = EPOCH_ALL(EPOCH_ALL.Usable & isfinite(EPOCH_ALL.AlphaPower), :);

if isempty(Euse)
    warning('No usable resting epochs. No baseline file created.');
    return;
end

Gseg = findgroups(Euse.Subject, Euse.SubjectNum, Euse.RestSegment);
seg_subject = splitapply(@(x)x(1), Euse.Subject, Gseg);
seg_subject_num = splitapply(@(x)x(1), Euse.SubjectNum, Gseg);
seg_id = splitapply(@(x)x(1), Euse.RestSegment, Gseg);
seg_n_epoch = splitapply(@numel, Euse.AlphaPower, Gseg);
seg_alpha = splitapply(@mean_omitnan, Euse.AlphaPower, Gseg);
seg_alpha_sd = splitapply(@std_omitnan, Euse.AlphaPower, Gseg);
seg_alpha_raw = splitapply(@mean_omitnan, Euse.AlphaPower_raw, Gseg);

SEG = table();
SEG.Subject = seg_subject;
SEG.SubjectNum = seg_subject_num;
SEG.RestSegment = seg_id;
SEG.N_Epochs = seg_n_epoch;
SEG.AlphaPower_segment = seg_alpha;
SEG.AlphaPower_segment_SD = seg_alpha_sd;
SEG.AlphaPower_raw_segment = seg_alpha_raw;
SEG.AlphaMetric = repmat(string(metric_name(use_log10_power)), height(SEG), 1);

writetable(SEG, out_segment_csv);
fprintf('Segment-level output saved to:\n%s\n', out_segment_csv);

%% ===================== Subject-level baseline summary =====================
Gsub = findgroups(Euse.Subject, Euse.SubjectNum);
base_subject = splitapply(@(x)x(1), Euse.Subject, Gsub);
base_subject_num = splitapply(@(x)x(1), Euse.SubjectNum, Gsub);
base_n_epoch = splitapply(@numel, Euse.AlphaPower, Gsub);
base_n_segment = splitapply(@(x)numel(unique(x)), Euse.RestSegment, Gsub);
base_alpha = splitapply(@mean_omitnan, Euse.AlphaPower, Gsub);
base_alpha_sd = splitapply(@std_omitnan, Euse.AlphaPower, Gsub);
base_alpha_raw = splitapply(@mean_omitnan, Euse.AlphaPower_raw, Gsub);

BASE = table();
BASE.Subject = base_subject;
BASE.SubjectNum = base_subject_num;
BASE.N_RestSegments = base_n_segment;
BASE.N_Epochs = base_n_epoch;
BASE.AlphaPower_baseline = base_alpha;
BASE.AlphaPower_baseline_SD = base_alpha_sd;
BASE.AlphaPower_raw_baseline = base_alpha_raw;
BASE.AlphaMetric = repmat(string(metric_name(use_log10_power)), height(BASE), 1);

writetable(BASE, out_baseline_csv);
fprintf('Subject-level baseline output saved to:\n%s\n', out_baseline_csv);

fprintf('\nAll done. Use alpha_rest_baseline_2s.csv in Alpha_VOFF_LMM_TypeIII.m.\n');

%% ===================== Helper functions =====================
function files = collect_files(datadir, patterns)
    files = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
    seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for p = 1:numel(patterns)
        tmp = dir(fullfile(datadir, patterns{p}));
        for i = 1:numel(tmp)
            if tmp(i).isdir
                continue;
            end
            key = lower(tmp(i).name);
            if ~isKey(seen, key)
                files(end+1) = tmp(i); %#ok<AGROW>
                seen(key) = true;
            end
        end
    end
end

function s = canonical_label(x)
    if isstring(x), x = char(x); end
    if ~ischar(x), x = char(string(x)); end
    s = upper(regexprep(x, '[\s_\-]', ''));
end

function s = canonical_event_type(x)
    if ischar(x) || isstring(x)
        s = char(x);
    elseif isnumeric(x)
        s = num2str(x);
    else
        s = '';
    end
    s = upper(regexprep(strtrim(s), '\s+', ''));
end

function [alpha_raw, alpha_out] = compute_alpha_for_channels(seg, fs, alpha_band, use_log10_power, welch_window_sec, overlap_frac, nfft_min)
    % seg: channels x samples
    nchan = size(seg, 1);
    nsamp = size(seg, 2);
    alpha_raw = nan(nchan, 1);
    alpha_out = nan(nchan, 1);

    if isempty(welch_window_sec)
        win_len = nsamp;
    else
        win_len = min(nsamp, max(1, round(welch_window_sec * fs)));
    end

    win = hamming(win_len);
    noverlap = round(overlap_frac * win_len);
    noverlap = min(noverlap, win_len - 1);
    nfft = max(nfft_min, 2^nextpow2(win_len));

    for ch = 1:nchan
        x = double(seg(ch, :));
        x = x - mean(x, 'omitnan');
        x = detrend(x);

        if all(~isfinite(x)) || std(x, 'omitnan') == 0
            continue;
        end

        [pxx, f] = pwelch(x, win, noverlap, nfft, fs);
        idx = f >= alpha_band(1) & f <= alpha_band(2);

        if ~any(idx)
            continue;
        end

        ap = mean(pxx(idx), 'omitnan');
        alpha_raw(ch) = ap;

        if use_log10_power
            alpha_out(ch) = log10(max(ap, realmin));
        else
            alpha_out(ch) = ap;
        end
    end
end

function y = mean_omitnan(x)
    y = mean(x, 'omitnan');
end

function y = std_omitnan(x)
    y = std(x, 'omitnan');
end

function s = metric_name(use_log10_power)
    if use_log10_power
        s = 'log10_mean_PSD_8_12Hz';
    else
        s = 'mean_PSD_8_12Hz';
    end
end

function n = extract_subject_number(s)
    % Extract the first integer from a file/subject string.
    % Examples: '01_REST_ica_clean' -> 1; 'sub-033_rest' -> 33.
    s = char(string(s));
    tok = regexp(s, '\d+', 'match', 'once');
    if isempty(tok)
        n = NaN;
    else
        n = str2double(tok);
    end
end
