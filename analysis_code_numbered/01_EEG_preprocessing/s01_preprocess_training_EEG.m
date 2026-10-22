%% Fully automatic pipeline (TRAINING ONLY; EEG-only analysis, auxiliary channels kept):
% - Remove TEST continuous segments marked by S  88 (until next S  111)
% - Keep auxiliary channels in saved datasets: ECG, EDA, ACCx, ACCy, ACCz
% - Exclude auxiliary channels from artifact rejection, PSD bad-channel detection, filtering, ICA, and ocular-IC detection
% - PSD-only bad channels on EEG channels only (EOG whitelist protected)
% - FieldTrip filtering on EEG channels only
% - ICA on EEG channels only; weights copied back to the full dataset
% - Auto remove ocular ICs (ICLabel optional + VEOG/HEOG corr), MAX 5 ICs
% - Do NOT resample the final saved datasets
clear; clc;

%% ===================== User settings =====================
datadir = 'D:\online_learning\training_psychophysiology';
filt = '*.vhdr';

% ---- Auxiliary channels: keep in output, exclude from EEG analysis ----
aux_channels = {'ECG','EDA','ACCx','ACCy','ACCz'};

% ---- Triggers ----
trig_train = 'S111';
trig_test  = 'S 88';

% ---- PSD-only bad channel params (EEGLAB spectopo) ----
psd_freqrange = [1 80];
psd_hp_check  = 0.5;      % set [] to disable
psd_kMAD      = 5;
psd_rule_any  = true;     % true:any metric abnormal; false:>=2 metrics abnormal

band_broad = [1 80];
band_emg   = [30 80];
band_erp   = [1 15];
band_drift = [0.5 2];
band_ref   = [2 8];

% ---- EOG peri-orbital whitelist (do NOT remove by PSD) ----
eog_whitelist = {'E8','E9','E126','E127','E1','E32','E17','E25'};

% ---- Whitelist "hard-bad" upgrade ----
flatline_sec = 5;         % (1) flatline >= 5 sec
rep_thr      = 0.20;      % (2) clipping proxy: mean(diff(x)==0) >= 0.20
win_rms_sec  = 1;         % RMS window length (sec)
sustain_frac_thr = 0.50;  % sustained fraction threshold

% ---- Filtering (FieldTrip) ----
mains_is_50hz = true; % set false if 60 Hz
lp_freq = 30;
hp_freq = 0.1;

% ---- ICA ----
ica_hp = 2;        % Hz high-pass for ICA
ica_rs = 100;      % Hz resample for ICA

% ---- Auto ocular IC removal ----
use_iclabel  = false;     % requires ICLabel plugin
thr_eye_prob = 0.90;      % ICLabel Eye probability threshold
thr_eye_corr = 0.40;      % corr(IC, VEOG/HEOG) threshold
max_rm_ic    = 5;         % <<< MAX delete ICs

%% ===================== Start EEGLAB =====================
eeglab; close(gcf);

%% ===================== File loop =====================
cd(datadir);
outputdir = datadir;
files = dir(fullfile(datadir, filt));

badchan = cell(numel(files),1);

for curfile = 1:length(files)
    file = files(curfile).name;
    EEG  = pop_loadbv(datadir, file);
    [~, nam, ext] = fileparts(file);

    filenam = nam; % no pre/post here

    fprintf('\n==============================\n');
    fprintf('Working on %s\n', [filenam ext]);
    fprintf('==============================\n');

    %% ----------------------------------------------------------
    % 0) Remove TEST segments: delete continuous data from S  88 to next S  111
    %% ----------------------------------------------------------
    EEG = eeg_checkset(EEG, 'eventconsistency');

    % Collect event types as strings
    evTypes = cell(numel(EEG.event),1);
    for k = 1:numel(EEG.event)
        if ischar(EEG.event(k).type)
            evTypes{k} = EEG.event(k).type;
        elseif isnumeric(EEG.event(k).type)
            evTypes{k} = num2str(EEG.event(k).type);
        else
            evTypes{k} = '';
        end
    end

    idx_test  = find(strcmp(evTypes, trig_test));
    idx_train = find(strcmp(evTypes, trig_train));

    rm_ranges = []; % [start end] in points
    if ~isempty(idx_test)
        % Sort by latency (just in case)
        test_lat = arrayfun(@(e)e.latency, EEG.event(idx_test));
        [~, ord] = sort(test_lat);
        idx_test = idx_test(ord);

        for ii = 1:numel(idx_test)
            lat_test = EEG.event(idx_test(ii)).latency;

            % Find next training event AFTER this test event
            next_train_candidates = idx_train(arrayfun(@(e)e.latency, EEG.event(idx_train)) > lat_test);

            if ~isempty(next_train_candidates)
                lat_next_train = min(arrayfun(@(e)e.latency, EEG.event(next_train_candidates)));
                seg_start = max(1, floor(lat_test));
                seg_end   = min(EEG.pnts, floor(lat_next_train)-1);
            else
                % No next training -> remove to end
                seg_start = max(1, floor(lat_test));
                seg_end   = EEG.pnts;
            end

            if seg_end > seg_start
                rm_ranges(end+1,:) = [seg_start seg_end]; %#ok<SAGROW>
            end
        end

        % Merge overlapping ranges
        if ~isempty(rm_ranges)
            rm_ranges = sortrows(rm_ranges, 1);
            merged = rm_ranges(1,:);
            for r = 2:size(rm_ranges,1)
                if rm_ranges(r,1) <= merged(end,2) + 1
                    merged(end,2) = max(merged(end,2), rm_ranges(r,2));
                else
                    merged = [merged; rm_ranges(r,:)]; %#ok<AGROW>
                end
            end
            rm_ranges = merged;

            fprintf('[Cut TEST] Removing %d segment(s) based on %s -> next %s\n', ...
                size(rm_ranges,1), trig_test, trig_train);
            disp(rm_ranges);

            EEG = pop_select(EEG, 'nopoint', rm_ranges);
            EEG = eeg_checkset(EEG, 'eventconsistency');
        end
    else
        fprintf('[Cut TEST] No %s found. Skip cutting.\n', trig_test);
    end

    %% Identify EEG vs auxiliary channels
    % IMPORTANT:
    % Auxiliary channels are kept in the dataset, but are excluded from all EEG analyses.
    [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels);

    fprintf('[Channels] EEG channels used for analysis: %d | Aux channels kept but excluded: %d\n', ...
        numel(eeg_idx), numel(aux_idx));
    if ~isempty(aux_idx)
        fprintf('[Channels] Aux kept: %s\n', strjoin({EEG.chanlocs(aux_idx).labels}, ', '));
    end

    %% Muscle artefact rejection using EEG channels only
    [~, rejections, ~] = pop_rejcont(EEG, ...
        'threshold', 10, ...
        'elecrange', eeg_idx, ...
        'onlyreturnselection', 'on', ...
        'taper', 'hamming', ...
        'freqlimit', [30 125]);

    % Remove bad continuous time segments from all channels, so EEG and aux remain synchronized
    EEG = pop_select(EEG, 'nopoint', rejections);
    EEG = eeg_checkset(EEG, 'eventconsistency');

    %% Channel label prefix for EEG channels only
    % Do not rename ECG/EDA/ACC channels.
    [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels);
    for i = eeg_idx(:)'
        oldlab = strtrim(EEG.chanlocs(i).labels);
        if isempty(regexp(oldlab, '^E\d+$', 'once')) && ~isempty(regexp(oldlab, '^\d+$', 'once'))
            EEG.chanlocs(i).labels = ['E' oldlab];
        end
    end

    %% Load chanlocs for EEG channels; auxiliary channels are retained without being used for EEG analyses
    EEG = pop_chanedit(EEG, ...
        {'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc'}, ...
        'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BESA\\egi128_GSN.sfp');

    [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels);

    %% ----------------------------------------------------------
    % 1) PSD-only bad channels on EEG channels only (with EOG whitelist protection)
    %% ----------------------------------------------------------
    orig_chanlocs = EEG.chanlocs;

    [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels);
    eeg_labels = {EEG.chanlocs(eeg_idx).labels};
    eog_local_idx = find(ismember(eeg_labels, eog_whitelist));

    EEG_check = EEG;
    if ~isempty(psd_hp_check)
        try
            % High-pass only EEG channels for PSD checking; aux channels are untouched
            EEG_check_eeg = pop_select(EEG_check, 'channel', eeg_idx);
            EEG_check_eeg = pop_eegfiltnew(EEG_check_eeg, psd_hp_check, []);
            EEG_check.data(eeg_idx,:) = EEG_check_eeg.data;
        catch
            warning('pop_eegfiltnew failed; PSD check continues without high-pass.');
        end
    end

    bad_final_idx = [];
    try
        % spectopo is run on EEG channels only
        [spectra_db, freqs] = spectopo(EEG_check.data(eeg_idx,:), 0, EEG_check.srate, ...
            'plot','off', 'freqrange', psd_freqrange);
    catch ME
        warning('spectopo failed (%s). Skip PSD bad channel detection for this file.', ME.message);
        spectra_db = [];
        freqs = [];
    end

    bad_psd_global_idx = [];
    if ~isempty(spectra_db)
        P = 10.^(spectra_db/10);
        bandpow = @(P,f,fr) trapz(f(f>=fr(1) & f<=fr(2)), P(:, f>=fr(1) & f<=fr(2)), 2);

        P_broad = bandpow(P, freqs, band_broad);
        P_emg   = bandpow(P, freqs, band_emg);
        P_erp   = bandpow(P, freqs, band_erp);
        P_drift = bandpow(P, freqs, band_drift);
        P_ref   = bandpow(P, freqs, band_ref);

        emg_ratio   = P_emg   ./ max(P_erp, 1e-12);
        drift_ratio = P_drift ./ max(P_ref, 1e-12);

        thr = @(x,k) median(x) + k*mad(x,1);
        thr_broad = thr(P_broad, psd_kMAD);
        thr_emg   = thr(emg_ratio, psd_kMAD);
        thr_drift = thr(drift_ratio, psd_kMAD);

        flag_broad = P_broad > thr_broad;
        flag_emg   = emg_ratio > thr_emg;
        flag_drift = drift_ratio > thr_drift;
        flag_count = double(flag_broad) + double(flag_emg) + double(flag_drift);

        if psd_rule_any
            bad_psd_local_idx = find(flag_count >= 1);
        else
            bad_psd_local_idx = find(flag_count >= 2);
        end

        % Protect EOG whitelist: do not remove by PSD
        bad_psd_local_idx = setdiff(bad_psd_local_idx, eog_local_idx);
        bad_psd_global_idx = eeg_idx(bad_psd_local_idx);

        % Hard-bad check for whitelisted EOG channels only
        bad_eog_hard_global_idx = [];

        if ~isempty(eog_local_idx)
            broad_q999 = quantile(P_broad, 0.999);

            win_rms  = max(1, round(win_rms_sec * EEG_check.srate));
            step_rms = max(1, round(win_rms/2));

            rms_all = [];
            for ch0 = eeg_idx(:)'
                x0 = double(EEG_check.data(ch0,:));
                for s0 = 1:step_rms:(length(x0)-win_rms+1)
                    seg0 = x0(s0:s0+win_rms-1);
                    rms_all(end+1,1) = sqrt(mean(seg0.^2)); %#ok<SAGROW>
                end
            end
            rms_q99 = quantile(rms_all, 0.99);

            win_flat  = max(1, round(flatline_sec * EEG_check.srate));
            step_flat = max(1, round(win_flat/5));
            tol_flat  = 1e-6;

            for ii = 1:numel(eog_local_idx)
                local_ch = eog_local_idx(ii);
                ch = eeg_idx(local_ch);  % global channel index in full EEG
                x  = double(EEG_check.data(ch,:));

                % (1) flatline
                is_flat = false;
                for s0 = 1:step_flat:(length(x)-win_flat+1)
                    seg0 = x(s0:s0+win_flat-1);
                    if (max(seg0) - min(seg0)) < tol_flat
                        is_flat = true; break;
                    end
                end

                % (2) clipping proxy
                rep_frac = mean(diff(x) == 0);
                is_clip  = (rep_frac >= rep_thr);

                % (3) extreme broadband + sustained
                is_extreme_broad = (P_broad(local_ch) > broad_q999);

                is_sustained = false;
                if is_extreme_broad
                    nwin = 0; nhigh = 0;
                    for s0 = 1:step_rms:(length(x)-win_rms+1)
                        seg0 = x(s0:s0+win_rms-1);
                        r = sqrt(mean(seg0.^2));
                        nwin = nwin + 1;
                        if r > rms_q99
                            nhigh = nhigh + 1;
                        end
                    end
                    if nwin > 0
                        frac_high = nhigh / nwin;
                        is_sustained = (frac_high >= sustain_frac_thr);
                    end
                end
                is_noise_hard = is_extreme_broad && is_sustained;

                if is_flat || is_clip || is_noise_hard
                    bad_eog_hard_global_idx(end+1,1) = ch; %#ok<SAGROW>
                    fprintf('EOG-hard-bad: %s | flat=%d clip=%d noise=%d (rep=%.2f, P_broad=%.3g)\n', ...
                        EEG_check.chanlocs(ch).labels, is_flat, is_clip, is_noise_hard, rep_frac, P_broad(local_ch));
                end
            end
        end

        bad_final_idx = unique([bad_psd_global_idx(:); bad_eog_hard_global_idx(:)]);
    else
        bad_final_idx = [];
    end

    % Safety: never remove auxiliary channels here
    bad_final_idx = setdiff(bad_final_idx, aux_idx);

    fprintf('PSD bad EEG channels removed/interpolated: %d\n', numel(bad_final_idx));
    if ~isempty(bad_final_idx)
        disp({EEG.chanlocs(bad_final_idx).labels});
    end
    badchan{curfile} = {EEG.chanlocs(bad_final_idx).labels};

    if ~isempty(bad_final_idx)
        EEG = pop_select(EEG, 'nochannel', bad_final_idx);
        EEG = eeg_checkset(EEG);

        % Interpolate removed EEG channels back into the full dataset.
        % Auxiliary channels are still present and are not interpolated or analyzed.
        EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
        EEG = eeg_checkset(EEG);
    end

    [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels);

    %% ----------------------------------------------------------
    % 2) Filtering + notch (FieldTrip), EEG channels only
    %% ----------------------------------------------------------
    EEG_eegonly = pop_select(EEG, 'channel', eeg_idx);
    ft_data = eeglab2fieldtrip(EEG_eegonly, 'preprocessing', 'none');

    cfg = [];
    cfg.lpfilter   = 'yes';
    cfg.hpfilter   = 'yes';
    cfg.lpfreq     = lp_freq;
    cfg.hpfreq     = hp_freq;
    cfg.lpfiltdir  = 'twopass';
    cfg.hpfiltdir  = 'twopass';
    cfg.hpfiltord  = 4;
    cfg.detrend    = 'yes';
    cfg.bsfilter   = 'yes';

    if mains_is_50hz
        cfg.bsfreq = [49 51; 99 101];
    else
        cfg.bsfreq = [59 61; 119 121];
    end

    tempdata = ft_preprocessing(cfg, ft_data);
    EEG.data(eeg_idx,:) = tempdata.trial{1,1};
    EEG = eeg_checkset(EEG);

    %% ----------------------------------------------------------
    % 2.5) Epoch ONLY training trigger (after test cut)
    %     <<< Epoch triggers before ICA were already relabelled to the training triggers >>>
    %% ----------------------------------------------------------
   % EEG = pop_epoch(EEG, {trig_train}, [-0.2 1]);
    EEG = eeg_checkset(EEG);

    %% ----------------------------------------------------------
    % 3) ICA (HP + RS) then copy weights back
    %% ----------------------------------------------------------
    [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels);

    % ICA is trained on EEG channels only.
    EEG_forICA = pop_select(EEG, 'channel', eeg_idx);
    EEG_forICA = pop_eegfiltnew(EEG_forICA, 'locutoff', ica_hp);
    EEG_forICA = pop_resample(EEG_forICA, ica_rs); % resample only the temporary ICA copy
    EEG_forICA = pop_runica(EEG_forICA, 'icatype', 'runica', 'extended', 1, 'interrupt', 'off');

    % Copy ICA weights back to the full dataset.
    % icachansind must point to the EEG channels in the full dataset.
    EEG.icaweights  = EEG_forICA.icaweights;
    EEG.icasphere   = EEG_forICA.icasphere;
    EEG.icawinv     = EEG_forICA.icawinv;
    EEG.icachansind = eeg_idx;
    EEG = eeg_checkset(EEG);

    %% Save ICA-weighted dataset without final resampling
    pop_saveset(EEG, 'filename', [filenam '_ica.set'], 'filepath', outputdir);

    %% ----------------------------------------------------------
    % 4) Auto remove ocular ICs (ICLabel + VEOG/HEOG corr), MAX 5
    %% ----------------------------------------------------------
    [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels);

    idx_veog = find(ismember({EEG.chanlocs.labels}, {'E8','E9','E126','E127'}));
    idx_heog = find(ismember({EEG.chanlocs.labels}, {'E1','E32','E17','E25'}));

    have_veog = (numel(idx_veog) == 4);
    have_heog = (numel(idx_heog) == 4);

    VEOG = [];
    HEOG = [];
    if have_veog
        VEOG = mean(EEG.data(idx_veog(1:2),:),1) - mean(EEG.data(idx_veog(3:4),:),1);
    end
    if have_heog
        HEOG = mean(EEG.data(idx_heog(1:2),:),1) - mean(EEG.data(idx_heog(3:4),:),1);
    end

    W = EEG.icaweights * EEG.icasphere;
    icaact = W * EEG.data(EEG.icachansind, :);

    rm_iclabel = [];
    if use_iclabel
        try
            % Run ICLabel on an EEG-only temporary dataset to avoid auxiliary-channel issues.
            EEG_iclabel = pop_select(EEG, 'channel', EEG.icachansind);
            EEG_iclabel.icaweights  = EEG.icaweights;
            EEG_iclabel.icasphere   = EEG.icasphere;
            EEG_iclabel.icawinv     = EEG.icawinv;
            EEG_iclabel.icachansind = 1:numel(EEG.icachansind);
            EEG_iclabel = eeg_checkset(EEG_iclabel);
            EEG_iclabel = pop_iclabel(EEG_iclabel, 'default');
            classif = EEG_iclabel.etc.ic_classification.ICLabel.classifications;
            p_eye = classif(:,3);
            rm_iclabel = find(p_eye >= thr_eye_prob);
        catch ME
            warning('ICLabel not available or failed (%s). Skip ICLabel step.', ME.message);
            rm_iclabel = [];
        end
    end

    rm_corr = [];
    rmax = zeros(size(icaact,1),1);
    if ~isempty(VEOG) || ~isempty(HEOG)
        for ic = 1:size(icaact,1)
            rlist = [];
            if ~isempty(VEOG), rlist(end+1) = corr(icaact(ic,:)', VEOG'); end %#ok<AGROW>
            if ~isempty(HEOG), rlist(end+1) = corr(icaact(ic,:)', HEOG'); end %#ok<AGROW>
            if ~isempty(rlist)
                rmax(ic) = max(abs(rlist));
            else
                rmax(ic) = 0;
            end
        end
        rm_corr = find(rmax >= thr_eye_corr);
    end

    rmIC_all = unique([rm_iclabel(:); rm_corr(:)]);

    % Enforce MAX 5 deletion: prioritize by rmax (desc). If rmax not available -> keep first.
    rmIC = rmIC_all;
    if ~isempty(rmIC) && numel(rmIC) > max_rm_ic
        [~,ord] = sort(rmax(rmIC), 'desc');
        rmIC = rmIC(ord(1:max_rm_ic));
    end

    fprintf('Auto ocular IC removal: ICLabel=%d, Corr=%d, Union=%d, Keep(MAX%d)=%d\n', ...
        numel(rm_iclabel), numel(rm_corr), numel(rmIC_all), max_rm_ic, numel(rmIC));

    % Save removed IC topomap BEFORE pop_subcomp
    icFigDir = outputdir;
    if ~isempty(rmIC)
        figure;
        try
            pop_topoplot(EEG, 0, rmIC, sprintf('Removed ICs: %s', filenam), 0, 'electrodes','on');
            outPng = fullfile(icFigDir, sprintf('%s_removedIC_topo.png', filenam));
            print(gcf, outPng, '-dpng', '-r200');
        catch ME
            warning('Failed to save topoplot for removed ICs: %s', ME.message);
        end
        close(gcf);
    end

    if ~isempty(rmIC)
        EEG = pop_subcomp(EEG, rmIC, 0);
        EEG = eeg_checkset(EEG);
    end

    %% Save final
    EEG = pop_saveset(EEG, 'filename', [filenam '_ica_clean.set'], 'filepath', outputdir);

end

save(fullfile(datadir,'badchan_psd_auto.mat'),'badchan');
fprintf('\nAll done. Output saved to: %s\n', datadir);

%% ===================== Local helper functions =====================
function [aux_idx, eeg_idx] = get_eeg_aux_indices(EEG, aux_channels)
% Return global channel indices for auxiliary and EEG channels.
% Matching ignores spaces and is case-insensitive.
    labels = {EEG.chanlocs.labels};
    labels_norm = cellfun(@normalize_chan_label, labels, 'UniformOutput', false);
    aux_norm = cellfun(@normalize_chan_label, aux_channels, 'UniformOutput', false);

    aux_idx = find(ismember(labels_norm, aux_norm));
    eeg_idx = setdiff(1:EEG.nbchan, aux_idx);

    if isempty(eeg_idx)
        error('No EEG channels remain after excluding auxiliary channels.');
    end
end

function out = normalize_chan_label(in)
% Normalize channel labels for robust matching.
    if isstring(in)
        in = char(in);
    end
    out = lower(strrep(strtrim(in), ' ', ''));
end
