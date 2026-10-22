%% Fully automatic pipeline (RESTING STATE: KEEP 2x EYES-OPEN 100s, NO EPOCH):
% - Remove EYES-CLOSED segments: delete continuous data from S222 to next S111
% - Keep ONLY TWO segments: 100s after S111 (first two S111 segments)
% - PSD-only bad channels (EOG whitelist protected)
% - FieldTrip filtering
% - ICA on continuous data
% - Auto remove ocular ICs (ICLabel optional + VEOG/HEOG corr), MAX 5 ICs

clear; clc;

%% ===================== User settings =====================
resamp_srate = 250;
datadir = 'D:\online_learning\resting_state';
filt = '*.vhdr';

% ---- Triggers ----
trig_open  = 'S111';  % eyes open
trig_close = 'S222';  % eyes closed
open_keep_sec = 100;  % keep 100s after each S111
open_keep_n   = 2;    % keep only TWO open segments

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
outputdir = datadir;
files = dir(fullfile(datadir, filt));
badchan = cell(numel(files),1);
%%
for curfile = 22%1:length(files)
    file = files(curfile).name;
    EEG  = pop_loadbv(datadir, file);
    [~, nam, ext] = fileparts(file);
    filenam = nam;

    fprintf('\n==============================\n');
    fprintf('Working on %s\n', [filenam ext]);
    fprintf('==============================\n');

    %% ----------------------------------------------------------
    % 0) Cut eyes-closed segments: S222 -> next S111
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

    idx_open  = find(strcmp(evTypes, trig_open));
    idx_close = find(strcmp(evTypes, trig_close));

    rm_ranges = []; % [start end] in points (to delete)

    if ~isempty(idx_close)
        % sort close events by latency
        close_lat = arrayfun(@(e)e.latency, EEG.event(idx_close));
        [~, ord] = sort(close_lat);
        idx_close = idx_close(ord);

        open_lat_all = arrayfun(@(e)e.latency, EEG.event(idx_open));

        for ii = 1:numel(idx_close)
            lat_close = EEG.event(idx_close(ii)).latency;

            % next open AFTER this close
            next_open_lat = open_lat_all(open_lat_all > lat_close);
            if ~isempty(next_open_lat)
                lat_next_open = min(next_open_lat);
                seg_start = max(1, floor(lat_close));
                seg_end   = min(EEG.pnts, floor(lat_next_open)-1);
            else
                % no next open -> delete to end
                seg_start = max(1, floor(lat_close));
                seg_end   = EEG.pnts;
            end

            if seg_end > seg_start
                rm_ranges(end+1,:) = [seg_start seg_end]; %#ok<SAGROW>
            end
        end

        % merge overlaps
        rm_ranges = sortrows(rm_ranges, 1);
        merged = rm_ranges(1,:);
        for r = 2:size(rm_ranges,1)
            if rm_ranges(r,1) <= merged(end,2)+1
                merged(end,2) = max(merged(end,2), rm_ranges(r,2));
            else
                merged = [merged; rm_ranges(r,:)]; %#ok<AGROW>
            end
        end
        rm_ranges = merged;

        fprintf('[Cut CLOSED] Removing %d segment(s): %s -> next %s\n', size(rm_ranges,1), trig_close, trig_open);
        disp(rm_ranges);

        EEG = pop_select(EEG, 'nopoint', rm_ranges);
        EEG = eeg_checkset(EEG, 'eventconsistency');
    else
        fprintf('[Cut CLOSED] No %s found. Skip.\n', trig_close);
    end

    %% ----------------------------------------------------------
    % 0.5) Keep ONLY TWO open segments: 100s after S111 (first two)
    %% ----------------------------------------------------------
    EEG = eeg_checkset(EEG, 'eventconsistency');

    % refresh evTypes after cutting
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
    idx_open = find(strcmp(evTypes, trig_open));

    if isempty(idx_open)
        warning('[Keep OPEN] No %s found after cutting. Skip this file.', trig_open);
        continue;
    end

    open_lat = arrayfun(@(e)e.latency, EEG.event(idx_open));
    [open_lat, ord] = sort(open_lat);
    idx_open = idx_open(ord);

    n_keep = min(open_keep_n, numel(idx_open));
    seg_len_pnts = round(open_keep_sec * EEG.srate);

    keep_ranges = zeros(n_keep, 2);
    for ii = 1:n_keep
        lat0 = EEG.event(idx_open(ii)).latency;
        s0 = max(1, floor(lat0));
        e0 = min(EEG.pnts, s0 + seg_len_pnts - 1);
        keep_ranges(ii,:) = [s0 e0];
    end

    % merge (in case overlap)
    keep_ranges = sortrows(keep_ranges, 1);
    merged = keep_ranges(1,:);
    for r = 2:size(keep_ranges,1)
        if keep_ranges(r,1) <= merged(end,2)+1
            merged(end,2) = max(merged(end,2), keep_ranges(r,2));
        else
            merged = [merged; keep_ranges(r,:)]; %#ok<AGROW>
        end
    end
    keep_ranges = merged;

    fprintf('[Keep OPEN] Keeping %d segment(s): first %d x %ds after %s\n', ...
        size(keep_ranges,1), n_keep, open_keep_sec, trig_open);
    disp(keep_ranges);

    EEG = pop_select(EEG, 'point', keep_ranges);
    EEG = eeg_checkset(EEG, 'eventconsistency');

    %% Remove non-EEG channels
    EEG = pop_select(EEG, 'nochannel', {'ECG','EDA','ACCx','ACCy','ACCz'});

    %% Muscle artefact rejection (your original)
    [~, rejections, ~] = pop_rejcont(EEG, ...
        'threshold', 10, ...
        'elecrange', 1:EEG.nbchan, ...
        'onlyreturnselection', 'on', ...
        'taper', 'hamming', ...
        'freqlimit', [30 125]);
    EEG = pop_select(EEG, 'nopoint', rejections);

    %% Channel label prefix (your original)
    for i = 1:length(EEG.chanlocs)
        EEG.chanlocs(i).labels = ['E' EEG.chanlocs(i).labels];
    end

    %% Load chanlocs (your original)
    EEG = pop_chanedit(EEG, ...
        {'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc'}, ...
        'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BESA\\egi128_GSN.sfp');

    %% ----------------------------------------------------------
    % 1) PSD-only bad channels (with EOG whitelist protection)
    %% ----------------------------------------------------------
    orig_chanlocs = EEG.chanlocs;
    eog_idx = find(ismember({EEG.chanlocs.labels}, eog_whitelist));

    EEG_check = EEG;
    if ~isempty(psd_hp_check)
        try
            EEG_check = pop_eegfiltnew(EEG_check, psd_hp_check, []);
        catch
            warning('pop_eegfiltnew failed; PSD check continues without high-pass.');
        end
    end

    bad_final_idx = [];
    try
        [spectra_db, freqs] = spectopo(EEG_check.data, 0, EEG_check.srate, ...
            'plot','off', 'freqrange', psd_freqrange);
    catch ME
        warning('spectopo failed (%s). Skip PSD bad channel detection for this file.', ME.message);
        spectra_db = [];
        freqs = [];
    end

    bad_psd_idx = [];
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
            bad_psd_idx = find(flag_count >= 1);
        else
            bad_psd_idx = find(flag_count >= 2);
        end

        % Protect EOG whitelist: do not remove by PSD
        bad_psd_idx = setdiff(bad_psd_idx, eog_idx);

        % Hard-bad check for whitelisted channels
        bad_eog_hard_idx = [];

        if ~isempty(eog_idx)
            broad_q999 = quantile(P_broad, 0.999);

            win_rms  = max(1, round(win_rms_sec * EEG_check.srate));
            step_rms = max(1, round(win_rms/2));

            rms_all = [];
            for ch0 = 1:EEG_check.nbchan
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

            for ii = 1:numel(eog_idx)
                ch = eog_idx(ii);
                x  = double(EEG_check.data(ch,:));

                % (1) flatline
                is_flat = false;
                for s = 1:step_flat:(length(x)-win_flat+1)
                    seg = x(s:s+win_flat-1);
                    if (max(seg) - min(seg)) < tol_flat
                        is_flat = true; break;
                    end
                end

                % (2) clipping proxy
                rep_frac = mean(diff(x) == 0);
                is_clip  = (rep_frac >= rep_thr);

                % (3) extreme broadband + sustained
                is_extreme_broad = (P_broad(ch) > broad_q999);

                is_sustained = false;
                if is_extreme_broad
                    nwin = 0; nhigh = 0;
                    for s = 1:step_rms:(length(x)-win_rms+1)
                        seg = x(s:s+win_rms-1);
                        r = sqrt(mean(seg.^2));
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
                    bad_eog_hard_idx(end+1,1) = ch; %#ok<SAGROW>
                    fprintf('EOG-hard-bad: %s | flat=%d clip=%d noise=%d (rep=%.2f, P_broad=%.3g)\n', ...
                        EEG_check.chanlocs(ch).labels, is_flat, is_clip, is_noise_hard, rep_frac, P_broad(ch));
                end
            end
        end

        bad_final_idx = unique([bad_psd_idx(:); bad_eog_hard_idx(:)]);
    else
        bad_final_idx = [];
    end

    fprintf('PSD bad channels removed: %d\n', numel(bad_final_idx));
    if ~isempty(bad_final_idx), disp({EEG.chanlocs(bad_final_idx).labels}); end
    badchan{curfile} = {EEG.chanlocs(bad_final_idx).labels};

    if ~isempty(bad_final_idx)
        EEG = pop_select(EEG, 'nochannel', bad_final_idx);
        EEG = eeg_checkset(EEG);
        EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
        EEG = eeg_checkset(EEG);
    end

    %% ----------------------------------------------------------
    % 2) Filtering + notch (FieldTrip)
    %% ----------------------------------------------------------
    ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

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
    EEG.data = tempdata.trial{1,1};
    EEG = eeg_checkset(EEG);

    %% ----------------------------------------------------------
    % 3) ICA (HP + RS) then copy weights back (continuous)
    %% ----------------------------------------------------------
    EEG_forICA = pop_eegfiltnew(EEG, 'locutoff', ica_hp);
    EEG_forICA = pop_resample(EEG_forICA, ica_rs);
    EEG_forICA = pop_runica(EEG_forICA, 'icatype', 'runica', 'extended', 1, 'interrupt', 'off');

    EEG.icaweights  = EEG_forICA.icaweights;
    EEG.icasphere   = EEG_forICA.icasphere;
    EEG.icawinv     = EEG_forICA.icawinv;
    EEG.icachansind = EEG_forICA.icachansind;
    EEG = eeg_checkset(EEG);

    %% Resample final
    EEG = pop_resample(EEG, resamp_srate);
    EEG = pop_saveset(EEG, 'filename', [filenam '_ica.set'], 'filepath', outputdir);

    %% ----------------------------------------------------------
    % 4) Auto remove ocular ICs (ICLabel + VEOG/HEOG corr), MAX 5
    %% ----------------------------------------------------------
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
            EEG = pop_iclabel(EEG, 'default');
            classif = EEG.etc.ic_classification.ICLabel.classifications;
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
    rmIC = rmIC_all;

    if ~isempty(rmIC) && numel(rmIC) > max_rm_ic
        [~,ord] = sort(rmax(rmIC), 'desc');
        rmIC = rmIC(ord(1:max_rm_ic));
    end

    fprintf('Auto ocular IC removal: ICLabel=%d, Corr=%d, Union=%d, Keep(MAX%d)=%d\n', ...
        numel(rm_iclabel), numel(rm_corr), numel(rmIC_all), max_rm_ic, numel(rmIC));

    if ~isempty(rmIC)
        figure;
        try
            pop_topoplot(EEG, 0, rmIC, sprintf('Removed ICs: %s', filenam), 0, 'electrodes','on');
            outPng = fullfile(outputdir, sprintf('%s_removedIC_topo.png', filenam));
            exportgraphics(gcf, outPng, 'Resolution', 200);
        catch ME
            warning('Failed to save topoplot for removed ICs: %s', ME.message);
        end
        close(gcf);
    end

    if ~isempty(rmIC)
        EEG = pop_subcomp(EEG, rmIC, 0);
        EEG = eeg_checkset(EEG);
    end

    EEG = pop_saveset(EEG, 'filename', [filenam '_ica_clean.set'], 'filepath', outputdir);
end

save(fullfile(datadir,'badchan_psd_auto.mat'),'badchan');
fprintf('\nAll done. Output saved to: %s\n', datadir);