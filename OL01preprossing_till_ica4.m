% Fully automatic pipeline:
% PSD-only bad channels (with EOG whitelist) + FT filtering + ICA + auto remove ocular ICs (ICLabel + VEOG/HEOG corr)
clear; clc;

%% Settings
resamp_srate = 250;
datadir = {'D:\online_learning\equiprobability_EEG\pre'; ...
           'D:\online_learning\equiprobability_EEG\post'};
filt = '*.vhdr';

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

% ---- Whitelist "hard-bad" upgrade: remove whitelist ONLY if any strong criterion holds ----
flatline_sec = 5;         % (1) flatline >= 5 sec
rep_thr      = 0.20;      % (2) clipping proxy: mean(diff(x)==0) >= 0.20
win_rms_sec  = 1;         % RMS window length (sec) for sustained check
sustain_frac_thr = 0.50;  % (3) sustained: fraction of high-RMS windows >= 0.50

% ---- Filtering (FieldTrip) ----
mains_is_50hz = true; % set false if 60 Hz
lp_freq = 110;
hp_freq = 0.1;

% ---- ICA ----
ica_hp = 1;        % Hz high-pass for ICA
ica_rs = 100;      % Hz resample for ICA

% ---- Auto ocular IC removal ----
use_iclabel = true;      % requires ICLabel plugin
thr_eye_prob = 0.90;     % ICLabel Eye probability threshold
thr_eye_corr = 0.40;     % corr(IC, VEOG/HEOG) threshold

%% Start EEGLAB
eeglab;
close(gcf);

%% Loop
for md = 1:2
    cd(datadir{md});
    outputdir = datadir{md};
    files = dir(filt);

    for curfile = 25%1:length(files)
        file = files(curfile).name;
        EEG = pop_loadbv(pwd, file);
        [~, nam, ext] = fileparts(file);

        if md == 1
            filenam = [nam(1:4), 'pre'];
        else
            filenam = [nam(1:4), 'post'];
        end

        fprintf('\n==============================\n');
        fprintf('Working on %s\n', [filenam ext]);
        fprintf('==============================\n');

        %% Resample
        EEG = pop_resample(EEG, resamp_srate);

        %% Remove non-EEG channels
        EEG = pop_select(EEG, 'nochannel', {'ECG','EDA','ACCx','ACCy','ACCz'});

        %% Muscle artefact rejection (your original)
        [~, rejections, ~] = pop_rejcont(EEG, ...
            'threshold', 10, ...
            'elecrange', 1:EEG.nbchan, ...
            'onlyreturnselection', 'on', ...
            'taper', 'hamming', ...
            'freqlimit', [110 125]);

        xmax(md,curfile,1) = EEG.xmax;
        EEG = pop_select(EEG, 'nopoint', rejections);
        xmax(md,curfile,2) = EEG.xmax;
        xmax(md,curfile,3) = xmax(md,curfile,1) - xmax(md,curfile,2);

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
bad_psd_idx
            % Protect EOG whitelist: do not remove by PSD
            bad_psd_idx = setdiff(bad_psd_idx, eog_idx);
bad_psd_idx
            % ==========================================================
            % Hard-bad check for whitelisted channels: (1) flatline OR (2) clipping OR (3) extreme broadband + sustained
            % ==========================================================
            bad_eog_hard_idx = [];

            if ~isempty(eog_idx)
                % (3a) extreme broadband threshold: q99.9 across ALL channels
                broad_q999 = quantile(P_broad, 0.999);

                % (3b) sustained high RMS threshold: q99 across ALL channels' RMS windows
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

                % (1) flatline
                win_flat  = max(1, round(flatline_sec * EEG_check.srate));
                step_flat = max(1, round(win_flat/5));
                tol_flat  = 1e-6;

                for ii = 1:numel(eog_idx)
                    ch = eog_idx(ii);
                    x  = double(EEG_check.data(ch,:));

                    % ---- (1) flatline >= flatline_sec ----
                    is_flat = false;
                    for s = 1:step_flat:(length(x)-win_flat+1)
                        seg = x(s:s+win_flat-1);
                        if (max(seg) - min(seg)) < tol_flat
                            is_flat = true; break;
                        end
                    end

                    % ---- (2) clipping/saturation proxy: many repeated samples ----
                    rep_frac = mean(diff(x) == 0);
                    is_clip  = (rep_frac >= rep_thr);

                    % ---- (3) extreme broadband + sustained ----
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

                    % ---- Allow removal if ANY strong criterion holds ----
                    if is_flat || is_clip || is_noise_hard
                        bad_eog_hard_idx(end+1,1) = ch; %#ok<SAGROW>
                        fprintf('EOG-hard-bad: %s | flat=%d clip=%d noise=%d (rep=%.2f, P_broad=%.3g)\n', ...
                            EEG_check.chanlocs(ch).labels, is_flat, is_clip, is_noise_hard, rep_frac, P_broad(ch));
                    end
                end
            end

            % Final bad channels = PSD-bad (non-whitelist) + hard-bad(whitelist)
            bad_final_idx = unique([bad_psd_idx(:); bad_eog_hard_idx(:)]);
        else
            bad_final_idx = [];
        end
        bad_psd_idx
    if curfile==25
        bad_final_idx = unique([bad_final_idx(:);8; 126;127]);
    end
        fprintf('PSD bad channels removed: %d\n', numel(bad_final_idx));
        if ~isempty(bad_final_idx)
            disp({EEG.chanlocs(bad_final_idx).labels});
        end

        if ~isempty(bad_final_idx)
            EEG = pop_select(EEG, 'nochannel', bad_final_idx);
            EEG = eeg_checkset(EEG);
            EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
            EEG = eeg_checkset(EEG);
        end

        %% ----------------------------------------------------------
        % 2) Filtering + notch (your FieldTrip block)
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

        %% ----------------------------------------------------------
        % 3) ICA (1 Hz HP + 100 Hz RS) then copy weights back
        %% ----------------------------------------------------------
        EEG_forICA = pop_eegfiltnew(EEG, 'locutoff', ica_hp);
        EEG_forICA = pop_resample(EEG_forICA, ica_rs);
        EEG_forICA = pop_runica(EEG_forICA, 'icatype', 'runica', 'extended', 1, 'interrupt', 'off');

        EEG.icaweights = EEG_forICA.icaweights;
        EEG.icasphere  = EEG_forICA.icasphere;
        EEG.icawinv    = EEG_forICA.icawinv;    % keep consistent
        EEG.icachansind = EEG_forICA.icachansind;

        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG, 'filename', [filenam '_ica.set'], 'filepath', outputdir);

        %% ----------------------------------------------------------
        % 4) Auto remove ocular ICs (ICLabel + VEOG/HEOG corr)
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
                warning('ICLabel not available or failed (%s). Will skip ICLabel step.', ME.message);
                rm_iclabel = [];
            end
        end

        rm_corr = [];
        if ~isempty(VEOG) || ~isempty(HEOG)
            rmax = zeros(size(icaact,1),1);
            for ic = 1:size(icaact,1)
                rlist = [];
                if ~isempty(VEOG)
                    rlist(end+1) = corr(icaact(ic,:)', VEOG');
                end
                if ~isempty(HEOG)
                    rlist(end+1) = corr(icaact(ic,:)', HEOG');
                end
                rmax(ic) = max(abs(rlist));
            end
            rm_corr = find(rmax >= thr_eye_corr);
        end

        rmIC = unique([rm_iclabel(:); rm_corr(:)]);

        fprintf('Auto ocular IC removal: ICLabel=%d, Corr=%d, Union=%d\n', ...
            numel(rm_iclabel), numel(rm_corr), numel(rmIC));

        if ~isempty(rmIC)
            EEG = pop_subcomp(EEG, rmIC, 0);
            EEG = eeg_checkset(EEG);
        end

        %% Save final
        EEG = pop_saveset(EEG, 'filename', [filenam '_ica_clean.set'], 'filepath', outputdir);

    end

    fprintf('Done folder %d/%d\n', md, 2);
end

% save(fullfile(pwd,'_badchan_psd_auto.mat'),'badchan');
