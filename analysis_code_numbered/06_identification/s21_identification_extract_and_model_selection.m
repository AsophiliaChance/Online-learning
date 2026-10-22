%% ============================================================
% data_extraction_ACC_RT_fixedFE_randomCompare_v12.m
%
% Purpose
%   Extract identification-test ACC/RT data from pre/post Excel files and
%   fit final ACC and RT models with:
%       fixed effects: group*session + tone + speaker
%       random effects: selected by RI vs RS comparison
%
% Main changes relative to earlier exploratory script
%   1) All fitglme/fitlme calls use DummyVarCoding = 'effects', so
%      anova(model) can be interpreted as Type-III / marginal tests.
%   2) No p-value-based fixed-effect selection. Tone and speaker are always
%      retained as planned control factors.
%   3) RT unit is handled explicitly/automatically and converted to ms before
%      thresholding. RT analysis uses log(RT_ms).
%   4) Random-effects structure is still compared:
%        RI: (1|subj) + (1|item)
%        RS: (1 + session|subj) + (1|item)
%      The simpler RI model is retained unless RS converges and improves fit.
%
% Input folders
%   D:\online_learning\indentification\pre
%   D:\online_learning\indentification\post
%
% Outputs
%   D:\online_learning\indentification\ACC_RT_fixedFE_randomCompare_v12
% ============================================================

clear; clc;

%% ===================== USER SETTINGS =====================
rootDir = 'D:\online_learning\indentification';
preDir  = fullfile(rootDir, 'pre');
postDir = fullfile(rootDir, 'post');

outDir = fullfile(rootDir, 'ACC_RT_fixedFE_randomCompare_v12');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Group definition: IDs in group1_ids are coded as G1; others are G2.
group1_ids = [1:6, 14:20, 28:29, 33];
% group1_ids = [1, 6, 9, 11, 13, 14, 15, 17, 18, 22, 24, 27, 28, 29, 33];

% RT unit handling.
%   'auto' : median RT > 10 is treated as ms; otherwise seconds.
%   's'    : raw RT is seconds.
%   'ms'   : raw RT is milliseconds.
RT_UNIT = 'auto';
RT_MIN_MS = 200;
RT_MAX_MS = 3000;

% Optional within-subject/session trimming on logRT.
DO_MAD_TRIM = true;
MAD_K = 3;

% Random-effect model selection threshold.
ALPHA_RANDOM = 0.05;
ALLOW_FALLBACK_IF_RS_FAIL = true;

% If true, also save model objects to .mat.
SAVE_MODELS = true;
%% ==========================================================

assert(isfolder(preDir),  'pre folder not found: %s', preDir);
assert(isfolder(postDir), 'post folder not found: %s', postDir);

folders  = {preDir, postDir};
sessions = {'pre','post'};

reportFile = fullfile(outDir, 'ACC_RT_fixedFE_randomCompare_v12_report.txt');
diary(reportFile);

fprintf('============================================================\n');
fprintf('ACC/RT fixed-FE final models with random-effect comparison\n');
fprintf('Generated: %s\n', datestr(now));
fprintf('Output directory: %s\n', outDir);
fprintf('============================================================\n\n');

%% ===================== LOAD AND BUILD LONG TABLE =====================
All = table();

tone_template = {};
tone_template_src = '';
speaker_template = {};
speaker_template_src = '';

for ff = 1:numel(folders)
    thisDir = folders{ff};
    sessionLabel = sessions{ff};

    files = dir(fullfile(thisDir, '*.xlsx'));
    files = files(~startsWith({files.name}, '~$'));

    fprintf('[LOAD] session=%s, folder=%s, files=%d\n', sessionLabel, thisDir, numel(files));

    for k = 1:numel(files)
        fname = files(k).name;
        fpath = fullfile(files(k).folder, fname);

        if length(fname) < 4
            warning('Skip file with short filename: %s', fname);
            continue;
        end

        % Original rule: participant ID is characters 3-4 in filename.
        id_str = fname(3:4);
        id_num = str2double(id_str);
        if isnan(id_num)
            warning('Skip file; cannot parse participant ID from characters 3-4: %s', fname);
            continue;
        end

        if ismember(id_num, group1_ids)
            grp = 'G1';
        else
            grp = 'G2';
        end

        opts = detectImportOptions(fpath);
        T = readtable(fpath, opts);

        vnames = T.Properties.VariableNames;
        vlow = lower(vnames);

        % Keep the original flexible matching, but report selected variables.
        idxAcc = find(contains(vlow,'corr') & contains(vlow,'mean'), 1);
        if isempty(idxAcc), idxAcc = find(contains(vlow,'corr_mean'), 1); end

        idxRT = find(contains(vlow,'rt') & contains(vlow,'raw'), 1);
        if isempty(idxRT), idxRT = find(contains(vlow,'rt_raw'), 1); end

        idxSign = find(strcmp(vlow,'sign_123'), 1);
        if isempty(idxSign), idxSign = find(contains(vlow,'sign_123'), 1); end

        idxSyl = find(contains(vlow,'syllable'), 1);
        if isempty(idxSyl)
            idxSyl = 1;
            warning('[speaker] Syllable column not found in %s. Falling back to first column: %s', fname, vnames{idxSyl});
        end

        if isempty(idxAcc) || isempty(idxRT)
            warning('Skip file with missing ACC or RT column: %s', fname);
            fprintf('Available columns:\n'); disp(vnames');
            continue;
        end

        acc = local_to_double(T.(vnames{idxAcc})); acc = acc(:);
        rt  = local_to_double(T.(vnames{idxRT}));  rt  = rt(:);

        validMask = ~(isnan(acc) & isnan(rt));
        acc = acc(validMask);
        rt  = rt(validMask);

        if isempty(acc)
            warning('Skip file; no valid trial rows after filtering: %s', fname);
            continue;
        end

        % Keep previous script's template behavior. This assumes trial order is identical
        % across all files. If order differs, replace this with per-file extraction.
        if isempty(tone_template)
            if isempty(idxSign)
                error('Tone template not found. First valid file must include sign_123. File=%s', fname);
            end
            s = local_to_cellstr(T.(vnames{idxSign}));
            s = s(:);
            s = s(validMask);
            s = local_norm_sign_123(s);
            tone_template = s(:);
            tone_template_src = fname;

            fprintf('[tone] Template length=%d, source=%s\n', numel(tone_template), tone_template_src);
            fprintf('[tone] Counts: 1=%d, 2=%d, 3=%d, empty=%d\n', ...
                sum(strcmp(tone_template,'1')), sum(strcmp(tone_template,'2')), ...
                sum(strcmp(tone_template,'3')), sum(cellfun(@isempty,tone_template)));
        end

        if isempty(speaker_template)
            syl = local_to_cellstr(T.(vnames{idxSyl}));
            syl = syl(:);
            syl = syl(validMask);
            spk = local_norm_speaker_from_syllable(syl);
            speaker_template = spk(:);
            speaker_template_src = fname;

            fprintf('[speaker] Template length=%d, source=%s\n', numel(speaker_template), speaker_template_src);
            fprintf('[speaker] Unique speaker IDs: %s\n', strjoin(unique(string(speaker_template))', ', '));
        end

        n = numel(acc);
        tmp = table();
        tmp.subj    = repmat(id_num, n, 1);
        tmp.group   = repmat({grp}, n, 1);
        tmp.session = repmat({sessionLabel}, n, 1);
        tmp.trial   = (1:n)';
        tmp.item    = tmp.trial;  % trial number used as item index, as in previous script
        tmp.acc     = acc;
        tmp.rt_raw  = rt;

        if numel(tone_template) < n
            toneVec = [tone_template; repmat({''}, n-numel(tone_template), 1)];
        else
            toneVec = tone_template(1:n);
        end
        tmp.tone = toneVec;

        if numel(speaker_template) < n
            spkVec = [speaker_template; repmat({'1'}, n-numel(speaker_template), 1)];
        else
            spkVec = speaker_template(1:n);
        end
        tmp.speaker = spkVec;

        All = [All; tmp]; %#ok<AGROW>
    end
end

if isempty(All)
    error('No valid data loaded.');
end

%% ===================== FORMAT TABLE =====================
All.subj    = categorical(All.subj);
All.group   = categorical(All.group, {'G2','G1'});
All.session = categorical(All.session, {'pre','post'}, 'Ordinal', false);
All.tone    = categorical(All.tone, {'1','2','3'});
All.speaker = categorical(All.speaker);
All.item    = categorical(All.item);
All.trial   = categorical(All.trial);
All.acc     = double(All.acc);
All.rt_raw  = double(All.rt_raw);

All = All(~isundefined(All.group) & ~isundefined(All.session), :);

writetable(All, fullfile(outDir, 'LMM_long_table_withToneSpeaker_itemTrial_v12.csv'));
save(fullfile(outDir, 'LMM_long_table_withToneSpeaker_itemTrial_v12.mat'), 'All');
save(fullfile(outDir, 'tone_speaker_templates_v12.mat'), ...
    'tone_template','tone_template_src','speaker_template','speaker_template_src');

fprintf('\nLoaded rows: %d\n', height(All));
fprintf('Subjects: %d\n', numel(categories(All.subj)));
fprintf('Tone undefined: %d / %d\n', sum(isundefined(All.tone)), height(All));
fprintf('Speaker levels: %s\n', strjoin(cellstr(categories(All.speaker))', ', '));

disp('Counts by group x session:');
disp(groupsummary(All, {'group','session'}));

%% =================================================================
% 1) ACC GLMM
%% =================================================================
fprintf('\n============================================================\n');
fprintf('ACC GLMM\n');
fprintf('============================================================\n');

accTbl = All(~isnan(All.acc) & ~isundefined(All.tone) & ~isundefined(All.speaker), :);
accTbl.acc = double(accTbl.acc);

% Ensure binary ACC values.
badAcc = ~(accTbl.acc == 0 | accTbl.acc == 1);
if any(badAcc)
    warning('[ACC] Non-binary ACC values found. Rounding to nearest 0/1 for GLMM.');
    accTbl.acc = round(accTbl.acc);
    accTbl.acc(accTbl.acc < 0) = 0;
    accTbl.acc(accTbl.acc > 1) = 1;
end

accTbl.group   = removecats(accTbl.group);
accTbl.session = removecats(accTbl.session);
accTbl.tone    = removecats(accTbl.tone);
accTbl.speaker = removecats(accTbl.speaker);
accTbl.subj    = removecats(accTbl.subj);
accTbl.item    = removecats(accTbl.item);

fprintf('[ACC] Usable rows: %d\n', height(accTbl));
fprintf('[ACC] Subjects: %d\n', numel(categories(accTbl.subj)));
disp('[ACC] Observed mean by group x session:');
disp(groupsummary(accTbl, {'group','session'}, 'mean', 'acc'));
disp('[ACC] Observed mean by tone:');
disp(groupsummary(accTbl, 'tone', 'mean', 'acc'));
disp('[ACC] Observed mean by speaker:');
disp(groupsummary(accTbl, 'speaker', 'mean', 'acc'));

% Planned fixed effects: keep tone and speaker; no fixed-effect selection.
fe_acc = 'group*session + tone + speaker';
f_acc_RI = ['acc ~ ' fe_acc ' + (1|subj) + (1|item)'];
f_acc_RS = ['acc ~ ' fe_acc ' + (1 + session|subj) + (1|item)'];

[glmm_acc_final, acc_random_structure, acc_compare_random, acc_final_formula] = ...
    local_fit_ACC_random_compare(accTbl, f_acc_RI, f_acc_RS, ALPHA_RANDOM, ALLOW_FALLBACK_IF_RS_FAIL);

fprintf('\n===== FINAL ACC MODEL =====\n');
fprintf('Chosen random structure: %s\n', acc_random_structure);
fprintf('Final formula: %s\n', acc_final_formula);
disp(glmm_acc_final);

fprintf('\n[ACC] Type-III / marginal Wald tests with effects coding:\n');
acc_anova = anova(glmm_acc_final);
disp(acc_anova);

fprintf('\n[ACC] Fixed-effect coefficients:\n');
disp(glmm_acc_final.Coefficients);

local_safe_writetable(acc_anova, fullfile(outDir, 'ACC_final_ANOVA_TypeIII_marginal_Wald.csv'));
local_safe_writetable(glmm_acc_final.Coefficients, fullfile(outDir, 'ACC_final_fixed_coefficients.csv'));
if ~isempty(acc_compare_random)
    local_safe_writetable(acc_compare_random, fullfile(outDir, 'ACC_random_RI_vs_RS_compare.csv'));
end

%% =================================================================
% 2) RT LMM
%% =================================================================
fprintf('\n============================================================\n');
fprintf('RT LMM\n');
fprintf('============================================================\n');

lmm_rt_final_ML = [];
lmm_rt_final_REML = [];
rt_compare_random = [];
rt_random_structure = '';
rt_final_formula = '';

rtTbl = All(All.acc == 1 & ~isnan(All.rt_raw) & All.rt_raw > 0 & ...
            ~isundefined(All.tone) & ~isundefined(All.speaker), :);
fprintf('[RT] Correct + valid RT + defined tone/speaker rows: %d\n', height(rtTbl));

if height(rtTbl) == 0
    warning('[RT] No usable RT rows. Skipping RT LMM.');
else
    % Detect/convert RT unit to milliseconds.
    [rt_ms, unitDetected] = local_convert_rt_to_ms(rtTbl.rt_raw, RT_UNIT);
    rtTbl.rt_ms = rt_ms;
    fprintf('[RT] RT_UNIT setting: %s; detected/used unit: %s\n', RT_UNIT, unitDetected);
    fprintf('[RT] RT_ms summary before threshold: median=%.2f, min=%.2f, max=%.2f\n', ...
        median(rtTbl.rt_ms, 'omitnan'), min(rtTbl.rt_ms), max(rtTbl.rt_ms));

    rtTbl = rtTbl(rtTbl.rt_ms >= RT_MIN_MS & rtTbl.rt_ms <= RT_MAX_MS, :);
    fprintf('[RT] After threshold [%g, %g] ms: %d rows\n', RT_MIN_MS, RT_MAX_MS, height(rtTbl));

    if height(rtTbl) == 0
        warning('[RT] No rows left after RT threshold. Skipping RT LMM.');
    else
        rtTbl.logRT = log(rtTbl.rt_ms);

        if DO_MAD_TRIM
            G = findgroups(rtTbl.subj, rtTbl.session);
            keep = true(height(rtTbl), 1);
            for g = 1:max(double(G))
                idx = (G == g);
                x = rtTbl.logRT(idx);
                medx = median(x, 'omitnan');
                madx = mad(x, 1);
                if madx == 0 || isnan(madx)
                    continue;
                end
                lo = medx - MAD_K * madx;
                hi = medx + MAD_K * madx;
                keep(idx) = (x >= lo & x <= hi);
            end
            rtTbl = rtTbl(keep, :);
            fprintf('[RT] After MAD trimming within subject x session (K=%g): %d rows\n', MAD_K, height(rtTbl));
        end

        rtTbl.group   = removecats(rtTbl.group);
        rtTbl.session = removecats(rtTbl.session);
        rtTbl.tone    = removecats(rtTbl.tone);
        rtTbl.speaker = removecats(rtTbl.speaker);
        rtTbl.subj    = removecats(rtTbl.subj);
        rtTbl.item    = removecats(rtTbl.item);

        writetable(rtTbl, fullfile(outDir, 'RT_table_used_itemTrial_v12.csv'));
        save(fullfile(outDir, 'RT_table_used_itemTrial_v12.mat'), ...
            'rtTbl','RT_UNIT','unitDetected','RT_MIN_MS','RT_MAX_MS','DO_MAD_TRIM','MAD_K');

        fprintf('[RT] Final rows: %d\n', height(rtTbl));
        fprintf('[RT] Subjects: %d\n', numel(categories(rtTbl.subj)));
        disp('[RT] Observed RT_ms by group x session:');
        disp(groupsummary(rtTbl, {'group','session'}, {'mean','std'}, 'rt_ms'));

        % Planned fixed effects: keep tone and speaker; no fixed-effect selection.
        fe_rt = 'group*session + tone + speaker';
        f_rt_RI = ['logRT ~ ' fe_rt ' + (1|subj) + (1|item)'];
        f_rt_RS = ['logRT ~ ' fe_rt ' + (1 + session|subj) + (1|item)'];

        [lmm_rt_final_ML, lmm_rt_final_REML, rt_random_structure, rt_compare_random, rt_final_formula] = ...
            local_fit_RT_random_compare(rtTbl, f_rt_RI, f_rt_RS, ALPHA_RANDOM, ALLOW_FALLBACK_IF_RS_FAIL);

        fprintf('\n===== FINAL RT MODEL =====\n');
        fprintf('Chosen random structure: %s\n', rt_random_structure);
        fprintf('Final formula: %s\n', rt_final_formula);

        fprintf('\n[RT] Final model summary, ML:\n');
        disp(lmm_rt_final_ML);
        fprintf('\n[RT] Type-III / marginal tests on ML model, effects coding:\n');
        rt_anova_ML = local_lme_anova_satt(lmm_rt_final_ML);
        disp(rt_anova_ML);

        fprintf('\n[RT] Final model summary, REML refit for reporting:\n');
        disp(lmm_rt_final_REML);
        fprintf('\n[RT] Type-III / marginal tests on REML model, effects coding:\n');
        rt_anova_REML = local_lme_anova_satt(lmm_rt_final_REML);
        disp(rt_anova_REML);

        fprintf('\n[RT] Fixed-effect coefficients, REML:\n');
        disp(lmm_rt_final_REML.Coefficients);

        local_safe_writetable(rt_anova_ML, fullfile(outDir, 'RT_final_ANOVA_ML_TypeIII_marginal.csv'));
        local_safe_writetable(rt_anova_REML, fullfile(outDir, 'RT_final_ANOVA_REML_TypeIII_marginal.csv'));
        local_safe_writetable(lmm_rt_final_REML.Coefficients, fullfile(outDir, 'RT_final_REML_fixed_coefficients.csv'));
        if ~isempty(rt_compare_random)
            local_safe_writetable(rt_compare_random, fullfile(outDir, 'RT_random_RI_vs_RS_compare.csv'));
        end
    end
end

%% ===================== SAVE MODEL OBJECTS =====================
if SAVE_MODELS
    save(fullfile(outDir, 'models_ACC_RT_fixedFE_randomCompare_v12.mat'), ...
        'glmm_acc_final','acc_random_structure','acc_compare_random','acc_final_formula', ...
        'lmm_rt_final_ML','lmm_rt_final_REML','rt_random_structure','rt_compare_random','rt_final_formula');
end

fprintf('\n============================================================\n');
fprintf('Done. Report saved to:\n%s\n', reportFile);
fprintf('============================================================\n');

diary off;

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function [mdl_final, chosen_random, cmp_random, final_formula] = local_fit_ACC_random_compare(T, f_RI, f_RS, alpha_random, allow_fallback)
    fprintf('\n[ACC] Random-effects comparison\n');
    fprintf('RI formula: %s\n', f_RI);
    fprintf('RS formula: %s\n', f_RS);

    mdl_RI = fitglme(T, f_RI, ...
        'Distribution','Binomial', ...
        'Link','logit', ...
        'FitMethod','Laplace', ...
        'DummyVarCoding','effects');

    fprintf('[ACC_RI] '); local_print_fitstats(mdl_RI);

    mdl_RS = [];
    rs_ok = false;
    try
        mdl_RS = fitglme(T, f_RS, ...
            'Distribution','Binomial', ...
            'Link','logit', ...
            'FitMethod','Laplace', ...
            'DummyVarCoding','effects');
        rs_ok = true;
        fprintf('[ACC_RS] '); local_print_fitstats(mdl_RS);
    catch ME
        warning('[ACC] RS model failed: %s', ME.message);
        if ~allow_fallback
            rethrow(ME);
        end
    end

    chosen_random = 'RI';
    mdl_final = mdl_RI;
    final_formula = f_RI;
    cmp_random = [];

    if rs_ok
        try
            cmp_random = compare(mdl_RI, mdl_RS);
            disp('[ACC] RI vs RS comparison:');
            disp(cmp_random);
            p = local_extract_compare_p(cmp_random);
            fprintf('[ACC] Extracted RI-vs-RS p = %.12g\n', p);
            if ~isnan(p) && p < alpha_random
                chosen_random = 'RS';
                mdl_final = mdl_RS;
                final_formula = f_RS;
                fprintf('[ACC] Decision: use RS model.\n');
            else
                fprintf('[ACC] Decision: keep RI model.\n');
            end
        catch ME
            warning('[ACC] compare(RI, RS) failed: %s', ME.message);
            fprintf('[ACC] Decision: keep RI model.\n');
        end
    else
        fprintf('[ACC] RS unavailable. Decision: keep RI model.\n');
    end
end

function [mdl_final_ML, mdl_final_REML, chosen_random, cmp_random, final_formula] = local_fit_RT_random_compare(T, f_RI, f_RS, alpha_random, allow_fallback)
    fprintf('\n[RT] Random-effects comparison using ML\n');
    fprintf('RI formula: %s\n', f_RI);
    fprintf('RS formula: %s\n', f_RS);

    mdl_RI_ML = fitlme(T, f_RI, ...
        'FitMethod','ML', ...
        'DummyVarCoding','effects');

    fprintf('[RT_RI_ML] '); local_print_fitstats(mdl_RI_ML);

    mdl_RS_ML = [];
    rs_ok = false;
    try
        mdl_RS_ML = fitlme(T, f_RS, ...
            'FitMethod','ML', ...
            'DummyVarCoding','effects');
        rs_ok = true;
        fprintf('[RT_RS_ML] '); local_print_fitstats(mdl_RS_ML);
    catch ME
        warning('[RT] RS model failed: %s', ME.message);
        if ~allow_fallback
            rethrow(ME);
        end
    end

    chosen_random = 'RI';
    mdl_final_ML = mdl_RI_ML;
    final_formula = f_RI;
    cmp_random = [];

    if rs_ok
        try
            cmp_random = compare(mdl_RI_ML, mdl_RS_ML);
            disp('[RT] RI vs RS comparison:');
            disp(cmp_random);
            p = local_extract_compare_p(cmp_random);
            fprintf('[RT] Extracted RI-vs-RS p = %.12g\n', p);
            if ~isnan(p) && p < alpha_random
                chosen_random = 'RS';
                mdl_final_ML = mdl_RS_ML;
                final_formula = f_RS;
                fprintf('[RT] Decision: use RS model.\n');
            else
                fprintf('[RT] Decision: keep RI model.\n');
            end
        catch ME
            warning('[RT] compare(RI, RS) failed: %s', ME.message);
            fprintf('[RT] Decision: keep RI model.\n');
        end
    else
        fprintf('[RT] RS unavailable. Decision: keep RI model.\n');
    end

    % Refit final selected model using REML for reporting.
    mdl_final_REML = fitlme(T, final_formula, ...
        'FitMethod','REML', ...
        'DummyVarCoding','effects');
end

function [rt_ms, unitDetected] = local_convert_rt_to_ms(rt_raw, unitSetting)
    x = double(rt_raw);
    switch lower(strtrim(unitSetting))
        case {'ms','millisecond','milliseconds'}
            rt_ms = x;
            unitDetected = 'ms';
        case {'s','sec','second','seconds'}
            rt_ms = x * 1000;
            unitDetected = 's';
        case 'auto'
            medx = median(x, 'omitnan');
            if medx > 10
                rt_ms = x;
                unitDetected = 'ms(auto)';
            else
                rt_ms = x * 1000;
                unitDetected = 's(auto)';
            end
        otherwise
            error('Unknown RT_UNIT: %s. Use auto, s, or ms.', unitSetting);
    end
end

function A = local_lme_anova_satt(mdl)
    try
        A = anova(mdl, 'DFMethod','Satterthwaite');
    catch
        A = anova(mdl);
    end
end

function x = local_to_double(x)
    if islogical(x), x = double(x); end
    if iscell(x),    x = str2double(string(x)); end
    if isstring(x) || ischar(x), x = str2double(string(x)); end
    if isnumeric(x), x = double(x); end
end

function c = local_to_cellstr(x)
    if iscell(x)
        c = cellstr(string(x));
    elseif isnumeric(x) || islogical(x)
        c = cellstr(string(double(x)));
    else
        c = cellstr(string(x));
    end
end

function s_out = local_norm_sign_123(s_in)
    s_out = s_in;
    for ii = 1:numel(s_out)
        if isempty(s_out{ii})
            s_out{ii} = '';
            continue;
        end
        ss = strtrim(s_out{ii});

        tok = regexp(ss, '^\[\s*''''\s*([123])\s*''''\s*\]$', 'tokens', 'once');
        if ~isempty(tok)
            s_out{ii} = tok{1};
            continue;
        end

        tok = regexp(ss, '([123])', 'tokens', 'once');
        if ~isempty(tok)
            s_out{ii} = tok{1};
        else
            s_out{ii} = '';
        end
    end
end

function spk_out = local_norm_speaker_from_syllable(syl_in)
    spk_out = cell(size(syl_in));
    for ii = 1:numel(syl_in)
        if isempty(syl_in{ii})
            spk_out{ii} = '1';
            continue;
        end

        ss = strtrim(syl_in{ii});
        ss = strrep(ss, '\\', '/');

        parts = regexp(ss, '/+', 'split');
        base = parts{end};
        base = regexprep(base, '\.[A-Za-z0-9]+$', '');

        tok = regexp(base, '_(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            spk_out{ii} = tok{1};
        else
            spk_out{ii} = '1';
        end
    end
end

function local_print_fitstats(mdl)
    try
        aic = mdl.ModelCriterion.AIC;
        bic = mdl.ModelCriterion.BIC;
    catch
        aic = NaN;
        bic = NaN;
    end
    try
        ll = mdl.LogLikelihood;
    catch
        ll = NaN;
    end
    fprintf('LogLik=%10.3f   AIC=%10.3f   BIC=%10.3f\n', ll, aic, bic);
end

function p = local_extract_compare_p(C)
    p = NaN;
    try
        if istable(C)
            vn = C.Properties.VariableNames;
            idx = find(strcmpi(vn, 'pValue') | strcmpi(vn, 'pvalue'), 1);
            if ~isempty(idx) && height(C) >= 2
                p = C{2, idx};
                return;
            end
        end
    catch
    end

    try
        p = C{2,end};
    catch
        try
            p = C(end,end);
        catch
            p = NaN;
        end
    end

    if iscell(p), p = p{1}; end
    if isstring(p) || ischar(p), p = str2double(string(p)); end
    if ~isnumeric(p) || isempty(p), p = NaN; end
end

function local_safe_writetable(X, outPath)
    try
        if istable(X)
            writetable(X, outPath);
            return;
        end
    catch
    end

    try
        T = dataset2table(X); %#ok<DTSTBL>
        writetable(T, outPath);
        return;
    catch
    end

    try
        T = struct2table(X);
        writetable(T, outPath);
        return;
    catch
    end

    % Last-resort text dump if conversion to table fails.
    [folder, base, ~] = fileparts(outPath);
    txtPath = fullfile(folder, [base '.txt']);
    fid = fopen(txtPath, 'w');
    if fid > 0
        fprintf(fid, '%s\n', evalc('disp(X)'));
        fclose(fid);
        warning('Could not write table to CSV. Wrote text dump instead: %s', txtPath);
    else
        warning('Could not write output: %s', outPath);
    end
end
