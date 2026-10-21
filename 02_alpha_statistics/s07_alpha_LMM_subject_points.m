%% Alpha VOFF LMM with random-effects comparison and Type-III fixed-effect tests
% Input:
%   1) alpha_VOFF_block_level_0p3_2p3.csv
%      Required information: Subject, Block, alpha-power variable
%   2) alpha_rest_baseline_2s.csv  (optional but recommended)
%      Required information: Subject, AlphaPower_baseline
%
% Model:
%   With baseline:
%     AlphaPower ~ Group*Block + AlphaPower_baseline + (1|Subject)
%     AlphaPower ~ Group*Block + AlphaPower_baseline + (1+Block|Subject)
%   Without baseline:
%     AlphaPower ~ Group*Block + (1|Subject)
%     AlphaPower ~ Group*Block + (1+Block|Subject)
%
% Random-effects comparison:
%   RI vs RS fitted with ML, same fixed-effects structure.
%   If RS improves model fit and converges, retain RS; otherwise retain RI.
%
% Final model:
%   Refit selected model with REML.
%   Fixed effects assessed using Type-III / marginal ANOVA with Satterthwaite DF.

clear; clc;

%% ===================== User settings =====================
base_dir = 'D:\online_learning\training_psychophysiology';

alpha_block_csv = fullfile(base_dir, 'alpha_VOFF_block_level_0p3_2p3.csv');
rest_base_csv   = fullfile(base_dir, 'alpha_rest_baseline_2s.csv');

out_dir = fullfile(base_dir, 'alpha_VOFF_LMM_TypeIII_randomEffects');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

summary_txt = fullfile(out_dir, 'alpha_VOFF_LMM_TypeIII_randomEffects_summary.txt');

% Group coding. Edit these IDs if your grouping is different.
group1_ids = [1:6, 14:20, 28:29, 33];
group1_label = 'Offline';
group2_label = 'Online';

% Reference category order. The first category is the reference-like level for coefficients.
% Effects coding is used for Type-III/marginal tests.
group_category_order = {group2_label, group1_label};
block_category_order = {'B1','B2','B3'};

% Candidate outcome names in the block-level alpha CSV.
outcome_candidates = { ...
    'AlphaPower', ...
    'AlphaPower_log10', ...
    'AlphaPower_ROI', ...
    'AlphaPower_mean', ...
    'AlphaPower_block_mean', ...
    'AlphaPower_BlockMean', ...
    'MeanAlphaPower', ...
    'ROI_AlphaPower', ...
    'Alpha_log10'};

baseline_candidates = { ...
    'AlphaPower_baseline', ...
    'AlphaBaseline', ...
    'AlphaPower_rest', ...
    'RestAlphaPower', ...
    'AlphaPower'};

use_satterthwaite = true;
alpha_level = 0.05;

% Use random-slope model only when the likelihood-ratio test is significant.
random_slope_p_threshold = 0.05;

% Plotting
make_plots = true;

% Raw alpha-power figure settings
% The current alpha pipeline uses log10-transformed band power, so the
% default unit is log10(microvolt squared). If AlphaPower is a mean PSD
% rather than integrated band power, change this to:
%   '\alpha power (log_{10}(\muV^2/Hz))'
alpha_y_axis_label = '\alpha power (log_{10}(\muV^2))';

%% ===================== Read block-level alpha table =====================
if ~exist(alpha_block_csv, 'file')
    error('Cannot find alpha block-level CSV:\n%s', alpha_block_csv);
end

Traw = readtable_compat(alpha_block_csv);
Traw.Properties.VariableNames = matlab.lang.makeValidName(Traw.Properties.VariableNames);

% Identify variables.
subj_var = find_var(Traw, {'Subject','SubjID','Participant','ID','File','Filename'});
block_var = find_var(Traw, {'Block','TrainingBlock','block'});
outcome_var = find_var(Traw, outcome_candidates);

if isempty(subj_var)
    error('No Subject-like variable found in %s.', alpha_block_csv);
end
if isempty(block_var)
    error('No Block-like variable found in %s.', alpha_block_csv);
end
if isempty(outcome_var)
    outcome_var = infer_alpha_outcome_var(Traw);
end
if isempty(outcome_var)
    error('Could not identify alpha outcome variable. Please rename your outcome column to AlphaPower.');
end

T = table();
T.SubjectRaw = string(Traw.(subj_var));
T.SubjectID  = extract_subject_id(T.SubjectRaw);
T.Subject    = categorical(string(T.SubjectID));
T.Block      = normalize_block(Traw.(block_var), block_category_order);
T.AlphaPower = double(Traw.(outcome_var));

% Remove rows without key information.
valid = ~isnan(T.SubjectID) & ~isundefined(T.Subject) & ~isundefined(T.Block) & isfinite(T.AlphaPower);
T = T(valid, :);

% Rebuild Group.
grp = strings(height(T), 1);
for i = 1:height(T)
    if ismember(T.SubjectID(i), group1_ids)
        grp(i) = group1_label;
    else
        grp(i) = group2_label;
    end
end
T.Group = categorical(grp, group_category_order, 'Ordinal', false);

%% ===================== Merge resting alpha baseline if available =====================
has_baseline = false;
if exist(rest_base_csv, 'file')
    Braw = readtable_compat(rest_base_csv);
    Braw.Properties.VariableNames = matlab.lang.makeValidName(Braw.Properties.VariableNames);

    b_subj_var = find_var(Braw, {'Subject','SubjID','Participant','ID','File','Filename'});
    b_alpha_var = find_var(Braw, baseline_candidates);

    if ~isempty(b_subj_var) && ~isempty(b_alpha_var)
        B = table();
        B.SubjectRaw = string(Braw.(b_subj_var));
        B.SubjectID = extract_subject_id(B.SubjectRaw);
        B.AlphaPower_baseline = double(Braw.(b_alpha_var));
        B = B(~isnan(B.SubjectID) & isfinite(B.AlphaPower_baseline), :);

        % Average if multiple baseline rows per subject.
        if ~isempty(B)
            Bmean = groupsummary(B, 'SubjectID', 'mean', 'AlphaPower_baseline');
            Bmean.Properties.VariableNames{'mean_AlphaPower_baseline'} = 'AlphaPower_baseline';
            Bmean = Bmean(:, {'SubjectID','AlphaPower_baseline'});

            T = outerjoin(T, Bmean, 'Keys', 'SubjectID', 'MergeKeys', true, 'Type', 'left');
            has_baseline = any(isfinite(T.AlphaPower_baseline));
        end
    else
        warning('Baseline CSV found but Subject or baseline alpha variable could not be identified. Continue without baseline.');
    end
else
    warning('Baseline CSV not found. Continue without AlphaPower_baseline covariate:\n%s', rest_base_csv);
end

% Keep only complete rows for the final model.
if has_baseline
    Tuse = T(isfinite(T.AlphaPower) & isfinite(T.AlphaPower_baseline) & ...
        ~isundefined(T.Group) & ~isundefined(T.Block) & ~isundefined(T.Subject), :);
else
    Tuse = T(isfinite(T.AlphaPower) & ...
        ~isundefined(T.Group) & ~isundefined(T.Block) & ~isundefined(T.Subject), :);
end

% Enforce category orders after subsetting.
Tuse.Group = categorical(string(Tuse.Group), group_category_order, 'Ordinal', false);
Tuse.Block = categorical(string(Tuse.Block), block_category_order, 'Ordinal', false);
Tuse.Subject = categorical(string(Tuse.Subject));

if height(Tuse) == 0
    error('No usable rows after missing-value handling.');
end

%% ===================== Open summary file =====================
fid = fopen(summary_txt, 'w');
if fid < 0
    error('Cannot open summary file for writing:\n%s', summary_txt);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf_both(fid, 'ALPHA VOFF LMM TYPE-III SUMMARY\n');
fprintf_both(fid, 'Generated: %s\n', datestr(now));
fprintf_both(fid, 'Input block CSV: %s\n', alpha_block_csv);
if exist(rest_base_csv, 'file')
    fprintf_both(fid, 'Input baseline CSV: %s\n', rest_base_csv);
else
    fprintf_both(fid, 'Input baseline CSV: NOT FOUND\n');
end
fprintf_both(fid, '\nOutcome variable in block CSV: %s\n', outcome_var);
fprintf_both(fid, 'Baseline covariate used: %d\n', has_baseline);
fprintf_both(fid, 'Rows used: %d\n', height(Tuse));
fprintf_both(fid, 'Subjects used: %d\n', numel(categories(Tuse.Subject)));

fprintf_both(fid, '\nCounts by Group x Block:\n');
disp_to_file(fid, groupsummary(Tuse, {'Group','Block'}));

fprintf_both(fid, '\nDescriptive statistics by Group x Block:\n');
D = groupsummary(Tuse, {'Group','Block'}, {'mean','std','numel'}, 'AlphaPower');
disp_to_file(fid, D);
writetable(D, fullfile(out_dir, 'alpha_descriptives_by_Group_Block.csv'));

if has_baseline
    fprintf_both(fid, '\nResting alpha baseline by Group:\n');
    Bdesc = groupsummary(Tuse, 'Group', {'mean','std','numel'}, 'AlphaPower_baseline');
    disp_to_file(fid, Bdesc);
    writetable(Bdesc, fullfile(out_dir, 'alpha_baseline_descriptives_by_Group.csv'));
end

%% ===================== Define models =====================
if has_baseline
    formula_RI = 'AlphaPower ~ Group*Block + AlphaPower_baseline + (1|Subject)';
    formula_RS = 'AlphaPower ~ Group*Block + AlphaPower_baseline + (1+Block|Subject)';
else
    formula_RI = 'AlphaPower ~ Group*Block + (1|Subject)';
    formula_RS = 'AlphaPower ~ Group*Block + (1+Block|Subject)';
end

fprintf_both(fid, '\n====================================================\n');
fprintf_both(fid, 'RANDOM-EFFECTS STRUCTURE COMPARISON\n');
fprintf_both(fid, '====================================================\n');
fprintf_both(fid, 'Random-intercept formula:\n%s\n', formula_RI);
fprintf_both(fid, 'Random-slope formula:\n%s\n', formula_RS);

%% ===================== Fit ML models for random-effects comparison =====================
use_RS = false;
RS_fit_ok = false;
random_p = NaN;
final_formula = formula_RI;
random_decision = 'Keep random-intercept model.';

try
    lme_RI_ML = fitlme(Tuse, formula_RI, ...
        'FitMethod', 'ML', ...
        'DummyVarCoding', 'effects');

    fprintf_both(fid, '\nRI model fit statistics, ML:\n');
    disp_to_file(fid, lme_RI_ML.ModelCriterion);

catch ME
    error('Random-intercept ML model failed: %s', ME.message);
end

try
    lme_RS_ML = fitlme(Tuse, formula_RS, ...
        'FitMethod', 'ML', ...
        'DummyVarCoding', 'effects');
    RS_fit_ok = true;

    fprintf_both(fid, '\nRS model fit statistics, ML:\n');
    disp_to_file(fid, lme_RS_ML.ModelCriterion);

    cmp_random = compare(lme_RI_ML, lme_RS_ML);
    fprintf_both(fid, '\nRandom-effects comparison RI vs RS, ML:\n');
    disp_to_file(fid, cmp_random);

    random_p = get_compare_pvalue(cmp_random);
    fprintf_both(fid, '\nExtracted random-effects comparison p-value = %.6g\n', random_p);

    if isfinite(random_p) && random_p < random_slope_p_threshold
        use_RS = true;
        final_formula = formula_RS;
        random_decision = 'Use random-slope model.';
    else
        use_RS = false;
        final_formula = formula_RI;
        random_decision = 'Keep random-intercept model.';
    end

catch ME
    RS_fit_ok = false;
    use_RS = false;
    final_formula = formula_RI;
    random_decision = sprintf('Random-slope model failed; keep random-intercept model. Error: %s', ME.message);
    fprintf_both(fid, '\nRS model failed or comparison failed:\n%s\n', ME.message);
end

fprintf_both(fid, '\nDecision: %s\n', random_decision);
fprintf_both(fid, 'Final selected formula:\n%s\n', final_formula);

%% ===================== Fit final REML model =====================
fprintf_both(fid, '\n====================================================\n');
fprintf_both(fid, 'FINAL MODEL, REML\n');
fprintf_both(fid, '====================================================\n');

lme_FINAL = fitlme(Tuse, final_formula, ...
    'FitMethod', 'REML', ...
    'DummyVarCoding', 'effects');

fprintf_both(fid, '\nFinal model summary, REML:\n');
disp_to_file(fid, lme_FINAL);

fprintf_both(fid, '\nFinal model fit statistics, REML:\n');
disp_to_file(fid, lme_FINAL.ModelCriterion);

%% ===================== Type-III / marginal ANOVA fixed-effect tests =====================
fprintf_both(fid, '\nType-III / marginal ANOVA fixed-effect tests on final model:\n');
try
    if use_satterthwaite
        a_final = anova(lme_FINAL, 'DFMethod', 'Satterthwaite');
        fprintf_both(fid, 'DFMethod: Satterthwaite\n');
    else
        a_final = anova(lme_FINAL);
        fprintf_both(fid, 'DFMethod: default\n');
    end
catch ME
    warning('Satterthwaite/default anova failed with message: %s. Trying default anova(lme_FINAL).', ME.message);
    a_final = anova(lme_FINAL);
    fprintf_both(fid, 'DFMethod: default fallback because requested ANOVA failed.\n');
end

disp_to_file(fid, a_final);
writetable(table_compat(a_final), fullfile(out_dir, 'alpha_final_TypeIII_ANOVA.csv'));

fprintf_both(fid, '\nFixed-effect coefficients from final model:\n');
disp_to_file(fid, lme_FINAL.Coefficients);
writetable(table_compat(lme_FINAL.Coefficients), fullfile(out_dir, 'alpha_final_fixed_coefficients.csv'));

%% ===================== Model-adjusted means =====================
fprintf_both(fid, '\n====================================================\n');
fprintf_both(fid, 'MODEL-ADJUSTED MEANS BY GROUP x BLOCK\n');
fprintf_both(fid, '====================================================\n');

try
    adjT = make_adjusted_means_table(Tuse, lme_FINAL, has_baseline, group_category_order, block_category_order, alpha_level);
    disp_to_file(fid, adjT);
    writetable(adjT, fullfile(out_dir, 'alpha_model_adjusted_means.csv'));
catch ME
    fprintf_both(fid, 'Model-adjusted means failed. Error:\n%s\n', ME.message);
    adjT = table();
end

%% ===================== Save analysis-ready table =====================
writetable(Tuse, fullfile(out_dir, 'alpha_LMM_analysis_table.csv'));

%% ===================== Plots =====================
if make_plots
    try
        plot_raw_means( ...
            Tuse, ...
            fullfile(out_dir, 'alpha_raw_mean_sem_by_Group_Block_600dpi.tif'), ...
            alpha_y_axis_label);
    catch ME
        fprintf_both(fid, 'Raw plot failed: %s\n', ME.message);
    end

    if exist('adjT', 'var') && ~isempty(adjT)
        try
            plot_adjusted_means(adjT, fullfile(out_dir, 'alpha_model_adjusted_means.png'));
        catch ME
            fprintf_both(fid, 'Adjusted-means plot failed: %s\n', ME.message);
        end
    end
end

fprintf_both(fid, '\nAll done. Output directory:\n%s\n', out_dir);

%% ===================== Local functions =====================
function T = readtable_compat(file)
    try
        T = readtable(file, 'VariableNamingRule', 'preserve');
    catch
        T = readtable(file);
    end
end

function v = find_var(T, candidates)
    v = '';
    names = T.Properties.VariableNames;
    names_low = lower(names);
    for i = 1:numel(candidates)
        c = matlab.lang.makeValidName(candidates{i});
        idx = find(strcmpi(names, c), 1);
        if ~isempty(idx)
            v = names{idx};
            return;
        end
        idx = find(strcmp(names_low, lower(c)), 1);
        if ~isempty(idx)
            v = names{idx};
            return;
        end
    end
end

function v = infer_alpha_outcome_var(T)
    v = '';
    names = T.Properties.VariableNames;
    for i = 1:numel(names)
        nm = lower(names{i});
        if contains(nm, 'alpha') && ~contains(nm, 'baseline') && ...
                ~contains(nm, 'std') && ~contains(nm, 'sem') && ...
                ~contains(nm, 'ntrial') && ~contains(nm, 'count') && ...
                isnumeric(T.(names{i}))
            v = names{i};
            return;
        end
    end
end

function id = extract_subject_id(x)
    x = string(x);
    id = nan(numel(x), 1);
    for i = 1:numel(x)
        tok = regexp(x(i), '\d+', 'match');
        if ~isempty(tok)
            % Use the first number in the subject/file string.
            id(i) = str2double(tok{1});
        end
    end
end

function B = normalize_block(x, block_order)
    if isnumeric(x)
        s = strings(numel(x), 1);
        for i = 1:numel(x)
            if isfinite(x(i))
                s(i) = "B" + string(round(x(i)));
            else
                s(i) = "";
            end
        end
    else
        s = string(x);
        s = strtrim(s);
        for i = 1:numel(s)
            if strlength(s(i)) == 0
                continue;
            end
            % Convert '1' -> 'B1', 'Block1' -> 'B1', 'block_2' -> 'B2'.
            tok = regexp(s(i), '\d+', 'match');
            if ~isempty(tok)
                s(i) = "B" + string(str2double(tok{1}));
            else
                s(i) = upper(s(i));
            end
        end
    end
    B = categorical(s, block_order, 'Ordinal', true);
end

function fprintf_both(fid, varargin)
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end

function disp_to_file(fid, obj)
    txt = evalc('disp(obj)');
    fprintf('%s\n', txt);
    fprintf(fid, '%s\n', txt);
end

function p = get_compare_pvalue(cmp)
    p = NaN;
    try
        if istable(cmp)
            if any(strcmp(cmp.Properties.VariableNames, 'pValue'))
                p = cmp.pValue(end);
                return;
            end
        end
    catch
    end

    try
        % Works for MATLAB dataset arrays and many table-like outputs.
        pvals = cmp.pValue;
        p = pvals(end);
        return;
    catch
    end

    try
        % Fallback: convert display text and extract the last numeric value on the last model-comparison row.
        txt = evalc('disp(cmp)');
        lines = regexp(txt, '\n', 'split');
        nums_last = [];
        for i = numel(lines):-1:1
            nums = regexp(lines{i}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
            if numel(nums) >= 2
                nums_last = nums;
                break;
            end
        end
        if ~isempty(nums_last)
            p = str2double(nums_last{end});
        end
    catch
        p = NaN;
    end
end

function Tout = table_compat(x)
    if istable(x)
        Tout = x;
    else
        try
            Tout = dataset2table(x);
        catch
            try
                Tout = struct2table(x);
            catch
                error('Cannot convert object to table for CSV writing.');
            end
        end
    end
end

function adjT = make_adjusted_means_table(Tuse, lme, has_baseline, group_order, block_order, alpha_level) %#ok<INUSD>
    % Build a prediction table whose categorical variables exactly match the
    % model table. This avoids MATLAB errors such as:
    % "ISMEMBER categorical inputs must both be ordinal or both be unordered."

    group_cats = categories(Tuse.Group);
    block_cats = categories(Tuse.Block);
    subj_cats  = categories(Tuse.Subject);

    [G, B] = ndgrid(string(group_cats), string(block_cats));
    n = numel(G);

    adjT = table();

    % IMPORTANT: match both categories and ordinal status to Tuse.
    adjT.Group = categorical(reshape(G, n, 1), group_cats, ...
        'Ordinal', isordinal(Tuse.Group));
    adjT.Block = categorical(reshape(B, n, 1), block_cats, ...
        'Ordinal', isordinal(Tuse.Block));

    % For population-level prediction, Conditional=false ignores the subject
    % random effect, but the variable must still exist and have matching categories.
    adjT.Subject = categorical(repmat(string(subj_cats{1}), n, 1), subj_cats, ...
        'Ordinal', isordinal(Tuse.Subject));

    if has_baseline
        adjT.AlphaPower_baseline = repmat(mean(Tuse.AlphaPower_baseline, 'omitnan'), n, 1);
    end

    try
        [ypred, yci] = predict(lme, adjT, 'Conditional', false, 'Alpha', alpha_level);
        adjT.AdjustedMean = ypred;
        adjT.CI_Lower = yci(:,1);
        adjT.CI_Upper = yci(:,2);
    catch
        ypred = predict(lme, adjT, 'Conditional', false);
        adjT.AdjustedMean = ypred;
        adjT.CI_Lower = nan(n,1);
        adjT.CI_Upper = nan(n,1);
    end
end

function plot_raw_means(Tuse, out_tif, y_axis_label)
    % Publication-style raw alpha-power plot matched to the production plot:
    %   - one point per participant and block
    %   - no participant trajectory lines
    %   - group mean connected across blocks
    %   - mean +/- between-subject SEM
    %   - Nature-style colors
    %   - no grid
    %   - 600 dpi TIFF output

    %% ---------- Style ----------
    font_name = 'Times New Roman';
    font_size = 14;
    panel_label_size = 16;

    fig_width_cm = 11.67;
    fig_height_cm = 8.33;

    color_offline = [60, 84, 136] ./ 255;  % #3C5488
    color_online  = [230, 75, 53] ./ 255;  % #E64B35

    group_order = {'Offline', 'Online'};
    group_colors = [color_offline; color_online];

    block_order = {'B1','B2','B3'};
    block_labels = {'Block 1','Block 2','Block 3'};

    group_x_offset = 0.10;
    jitter_width = 0.05;

    subject_marker_size = 23;
    mean_marker_size = 4;
    mean_line_width = 1.8;

    %% ---------- Aggregate to Subject x Group x Block ----------
    % This avoids treating duplicate rows from the same participant/block as
    % independent observations.
    [Gsubj, subj_col, group_col, block_col] = findgroups( ...
        string(Tuse.Subject), ...
        string(Tuse.Group), ...
        string(Tuse.Block));

    alpha_subject = splitapply(@mean, Tuse.AlphaPower, Gsubj);

    subjectPlot = table( ...
        subj_col, group_col, block_col, alpha_subject, ...
        'VariableNames', {'Subject','Group','Block','AlphaPower'});

    block_num = nan(height(subjectPlot),1);
    for b = 1:numel(block_order)
        block_num(strcmpi(subjectPlot.Block, block_order{b})) = b;
    end
    subjectPlot.BlockNum = block_num;
    subjectPlot = subjectPlot(isfinite(subjectPlot.AlphaPower) & ...
        isfinite(subjectPlot.BlockNum), :);

    %% ---------- Group mean and between-subject SEM ----------
    [Gsum, summary_group, summary_block, summary_block_num] = findgroups( ...
        subjectPlot.Group, subjectPlot.Block, subjectPlot.BlockNum);

    mean_alpha = splitapply(@mean, subjectPlot.AlphaPower, Gsum);
    sd_alpha   = splitapply(@std,  subjectPlot.AlphaPower, Gsum);
    n_subjects = splitapply(@numel, subjectPlot.AlphaPower, Gsum);
    sem_alpha  = sd_alpha ./ sqrt(n_subjects);

    summaryPlot = table( ...
        summary_group, summary_block, summary_block_num, ...
        mean_alpha, sd_alpha, sem_alpha, n_subjects, ...
        'VariableNames', ...
        {'Group','Block','BlockNum','MeanAlpha','SD','SEM','N'});

    writetable(subjectPlot, ...
        fullfile(fileparts(out_tif), 'alpha_raw_subject_points.csv'));
    writetable(summaryPlot, ...
        fullfile(fileparts(out_tif), 'alpha_raw_group_mean_SEM.csv'));

    %% ---------- Dynamic y limits ----------
    y_all = [ ...
        subjectPlot.AlphaPower; ...
        summaryPlot.MeanAlpha - summaryPlot.SEM; ...
        summaryPlot.MeanAlpha + summaryPlot.SEM];

    y_min_data = min(y_all, [], 'omitnan');
    y_max_data = max(y_all, [], 'omitnan');
    y_range = y_max_data - y_min_data;

    if ~isfinite(y_range) || y_range <= 0
        y_range = max(abs(y_max_data), 1);
    end

    y_min = y_min_data - 0.10 .* y_range;
    y_max = y_max_data + 0.30 .* y_range;

    %% ---------- Figure ----------
    fig = figure('Color','w');
    set(fig, ...
        'Units','centimeters', ...
        'Position',[2, 2, fig_width_cm, fig_height_cm]);

    ax = axes(fig);
    hold(ax,'on');

    legend_handles = gobjects(numel(group_order),1);
    group_shifts = [-group_x_offset, group_x_offset];

    for g = 1:numel(group_order)
        this_group = string(group_order{g});
        this_color = group_colors(g,:);
        this_shift = group_shifts(g);

        %% Individual participant points
        group_data = subjectPlot(strcmpi(subjectPlot.Group, this_group), :);
        group_data = sortrows(group_data, {'Subject','BlockNum'});

        subjects = unique(group_data.Subject, 'stable');
        n = numel(subjects);

        if n <= 1
            subject_jitter = 0;
        else
            subject_jitter = linspace(-jitter_width, jitter_width, n);
        end

        light_color = 0.65 .* [1 1 1] + 0.35 .* this_color;

        for s = 1:n
            this_subject = subjects(s);
            tmp = group_data(strcmpi(group_data.Subject, this_subject), :);
            tmp = sortrows(tmp, 'BlockNum');

            scatter(ax, ...
                tmp.BlockNum + this_shift + subject_jitter(s), ...
                tmp.AlphaPower, ...
                subject_marker_size, ...
                'Marker','o', ...
                'MarkerFaceColor',light_color, ...
                'MarkerEdgeColor',this_color, ...
                'LineWidth',0.6, ...
                'HandleVisibility','off');
        end

        %% Group mean +/- SEM
        mean_data = summaryPlot(strcmpi(summaryPlot.Group, this_group), :);
        mean_data = sortrows(mean_data, 'BlockNum');

        x_mean = mean_data.BlockNum + this_shift;

        legend_handles(g) = errorbar(ax, ...
            x_mean, ...
            mean_data.MeanAlpha, ...
            mean_data.SEM, ...
            '-o', ...
            'Color',this_color, ...
            'MarkerFaceColor',this_color, ...
            'MarkerEdgeColor',this_color, ...
            'MarkerSize',mean_marker_size, ...
            'LineWidth',mean_line_width, ...
            'CapSize',10, ...
            'DisplayName',char(this_group));
    end

    %% ---------- Significance: Offline Block 1 vs Block 3 ----------
    % Place the bracket close to the x-axis. The bracket opens upward,
    % and the asterisk sits tightly above the bracket.
    offline_data = summaryPlot(strcmpi(summaryPlot.Group, 'Offline'), :);
    offline_data = sortrows(offline_data, 'BlockNum');

    if all(ismember([1,3], offline_data.BlockNum))
        x1 = 1 - group_x_offset;
        x2 = 3 - group_x_offset;

        bracket_y = y_min + 0.05 .* (y_max - y_min);
        bracket_h = 0.01 .* (y_max - y_min);

        % Upward-opening bracket: ?
        plot(ax, ...
            [x1, x1, x2, x2], ...
            [bracket_y + bracket_h, bracket_y, bracket_y, bracket_y + bracket_h], ...
            'Color',color_offline, ...
            'LineWidth',1.3, ...
            'HandleVisibility','off');

        % Align the asterisk with the Block 2 x-tick and keep it tight to the bracket.
        text(ax, ...
            2, ...
            bracket_y + bracket_h - 0.09 .* (y_max-y_min), ...
            '*', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontName',font_name, ...
            'FontSize',16, ...
            'FontWeight','bold', ...
            'Color',color_offline, ...
            'Clipping','off');
    end

    %% ---------- Axes ----------
    set(ax, ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'LineWidth',1.0, ...
        'Box','off', ...
        'TickDir','out', ...
        'XColor','k', ...
        'YColor','k');

    xlim(ax,[0.6, 3.4]);
    ylim(ax,[y_min, y_max]);

    xticks(ax,1:3);
    xticklabels(ax,block_labels);

    xlabel(ax,'');
    ylabel(ax, ...
        y_axis_label, ...
        'Interpreter','tex', ...
        'FontName',font_name, ...
        'FontSize',font_size);

    title(ax, ...
        '\alpha power during training', ...
        'Interpreter','tex', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','normal');

    % Panel label b at the top-left corner of the full figure.
    % Use a figure-level annotation so the label is always visible.
    annotation(fig, 'textbox', [0.000, 0.999, 0.03, 0.03], ...
        'String', 'b', ...
        'LineStyle', 'none', ...
        'FitBoxToText', 'on', ...
        'FontName', font_name, ...
        'FontSize', panel_label_size, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    % Remove every grid.
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';

    lgd = legend(ax, ...
        legend_handles, ...
        group_order, ...
        'Location','northwest', ...
        'Box','off', ...
        'FontName',font_name, ...
        'FontSize',font_size);

    % Move the legend slightly downward.
    lgd.Units = 'normalized';
    p = lgd.Position;
    p(2) = min(0.98-p(4), p(2)-0.01);
    lgd.Position = p;

    hold(ax,'off');

    %% ---------- Save only 600 dpi TIFF ----------
    try
        exportgraphics(fig, out_tif, 'Resolution',600);
    catch
        print(fig, out_tif, '-dtiff', '-r600');
    end

    close(fig);
end

function plot_adjusted_means(adjT, out_png)
    % Robust adjusted-means plot using string comparisons to avoid ordered /
    % unordered categorical mismatches.

    groups = string(categories(adjT.Group));
    blocks = string(categories(adjT.Block));

    fig = figure('Color','w'); hold on;
    x = 1:numel(blocks);

    for g = 1:numel(groups)
        y  = nan(numel(blocks),1);
        lo = nan(numel(blocks),1);
        hi = nan(numel(blocks),1);

        for b = 1:numel(blocks)
            idx = string(adjT.Group) == groups(g) & string(adjT.Block) == blocks(b);
            if any(idx)
                y(b)  = adjT.AdjustedMean(idx);
                lo(b) = adjT.CI_Lower(idx);
                hi(b) = adjT.CI_Upper(idx);
            end
        end

        if all(isfinite(lo)) && all(isfinite(hi))
            err = [y - lo, hi - y];
            errorbar(x, y, err(:,1), err(:,2), '-o', 'LineWidth', 1.5);
        else
            plot(x, y, '-o', 'LineWidth', 1.5);
        end
    end

    xlim([0.8 numel(blocks)+0.2]);
    set(gca, 'XTick', x, 'XTickLabel', cellstr(blocks));
    xlabel('Block'); ylabel('Model-adjusted alpha power');
    title('Model-adjusted alpha power');
    legend(cellstr(groups), 'Location','best');
    box off;
    print(fig, out_png, '-dpng', '-r200');
    close(fig);
end
