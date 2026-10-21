%% ============================================================
% EDA LMM from NeuroKit2 long table
% Random-effects structure comparison + Type-III fixed-effect ANOVA
%
% Input priority:
%   1) physio_training_long_neurokit_scr_cleaned.csv
%   2) physio_training_long_neurokit.csv
%
% Outcomes:
%   SCL_mean, SCR_rate
%
% Model per outcome:
%   outcome ~ Group*Block + corresponding_baseline + random effect
%
% Random-effects structure:
%   RI: outcome ~ Group*Block + baseline + (1|Subject)
%   RS: outcome ~ Group*Block + baseline + (1+Block|Subject)
%
% Procedure:
%   1) Compare RI vs RS using ML
%   2) Select RS only if it converges and improves fit, p < 0.05
%   3) Refit selected model using REML
%   4) Assess fixed effects using Type-III / marginal ANOVA
%      with effects coding and Satterthwaite df if available
% ============================================================

clear; clc;

%% ---------------- PATHS ----------------
base_dir = 'D:\online_learning\training_psychophysiology\neurokit_out';

clean_csv = fullfile(base_dir, 'physio_training_long_neurokit_scr_cleaned.csv');
raw_csv   = fullfile(base_dir, 'physio_training_long_neurokit.csv');

out_dir = fullfile(base_dir, 'lmm_stats_from_neurokit_logSCRrate');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

out_txt = fullfile(out_dir, 'EDA_LMM_TypeIII_summary.txt');
fig_hr  = fullfile(out_dir, 'SCL_mean_curve_from_NeuroKit_TypeIII.png');
fig_rm  = fullfile(out_dir, 'SCR_rate_curve_from_NeuroKit_TypeIII.png');

%% ---------------- LOAD ----------------
if isfile(clean_csv)
    in_csv = clean_csv;
else
    in_csv = raw_csv;
end

assert(isfile(in_csv), 'Input CSV not found: %s', in_csv);
T = readtable(in_csv);

needCols = {'Subject','Group','Block','SCL_mean','SCR_rate','SCL_baseline','SCR_baseline'};
for i = 1:numel(needCols)
    assert(any(strcmp(T.Properties.VariableNames, needCols{i})), ...
        'Missing required column: %s', needCols{i});
end

% Preserve raw values; ln transforms require strictly positive values.
T.SCR_rate_raw = T.SCR_rate;
T.SCR_baseline_raw = T.SCR_baseline;
assert(all(T.SCR_rate(isfinite(T.SCR_rate)) > 0), 'Nonpositive SCR_rate');
assert(all(T.SCR_baseline(isfinite(T.SCR_baseline)) > 0), 'Nonpositive baseline SCR_rate');
T.SCR_rate = log(T.SCR_rate);
T.SCR_baseline = log(T.SCR_baseline);
%% ---------------- BASIC FORMATTING ----------------
T.Subject = categorical(string(T.Subject));
T.Group   = normalize_group_column(T.Group);
T.Block   = categorical(string(T.Block), {'B1','B2','B3'}, 'Ordinal', true);

% Remove unused categories
T.Subject = removecats(T.Subject);
T.Group   = removecats(T.Group);
T.Block   = removecats(T.Block);

%% ---------------- WRITE HEADER ----------------
fid = fopen(out_txt, 'w');
assert(fid > 0, 'Cannot open output file: %s', out_txt);

fprintf(fid, 'EDA LMM SUMMARY FROM NEUROKIT2 LONG TABLE\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Input CSV: %s\n\n', in_csv);

fprintf(fid, 'SCL_mean: original scale; SCR_rate and SCR baseline: natural logarithm.\n');
fprintf(fid, 'Total rows loaded: %d\n', height(T));
fprintf(fid, 'Unique subjects loaded: %d\n\n', numel(categories(T.Subject)));

fprintf(fid, 'Counts by Group x Block before outcome-specific missing handling:\n');
fprintf(fid, '%s\n', evalc('disp(groupsummary(T, {''Group'',''Block''}))'));

%% ---------------- BASELINE CHECKS ----------------
fprintf(fid, '\n====================================================\n');
fprintf(fid, 'BASELINE GROUP CHECKS\n');
fprintf(fid, '====================================================\n');

baseline_ttest(T, 'SCL_baseline', 'REST SCL_mean baseline', fid);
baseline_ttest(T, 'SCR_baseline', 'REST SCR_rate baseline', fid);

%% ---------------- LMM: SCL_mean ----------------
run_one_lmm_typeIII(T, ...
    'SCL_mean', ...
    'SCL_baseline', ...
    'SCL_mean', ...
    'recorded units', ...
    fid);

%% ---------------- LMM: SCR_rate ----------------
run_one_lmm_typeIII(T, ...
    'SCR_rate', ...
    'SCR_baseline', ...
    'SCR_rate', ...
    'ln(responses/min)', ...
    fid);

fclose(fid);
fprintf('Wrote EDA summary:\n%s\n', out_txt);

%% ---------------- PLOTS ----------------
plot_mean_sem_curve(T, 'SCL_mean', 'Mean SCL (recorded units)', fig_hr);
plot_mean_sem_curve(T, 'SCR_rate', 'ln SCR rate (responses/min)', fig_rm);

fprintf('Saved figures:\n  %s\n  %s\n', fig_hr, fig_rm);
fprintf('DONE.\n');


%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function G = normalize_group_column(x)
    xs = string(x);
    xs = strtrim(xs);

    % Accept both Group1/Group2 and G1/G2 labels
    xs(strcmpi(xs, 'G1')) = "Group1";
    xs(strcmpi(xs, 'G2')) = "Group2";
    xs(strcmpi(xs, '1'))  = "Group1";
    xs(strcmpi(xs, '2'))  = "Group2";

    G = categorical(xs, {'Group1','Group2'});
end

function baseline_ttest(T, baseline_var, label, fid)
    Tsub = unique(T(:, {'Subject','Group', baseline_var}), 'rows');
    y = Tsub.(baseline_var);

    g1 = y(Tsub.Group == 'Group1');
    g2 = y(Tsub.Group == 'Group2');

    g1 = g1(~isnan(g1));
    g2 = g2(~isnan(g2));

    fprintf(fid, '\n--- %s ---\n', label);
    fprintf(fid, 'Group1 N=%d, mean=%.4f, SD=%.4f\n', numel(g1), mean(g1,'omitnan'), std(g1,'omitnan'));
    fprintf(fid, 'Group2 N=%d, mean=%.4f, SD=%.4f\n', numel(g2), mean(g2,'omitnan'), std(g2,'omitnan'));

    if numel(g1) >= 2 && numel(g2) >= 2
        [~, p, ~, stats] = ttest2(g1, g2);
        fprintf(fid, 'Independent-samples t-test: t(%g)=%.4f, p=%.6f\n', stats.df, stats.tstat, p);
    else
        fprintf(fid, 'Independent-samples t-test not run: insufficient data.\n');
    end
end

function run_one_lmm_typeIII(T, outcome_var, baseline_var, label, unit_label, fid)
    fprintf(fid, '\n====================================================\n');
    fprintf(fid, '%s MODEL\n', label);
    fprintf(fid, '====================================================\n');

    % Outcome-specific missing-data handling
    Tuse = T(~isnan(T.(outcome_var)) & ~isnan(T.(baseline_var)), :);
    Tuse.Subject = removecats(Tuse.Subject);
    Tuse.Group   = removecats(Tuse.Group);
    Tuse.Block   = removecats(Tuse.Block);

    fprintf(fid, 'Rows used: %d\n', height(Tuse));
    fprintf(fid, 'Subjects used: %d\n', numel(categories(Tuse.Subject)));
    fprintf(fid, 'Counts by Group x Block:\n');
    fprintf(fid, '%s\n', evalc('disp(groupsummary(Tuse, {''Group'',''Block''}))'));

    if height(Tuse) < 6 || numel(categories(Tuse.Subject)) < 3
        fprintf(fid, 'Model skipped: insufficient data.\n');
        return;
    end

    formula_RI = sprintf('%s ~ Group*Block + %s + (1|Subject)', outcome_var, baseline_var);
    formula_RS = sprintf('%s ~ Group*Block + %s + (1+Block|Subject)', outcome_var, baseline_var);

    fprintf(fid, '\nRandom-intercept formula:\n%s\n', formula_RI);
    fprintf(fid, 'Random-slope formula:\n%s\n', formula_RS);

    % ---------------- Random-effects comparison using ML ----------------
    use_RS = false;
    final_formula = formula_RI;

    try
        lme_RI_ML = fitlme(Tuse, formula_RI, ...
            'FitMethod', 'ML', ...
            'DummyVarCoding', 'effects');

        fprintf(fid, '\nRI model fit statistics, ML:\n');
        fprintf(fid, '%s\n', evalc('disp(lme_RI_ML.ModelCriterion)'));

        try
            lme_RS_ML = fitlme(Tuse, formula_RS, ...
                'FitMethod', 'ML', ...
                'DummyVarCoding', 'effects');

            fprintf(fid, 'RS model fit statistics, ML:\n');
            fprintf(fid, '%s\n', evalc('disp(lme_RS_ML.ModelCriterion)'));

            cmp = compare(lme_RI_ML, lme_RS_ML);
            fprintf(fid, '\nRandom-effects comparison RI vs RS, ML:\n');
            fprintf(fid, '%s\n', evalc('disp(cmp)'));

            p_cmp = extract_compare_p(cmp);
            fprintf(fid, 'Extracted p-value for RI vs RS comparison = %.12g\n', p_cmp);

            if ~isnan(p_cmp) && p_cmp < 0.05
                use_RS = true;
                final_formula = formula_RS;
                fprintf(fid, 'Decision: Use random-slope model.\n');
            else
                fprintf(fid, 'Decision: Keep random-intercept model.\n');
            end

        catch ME_RS
            fprintf(fid, '\nRS model failed or did not converge reliably:\n%s\n', ME_RS.message);
            fprintf(fid, 'Decision: Keep random-intercept model.\n');
        end

    catch ME_RI
        fprintf(fid, 'RI model failed; cannot fit %s model. Error:\n%s\n', label, ME_RI.message);
        return;
    end

    % ---------------- Final model using REML ----------------
    fprintf(fid, '\nFinal selected formula:\n%s\n', final_formula);

    try
        lme_FINAL = fitlme(Tuse, final_formula, ...
            'FitMethod', 'REML', ...
            'DummyVarCoding', 'effects');

        fprintf(fid, '\nFinal model summary, REML:\n');
        fprintf(fid, '%s\n', evalc('disp(lme_FINAL)'));

        fprintf(fid, '\nFinal model fit statistics, REML:\n');
        fprintf(fid, '%s\n', evalc('disp(lme_FINAL.ModelCriterion)'));

        % Conditional residual diagnostics for the FINAL selected model.
        output_folder = fileparts(fopen(fid));
        r = residuals(lme_FINAL, 'Conditional', true);
        yh = fitted(lme_FINAL, 'Conditional', true);
        diagnostic_table = Tuse(:, {'Subject','Group','Block'});
        diagnostic_table.Fitted = yh;
        diagnostic_table.Residual = r;
        writetable(diagnostic_table, fullfile(output_folder, [outcome_var '_final_residuals.csv']));
        fprintf(fid, 'Conditional residual skewness: %.6f\n', skewness(r,0));
        f = figure('Visible','off','Color','w','Position',[100 100 1200 350]);
        subplot(1,3,1); qqplot(r); title('Conditional residual Q-Q');
        subplot(1,3,2); scatter(yh,r); yline(0,'--'); xlabel('Fitted'); ylabel('Residual');
        subplot(1,3,3); scatter(yh,sqrt(abs(r/std(r)))); xlabel('Fitted'); ylabel('sqrt(|standardized residual|)');
        saveas(f, fullfile(output_folder, [outcome_var '_final_diagnostics.png'])); close(f);
        % Type-III / marginal fixed-effect tests
        fprintf(fid, '\nType-III / marginal ANOVA fixed-effect tests on final model:\n');
        try
            a_final = anova(lme_FINAL, 'DFMethod', 'Satterthwaite');
            fprintf(fid, 'DFMethod: Satterthwaite\n');
        catch ME_A
            fprintf(fid, 'Satterthwaite DF failed: %s\n', ME_A.message);
            fprintf(fid, 'Falling back to default anova(lme_FINAL).\n');
            a_final = anova(lme_FINAL);
        end
        fprintf(fid, '%s\n', evalc('disp(a_final)'));

        fprintf(fid, '\nFixed-effect coefficients from final model:\n');
        fprintf(fid, '%s\n', evalc('disp(lme_FINAL.Coefficients)'));

        % Model-adjusted means by Group x Block, baseline fixed at sample mean
        fprintf(fid, '\nModel-adjusted means by Group x Block, baseline fixed at sample mean:\n');
        pred_tbl = adjusted_means_table(lme_FINAL, Tuse, outcome_var, baseline_var);
        fprintf(fid, '%s\n', evalc('disp(pred_tbl)'));

    catch ME_FINAL
        fprintf(fid, 'Final REML model failed for %s. Error:\n%s\n', label, ME_FINAL.message);
    end
end

function p = extract_compare_p(cmp)
    p = NaN;
    if isa(cmp, 'dataset')
        names = cmp.Properties.VarNames;
    else
        names = cmp.Properties.VariableNames;
    end

    possible = {'pValue','pvalue','p','PValue'};
    for i = 1:numel(possible)
        hit = strcmp(names, possible{i});
        if any(hit)
            vals = cmp.(possible{i});
            if numel(vals) >= 2
                p = vals(end);
            elseif numel(vals) == 1
                p = vals(1);
            end
            return;
        end
    end

    % Fallback: search for any variable containing p
    for i = 1:numel(names)
        if contains(lower(names{i}), 'p')
            vals = cmp.(names{i});
            if isnumeric(vals)
                p = vals(end);
                return;
            end
        end
    end
end

function pred_tbl = adjusted_means_table(lme, Tuse, outcome_var, baseline_var)
    groups = categories(Tuse.Group);
    blocks = categories(Tuse.Block);
    b0 = mean(Tuse.(baseline_var), 'omitnan');

    rows = {};
    rr = 0;

    for gi = 1:numel(groups)
        for bi = 1:numel(blocks)
            rr = rr + 1;
            rows{rr,1} = groups{gi}; %#ok<AGROW>
            rows{rr,2} = blocks{bi}; %#ok<AGROW>
        end
    end

    pred_tbl = cell2table(rows, 'VariableNames', {'Group','Block'});
    pred_tbl.Group = categorical(string(pred_tbl.Group), categories(Tuse.Group));
    pred_tbl.Block = categorical(string(pred_tbl.Block), categories(Tuse.Block), 'Ordinal', true);
    pred_tbl.(baseline_var) = repmat(b0, height(pred_tbl), 1);
    subject_levels = categories(Tuse.Subject);
    pred_tbl.Subject = categorical(repmat(string(subject_levels{1}), height(pred_tbl), 1), subject_levels);

    yhat = predict(lme, pred_tbl, 'Conditional', false);
    pred_tbl.Predicted = yhat;
    if strcmp(outcome_var, 'SCR_rate')
        pred_tbl.BackTransformed_responses_per_min = exp(yhat); % geometric-scale prediction, not arithmetic mean
    end
    pred_tbl.Outcome = repmat(string(outcome_var), height(pred_tbl), 1);
end

function plot_mean_sem_curve(T, outcome_var, y_label, out_png)
    Tuse = T(~isnan(T.(outcome_var)), :);
    Tuse.Group = removecats(Tuse.Group);
    Tuse.Block = removecats(Tuse.Block);

    groups = categories(Tuse.Group);
    blocks = categories(Tuse.Block);

    means = nan(numel(groups), numel(blocks));
    sems  = nan(numel(groups), numel(blocks));

    for gi = 1:numel(groups)
        for bi = 1:numel(blocks)
            y = Tuse.(outcome_var)(Tuse.Group == groups{gi} & Tuse.Block == blocks{bi});
            means(gi,bi) = mean(y, 'omitnan');
            sems(gi,bi)  = std(y, 'omitnan') / sqrt(sum(~isnan(y)));
        end
    end

    f = figure('Color','w','Position',[200 200 700 450]); hold on;
    x = 1:numel(blocks);
    for gi = 1:numel(groups)
        errorbar(x, means(gi,:), sems(gi,:), '-o', 'LineWidth', 1.5);
    end
    xticks(x); xticklabels(blocks);
    xlabel('Training block');
    ylabel(y_label);
    legend(groups, 'Location','best');
    title(sprintf('%s across training blocks (mean \\pm SEM)', outcome_var));
    grid on;
    saveas(f, out_png);
    close(f);
end
