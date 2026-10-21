%% ============================================================
% Production accuracy: simple-effects analysis
%
% Outcome:
%   Empirical-logit-transformed production accuracy per
%   subject x session x tone
%
% Full model:
%   logit_acc ~ group*session + tone + (1|subj)
%   logit_acc ~ group*session + tone + (1+session_num|subj)
%
% Follow-up tests after a significant group x session interaction:
%   A) Offline vs Online within each session
%   B) Omnibus Session effect within each group
%   C) All six pairwise Session contrasts within each group
%
% Multiple-comparison correction:
%   - Holm correction across the four within-session group contrasts
%   - Holm correction across the two omnibus within-group Session tests
%   - Holm correction across all twelve within-group pairwise Session contrasts
%
% Important:
%   - All contrasts are tested in the same full model.
%   - Model-predicted production accuracy is calculated exactly as in the
%     main interaction-analysis script: predictions are back-transformed
%     separately for tone2, tone3, and tone4, then averaged on the
%     probability scale.
%   - Simple-effect contrasts remain tested on the empirical-logit scale.
%   - Confidence intervals for the averaged predicted probabilities are
%     obtained by parametric simulation of the fixed effects.
%   - Omnibus tests do not have a single effect estimate or 95% CI because
%     they test several contrasts jointly.
%% ============================================================

clear; clc;

%% ---------------- PATHS ----------------
base_dir  = 'D:\online_learning\production';
inputFile = fullfile(base_dir, 'tone_count_summary.xlsx');

outDir = fullfile(base_dir, 'production_simple_effects_revised');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

outTxt = fullfile(outDir, 'production_simple_effects_revised_report.txt');
assert(isfile(inputFile), 'Missing file: %s', inputFile);

diary(outTxt);

disp('============================================================');
disp('Production simple-effects analysis');
disp(['Generated: ', datestr(now)]);
disp(['Input file: ', inputFile]);
disp(['Output directory: ', outDir]);
disp('============================================================');

%% ---------------- USER SETTINGS ----------------
group1_ids = [1:6, 14:20, 28:29, 33];

group1_label = 'Offline';
group2_label = 'Online';

group_order = {group2_label, group1_label};
test_sessions_keep = {'test1','test2','test3','test4'};
toneVars = {'tone2','tone3','tone4'};

BINOM_N = 15;
alpha_level = 0.05;
random_slope_p_threshold = 0.05;
use_random_slope_candidate = true;

%% ---------------- LOAD DATA ----------------
T = readtable(inputFile, 'TextType', 'string');

reqVars = [{'subject','session'}, toneVars];
for k = 1:numel(reqVars)
    assert(ismember(reqVars{k}, T.Properties.VariableNames), ...
        'Missing required column: %s', reqVars{k});
end

T.subject = string(T.subject);
T.session = lower(strtrim(string(T.session)));

T = T(ismember(T.session, test_sessions_keep), :);
assert(~isempty(T), 'No rows left after filtering test1-test4.');

%% ---------------- SUBJECT ID ----------------
subj_num = nan(height(T),1);

for i = 1:height(T)
    s = strtrim(T.subject(i));
    x = str2double(s);

    if isnan(x)
        tok = regexp(char(s), '\d+', 'match', 'once');
        if ~isempty(tok)
            x = str2double(tok);
        end
    end

    subj_num(i) = x;
end

assert(all(~isnan(subj_num)), ...
    'Some subject IDs could not be converted to numeric.');

T.subj_num = subj_num;
T.subj = categorical(T.subject);

%% ---------------- GROUP ----------------
group = strings(height(T),1);
group(ismember(T.subj_num, group1_ids)) = string(group1_label);
group(group == "") = string(group2_label);

T.group = categorical(group, group_order, 'Ordinal', false);

%% ---------------- SESSION ----------------
T.session = categorical( ...
    T.session, test_sessions_keep, 'Ordinal', false);

session_num = nan(height(T),1);
for i = 1:numel(test_sessions_keep)
    session_num(T.session == test_sessions_keep{i}) = i;
end
T.session_num = session_num;

%% ============================================================
% BUILD LONG-FORM DATA
%% ============================================================
nRows = height(T) * numel(toneVars);

subj_col        = strings(nRows,1);
group_col       = strings(nRows,1);
session_col     = strings(nRows,1);
session_num_col = zeros(nRows,1);
tone_col        = strings(nRows,1);
correct_col     = nan(nRows,1);
nTrials_col     = repmat(BINOM_N, nRows, 1);

row_idx = 0;

for i = 1:height(T)
    for tt = 1:numel(toneVars)
        row_idx = row_idx + 1;
        toneName = toneVars{tt};
        this_count = T.(toneName)(i);

        if iscell(this_count)
            this_count = this_count{1};
        end

        if isstring(this_count) || ischar(this_count)
            this_count = str2double(string(this_count));
        end

        subj_col(row_idx)        = string(T.subj(i));
        group_col(row_idx)       = string(T.group(i));
        session_col(row_idx)     = string(T.session(i));
        session_num_col(row_idx) = T.session_num(i);
        tone_col(row_idx)        = string(toneName);
        correct_col(row_idx)     = double(this_count);
    end
end

D = table;
D.subj = categorical(subj_col);
D.group = categorical(group_col, group_order, 'Ordinal', false);
D.session = categorical(session_col, test_sessions_keep, 'Ordinal', false);
D.session_num = session_num_col;
D.tone = categorical(tone_col, toneVars, 'Ordinal', false);
D.correct_raw = correct_col;
D.nTrials = nTrials_col;

D.correct = round(D.correct_raw);

bad = isnan(D.correct) | isnan(D.nTrials) | ...
      D.correct < 0 | D.nTrials <= 0 | ...
      D.correct > D.nTrials | ...
      abs(D.correct_raw - round(D.correct_raw)) > 1e-8;

if any(bad)
    badFile = fullfile(outDir, 'bad_count_rows.csv');
    safe_writetable(D(bad,:), badFile);
    error('Invalid count rows found. Saved to: %s', badFile);
end

D.acc_raw = D.correct ./ D.nTrials;
D.acc_adj = (D.correct + 0.5) ./ (D.nTrials + 1);
D.logit_acc = log(D.acc_adj ./ (1 - D.acc_adj));

safe_writetable(D, ...
    fullfile(outDir, 'production_simple_effects_long_table.csv'));

fprintf('Rows used: %d\n', height(D));
fprintf('Subjects used: %d\n', numel(categories(D.subj)));

disp('Counts by group x session:');
disp(groupsummary(D, {'group','session'}));

%% ============================================================
% RANDOM-EFFECTS COMPARISON
%% ============================================================
formula_RI = ...
    'logit_acc ~ group*session + tone + (1|subj)';

formula_RS = ...
    'logit_acc ~ group*session + tone + (1+session_num|subj)';

disp(' ');
disp('============================================================');
disp('RANDOM-EFFECTS COMPARISON');
disp('============================================================');

lme_RI_ML = fitlme(D, formula_RI, ...
    'FitMethod', 'ML', ...
    'DummyVarCoding', 'effects');

use_RS = false;
final_formula = formula_RI;
p_random = NaN;

if use_random_slope_candidate
    try
        lme_RS_ML = fitlme(D, formula_RS, ...
            'FitMethod', 'ML', ...
            'DummyVarCoding', 'effects');

        cmp_rand = compare(lme_RI_ML, lme_RS_ML);
        disp(cmp_rand);

        p_random = get_compare_pvalue(cmp_rand);
        fprintf('Extracted RI vs RS p-value = %.6g\n', p_random);

        safe_writetable(cmp_rand, ...
            fullfile(outDir, ...
            'random_effects_comparison_RI_vs_RS.csv'));

        if isfinite(p_random) && p_random < random_slope_p_threshold
            use_RS = true;
            final_formula = formula_RS;
        end

    catch ME
        fprintf('Random-slope model failed: %s\n', ME.message);
    end
end

if use_RS
    disp('Selected model: random intercept + linear session slope.');
else
    disp('Selected model: random intercept only.');
end

disp(['Final formula: ', final_formula]);

%% ============================================================
% FINAL MODEL
%% ============================================================
finalModel = fitlme(D, final_formula, ...
    'FitMethod', 'REML', ...
    'DummyVarCoding', 'effects');

disp(' ');
disp('============================================================');
disp('FINAL REML MODEL');
disp('============================================================');
disp(finalModel);

safe_writetable(finalModel.Coefficients, ...
    fullfile(outDir, ...
    'production_simple_effects_final_coefficients.csv'));

try
    anovaTbl = anova(finalModel, ...
        'DFMethod', 'Satterthwaite');
catch
    anovaTbl = anova(finalModel);
end

safe_writetable(anovaTbl, ...
    fullfile(outDir, ...
    'production_simple_effects_final_ANOVA.csv'));

coef_names = string(finalModel.CoefficientNames(:));
beta = fixedEffects(finalModel);
Vbeta = finalModel.CoefficientCovariance;

%% ============================================================
% MODEL-ADJUSTED MEANS AND 95% CIs
%% ============================================================
adjT = make_adjusted_means_table( ...
    D, finalModel, ...
    group_order, test_sessions_keep, toneVars, alpha_level);

safe_writetable(adjT, ...
    fullfile(outDir, ...
    'production_adjusted_means_group_session.csv'));

disp(' ');
disp('============================================================');
disp('MODEL-PREDICTED PRODUCTION ACCURACY');
disp('Tone-specific probabilities were predicted first and then averaged.');
disp('This matches the main interaction-analysis script.');
disp('CIs for averaged probabilities were obtained by parametric simulation.');
disp('============================================================');
disp(adjT);

%% ============================================================
% A. GROUP SIMPLE EFFECT WITHIN EACH SESSION
%% ============================================================
disp(' ');
disp('============================================================');
disp('GROUP SIMPLE EFFECTS WITHIN SESSION');
disp('Direction: Offline - Online');
disp('============================================================');

group_results = table();

for s = 1:numel(test_sessions_keep)
    session_label = test_sessions_keep{s};

    x_offline = marginal_design_row( ...
        coef_names, group_order, test_sessions_keep, toneVars, ...
        group1_label, session_label);

    x_online = marginal_design_row( ...
        coef_names, group_order, test_sessions_keep, toneVars, ...
        group2_label, session_label);

    H = x_offline - x_online;

    R = single_contrast_result( ...
        finalModel, H, beta, Vbeta, ...
        sprintf('Offline - Online at %s', session_label), ...
        alpha_level);

    R.Session = string(session_label);
    group_results = [group_results; R]; %#ok<AGROW>
end

group_results.p_Holm = holm_adjust(group_results.p_raw);
group_results.Significant_Holm = ...
    group_results.p_Holm < alpha_level;

disp(group_results);

safe_writetable(group_results, ...
    fullfile(outDir, ...
    'production_group_effects_within_session.csv'));

%% ============================================================
% B. OMNIBUS SESSION EFFECT WITHIN EACH GROUP
%% ============================================================
disp(' ');
disp('============================================================');
disp('OMNIBUS SESSION SIMPLE EFFECTS WITHIN GROUP');
disp('Joint contrasts: test2-test1, test3-test1, test4-test1');
disp('No single estimate or CI is reported for this multi-df test.');
disp('============================================================');

omnibus_results = table();

for g = 1:numel(group_order)
    group_label = group_order{g};

    x1 = marginal_design_row( ...
        coef_names, group_order, test_sessions_keep, toneVars, ...
        group_label, 'test1');

    x2 = marginal_design_row( ...
        coef_names, group_order, test_sessions_keep, toneVars, ...
        group_label, 'test2');

    x3 = marginal_design_row( ...
        coef_names, group_order, test_sessions_keep, toneVars, ...
        group_label, 'test3');

    x4 = marginal_design_row( ...
        coef_names, group_order, test_sessions_keep, toneVars, ...
        group_label, 'test4');

    H_joint = [x2-x1; x3-x1; x4-x1];

    [p,F,df1,df2] = run_coef_test(finalModel, H_joint);

    row = table( ...
        string(group_label), F, df1, df2, p, ...
        'VariableNames', ...
        {'Group','F','DF1','DF2','p_raw'});

    omnibus_results = [omnibus_results; row]; %#ok<AGROW>
end

omnibus_results.p_Holm = holm_adjust(omnibus_results.p_raw);
omnibus_results.Significant_Holm = ...
    omnibus_results.p_Holm < alpha_level;

disp(omnibus_results);

safe_writetable(omnibus_results, ...
    fullfile(outDir, ...
    'production_omnibus_session_effects_within_group.csv'));

%% ============================================================
% C. ALL PAIRWISE SESSION EFFECTS WITHIN EACH GROUP
%% ============================================================
disp(' ');
disp('============================================================');
disp('PAIRWISE SESSION SIMPLE EFFECTS WITHIN GROUP');
disp('Holm correction is applied across all 12 contrasts.');
disp('============================================================');

pair_defs = { ...
    'test2','test1','test2 - test1'; ...
    'test3','test1','test3 - test1'; ...
    'test4','test1','test4 - test1'; ...
    'test3','test2','test3 - test2'; ...
    'test4','test2','test4 - test2'; ...
    'test4','test3','test4 - test3'};

pairwise_results = table();

for g = 1:numel(group_order)
    group_label = group_order{g};

    for k = 1:size(pair_defs,1)
        session_a = pair_defs{k,1};
        session_b = pair_defs{k,2};
        label = pair_defs{k,3};

        x_a = marginal_design_row( ...
            coef_names, group_order, test_sessions_keep, toneVars, ...
            group_label, session_a);

        x_b = marginal_design_row( ...
            coef_names, group_order, test_sessions_keep, toneVars, ...
            group_label, session_b);

        H = x_a - x_b;

        R = single_contrast_result( ...
            finalModel, H, beta, Vbeta, ...
            sprintf('%s within %s', label, group_label), ...
            alpha_level);

        R.Group = string(group_label);
        R.SessionContrast = string(label);

        pairwise_results = [pairwise_results; R]; %#ok<AGROW>
    end
end

pairwise_results.p_Holm = holm_adjust(pairwise_results.p_raw);
pairwise_results.Significant_Holm = ...
    pairwise_results.p_Holm < alpha_level;

disp(pairwise_results);

safe_writetable(pairwise_results, ...
    fullfile(outDir, ...
    'production_pairwise_session_effects_within_group.csv'));

%% ============================================================
% SAVE ANALYSIS TABLE
%% ============================================================
safe_writetable(D, ...
    fullfile(outDir, ...
    'production_simple_effects_analysis_table.csv'));

disp(' ');
disp('============================================================');
disp('INTERPRETATION');
disp('============================================================');
disp('- Use p_Holm for inference within each contrast family.');
disp('- Positive Offline-Online estimates indicate higher Offline accuracy.');
disp('- Positive later-minus-earlier estimates indicate improvement.');
disp('- Single contrasts have estimates, SEs, 95% CIs, and p values.');
disp('- Omnibus multi-df tests have F, dfs, and p values, but no single CI.');

disp(' ');
disp('============================================================');
disp('Done.');
disp(['Report saved to: ', outTxt]);
disp('============================================================');

diary off;

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function p = get_compare_pvalue(cmp)
    p = NaN;

    try
        if istable(cmp) && ...
                ismember('pValue', cmp.Properties.VariableNames)
            p = cmp.pValue(end);
            return;
        end
    catch
    end

    try
        pvals = cmp.pValue;
        p = pvals(end);
    catch
    end
end

function x = marginal_design_row( ...
        coef_names, group_order, session_order, tone_order, ...
        group_value, session_value)

    x_all = zeros(numel(tone_order), numel(coef_names));

    for t = 1:numel(tone_order)
        x_all(t,:) = fixed_design_row( ...
            coef_names, ...
            group_order, session_order, tone_order, ...
            group_value, session_value, tone_order{t});
    end

    x = mean(x_all,1);
end

function x = fixed_design_row( ...
        coef_names, group_order, session_order, tone_order, ...
        group_value, session_value, tone_value)

    p = numel(coef_names);
    x = zeros(1,p);

    for j = 1:p
        nm = char(coef_names(j));

        if strcmp(nm,'(Intercept)')
            x(j) = 1;
            continue;
        end

        has_group = contains(nm,'group_');
        has_session = contains(nm,'session_');
        has_tone = contains(nm,'tone_');

        if has_group && has_session
            g_level = extract_factor_level(nm,'group_');
            s_level = extract_factor_level(nm,'session_');

            x(j) = ...
                effects_code_value( ...
                group_value,g_level,group_order) .* ...
                effects_code_value( ...
                session_value,s_level,session_order);

        elseif has_group
            g_level = extract_factor_level(nm,'group_');

            x(j) = effects_code_value( ...
                group_value,g_level,group_order);

        elseif has_session
            s_level = extract_factor_level(nm,'session_');

            x(j) = effects_code_value( ...
                session_value,s_level,session_order);

        elseif has_tone
            t_level = extract_factor_level(nm,'tone_');

            x(j) = effects_code_value( ...
                tone_value,t_level,tone_order);
        end
    end
end

function level = extract_factor_level(coef_name,prefix)
    expr = [regexptranslate('escape',prefix), '([^:]+)'];
    tok = regexp(coef_name,expr,'tokens','once');

    if isempty(tok)
        error('Could not parse coefficient name: %s',coef_name);
    end

    level = tok{1};
end

function value = effects_code_value( ...
        target_level,coefficient_level,category_order)

    target_level = char(string(target_level));
    coefficient_level = char(string(coefficient_level));
    category_order = cellstr(string(category_order));

    last_level = category_order{end};

    if strcmp(target_level,coefficient_level)
        value = 1;
    elseif strcmp(target_level,last_level)
        value = -1;
    else
        value = 0;
    end
end

function R = single_contrast_result( ...
        lme,H,beta,Vbeta,label,alpha_level)

    estimate = H*beta;
    variance = H*Vbeta*H';
    variance = max(variance,0);
    SE = sqrt(variance);

    [p,F,df1,df2] = run_coef_test(lme,H);

    if SE > 0 && isfinite(df2)
        t_value = estimate/SE;
        tcrit = tinv(1-alpha_level/2,df2);
        CI_Lower = estimate-tcrit*SE;
        CI_Upper = estimate+tcrit*SE;
    else
        t_value = NaN;
        CI_Lower = NaN;
        CI_Upper = NaN;
    end

    R = table( ...
        string(label), ...
        estimate,SE,t_value,F,df1,df2, ...
        CI_Lower,CI_Upper,p, ...
        'VariableNames',{ ...
        'Contrast','Estimate','SE','t','F', ...
        'DF1','DF2','CI_Lower','CI_Upper','p_raw'});
end

function [p,F,df1,df2] = run_coef_test(lme,H)
    try
        [p,F,df1,df2] = coefTest( ...
            lme,H,0,'DFMethod','Satterthwaite');
    catch
        try
            [p,F,df1,df2] = coefTest(lme,H,0);
        catch
            p = coefTest(lme,H);
            F = NaN;
            df1 = rank(H);
            df2 = NaN;
        end
    end
end

function p_adj = holm_adjust(p_raw)
    p_raw = double(p_raw(:));
    p_adj = nan(size(p_raw));

    valid = isfinite(p_raw);
    pv = p_raw(valid);

    if isempty(pv)
        return;
    end

    [p_sorted,order] = sort(pv);
    m = numel(p_sorted);

    adjusted_sorted = zeros(m,1);

    for i = 1:m
        adjusted_sorted(i) = ...
            min(1,(m-i+1)*p_sorted(i));
    end

    adjusted_sorted = cummax(adjusted_sorted);

    adjusted = zeros(m,1);
    adjusted(order) = adjusted_sorted;

    p_adj(valid) = adjusted;
end

function adjT = make_adjusted_means_table( ...
        D,lme,group_order,session_order,tone_order,alpha_level)
    % Match the main interaction-analysis script exactly:
    %   1) predict each Group x Session x Tone combination on the logit scale;
    %   2) back-transform each prediction to probability;
    %   3) average the three probabilities across tones.
    %
    % The reported group x session probability is therefore the arithmetic
    % mean of tone-specific predicted probabilities.
    %
    % For the 95% CI of this averaged probability, use a parametric
    % simulation from the fixed-effect coefficient covariance matrix.

    nRows = numel(group_order)*numel(session_order);

    group_col = strings(nRows,1);
    session_col = strings(nRows,1);
    adjusted_probability = nan(nRows,1);
    ci_lower_probability = nan(nRows,1);
    ci_upper_probability = nan(nRows,1);

    beta = fixedEffects(lme);
    Vbeta = lme.CoefficientCovariance;
    coef_names = string(lme.CoefficientNames(:));

    % Reproducible simulation-based CI for the nonlinear average probability.
    rng(20260715);
    nSim = 20000;
    beta_sim = mvnrnd(beta(:)', Vbeta, nSim);

    row = 0;

    for g = 1:numel(group_order)
        for s = 1:numel(session_order)
            row = row+1;

            group_col(row) = string(group_order{g});
            session_col(row) = string(session_order{s});

            tone_prob = nan(numel(tone_order),1);
            tone_prob_sim = nan(nSim,numel(tone_order));

            for t = 1:numel(tone_order)
                newT = table();

                newT.group = categorical( ...
                    string(group_order{g}), ...
                    group_order,'Ordinal',false);

                newT.session = categorical( ...
                    string(session_order{s}), ...
                    session_order,'Ordinal',false);

                newT.tone = categorical( ...
                    string(tone_order{t}), ...
                    tone_order,'Ordinal',false);

                newT.session_num = s;
                newT.subj = D.subj(1);

                try
                    pred_logit = predict( ...
                        lme,newT,'Conditional',false);
                catch
                    pred_logit = predict(lme,newT);
                end

                tone_prob(t) = invlogit(pred_logit);

                x = fixed_design_row( ...
                    coef_names,group_order,session_order,tone_order, ...
                    group_order{g},session_order{s},tone_order{t});

                tone_prob_sim(:,t) = invlogit(beta_sim*x');
            end

            adjusted_probability(row) = mean(tone_prob,'omitnan');

            sim_mean_probability = mean(tone_prob_sim,2,'omitnan');
            ci = quantile(sim_mean_probability, ...
                [alpha_level/2, 1-alpha_level/2]);

            ci_lower_probability(row) = ci(1);
            ci_upper_probability(row) = ci(2);
        end
    end

    adjT = table;
    adjT.Group = group_col;
    adjT.Session = session_col;
    adjT.AdjustedProbability = adjusted_probability;
    adjT.CI_Lower_probability = ci_lower_probability;
    adjT.CI_Upper_probability = ci_upper_probability;
end

function y = invlogit(x)
    y = 1./(1+exp(-x));
end

function safe_writetable(X,outPath)
    try
        if istable(X)
            T = X;
        else
            T = to_table_for_write(X);
        end

        writetable(T,outPath);

    catch ME
        warning('safe_writetable failed for %s: %s', ...
            outPath,ME.message);

        try
            fid = fopen(strrep(outPath,'.csv','.txt'),'w');
            fprintf(fid,'%s\n',evalc('disp(X)'));
            fclose(fid);
        catch
        end
    end
end

function T = to_table_for_write(X)
    if istable(X)
        T = X;
        return;
    end

    if exist('dataset2table','file') == 2
        try
            T = dataset2table(X);
            return;
        catch
        end
    end

    try
        varNames = X.Properties.VarNames;
        vals = cell(1,numel(varNames));

        for ii = 1:numel(varNames)
            vals{ii} = X.(varNames{ii});
        end

        T = table(vals{:}, ...
            'VariableNames', ...
            matlab.lang.makeValidName(varNames));
        return;
    catch
    end

    try
        varNames = X.Properties.VariableNames;
        vals = cell(1,numel(varNames));

        for ii = 1:numel(varNames)
            vals{ii} = X.(varNames{ii});
        end

        T = table(vals{:}, ...
            'VariableNames', ...
            matlab.lang.makeValidName(varNames));
        return;
    catch
    end

    try
        T = struct2table(X);
        return;
    catch
    end

    T = table( ...
        string(evalc('disp(X)')), ...
        'VariableNames',{'DisplayText'});
end
