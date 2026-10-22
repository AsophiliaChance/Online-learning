%% Alpha VOFF simple-effects analysis for a significant Group x Block interaction
%
% Standalone script based on the main Alpha VOFF LMM analysis.
%
% Input:
%   1) alpha_VOFF_block_level_0p3_2p3.csv
%   2) alpha_rest_baseline_2s.csv (optional but recommended)
%
% Full model:
%   With baseline:
%     AlphaPower ~ Group*Block + AlphaPower_baseline + (1|Subject)
%     AlphaPower ~ Group*Block + AlphaPower_baseline + (1+Block|Subject)
%
%   Without baseline:
%     AlphaPower ~ Group*Block + (1|Subject)
%     AlphaPower ~ Group*Block + (1+Block|Subject)
%
% Follow-up tests:
%   A) Group1 vs Group2 within B1, B2, and B3
%   B) Omnibus Block effect within Group1 and Group2
%   C) Pairwise Block contrasts within each Group:
%        B2-B1, B3-B1, and B3-B2
%
% Multiple-comparison correction:
%   - Holm correction across the three within-Block group comparisons
%   - Holm correction across the two within-Group omnibus Block tests
%   - Holm correction across the six pairwise within-Group Block comparisons
%
% Notes:
%   - All follow-up contrasts are tested in the original full model.
%   - The data are not split and refitted separately by Group.
%   - Positive Group1-Group2 estimates indicate higher alpha power in Group1.
%   - Positive B2-B1 estimates indicate higher alpha power in B2.
%
% Output:
%   alpha_simple_effects_summary.txt
%   alpha_adjusted_means_Group_Block.csv
%   alpha_group_effects_within_Block.csv
%   alpha_omnibus_Block_effects_within_Group.csv
%   alpha_pairwise_Block_effects_within_Group.csv
%   alpha_simple_effects_final_model_ANOVA.csv
%   alpha_simple_effects_final_model_coefficients.csv
%   alpha_simple_effects_analysis_table.csv

clear; clc;

%% ===================== User settings =====================
base_dir = 'D:\online_learning\training_psychophysiology';

alpha_block_csv = fullfile(base_dir, ...
    'alpha_VOFF_block_level_0p3_2p3.csv');
rest_base_csv = fullfile(base_dir, ...
    'alpha_rest_baseline_2s.csv');

out_dir = fullfile(base_dir, ...
    'alpha_VOFF_simple_effects');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

summary_txt = fullfile(out_dir, ...
    'alpha_simple_effects_summary.txt');

% Group coding
group1_ids = [1:6, 14:20, 28:29, 33];
group1_label = 'Group1';
group2_label = 'Group2';

% Category order must match the main analysis
group_category_order = {group2_label, group1_label};
block_category_order = {'B1','B2','B3'};

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

alpha_level = 0.05;
random_slope_p_threshold = 0.05;
use_satterthwaite = true;

%% ===================== Read block-level alpha data =====================
if ~exist(alpha_block_csv, 'file')
    error('Cannot find alpha block-level CSV:\n%s', alpha_block_csv);
end

Traw = readtable_compat(alpha_block_csv);
Traw.Properties.VariableNames = ...
    matlab.lang.makeValidName(Traw.Properties.VariableNames);

subj_var = find_var(Traw, ...
    {'Subject','SubjID','Participant','ID','File','Filename'});
block_var = find_var(Traw, ...
    {'Block','TrainingBlock','block'});
outcome_var = find_var(Traw, outcome_candidates);

if isempty(subj_var)
    error('No Subject-like variable was found.');
end
if isempty(block_var)
    error('No Block-like variable was found.');
end
if isempty(outcome_var)
    outcome_var = infer_alpha_outcome_var(Traw);
end
if isempty(outcome_var)
    error(['Could not identify the alpha outcome variable. ' ...
        'Rename it to AlphaPower or add its name to outcome_candidates.']);
end

T = table();
T.SubjectRaw = string(Traw.(subj_var));
T.SubjectID = extract_subject_id(T.SubjectRaw);
T.Subject = categorical(string(T.SubjectID));
T.Block = normalize_block( ...
    Traw.(block_var), block_category_order);
T.AlphaPower = double(Traw.(outcome_var));

valid = ~isnan(T.SubjectID) & ...
    ~isundefined(T.Subject) & ...
    ~isundefined(T.Block) & ...
    isfinite(T.AlphaPower);
T = T(valid,:);

% Reconstruct Group
grp = strings(height(T),1);
for i = 1:height(T)
    if ismember(T.SubjectID(i), group1_ids)
        grp(i) = group1_label;
    else
        grp(i) = group2_label;
    end
end
T.Group = categorical( ...
    grp, group_category_order, 'Ordinal', false);

%% ===================== Merge resting alpha baseline =====================
has_baseline = false;

if exist(rest_base_csv, 'file')
    Braw = readtable_compat(rest_base_csv);
    Braw.Properties.VariableNames = ...
        matlab.lang.makeValidName(Braw.Properties.VariableNames);

    b_subj_var = find_var(Braw, ...
        {'Subject','SubjID','Participant','ID','File','Filename'});
    b_alpha_var = find_var(Braw, baseline_candidates);

    if ~isempty(b_subj_var) && ~isempty(b_alpha_var)
        B = table();
        B.SubjectRaw = string(Braw.(b_subj_var));
        B.SubjectID = extract_subject_id(B.SubjectRaw);
        B.AlphaPower_baseline = double(Braw.(b_alpha_var));

        B = B(~isnan(B.SubjectID) & ...
            isfinite(B.AlphaPower_baseline),:);

        if ~isempty(B)
            Bmean = groupsummary( ...
                B, 'SubjectID', 'mean', 'AlphaPower_baseline');

            mean_name = 'mean_AlphaPower_baseline';
            if ismember(mean_name, Bmean.Properties.VariableNames)
                Bmean.Properties.VariableNames{mean_name} = ...
                    'AlphaPower_baseline';
            end

            Bmean = Bmean(:, ...
                {'SubjectID','AlphaPower_baseline'});

            T = outerjoin( ...
                T, Bmean, ...
                'Keys','SubjectID', ...
                'MergeKeys',true, ...
                'Type','left');

            has_baseline = any( ...
                isfinite(T.AlphaPower_baseline));
        end
    else
        warning(['Baseline file was found, but its Subject or ' ...
            'alpha variable could not be identified.']);
    end
else
    warning('Baseline CSV not found. Fitting models without baseline.');
end

%% ===================== Keep complete model rows =====================
if has_baseline
    keep = isfinite(T.AlphaPower) & ...
        isfinite(T.AlphaPower_baseline) & ...
        ~isundefined(T.Group) & ...
        ~isundefined(T.Block) & ...
        ~isundefined(T.Subject);
else
    keep = isfinite(T.AlphaPower) & ...
        ~isundefined(T.Group) & ...
        ~isundefined(T.Block) & ...
        ~isundefined(T.Subject);
end

Tuse = T(keep,:);

Tuse.Group = categorical( ...
    string(Tuse.Group), ...
    group_category_order, ...
    'Ordinal',false);

Tuse.Block = categorical( ...
    string(Tuse.Block), ...
    block_category_order, ...
    'Ordinal',false);

Tuse.Subject = categorical(string(Tuse.Subject));

if isempty(Tuse)
    error('No usable observations remain.');
end

%% ===================== Open summary =====================
fid = fopen(summary_txt,'w');
if fid < 0
    error('Cannot open summary file:\n%s', summary_txt);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf_both(fid, ...
    'ALPHA VOFF SIMPLE-EFFECTS ANALYSIS\n');
fprintf_both(fid, ...
    'Generated: %s\n', datestr(now));
fprintf_both(fid, ...
    'Input block CSV: %s\n', alpha_block_csv);
fprintf_both(fid, ...
    'Input baseline CSV: %s\n', rest_base_csv);
fprintf_both(fid, ...
    'Outcome variable: %s\n', outcome_var);
fprintf_both(fid, ...
    'Baseline covariate used: %d\n', has_baseline);
fprintf_both(fid, ...
    'Rows used: %d\n', height(Tuse));
fprintf_both(fid, ...
    'Subjects used: %d\n', ...
    numel(categories(Tuse.Subject)));

fprintf_both(fid, ...
    '\nCounts by Group x Block:\n');
disp_to_file(fid, ...
    groupsummary(Tuse, {'Group','Block'}));

fprintf_both(fid, ...
    '\nRaw descriptive statistics by Group x Block:\n');
D = groupsummary( ...
    Tuse, {'Group','Block'}, ...
    {'mean','std','numel'}, 'AlphaPower');
disp_to_file(fid,D);
writetable(D, fullfile( ...
    out_dir, 'alpha_raw_descriptives_Group_Block.csv'));

%% ===================== Define candidate models =====================
if has_baseline
    formula_RI = [ ...
        'AlphaPower ~ Group*Block + ' ...
        'AlphaPower_baseline + (1|Subject)'];

    formula_RS = [ ...
        'AlphaPower ~ Group*Block + ' ...
        'AlphaPower_baseline + (1+Block|Subject)'];
else
    formula_RI = ...
        'AlphaPower ~ Group*Block + (1|Subject)';

    formula_RS = ...
        'AlphaPower ~ Group*Block + (1+Block|Subject)';
end

%% ===================== Random-effects comparison =====================
fprintf_both(fid, ...
    '\n====================================================\n');
fprintf_both(fid, ...
    'RANDOM-EFFECTS COMPARISON\n');
fprintf_both(fid, ...
    '====================================================\n');

lme_RI_ML = fitlme( ...
    Tuse, formula_RI, ...
    'FitMethod','ML', ...
    'DummyVarCoding','effects');

use_RS = false;
final_formula = formula_RI;
p_random = NaN;

try
    lme_RS_ML = fitlme( ...
        Tuse, formula_RS, ...
        'FitMethod','ML', ...
        'DummyVarCoding','effects');

    cmp_random = compare(lme_RI_ML,lme_RS_ML);
    disp_to_file(fid,cmp_random);

    p_random = get_compare_pvalue(cmp_random);
    fprintf_both(fid, ...
        'Extracted RI-vs-RS p = %.6g\n', p_random);

    if isfinite(p_random) && ...
            p_random < random_slope_p_threshold
        use_RS = true;
        final_formula = formula_RS;
    end
catch ME
    fprintf_both(fid, ...
        ['Random-slope model failed. ' ...
        'The random-intercept model was retained.\n%s\n'], ...
        ME.message);
end

if use_RS
    fprintf_both(fid, ...
        'Selected model: random intercept + Block slope.\n');
else
    fprintf_both(fid, ...
        'Selected model: random intercept only.\n');
end
fprintf_both(fid, ...
    'Final formula: %s\n', final_formula);

%% ===================== Fit final REML model =====================
lme = fitlme( ...
    Tuse, final_formula, ...
    'FitMethod','REML', ...
    'DummyVarCoding','effects');

fprintf_both(fid, ...
    '\n====================================================\n');
fprintf_both(fid, ...
    'FINAL REML MODEL\n');
fprintf_both(fid, ...
    '====================================================\n');
disp_to_file(fid,lme);

coef_names = string(lme.CoefficientNames(:));
beta = fixedEffects(lme);
Vbeta = lme.CoefficientCovariance;

fprintf_both(fid, ...
    '\nFixed-effect coefficient names:\n');
disp_to_file(fid,lme.Coefficients);

%% ===================== Type III ANOVA =====================
try
    if use_satterthwaite
        a_final = anova( ...
            lme,'DFMethod','Satterthwaite');
    else
        a_final = anova(lme);
    end
catch ME
    warning(['Requested ANOVA failed: %s. ' ...
        'Using default anova(lme).'], ME.message);
    a_final = anova(lme);
end

writetable( ...
    table_compat(a_final), ...
    fullfile(out_dir, ...
    'alpha_simple_effects_final_model_ANOVA.csv'));

writetable( ...
    table_compat(lme.Coefficients), ...
    fullfile(out_dir, ...
    'alpha_simple_effects_final_model_coefficients.csv'));

%% ===================== Covariate reference value =====================
if has_baseline
    baseline_value = mean( ...
        Tuse.AlphaPower_baseline,'omitnan');
else
    baseline_value = NaN;
end

%% ===================== Model-adjusted means =====================
adjT = make_adjusted_means_table( ...
    Tuse, lme, has_baseline, baseline_value, ...
    group_category_order, block_category_order, ...
    alpha_level);

writetable( ...
    adjT, ...
    fullfile(out_dir, ...
    'alpha_adjusted_means_Group_Block.csv'));

fprintf_both(fid, ...
    '\n====================================================\n');
fprintf_both(fid, ...
    'MODEL-ADJUSTED MEANS\n');
fprintf_both(fid, ...
    '====================================================\n');
disp_to_file(fid,adjT);

%% ===================== A. Group within each Block =====================
fprintf_both(fid, ...
    '\n====================================================\n');
fprintf_both(fid, ...
    'GROUP SIMPLE EFFECTS WITHIN BLOCK\n');
fprintf_both(fid, ...
    'Direction: Group1 - Group2\n');
fprintf_both(fid, ...
    '====================================================\n');

group_results = table();

for b = 1:numel(block_category_order)
    block_label = block_category_order{b};

    x_g1 = fixed_design_row( ...
        coef_names, ...
        group_category_order, ...
        block_category_order, ...
        group1_label, block_label, ...
        has_baseline, baseline_value);

    x_g2 = fixed_design_row( ...
        coef_names, ...
        group_category_order, ...
        block_category_order, ...
        group2_label, block_label, ...
        has_baseline, baseline_value);

    H = x_g1 - x_g2;

    R = single_contrast_result( ...
        lme,H,beta,Vbeta, ...
        sprintf('%s - %s at %s', ...
        group1_label,group2_label,block_label), ...
        alpha_level);

    R.Group = string( ...
        [group1_label '-' group2_label]);
    R.Block = string(block_label);

    group_results = [group_results; R]; %#ok<AGROW>
end

group_results.p_Holm = ...
    holm_adjust(group_results.p_raw);

group_results.Significant_Holm = ...
    group_results.p_Holm < alpha_level;

disp_to_file(fid,group_results);

writetable( ...
    group_results, ...
    fullfile(out_dir, ...
    'alpha_group_effects_within_Block.csv'));

%% ===================== B. Omnibus Block within Group =====================
fprintf_both(fid, ...
    '\n====================================================\n');
fprintf_both(fid, ...
    'OMNIBUS BLOCK SIMPLE EFFECTS WITHIN GROUP\n');
fprintf_both(fid, ...
    'Joint contrasts: B2-B1 and B3-B1\n');
fprintf_both(fid, ...
    '====================================================\n');

omnibus_results = table();

for g = 1:numel(group_category_order)
    group_label = group_category_order{g};

    x_b1 = fixed_design_row( ...
        coef_names, ...
        group_category_order, ...
        block_category_order, ...
        group_label,'B1', ...
        has_baseline,baseline_value);

    x_b2 = fixed_design_row( ...
        coef_names, ...
        group_category_order, ...
        block_category_order, ...
        group_label,'B2', ...
        has_baseline,baseline_value);

    x_b3 = fixed_design_row( ...
        coef_names, ...
        group_category_order, ...
        block_category_order, ...
        group_label,'B3', ...
        has_baseline,baseline_value);

    H_joint = [ ...
        x_b2 - x_b1; ...
        x_b3 - x_b1];

    [p,F,df1,df2] = run_coef_test(lme,H_joint);

    row = table( ...
        string(group_label),F,df1,df2,p, ...
        'VariableNames', ...
        {'Group','F','DF1','DF2','p_raw'});

    omnibus_results = [omnibus_results; row]; %#ok<AGROW>
end

omnibus_results.p_Holm = ...
    holm_adjust(omnibus_results.p_raw);

omnibus_results.Significant_Holm = ...
    omnibus_results.p_Holm < alpha_level;

disp_to_file(fid,omnibus_results);

writetable( ...
    omnibus_results, ...
    fullfile(out_dir, ...
    'alpha_omnibus_Block_effects_within_Group.csv'));

%% ===================== C. Pairwise Block within Group =====================
fprintf_both(fid, ...
    '\n====================================================\n');
fprintf_both(fid, ...
    'PAIRWISE BLOCK SIMPLE EFFECTS WITHIN GROUP\n');
fprintf_both(fid, ...
    'Directions: B2-B1, B3-B1, B3-B2\n');
fprintf_both(fid, ...
    '====================================================\n');

pair_defs = { ...
    'B2','B1','B2 - B1'; ...
    'B3','B1','B3 - B1'; ...
    'B3','B2','B3 - B2'};

pairwise_results = table();

for g = 1:numel(group_category_order)
    group_label = group_category_order{g};

    for k = 1:size(pair_defs,1)
        block_a = pair_defs{k,1};
        block_b = pair_defs{k,2};
        label = pair_defs{k,3};

        x_a = fixed_design_row( ...
            coef_names, ...
            group_category_order, ...
            block_category_order, ...
            group_label,block_a, ...
            has_baseline,baseline_value);

        x_b = fixed_design_row( ...
            coef_names, ...
            group_category_order, ...
            block_category_order, ...
            group_label,block_b, ...
            has_baseline,baseline_value);

        H = x_a - x_b;

        R = single_contrast_result( ...
            lme,H,beta,Vbeta, ...
            sprintf('%s within %s', ...
            label,group_label), ...
            alpha_level);

        R.Group = string(group_label);
        R.BlockContrast = string(label);

        pairwise_results = [pairwise_results; R]; %#ok<AGROW>
    end
end

pairwise_results.p_Holm = ...
    holm_adjust(pairwise_results.p_raw);

pairwise_results.Significant_Holm = ...
    pairwise_results.p_Holm < alpha_level;

disp_to_file(fid,pairwise_results);

writetable( ...
    pairwise_results, ...
    fullfile(out_dir, ...
    'alpha_pairwise_Block_effects_within_Group.csv'));

%% ===================== Save analysis table =====================
writetable( ...
    Tuse, ...
    fullfile(out_dir, ...
    'alpha_simple_effects_analysis_table.csv'));

fprintf_both(fid, ...
    '\n====================================================\n');
fprintf_both(fid, ...
    'INTERPRETATION\n');
fprintf_both(fid, ...
    '====================================================\n');

fprintf_both(fid, ...
    ['Use p_Holm for the corrected inference ' ...
    'within each contrast family.\n']);

fprintf_both(fid, ...
    ['Positive Group1-Group2 estimates indicate ' ...
    'higher alpha power in Group1.\n']);

fprintf_both(fid, ...
    ['Positive later-minus-earlier Block estimates ' ...
    'indicate higher alpha power in the later Block.\n']);

fprintf_both(fid, ...
    '\nAll done. Output directory:\n%s\n', out_dir);

%% ===================== Local functions =====================

function T = readtable_compat(file)
    try
        T = readtable( ...
            file,'VariableNamingRule','preserve');
    catch
        T = readtable(file);
    end
end

function v = find_var(T,candidates)
    v = '';
    names = T.Properties.VariableNames;
    names_low = lower(names);

    for i = 1:numel(candidates)
        c = matlab.lang.makeValidName(candidates{i});

        idx = find(strcmpi(names,c),1);
        if ~isempty(idx)
            v = names{idx};
            return;
        end

        idx = find(strcmp(names_low,lower(c)),1);
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

        if contains(nm,'alpha') && ...
                ~contains(nm,'baseline') && ...
                ~contains(nm,'std') && ...
                ~contains(nm,'sem') && ...
                ~contains(nm,'ntrial') && ...
                ~contains(nm,'count') && ...
                isnumeric(T.(names{i}))

            v = names{i};
            return;
        end
    end
end

function id = extract_subject_id(x)
    x = string(x);
    id = nan(numel(x),1);

    for i = 1:numel(x)
        tok = regexp(x(i),'\d+','match');

        if ~isempty(tok)
            id(i) = str2double(tok{1});
        end
    end
end

function B = normalize_block(x,block_order)
    if isnumeric(x)
        s = strings(numel(x),1);

        for i = 1:numel(x)
            if isfinite(x(i))
                s(i) = "B" + string(round(x(i)));
            else
                s(i) = "";
            end
        end
    else
        s = strtrim(string(x));

        for i = 1:numel(s)
            if strlength(s(i)) == 0
                continue;
            end

            tok = regexp(s(i),'\d+','match');

            if ~isempty(tok)
                s(i) = "B" + string(str2double(tok{1}));
            else
                s(i) = upper(s(i));
            end
        end
    end

    B = categorical( ...
        s,block_order,'Ordinal',false);
end

function fprintf_both(fid,varargin)
    fprintf(varargin{:});
    fprintf(fid,varargin{:});
end

function disp_to_file(fid,obj)
    txt = evalc('disp(obj)');
    fprintf('%s\n',txt);
    fprintf(fid,'%s\n',txt);
end

function p = get_compare_pvalue(cmp)
    p = NaN;

    try
        if istable(cmp) && ...
                any(strcmp( ...
                cmp.Properties.VariableNames,'pValue'))

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

function Tout = table_compat(x)
    if istable(x)
        Tout = x;
        return;
    end

    try
        Tout = dataset2table(x);
    catch
        try
            Tout = struct2table(x);
        catch
            error('Could not convert result to a MATLAB table.');
        end
    end
end

function x = fixed_design_row( ...
        coef_names,group_order,block_order, ...
        group_value,block_value, ...
        has_baseline,baseline_value)

    p = numel(coef_names);
    x = zeros(1,p);

    for j = 1:p
        nm = char(coef_names(j));

        if strcmp(nm,'(Intercept)')
            x(j) = 1;
            continue;
        end

        if has_baseline && ...
                strcmp(nm,'AlphaPower_baseline')
            x(j) = baseline_value;
            continue;
        end

        has_group = contains(nm,'Group_');
        has_block = contains(nm,'Block_');

        if has_group && has_block
            g_level = extract_factor_level(nm,'Group_');
            b_level = extract_factor_level(nm,'Block_');

            x(j) = ...
                effects_code_value( ...
                group_value,g_level,group_order) .* ...
                effects_code_value( ...
                block_value,b_level,block_order);

        elseif has_group
            g_level = extract_factor_level(nm,'Group_');

            x(j) = effects_code_value( ...
                group_value,g_level,group_order);

        elseif has_block
            b_level = extract_factor_level(nm,'Block_');

            x(j) = effects_code_value( ...
                block_value,b_level,block_order);
        end
    end
end

function level = extract_factor_level(coef_name,prefix)
    expr = [ ...
        regexptranslate('escape',prefix), ...
        '([^:]+)'];

    tok = regexp( ...
        coef_name,expr,'tokens','once');

    if isempty(tok)
        error( ...
            'Could not parse coefficient name: %s', ...
            coef_name);
    end

    level = tok{1};
end

function value = effects_code_value( ...
        target_level,coefficient_level,category_order)

    target_level = char(string(target_level));
    coefficient_level = ...
        char(string(coefficient_level));

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
        adjusted_sorted(i) = min( ...
            1,(m-i+1)*p_sorted(i));
    end

    adjusted_sorted = cummax(adjusted_sorted);

    adjusted = zeros(m,1);
    adjusted(order) = adjusted_sorted;

    p_adj(valid) = adjusted;
end

function adjT = make_adjusted_means_table( ...
        Tuse,lme,has_baseline,baseline_value, ...
        group_order,block_order,alpha_level)

    [G,B] = ndgrid( ...
        string(group_order), ...
        string(block_order));

    n = numel(G);
    adjT = table();

    adjT.Group = categorical( ...
        reshape(G,n,1), ...
        group_order,'Ordinal',false);

    adjT.Block = categorical( ...
        reshape(B,n,1), ...
        block_order,'Ordinal',false);

    subj_cats = categories(Tuse.Subject);

    adjT.Subject = categorical( ...
        repmat(string(subj_cats{1}),n,1), ...
        subj_cats,'Ordinal',false);

    if has_baseline
        adjT.AlphaPower_baseline = ...
            repmat(baseline_value,n,1);
    end

    try
        [ypred,yci] = predict( ...
            lme,adjT, ...
            'Conditional',false, ...
            'Alpha',alpha_level);

        adjT.AdjustedMean = ypred;
        adjT.CI_Lower = yci(:,1);
        adjT.CI_Upper = yci(:,2);
    catch
        ypred = predict( ...
            lme,adjT,'Conditional',false);

        adjT.AdjustedMean = ypred;
        adjT.CI_Lower = nan(n,1);
        adjT.CI_Upper = nan(n,1);
    end
end
