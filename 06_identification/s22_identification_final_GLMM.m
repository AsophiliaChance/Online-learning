%% ============================================================
% ACC_finalGLMM_noPlot_corrected.m
%
% Purpose:
%   Fit the final identification-accuracy GLMM and export model summaries
%   and model-predicted marginal probabilities without any plotting.
%
% Final model:
%   acc ~ group*session + tone + speaker + (1 + session|subj) + (1|item)
%
% Main corrections relative to the earlier plotting script:
%   1. Add DummyVarCoding = 'effects' for Type-III / marginal tests.
%   2. Do NOT manually construct the fixed-effect design matrix X.
%   3. Use predict(glme_FINAL, newT, 'Conditional', false) for population-
%      level model-predicted probabilities.
%   4. Remove all CI simulation and plotting code.
%
% Input:
%   D:\online_learning\indentification\LMM_long_table_withToneSpeaker_itemTrial_v12.mat
%   or
%   D:\online_learning\indentification\LMM_long_table_withToneSpeaker_itemTrial_v12.csv
%
% Required variables:
%   subj, group, session, tone, speaker, acc, item
%% ============================================================

clear; clc;

%% ---------------- PATHS ----------------
rootDir = 'D:\online_learning\indentification';

matFile = fullfile(rootDir, 'LMM_long_table_withToneSpeaker_itemTrial_v12.mat');
csvFile = fullfile(rootDir, 'LMM_long_table_withToneSpeaker_itemTrial_v12.csv');

outDir = fullfile(rootDir, 'ACC_finalGLMM_noPlot_corrected');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

out_report = fullfile(outDir, 'ACC_finalGLMM_noPlot_corrected_report.txt');
out_anova  = fullfile(outDir, 'ACC_finalGLMM_TypeIII_marginal_ANOVA.csv');
out_coef   = fullfile(outDir, 'ACC_finalGLMM_fixed_coefficients.csv');
out_full   = fullfile(outDir, 'ACC_predicted_probability_full_reference_grid.csv');
out_sess   = fullfile(outDir, 'ACC_predicted_probability_session_marginal.csv');
out_gs     = fullfile(outDir, 'ACC_predicted_probability_group_by_session.csv');

assert(isfolder(rootDir), 'rootDir not found: %s', rootDir);

%% ---------------- LOAD DATA ----------------
if exist(matFile, 'file')
    S = load(matFile);
    assert(isfield(S, 'All'), 'MAT file does not contain variable "All": %s', matFile);
    All = S.All;
elseif exist(csvFile, 'file')
    All = readtable(csvFile, 'TextType', 'string');
else
    error('Neither MAT nor CSV long table found in %s', rootDir);
end

reqVars = {'subj','group','session','tone','speaker','acc','item'};
for i = 1:numel(reqVars)
    assert(ismember(reqVars{i}, All.Properties.VariableNames), ...
        'Missing required variable: %s', reqVars{i});
end

%% ---------------- FORMAT DATA ----------------
All.subj    = categorical(All.subj);
All.group   = categorical(All.group);
All.session = categorical(All.session);
All.tone    = categorical(All.tone);
All.speaker = categorical(All.speaker);
All.item    = categorical(All.item);

% Keep pre/post ordering if both levels exist.
if all(ismember({'pre','post'}, cellstr(categories(All.session))))
    All.session = reordercats(All.session, {'pre','post'});
end

D = All(~isnan(All.acc) & ...
        ~isundefined(All.subj) & ...
        ~isundefined(All.group) & ...
        ~isundefined(All.session) & ...
        ~isundefined(All.tone) & ...
        ~isundefined(All.speaker) & ...
        ~isundefined(All.item), :);

D.acc = double(D.acc);

% Remove unused categorical levels after filtering.
D.group   = removecats(D.group);
D.session = removecats(D.session);
D.tone    = removecats(D.tone);
D.speaker = removecats(D.speaker);
D.subj    = removecats(D.subj);
D.item    = removecats(D.item);

% Basic response check.
badAcc = isnan(D.acc) | D.acc < 0 | D.acc > 1;
if any(badAcc)
    disp(D(badAcc, :));
    error('acc must be binary/proportion in [0, 1]. Invalid rows shown above.');
end

%% ---------------- REPORT START ----------------
fid = fopen(out_report, 'w');
assert(fid > 0, 'Cannot open report file: %s', out_report);
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, '============================================================\n');
fprintf(fid, 'Identification ACC final GLMM report, no plotting\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Root directory: %s\n', rootDir);
fprintf(fid, 'Output directory: %s\n', outDir);
fprintf(fid, '============================================================\n\n');

fprintf(fid, 'Data summary\n');
fprintf(fid, '------------\n');
fprintf(fid, 'Usable rows: %d\n', height(D));
fprintf(fid, 'Subjects: %d\n', numel(categories(D.subj)));
fprintf(fid, 'Items: %d\n', numel(categories(D.item)));
fprintf(fid, 'Groups: %s\n', strjoin(cellstr(categories(D.group)), ', '));
fprintf(fid, 'Sessions: %s\n', strjoin(cellstr(categories(D.session)), ', '));
fprintf(fid, 'Tones: %s\n', strjoin(cellstr(categories(D.tone)), ', '));
fprintf(fid, 'Speakers: %s\n\n', strjoin(cellstr(categories(D.speaker)), ', '));

fprintf(fid, 'Observed accuracy by session:\n%s\n', evalc('disp(groupsummary(D, ''session'', ''mean'', ''acc''))'));
fprintf(fid, 'Observed accuracy by group x session:\n%s\n', evalc('disp(groupsummary(D, {''group'',''session''}, ''mean'', ''acc''))'));

%% ============================================================
% FIT FINAL MODEL
%% ============================================================

final_formula = 'acc ~ group*session + tone + speaker + (1 + session|subj) + (1|item)';

fprintf(fid, '\n============================================================\n');
fprintf(fid, 'FINAL GLMM\n');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Formula:\n%s\n\n', final_formula);
fprintf(fid, 'Fit settings:\n');
fprintf(fid, 'Distribution = Binomial\n');
fprintf(fid, 'Link = logit\n');
fprintf(fid, 'FitMethod = Laplace\n');
fprintf(fid, 'DummyVarCoding = effects\n\n');

% Use effects coding so that anova(glme_FINAL) can be interpreted as
% Type-III / marginal tests for categorical predictors and interactions.
glme_FINAL = fitglme(D, final_formula, ...
    'Distribution', 'Binomial', ...
    'Link', 'logit', ...
    'FitMethod', 'Laplace', ...
    'DummyVarCoding', 'effects');

fprintf(fid, 'Model summary:\n%s\n', evalc('disp(glme_FINAL)'));

%% ---------------- ANOVA / WALD TESTS ----------------
fprintf(fid, '\n============================================================\n');
fprintf(fid, 'TYPE-III / MARGINAL ANOVA FIXED-EFFECT TESTS\n');
fprintf(fid, '============================================================\n');

try
    anovaTbl = anova(glme_FINAL);
    fprintf(fid, 'ANOVA output:\n%s\n', evalc('disp(anovaTbl)'));
    safe_writetable(anovaTbl, out_anova);
catch ME
    fprintf(fid, 'anova(glme_FINAL) failed:\n%s\n', ME.message);
end

%% ---------------- FIXED COEFFICIENTS ----------------
fprintf(fid, '\n============================================================\n');
fprintf(fid, 'FIXED-EFFECT COEFFICIENTS\n');
fprintf(fid, '============================================================\n');

coefTbl = glme_FINAL.Coefficients;
fprintf(fid, 'Coefficients:\n%s\n', evalc('disp(coefTbl)'));
safe_writetable(coefTbl, out_coef);

%% ============================================================
% MODEL-PREDICTED PROBABILITIES USING predict()
%% ============================================================

fprintf(fid, '\n============================================================\n');
fprintf(fid, 'MODEL-PREDICTED PROBABILITIES\n');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Predictions use predict(glme_FINAL, newT, ''Conditional'', false).\n');
fprintf(fid, 'Thus, predictions are population-level fixed-effect predictions.\n\n');

% Full reference grid: group x session x tone x speaker.
newT = build_reference_grid(D);

% Population-level predictions. Random effects are not included.
newT.pred_prob = predict(glme_FINAL, newT, 'Conditional', false);

safe_writetable(newT, out_full);
fprintf(fid, 'Saved full reference-grid predictions:\n%s\n\n', out_full);

% Session main effect: average equally across group, tone, and speaker.
[G_sess, sess_levels] = findgroups(newT.session);
pred_sess = splitapply(@(x) mean(x, 'omitnan'), newT.pred_prob, G_sess);
pred_session_tbl = table(sess_levels, pred_sess, ...
    'VariableNames', {'session','predicted_probability'});

safe_writetable(pred_session_tbl, out_sess);
fprintf(fid, 'Session marginal predicted probability, averaged equally across group, tone, and speaker:\n%s\n', ...
    evalc('disp(pred_session_tbl)'));

% Group x session interaction: average equally across tone and speaker.
[G_gs, g_levels, s_levels] = findgroups(newT.group, newT.session);
pred_gs = splitapply(@(x) mean(x, 'omitnan'), newT.pred_prob, G_gs);
pred_group_session_tbl = table(g_levels, s_levels, pred_gs, ...
    'VariableNames', {'group','session','predicted_probability'});

safe_writetable(pred_group_session_tbl, out_gs);
fprintf(fid, 'Group x session marginal predicted probability, averaged equally across tone and speaker:\n%s\n', ...
    evalc('disp(pred_group_session_tbl)'));

%% ---------------- DONE ----------------
fprintf(fid, '\n============================================================\n');
fprintf(fid, 'Done.\n');
fprintf(fid, 'Report saved to: %s\n', out_report);
fprintf(fid, '============================================================\n');

disp('Done.');
disp(['Report saved to: ', out_report]);
disp(['Output directory: ', outDir]);

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function newT = build_reference_grid(D)
    groupCats   = categories(D.group);
    sessionCats = categories(D.session);
    toneCats    = categories(D.tone);
    speakerCats = categories(D.speaker);

    nRows = numel(groupCats) * numel(sessionCats) * numel(toneCats) * numel(speakerCats);

    groupCol   = strings(nRows, 1);
    sessionCol = strings(nRows, 1);
    toneCol    = strings(nRows, 1);
    speakerCol = strings(nRows, 1);

    r = 0;
    for g = 1:numel(groupCats)
        for s = 1:numel(sessionCats)
            for t = 1:numel(toneCats)
                for sp = 1:numel(speakerCats)
                    r = r + 1;
                    groupCol(r)   = string(groupCats{g});
                    sessionCol(r) = string(sessionCats{s});
                    toneCol(r)    = string(toneCats{t});
                    speakerCol(r) = string(speakerCats{sp});
                end
            end
        end
    end

    newT = table;
    newT.group   = categorical(groupCol, string(groupCats));
    newT.session = categorical(sessionCol, string(sessionCats));
    newT.tone    = categorical(toneCol, string(toneCats));
    newT.speaker = categorical(speakerCol, string(speakerCats));

    % Grouping variables are required by the model formula even though
    % Conditional=false excludes their random effects from prediction.
    newT.subj = repmat(D.subj(1), height(newT), 1);
    newT.item = repmat(D.item(1), height(newT), 1);
end

function safe_writetable(X, outPath)
    % Robust table writer for MATLAB table/dataset/titled output objects.
    try
        if istable(X)
            writetable(X, outPath);
            return;
        end
    catch
    end

    try
        if isa(X, 'dataset')
            T = dataset2table(X);
            writetable(T, outPath);
            return;
        end
    catch
    end

    try
        T = struct2table(struct(X));
        writetable(T, outPath);
        return;
    catch
    end

    % Fallback: write display output as .txt if CSV writing fails.
    [p, n, ~] = fileparts(outPath);
    txtPath = fullfile(p, [n '.txt']);
    fid2 = fopen(txtPath, 'w');
    if fid2 > 0
        fprintf(fid2, '%s', evalc('disp(X)'));
        fclose(fid2);
    else
        warning('safe_writetable:Failed', 'Could not write %s or fallback txt.', outPath);
    end
end
