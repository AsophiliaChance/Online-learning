%% Direct differences in change: Offline minus Online, same final REML models.
% Alpha: all 3 Block pairs, Holm across 3 contrasts.
% Production: each later assessment minus pre-test, Holm across 3 contrasts.
% Estimates and unadjusted 95% CIs are effect sizes on the model scale.
% Alpha 10^estimate is a ratio of power-change ratios; exp(production estimate)
% is a ratio of empirical-logit odds-change ratios, not a GLMM odds ratio.
clear; clc;
base = 'D:\online_learning';
out = fullfile(base,'manuscript_revision','change_difference_results');
if ~exist(out,'dir'), mkdir(out); end
A = readtable(fullfile(base,'training_psychophysiology','alpha_VOFF_simple_effects','alpha_simple_effects_analysis_table.csv'));
A.Subject = categorical(string(A.Subject));
A.Group = categorical(string(A.Group), {'Group2','Group1'});
A.Block = categorical(string(A.Block), {'B1','B2','B3'});
ma = fitlme(A,'AlphaPower ~ Group*Block + AlphaPower_baseline + (1|Subject)', 'FitMethod','REML','DummyVarCoding','effects');
P = readtable(fullfile(base,'production','production_simple_effects_revised','production_simple_effects_analysis_table.csv'));
P.subj = categorical(string(P.subj));
P.group = categorical(string(P.group), {'Online','Offline'});
P.session = categorical(string(P.session), {'test1','test2','test3','test4'});
P.tone = categorical(string(P.tone), {'tone2','tone3','tone4'});
mp = fitlme(P,'logit_acc ~ group*session + tone + (1|subj)', 'FitMethod','REML','DummyVarCoding','effects');
models={ma,mp}; terms={'Group','group'}; times={'Block','session'};
groups={{'Group2','Group1'},{'Online','Offline'}};
levels={{'B1','B2','B3'},{'test1','test2','test3','test4'}};
pairs={ [2 1;3 1;3 2], [2 1;3 1;4 1] };
names={'alpha','production'};
for k=1:2
 m=models{k}; beta=fixedEffects(m); nm=m.CoefficientNames; V=m.CoefficientCovariance;
 result=table();
 for j=1:size(pairs{k},1)
  late=levels{k}{pairs{k}(j,1)};early=levels{k}{pairs{k}(j,2)};
  H = design(nm,terms{k},times{k},groups{k}{2},late,groups{k},levels{k}) ...
     -design(nm,terms{k},times{k},groups{k}{2},early,groups{k},levels{k}) ...
     -design(nm,terms{k},times{k},groups{k}{1},late,groups{k},levels{k}) ...
     +design(nm,terms{k},times{k},groups{k}{1},early,groups{k},levels{k});
  label=sprintf('(Offline %s-%s) - (Online %s-%s)',late,early,late,early);
  row=single_contrast_result(m,H,beta,V,label,.05);
  result=[result;row];
 end
 result.p_Holm=holm_adjust(result.p_raw);
 result.Significant_Holm=result.p_Holm<.05;
 if k==1
  result.ChangeRatio=10.^result.Estimate;
  result.Ratio_CI_Lower=10.^result.CI_Lower;
  result.Ratio_CI_Upper=10.^result.CI_Upper;
 else
  result.ChangeRatio=exp(result.Estimate);
  result.Ratio_CI_Lower=exp(result.CI_Lower);
  result.Ratio_CI_Upper=exp(result.CI_Upper);
 end
 writetable(result,fullfile(out,[names{k} '_change_difference.csv']));
 fprintf('%s\n',names{k});disp(result);
 fid=fopen(fullfile(out,[names{k} '_change_difference_report.txt']),'w');
 fprintf(fid,'%s\n',evalc('disp(result)'));
 fprintf(fid,'Direction: Offline change minus Online change. Holm correction across three contrasts. 95%% CIs unadjusted.\n');
 fprintf(fid,'%s\n',evalc('disp(anova(m,''DFMethod'',''Satterthwaite''))'));
 fclose(fid);
end
function x=design(nm,gprefix,tprefix,g,t,gorder,torder)
 x=zeros(1,numel(nm));
 for j=1:numel(nm)
  name=char(nm{j});gt=regexp(name,[gprefix '_([^:]+)'],'tokens','once');tt=regexp(name,[tprefix '_([^:]+)'],'tokens','once');
  if ~isempty(gt) && ~isempty(tt)
   x(j)=code(g,gt{1},gorder)*code(t,tt{1},torder);
  elseif ~isempty(gt),x(j)=code(g,gt{1},gorder);
  elseif ~isempty(tt),x(j)=code(t,tt{1},torder);
  end
 end
end
function v=code(level,coefficient,order)
 if strcmp(level,coefficient),v=1;elseif strcmp(level,order{end}),v=-1;else,v=0;end
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

