% Approximate partial eta squared is F*df1/(F*df1+df2); not total model variance explained.
%% Posterior alpha analysis: between-group changes and within-group B3/B1.
% Same model and validated contrasts as the teacher-revision analysis.
% Input data remain in D:/online_learning; outputs are new timestamped folders
% beside this script. No existing scripts, manuscript or results overwritten.
% Between-group family: three temporal changes, Holm3.
% Within-group family: Online/Offline B3-B1, Holm2.
% All confidence intervals are unadjusted pointwise 95% intervals.
clear; clc;
base='D:\online_learning';
parent=fullfile(fileparts(mfilename('fullpath')),'alpha_analysis_effect_sizes');
if ~exist(parent,'dir'),mkdir(parent);end
out=fullfile(parent,['results_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(out,'dir'),'Output already exists; rerun after one second.');
mkdir(out);
diary(fullfile(out,'analysis_log.txt'));
A=readtable(fullfile(base,'training_psychophysiology','alpha_VOFF_simple_effects','alpha_simple_effects_analysis_table.csv'));
A.Subject=categorical(string(A.Subject));
A.Group=categorical(string(A.Group),{'Group2','Group1'});
A.Block=categorical(string(A.Block),{'B1','B2','B3'});
assert(all(isfinite(A.AlphaPower)) && all(isfinite(A.AlphaPower_baseline)));
assert(~any(isundefined(A.Group)) && ~any(isundefined(A.Block)));

%% 1. Alpha: retain the original final model and between-group tests.
ma = fitlme(A,'AlphaPower ~ Group*Block + AlphaPower_baseline + (1|Subject)', ...
    'FitMethod','REML','DummyVarCoding','effects');
agroups = {'Group2','Group1'}; ablocks = {'B1','B2','B3'};
alphaDiff = between_changes(ma,'Group','Block',agroups,ablocks,[2 1;3 1;3 2]);
alphaDiff.p_Holm3 = holm_adjust(alphaDiff.p_raw);
alphaDiff.ApproxPartialEtaSquared=(alphaDiff.F.*alphaDiff.DF1)./(alphaDiff.F.*alphaDiff.DF1+alphaDiff.DF2);
alphaDiff.Significant_Holm3 = alphaDiff.p_Holm3 < .05;
alphaDiff.RatioOfPowerChangeRatios = 10.^alphaDiff.Estimate;
alphaDiff.Ratio_CI_Lower = 10.^alphaDiff.CI_Lower;
alphaDiff.Ratio_CI_Upper = 10.^alphaDiff.CI_Upper;
writetable(alphaDiff,fullfile(out,'alpha_between_group_changes.csv'));
check_change_predictions(ma,A,'Group','Block',agroups,alphaDiff,[2 1;3 1;3 2],ablocks);

%% 2. Alpha: clarify within-group B1-to-B3 direction and magnitude.
withinAlpha = table();
for g = 1:2
    H = design(ma.CoefficientNames,'Group','Block',agroups{g},'B3',agroups,ablocks) ...
      - design(ma.CoefficientNames,'Group','Block',agroups{g},'B1',agroups,ablocks);
    row = contrast_result(ma,H,[agroups{g} ': B3-B1']);
    if g==1, row.Group=string('Online'); else, row.Group=string('Offline'); end
    row.B3_B1_PowerRatio = 10.^row.Estimate;
    row.Ratio_CI_Lower = 10.^row.CI_Lower;
    row.Ratio_CI_Upper = 10.^row.CI_Upper;
    row.Direction = string('increase');
    if row.Estimate<0, row.Direction=string('decrease'); end
    withinAlpha = [withinAlpha;row]; %#ok<AGROW>
end
% These two descriptive decomposition intervals are pointwise. If using
% their p values as a separate family of within-group tests, use Holm2.
withinAlpha.p_Holm2 = holm_adjust(withinAlpha.p_raw);
withinAlpha.ApproxPartialEtaSquared=(withinAlpha.F.*withinAlpha.DF1)./(withinAlpha.F.*withinAlpha.DF1+withinAlpha.DF2);
assert(abs((withinAlpha.Estimate(2)-withinAlpha.Estimate(1))-alphaDiff.Estimate(2))<1e-8);
assert(abs(withinAlpha.B3_B1_PowerRatio(2)/withinAlpha.B3_B1_PowerRatio(1) ...
    - alphaDiff.RatioOfPowerChangeRatios(2))<1e-8);
writetable(withinAlpha,fullfile(out,'alpha_within_group_B3_B1_ratios.csv'));
write_model(ma,out,'alpha');

disp('BETWEEN-GROUP CHANGE CONTRASTS');disp(alphaDiff);
disp('WITHIN-GROUP B3/B1 POWER RATIOS');disp(withinAlpha);
fprintf('Output: %s\n',out);
diary off;

function R=between_changes(m,gprefix,tprefix,gorder,torder,pairs)
R=table();
for j=1:size(pairs,1)
    late=torder{pairs(j,1)}; early=torder{pairs(j,2)};
    H=design(m.CoefficientNames,gprefix,tprefix,gorder{2},late,gorder,torder) ...
     -design(m.CoefficientNames,gprefix,tprefix,gorder{2},early,gorder,torder) ...
     -design(m.CoefficientNames,gprefix,tprefix,gorder{1},late,gorder,torder) ...
     +design(m.CoefficientNames,gprefix,tprefix,gorder{1},early,gorder,torder);
    R=[R;contrast_result(m,H,sprintf('(Offline %s-%s) - (Online %s-%s)',late,early,late,early))]; %#ok<AGROW>
end
end

function R=contrast_result(m,H,label)
estimate=H*fixedEffects(m);
SE=sqrt(max(H*m.CoefficientCovariance*H',0));
[p,F,df1,df2]=coefTest(m,H,0,'DFMethod','Satterthwaite');
assert(isfinite(SE)&&SE>0&&isfinite(df2)&&df2>0);
tc=tinv(.975,df2);
R=table(string(label),estimate,SE,estimate/SE,F,df1,df2,estimate-tc*SE,estimate+tc*SE,p, ...
    'VariableNames',{'Contrast','Estimate','SE','t','F','DF1','DF2','CI_Lower','CI_Upper','p_raw'});
end

function x=design(nm,gprefix,tprefix,g,t,gorder,torder)
% Effects-coded contributions. Intercept/tone/baseline cancel in changes.
x=zeros(1,numel(nm));
for j=1:numel(nm)
    name=char(nm{j});
    gt=regexp(name,[gprefix '_([^:]+)'],'tokens','once');
    tt=regexp(name,[tprefix '_([^:]+)'],'tokens','once');
    if ~isempty(gt)&&~isempty(tt), x(j)=code(g,gt{1},gorder)*code(t,tt{1},torder);
    elseif ~isempty(gt), x(j)=code(g,gt{1},gorder);
    elseif ~isempty(tt), x(j)=code(t,tt{1},torder);
    end
end
end

function v=code(level,coefficient,order)
if strcmp(level,coefficient),v=1;elseif strcmp(level,order{end}),v=-1;else,v=0;end
end

function p_adj=holm_adjust(p)
[sorted,index]=sort(p(:));
adjusted=min(1,cummax((numel(p)-(1:numel(p))'+1).*sorted));
p_adj=zeros(size(p(:)));p_adj(index)=adjusted;
end

function check_change_predictions(m,T,gvar,tvar,groups,R,pairs,levels)
% Independent check using MATLAB predict, not the contrast-vector parser.
for j=1:size(pairs,1)
    Q=T(ones(4,1),:);
    Q.(gvar)=categorical({groups{2};groups{2};groups{1};groups{1}},groups);
    Q.(tvar)=categorical({levels{pairs(j,1)};levels{pairs(j,2)}; ...
        levels{pairs(j,1)};levels{pairs(j,2)}},levels);
    pred=predict(m,Q,'Conditional',false);
    assert(abs([1 -1 -1 1]*pred-R.Estimate(j))<1e-8,'Contrast/prediction mismatch.');
end
end

function write_model(m,out,name)
tests=anova(m,'DFMethod','Satterthwaite');
if isa(tests,'dataset'),tests=dataset2table(tests);end
tests.ApproxPartialEtaSquared=(tests.FStat.*tests.DF1)./(tests.FStat.*tests.DF1+tests.DF2);
writetable(tests,fullfile(out,[name '_TypeIII_ANOVA.csv']));
coefficients=m.Coefficients;
if isa(coefficients,'dataset'),coefficients=dataset2table(coefficients);end
writetable(coefficients,fullfile(out,[name '_coefficients.csv']));
rd=table(fitted(m),residuals(m),'VariableNames',{'Fitted','ConditionalResidual'});
writetable(rd,fullfile(out,[name '_residuals.csv']));
end

function check_original(R,path)
assert(exist(path,'file')==2,'Original results missing: %s',path);
old=readtable(path);
for field={'Estimate','SE','CI_Lower','CI_Upper','p_raw'}
    f=field{1};assert(max(abs(R.(f)-old.(f)))<1e-5,['Reproduction mismatch: ' f]);
end
end
