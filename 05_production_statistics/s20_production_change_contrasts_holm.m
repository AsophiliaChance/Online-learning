%% All six pairwise between-group differences in production change,
%% plus within-group Post-Pre changes (added 2026-10-18).
% Within-group family: Online/Offline Post-Pre (test4-test1), Holm2,
% following the within-group B3-B1 approach of alpha_change_analysis_with_effect_sizes.m.
% Reporting: Post-Pre in main text; other five contrasts in Supplementary.
% Model unchanged: all four sessions, empirical-logit outcome, tone covariate,
% participant random intercept, REML, effects coding, Satterthwaite df.
% Post-Pre is the main contrast (raw p); the other five form a Holm family.
% This reporting distinction does not establish original prespecification.
% The six contrasts are between-group change differences; the added
% within-group Post-Pre tests are reported as a separate Holm2 family.
% Confidence intervals are unadjusted pointwise 95% intervals.
% Does not modify any prior script or result directory.
clear; clc;
base='D:\online_learning';
% Results are saved inside the project folder F:\Documents\ChatGPT\online.
projectDir='F:\Documents\ChatGPT\online';
if ~exist(projectDir,'dir')
    projectDir=fileparts(mfilename('fullpath'));
end
assert(~isempty(projectDir) && exist(projectDir,'dir')==7,'Project folder not found.');
parent=fullfile(projectDir,'production_main_supplementary_results');
if ~exist(parent,'dir'),mkdir(parent);end
out=fullfile(parent,['results_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(out,'dir'),'Output exists; rerun after one second.');
mkdir(out);
diary(fullfile(out,'analysis_log.txt'));
P=readtable(fullfile(base,'production','production_simple_effects_revised','production_simple_effects_analysis_table.csv'));
P.subj=categorical(string(P.subj));
P.group=categorical(string(P.group),{'Online','Offline'});
levels={'test1','test2','test3','test4'};
P.session=categorical(string(P.session),levels);
P.tone=categorical(string(P.tone),{'tone2','tone3','tone4'});
assert(all(isfinite(P.logit_acc)) && ~any(isundefined(P.group)));
m=fitlme(P,'logit_acc ~ group*session + tone + (1|subj)', ...
    'FitMethod','REML','DummyVarCoding','effects');
pairs=fliplr(nchoosek(1:4,2));
R=between_changes(m,'group','session',{'Online','Offline'},levels,pairs);
assert(height(R)==6 && size(unique(sort(pairs,2),'rows'),1)==6);
mainIndex=pairs(:,1)==4 & pairs(:,2)==1;
R.p_Holm_supplementary=nan(6,1);
R.p_Holm_supplementary(~mainIndex)=holm_adjust(R.p_raw(~mainIndex));
R.p_reported=R.p_Holm_supplementary;
R.p_reported(mainIndex)=R.p_raw(mainIndex);
R.Significant_reported=R.p_reported<.05;
R.ApproxPartialEtaSquared=(R.F.*R.DF1)./(R.F.*R.DF1+R.DF2);
R.Correction=repmat(string('Holm across five supplementary comparisons'),6,1);
R.Correction(mainIndex)=string('Unadjusted main Post-Pre contrast');
R.ReportLocation=repmat(string('Supplementary'),6,1);
mainIndex=pairs(:,1)==4 & pairs(:,2)==1;
R.ReportLocation(mainIndex)=string('Main text');
assert(sum(mainIndex)==1);
check_change_predictions(m,P,'group','session',{'Online','Offline'},R,pairs,levels);
old=readtable(fullfile(base,'manuscript_revision','change_difference_results','production_change_difference.csv'));
for j=1:3
    idx=find(pairs(:,1)==j+1 & pairs(:,2)==1);
    for field={'Estimate','SE','CI_Lower','CI_Upper','p_raw'}
        f=field{1};assert(abs(R.(f)(idx)-old.(f)(j))<1e-5,'Original contrast not reproduced.');
    end
end
writetable(R,fullfile(out,'production_all_six_change_differences.csv'));
writetable(R(mainIndex,:),fullfile(out,'production_main_Post_minus_Pre.csv'));
writetable(R(~mainIndex,:),fullfile(out,'production_supplementary_other_five.csv'));

%% Within-group Post-Pre change in each group (same four-session model).
% Clarifies the direction and significance of change within each group.
% Separate family of two within-group tests, Holm2 (as for alpha B3-B1).
groupsW={'Online','Offline'};
withinProd=table();
for g=1:2
    H=design(m.CoefficientNames,'group','session',groupsW{g},'test4',groupsW,levels) ...
     -design(m.CoefficientNames,'group','session',groupsW{g},'test1',groupsW,levels);
    row=contrast_result(m,H,[groupsW{g} ': test4-test1']);
    row.Group=string(groupsW{g});
    row.Direction=string('increase');
    if row.Estimate<0, row.Direction=string('decrease'); end
    withinProd=[withinProd;row]; %#ok<AGROW>
end
withinProd.p_Holm2=holm_adjust(withinProd.p_raw);
withinProd.Significant_Holm2=withinProd.p_Holm2<.05;
withinProd.ApproxPartialEtaSquared=(withinProd.F.*withinProd.DF1)./(withinProd.F.*withinProd.DF1+withinProd.DF2);
% Consistency check: Offline change minus Online change equals the main Post-Pre contrast.
assert(abs((withinProd.Estimate(2)-withinProd.Estimate(1))-R.Estimate(mainIndex))<1e-8, ...
    'Within-group decomposition does not reproduce the main Post-Pre contrast.');
writetable(withinProd,fullfile(out,'production_within_group_Post_minus_Pre.csv'));
fprintf('WITHIN-GROUP POST-PRE (Holm across two groups)\n');disp(withinProd);

write_model(m,out,'production_four_session');
writetable(P,fullfile(out,'analysis_table.csv'));
fprintf('All pairwise contrasts: Offline change minus Online change.\n');disp(R);
fprintf('MAIN TEXT\n');disp(R(mainIndex,:));
fprintf('SUPPLEMENTARY\n');disp(R(~mainIndex,:));
fid=fopen(fullfile(out,'reporting_notes.txt'),'w');
fprintf(fid,'Main: Post-Pre; Supplementary: all remaining five temporal pairs.\n');
fprintf(fid,'Direction: Offline change minus Online change.\n');
fprintf(fid,'Effect units: empirical-logit; CIs: pointwise unadjusted 95%%.\n');
fprintf(fid,'Main Post-Pre is unadjusted. Holm correction covers ONLY the five supplementary contrasts.\n');
fprintf(fid,'All six contrasts use the SAME four-session model. No split-context models used here.\n');
fprintf(fid,'MATLAB predict cross-check and reproduction of original three pre-referenced contrasts passed.\n');
fprintf(fid,'Within-group Post-Pre (Online, Offline): separate family, Holm across two tests.\n');
fclose(fid);
diary off;
fprintf('Output: %s\n',out);

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
