%% ============================================================
% Questionnaire analysis (NO Open_feedback)
% Compatible with OLD MATLAB
%% ============================================================

clear; clc;

%% -------------------- Paths --------------------
infile = 'D:\online_learning\background_scales.xlsx';
outdir = 'D:\online_learning\questionnaire_stats';

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

outMerged   = fullfile(outdir, 'merged_data.csv');
outCont     = fullfile(outdir, 'continuous_stats.xlsx');
outBinary   = fullfile(outdir, 'strategy_stats.xlsx');
outTxt      = fullfile(outdir, 'report.txt');

%% -------------------- Variables --------------------
contVars = { ...
    'SPS','BAI', ...
    'Motivation','Attitude_to_Mandarin','Mood_during_learning', ...
    'Alertness','Enjoyment','Perceived_learning', ...
    'Pronunciation_difficulty','Test_stress', ...
    'Tutor_evaluation','Lack_of_interaction_impact'};

%% ============================================================
%% 1. READ BACKGROUND
%% ============================================================
bg = readtable(infile, 'Sheet', 'Background_information');
bg.Properties.VariableNames = matlab.lang.makeValidName(bg.Properties.VariableNames);

% ID → 数字编号
bg.SubjNum = nan(height(bg),1);
for i = 1:height(bg)
    idstr = char(string(bg.ID{i}));
    tok = regexp(idstr, '(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        bg.SubjNum(i) = str2double(tok{1});
    end
end

%% ============================================================
%% 2. READ TRAINING EXPERIENCE（核心稳定版）
%% ============================================================
[~,~,teCell] = xlsread(infile, 'training experience');

header = teCell(1,2:end);
items  = teCell(2:end,1);

subj = nan(1,length(header));
for j = 1:length(header)
    h = header{j};
    if isnumeric(h)
        subj(j) = h;
    else
        tok = regexp(char(h),'(\d+)','tokens','once');
        if ~isempty(tok)
            subj(j) = str2double(tok{1});
        end
    end
end

te = table();
te.SubjNum = subj(:);

for r = 1:length(items)
    name = matlab.lang.makeValidName(strtrim(char(items{r})));
    vals = teCell(r+1,2:end)';

    num = nan(length(vals),1);

    for k = 1:length(vals)
        v = vals{k};

        if isnumeric(v)
            num(k) = v;

        elseif ischar(v)
            tmp = str2double(strtrim(v));
            if ~isnan(tmp)
                num(k) = tmp;
            end

        elseif isstring(v)
            tmp = str2double(char(v));
            if ~isnan(tmp)
                num(k) = tmp;
            end
        end
    end

    te.(name) = num;
end

%% ============================================================
%% 3. MERGE
%% ============================================================
T = outerjoin(bg, te, 'Keys','SubjNum','MergeKeys',true);
T = sortrows(T,'SubjNum');

writetable(T, outMerged);

%% ============================================================
%% 4. CONTINUOUS ANALYSIS
%% ============================================================
res = table();

g = T.Group;
groups = unique(g(~isnan(g)));

for i = 1:length(contVars)
    v = contVars{i};

    if ~ismember(v, T.Properties.VariableNames)
        continue;
    end

    x = T.(v);

    x1 = x(g==groups(1));
    x2 = x(g==groups(2));

    x1 = x1(~isnan(x1));
    x2 = x2(~isnan(x2));

    n1 = length(x1); n2 = length(x2);

    m1 = mean(x1); m2 = mean(x2);
    sd1 = std(x1); sd2 = std(x2);

    % Welch t-test
    [~,p,~,st] = ttest2(x1,x2,'Vartype','unequal');

    % ranksum
    pr = ranksum(x1,x2);

    % effect size
    sp = sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2));
    d  = (m1-m2)/sp;

    new = table({v},n1,m1,sd1,n2,m2,sd2,st.tstat,p,pr,d,...
        'VariableNames',{'Var','n1','m1','sd1','n2','m2','sd2','t','p','p_rank','d'});

    res = [res; new];
end

writetable(res, outCont);

%% ============================================================
%% 5. STRATEGY_USE（二分类）
%% ============================================================
if ismember('Strategy_use',T.Properties.VariableNames)

    y = T.Strategy_use;

    y1 = y(g==groups(1));
    y2 = y(g==groups(2));

    y1 = y1(~isnan(y1));
    y2 = y2(~isnan(y2));

    yes1 = sum(y1==1); n1 = length(y1);
    yes2 = sum(y2==1); n2 = length(y2);

    tbl = [yes1 n1-yes1; yes2 n2-yes2];

    [~,pf] = fishertest(tbl);

    out = table(yes1,n1,yes2,n2,pf,...
        'VariableNames',{'yes1','n1','yes2','n2','Fisher_p'});

    writetable(out,outBinary);
end

%% ============================================================
%% 6. TEXT REPORT
%% ============================================================
fid = fopen(outTxt,'w');

fprintf(fid,'Group sizes:\n');
fprintf(fid,'Group1: %d\n',sum(g==groups(1)));
fprintf(fid,'Group2: %d\n\n',sum(g==groups(2)));

for i=1:height(res)
    fprintf(fid,'%s\n',res.Var{i});
    fprintf(fid,'  G1: M=%.2f SD=%.2f\n',res.m1(i),res.sd1(i));
    fprintf(fid,'  G2: M=%.2f SD=%.2f\n',res.m2(i),res.sd2(i));
    fprintf(fid,'  t=%.2f p=%.4f rank=%.4f d=%.2f\n\n',...
        res.t(i),res.p(i),res.p_rank(i),res.d(i));
end

fclose(fid);

disp('=== DONE ===');