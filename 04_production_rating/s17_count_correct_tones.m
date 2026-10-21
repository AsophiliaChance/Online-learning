clear; clc;

%% 文件
inputFile  = 'D:\online_learning\production\all_rater_finalresult.xlsx';
outputFile = 'D:\online_learning\production\tone_count_summary.xlsx';

%% 读取数据
T = readtable(inputFile,'TextType','string');

subj = string(T.subject_id);
sess = string(T.session_id);
tone = lower(string(T.final_code));

%% 分组 (被试 × session)
[G, subject, session] = findgroups(subj, sess);

%% 统计 tone 数量
tone2 = splitapply(@(x) sum(x=="tone2"), tone, G);
tone3 = splitapply(@(x) sum(x=="tone3"), tone, G);
tone4 = splitapply(@(x) sum(x=="tone4"), tone, G);

%% 根据 session 设置上限
limit = 15 * ones(length(session),1);   % 默认 test = 15

training_mask = startsWith(session,"training");
limit(training_mask) = 21;

%% 截断
tone2 = min(tone2, limit);
tone3 = min(tone3, limit);
tone4 = min(tone4, limit);

%% 输出
Result = table(subject,session,tone2,tone3,tone4);

writetable(Result,outputFile);

disp("统计完成，结果保存到：");
disp(outputFile);