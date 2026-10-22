clear; clc;

%% Files
inputFile  = 'D:\online_learning\production\all_rater_finalresult.xlsx';
outputFile = 'D:\online_learning\production\tone_count_summary.xlsx';

%% Read data
T = readtable(inputFile,'TextType','string');

subj = string(T.subject_id);
sess = string(T.session_id);
tone = lower(string(T.final_code));

%% Group by participant x session
[G, subject, session] = findgroups(subj, sess);

%% Count tones
tone2 = splitapply(@(x) sum(x=="tone2"), tone, G);
tone3 = splitapply(@(x) sum(x=="tone3"), tone, G);
tone4 = splitapply(@(x) sum(x=="tone4"), tone, G);

%% Set maximum per session
limit = 15 * ones(length(session),1);   % default: test = 15

training_mask = startsWith(session,"training");
limit(training_mask) = 21;

%% Cap counts
tone2 = min(tone2, limit);
tone3 = min(tone3, limit);
tone4 = min(tone4, limit);

%% Output
Result = table(subject,session,tone2,tone3,tone4);

writetable(Result,outputFile);

disp("Counting done. Results saved to:");
disp(outputFile);