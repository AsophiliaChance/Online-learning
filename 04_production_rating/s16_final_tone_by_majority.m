clear; clc;

%% 路径
inputFile  = 'D:\online_learning\production\all_rater.xlsx';
outputFile = 'D:\online_learning\production\all_rater_finalresult.xlsx';

%% 读取 Excel
T = readtable(inputFile, 'TextType', 'string');

% 第一列：刺激名
stim = string(T{:,1});

% 第2/3/4列：rater结果（如 tone2, tone3）
r2 = string(T{:,2});
r3 = string(T{:,3});
r4 = string(T{:,4});

N = length(stim);

%% 预分配
subject_id = strings(N,1);
session_id = strings(N,1);
final_code = strings(N,1);

%% 逐行处理
for i = 1:N
    
    %% ---------- 解析刺激名 ----------
    % 例子: 01_test1_11.wav
    name = erase(stim(i), ".wav");
    parts = split(name, "_");
    
    if numel(parts) >= 2
        subject_id(i) = parts(1);
        session_id(i) = parts(2);
    else
        subject_id(i) = "";
        session_id(i) = "";
    end

    %% ---------- 清理 rater 内容 ----------
    ratings = [r2(i), r3(i), r4(i)];
    ratings = strtrim(ratings);      % 去掉前后空格
    ratings = lower(ratings);        % 统一成小写
    
    % 可选：把缺失值统一处理
    ratings(ismissing(ratings)) = "";
    ratings(ratings == "nan") = "";
    ratings(ratings == "<missing>") = "";

    %% ---------- 多数投票 ----------
    % 若有两个及以上相同 -> 输出多数那个
    % 否则 -> ambiguous
    
    u = unique(ratings);
    
    % 如果你不想把空字符串参与投票，就去掉它
    u(u == "") = [];
    
    if isempty(u)
        final_code(i) = "ambiguous";
    else
        counts = zeros(length(u),1);
        for j = 1:length(u)
            counts(j) = sum(ratings == u(j));
        end
        
        if max(counts) >= 2
            temp = u(counts == max(counts));
            final_code(i) = temp(1);
        else
            final_code(i) = "ambiguous";
        end
    end
    
end

%% 输出结果
Tout = table(stim, subject_id, session_id, r2, r3, r4, final_code, ...
    'VariableNames', {'stimulus','subject_id','session_id','rater2','rater3','rater4','final_code'});

writetable(Tout, outputFile);

disp('完成，结果已保存到：');
disp(outputFile);