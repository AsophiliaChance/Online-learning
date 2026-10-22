clear; clc;

%% Paths
inputFile  = 'D:\online_learning\production\all_rater.xlsx';
outputFile = 'D:\online_learning\production\all_rater_finalresult.xlsx';

%% Read Excel
T = readtable(inputFile, 'TextType', 'string');

% Column 1: stimulus name
stim = string(T{:,1});

% Columns 2-4: rater judgements (e.g., tone2, tone3)
r2 = string(T{:,2});
r3 = string(T{:,3});
r4 = string(T{:,4});

N = length(stim);

%% Preallocate
subject_id = strings(N,1);
session_id = strings(N,1);
final_code = strings(N,1);

%% Process row by row
for i = 1:N
    
    %% ---------- Parse stimulus name ----------
    % Example: 01_test1_11.wav
    name = erase(stim(i), ".wav");
    parts = split(name, "_");
    
    if numel(parts) >= 2
        subject_id(i) = parts(1);
        session_id(i) = parts(2);
    else
        subject_id(i) = "";
        session_id(i) = "";
    end

    %% ---------- Clean rater entries ----------
    ratings = [r2(i), r3(i), r4(i)];
    ratings = strtrim(ratings);      % trim spaces
    ratings = lower(ratings);        % lower case
    
    % Optional: harmonise missing values
    ratings(ismissing(ratings)) = "";
    ratings(ratings == "nan") = "";
    ratings(ratings == "<missing>") = "";

    %% ---------- Majority vote ----------
    % If two or more raters agree -> majority category
    % Otherwise -> ambiguous
    
    u = unique(ratings);
    
    % Remove empty strings here if they should not count as votes
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

%% Write output
Tout = table(stim, subject_id, session_id, r2, r3, r4, final_code, ...
    'VariableNames', {'stimulus','subject_id','session_id','rater2','rater3','rater4','final_code'});

writetable(Tout, outputFile);

disp('Done. Results saved to:');
disp(outputFile);