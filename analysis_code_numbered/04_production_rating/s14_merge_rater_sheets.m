%% merge_rater1_audio_tonecat.m
% Read all Excel files under:
%   D:\online_learning\production\rater1
% and merge columns: audio, tone_cat
% into ONE big table, then save to CSV/XLSX.

clear; clc;

inDir = 'D:\online_learning\production\rater3';

% Output
outCsv  = fullfile(inDir, 'rater3_ALL_audio_tone_cat.csv');
outXlsx = fullfile(inDir, 'rater3_ALL_audio_tone_cat.xlsx');

% Excel patterns
patterns = {'*.xlsx','*.csv','*.xlsm'};

% Collect file list
files = [];
for i = 1:numel(patterns)
    files = [files; dir(fullfile(inDir, patterns{i}))]; %#ok<AGROW>
end
files = files(~startsWith({files.name}, '~$')); % skip temp Excel lock files

if isempty(files)
    error('No Excel files found in: %s', inDir);
end

BIG = table();

for k = 1:numel(files)
    fpath = fullfile(inDir, files(k).name);
    fprintf('Reading (%d/%d): %s\n', k, numel(files), files(k).name);

    % Use detectImportOptions for robustness
    try
        opts = detectImportOptions(fpath, 'NumHeaderLines', 0);
        T = readtable(fpath, opts);
    catch
        % Fallback for old MATLAB / weird files
        T = readtable(fpath);
    end

    % Normalize variable names (case-insensitive match)
    v = lower(string(T.Properties.VariableNames));
    idx_audio = find(v == "audio", 1);
    idx_tone  = find(v == "tone_cat", 1);

    if isempty(idx_audio) || isempty(idx_tone)
        warning('Skip %s: missing columns. Found vars: %s', ...
            files(k).name, strjoin(T.Properties.VariableNames, ', '));
        continue;
    end

    audio = T{:, idx_audio};
    tone  = T{:, idx_tone};

    % Convert to string for consistent merge
    audio = string(audio);
    tone  = string(tone);

    % Keep only non-empty rows (optional)
    keep = ~(audio == "" & tone == "");
    audio = audio(keep);
    tone  = tone(keep);

    Tadd = table( ...
        repmat(string(files(k).name), numel(audio), 1), ...
        audio(:), tone(:), ...
        'VariableNames', {'source_file','audio','tone_cat'} );

    BIG = [BIG; Tadd]; %#ok<AGROW>
end

if isempty(BIG)
    error('No usable rows merged. Check that Excel files contain columns: audio, tone_cat');
end

% Save
writetable(BIG, outCsv);
fprintf('Saved CSV: %s\n', outCsv);

try
    writetable(BIG, outXlsx);
    fprintf('Saved XLSX: %s\n', outXlsx);
catch ME
    warning('Failed to save XLSX (%s). CSV already saved.', ME.message);
end

fprintf('Done. Total rows: %d\n', height(BIG));