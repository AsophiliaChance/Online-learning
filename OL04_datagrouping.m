%% Split ERP files into 1way / 2way (ERP is 1x4 cell)
clear; clc;

rootDir = 'D:\online_learning\equiprobability_EEG\results';

% 1way ID rules
is1way_id = @(id) ( (id>=1  && id<=6)  || ...
                    (id>=14 && id<=20) || ...
                    (id>=28 && id<=29) || ...
                    (id==33) );

% Containers: each row = subject, each col = ERP component
Post1way = {};
Post2way = {};
Pre1way  = {};
Pre2way  = {};

fprintf('========== FILE GROUP ASSIGNMENT ==========\n');

%% ---------------- POST ----------------
postFiles = dir(fullfile(rootDir, '*post_erp.mat'));
[~, idx] = sort({postFiles.name});
postFiles = postFiles(idx);

for i = 1:numel(postFiles)
    fname = postFiles(i).name;
    fpath = fullfile(postFiles(i).folder, fname);

    % extract ID from 3rd-4th characters
    idStr = fname(3:4);
    id = str2double(idStr);

    if isnan(id)
        fprintf('[SKIP] %s (invalid ID: "%s")\n', fname, idStr);
        continue;
    end

    S = load(fpath, 'ERP');
    if ~isfield(S, 'ERP') || ~iscell(S.ERP) || numel(S.ERP) ~= 4
        fprintf('[SKIP] %s (ERP must be 1x4 cell)\n', fname);
        continue;
    end

    if is1way_id(id)
        Post1way(end+1, :) = S.ERP;
        fprintf('[POST][1way]  ID=%02d  %s\n', id, fname);
    else
        Post2way(end+1, :) = S.ERP;
        fprintf('[POST][2way]  ID=%02d  %s\n', id, fname);
    end
end

%% ---------------- PRE ----------------
preFiles = dir(fullfile(rootDir, '*pre__erp.mat'));
[~, idx] = sort({preFiles.name});
preFiles = preFiles(idx);

for i = 1:numel(preFiles)
    fname = preFiles(i).name;
    fpath = fullfile(preFiles(i).folder, fname);

    idStr = fname(3:4);
    id = str2double(idStr);

    if isnan(id)
        fprintf('[SKIP] %s (invalid ID: "%s")\n', fname, idStr);
        continue;
    end

    S = load(fpath, 'ERP');
    if ~isfield(S, 'ERP') || ~iscell(S.ERP) || numel(S.ERP) ~= 4
        fprintf('[SKIP] %s (ERP must be 1x4 cell)\n', fname);
        continue;
    end

    if is1way_id(id)
        Pre1way(end+1, :) = S.ERP;
        fprintf('[PRE ][1way]  ID=%02d  %s\n', id, fname);
    else
        Pre2way(end+1, :) = S.ERP;
        fprintf('[PRE ][2way]  ID=%02d  %s\n', id, fname);
    end
end

%% ---------------- SUMMARY ----------------
fprintf('\n========== SUMMARY ==========\n');
fprintf('POST 1way: %d subjects\n', size(Post1way,1));
fprintf('POST 2way: %d subjects\n', size(Post2way,1));
fprintf('PRE  1way: %d subjects\n', size(Pre1way,1));
fprintf('PRE  2way: %d subjects\n', size(Pre2way,1));

%% ---------------- SAVE ----------------
save(fullfile(rootDir, 'ERP_split_by_id.mat'), ...
    'Post1way','Post2way','Pre1way','Pre2way','-v7.3');
