%% ============================================================
% Plot ONLY the condition- and subject-averaged alpha topography
%
% Input:
%   alpha_condition_subject_average_channel.csv
%   Required columns: Channel, MeanAlphaPower
%
% Output:
%   alpha_topography_only_panel_c_600dpi.tif
%
% Changes relative to the full analysis script:
%   1) Does not load or process EEG datasets.
%   2) Does not calculate or save any CSV/MAT data.
%   3) ROI marker size is reduced from 6 to 3 (half size).
%   4) The colorbar is shifted to the left.
%   5) A lowercase panel label "c" is added at the upper-left corner.
%% ============================================================

clear; clc;

%% ---------------- PATHS ----------------
outDir = 'D:\online_learning\spectral_theta_alpha\alpha_topography_condition_subject_average';

inputCsv = fullfile(outDir, ...
    'alpha_condition_subject_average_channel.csv');

outTif = fullfile(outDir, ...
    'alpha_topography_only_panel_c_600dpi.tif');

% EGI HydroCel 129-channel location file; E1-E128 are used below.
lookupFile = ['S:\Program\matlab2019toolbox\eeglab_current\eeglab2024.0\', ...
              'plugins\dipfit\standard_BESA\GSN-HydroCel-129.sfp'];

%% ---------------- SETTINGS ----------------
expected_labels = arrayfun(@(x) sprintf('E%d', x), ...
    1:128, 'UniformOutput', false);

% Posterior alpha ROI.
alpha_roi = {'E70','E71','E72','E75','E76','E77','E67','E83'};

font_name = 'Times New Roman';
font_size = 12;
panel_font_size = 16;

fig_width_cm = 9.5;
fig_height_cm = 8.5;

% Original size was 6; use half of it here.
roi_marker_size = 3;

colorbar_label = '\alpha power (log_{10}(\muV^2/Hz))';

% Set to true only when a title is needed.
show_title = false;

%% ---------------- CHECK REQUIREMENTS ----------------
assert(exist(inputCsv, 'file') == 2, ...
    'Cannot find input CSV:\n%s', inputCsv);
assert(exist(lookupFile, 'file') == 2, ...
    'Cannot find channel-location file:\n%s', lookupFile);

if isempty(which('topoplot')) || isempty(which('readlocs'))
    error(['EEGLAB functions topoplot/readlocs are not on the MATLAB path.\n', ...
           'Add the EEGLAB folder to the MATLAB path before running this script.']);
end

%% ---------------- LOAD CHANNEL-LEVEL AVERAGE ----------------
T = readtable(inputCsv);

requiredVars = {'Channel', 'MeanAlphaPower'};
for iVar = 1:numel(requiredVars)
    if ~ismember(requiredVars{iVar}, T.Properties.VariableNames)
        error('Input CSV is missing the required column: %s', requiredVars{iVar});
    end
end

csv_labels = upper(strtrim(string(T.Channel)));
csv_values = double(T.MeanAlphaPower);

% Create one value for each EGI channel in the fixed E1-E128 order.
map_vector = nan(numel(expected_labels), 1);
for ch = 1:numel(expected_labels)
    idx = find(strcmpi(csv_labels, expected_labels{ch}), 1, 'first');
    if ~isempty(idx)
        map_vector(ch) = csv_values(idx);
    end
end

missing_channels = expected_labels(~isfinite(map_vector));
if ~isempty(missing_channels)
    warning('%d channels have missing/non-finite values and will be interpolated by topoplot.', ...
        numel(missing_channels));
end

finite_values = map_vector(isfinite(map_vector));
if isempty(finite_values)
    error('No finite MeanAlphaPower values were found in the input CSV.');
end

map_limits = [min(finite_values), max(finite_values)];
if diff(map_limits) <= 0
    expansion = 0.05 * max(abs(map_limits(1)), 1);
    map_limits = map_limits + [-expansion, expansion];
end

%% ---------------- BUILD E1-E128 CHANNEL LOCATIONS ----------------
lookupLocs = readlocs(lookupFile);
lookupLabels = upper(string({lookupLocs.labels}));

chanlocs128 = repmat(lookupLocs(1), 1, numel(expected_labels));
for ch = 1:numel(expected_labels)
    idx = find(strcmpi(lookupLabels, expected_labels{ch}), 1, 'first');
    if isempty(idx)
        error('Channel %s was not found in the lookup file.', expected_labels{ch});
    end
    chanlocs128(ch) = lookupLocs(idx);
end

roi_idx = find(ismember( ...
    upper(string(expected_labels)), upper(string(alpha_roi))));

%% ---------------- DRAW TOPOGRAPHY ONLY ----------------
fig = figure('Color', 'w');
set(fig, ...
    'Units', 'centimeters', ...
    'Position', [2, 2, fig_width_cm, fig_height_cm], ...
    'PaperPositionMode', 'auto');

% Leave enough space for the left-shifted colorbar.
ax = axes('Parent', fig, ...
    'Units', 'normalized', ...
    'Position', [0.07, 0.10, 0.66, 0.80]);

axes(ax); %#ok<LAXES>
topoplot(map_vector, chanlocs128, ...
    'maplimits', map_limits, ...
    'style', 'map', ...
    'electrodes', 'off', ...
    'numcontour', 6, ...
    'whitebk', 'on', ...
    'emarker2', {roi_idx, 'o', 'k', roi_marker_size, 1.2});

if show_title
    title('Condition- and subject-averaged \alpha power', ...
        'Interpreter', 'tex', ...
        'FontName', font_name, ...
        'FontSize', 14, ...
        'FontWeight', 'normal');
end

% Create the colorbar and then explicitly reset both positions, because
% MATLAB may resize the axes automatically when colorbar is created.
cb = colorbar('peer', ax);
set(ax, ...
    'Units', 'normalized', ...
    'Position', [0.0005, 0.10, 0.66, 0.80]);
set(cb, ...
    'Units', 'normalized', ...
    'Position', [0.555, 0.20, 0.028, 0.60], ...
    'FontName', font_name, ...
    'FontSize', font_size);

try
    cb.Label.String = colorbar_label;
    cb.Label.Interpreter = 'tex';
    cb.Label.FontName = font_name;
    cb.Label.FontSize = font_size;
catch
    % Compatibility with older MATLAB graphics handles.
    hcb_label = get(cb, 'YLabel');
    set(hcb_label, ...
        'String', colorbar_label, ...
        'Interpreter', 'tex', ...
        'FontName', font_name, ...
        'FontSize', font_size);
end

% Panel label at the upper-left corner of the whole figure.
annotation(fig, 'textbox', [0.015, 0.925, 0.08, 0.06], ...
    'String', 'c', ...
    'LineStyle', 'none', ...
    'Margin', 0, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontName', font_name, ...
    'FontSize', panel_font_size, ...
    'FontWeight', 'bold', ...
    'Color', 'k');

%% ---------------- SAVE FIGURE ----------------
try
    exportgraphics(fig, outTif, 'Resolution', 600);
catch
    print(fig, outTif, '-dtiff', '-r600');
end

fprintf('\nSaved topography-only figure:\n%s\n', outTif);
fprintf('ROI marker size: %g (half of the original size 6).\n', roi_marker_size);
fprintf('No EEG processing, CSV output, or MAT output was performed.\n');
