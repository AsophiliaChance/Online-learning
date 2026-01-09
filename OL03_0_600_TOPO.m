%% Grand average across subjects for 4 conditions (ERP is 1x4 cell)
clear; clc;

result_dir = 'D:\online_learning\equiprobability_EEG\results';

% 如果 FieldTrip 不在路径里，取消注释并填你的路径
% addpath('D:\MATLAB\fieldtrip-20240101');
ft_defaults;

files = dir(fullfile(result_dir, '*_erp.mat'));
if isempty(files)
    error('No *_erp.mat files found in: %s', result_dir);
end
fprintf('Found %d ERP files.\n', numel(files));

% subj{s,c} 存每个被试、每个条件的 timelock
subj = cell(numel(files), 4);
subj_names = cell(numel(files), 1);

for s = 1:numel(files)
    fpath = fullfile(files(s).folder, files(s).name);
    S = load(fpath, 'ERP');
    if ~isfield(S, 'ERP')
        error('File %s does not contain variable "ERP".', files(s).name);
    end

    ERP = S.ERP;
    subj_names{s} = files(s).name;

    % ---- ERP must be 1x4 cell ----
    if ~iscell(ERP) || numel(ERP) ~= 4
        error('ERP in %s is not 1x4 cell. Actual: %s (numel=%d)', ...
              files(s).name, class(ERP), numel(ERP));
    end

    for c = 1:4
        if isempty(ERP{c}) || ~isstruct(ERP{c}) || ~isfield(ERP{c}, 'avg')
            error('Invalid ERP{%d} in %s (expect FieldTrip timelock with .avg).', c, files(s).name);
        end
        subj{s,c} = ERP{c};
    end
end

%% Grand average per condition (across subjects)
GA = cell(1,4);

cfg = [];
cfg.keepindividual = 'no';   % 只输出整体平均（如果想保留被试维度：'yes'）

for c = 1:4
    inputs = subj(:,c);
    inputs = inputs(~cellfun('isempty', inputs));
    fprintf('Computing grand average: condition %d, n=%d\n', c, numel(inputs));
    GA{c} = ft_timelockgrandaverage(cfg, inputs{:});
end

%% Save
out_file = fullfile(result_dir, 'GA_allSubjects_4conds.mat');
%save(out_file, 'GA', 'subj_names', '-v7.3');
fprintf('Saved: %s\n', out_file);

%% （可选）把4个条件再平均成一个“总平均ERP”
 GA_all = ft_timelockgrandaverage(cfg, GA{:});
% save(fullfile(result_dir, 'GA_allSubjects_all4conds.mat'), 'GA_all', '-v7.3');
%% ===== 1) 4个条件再平均成“总平均ERP” =====
cfg = [];
cfg.keepindividual = 'no';
GA_all = ft_timelockgrandaverage(cfg, GA{:});

%% ===== 2) 准备 layout（尽量复用你给的 layout 写法） =====
example_tl = GA_all;  % 用总平均的 elec 即可（也可以换成任意一个条件/被试的 timelock）
cfg_layout = [];
cfg_layout.elec       = example_tl.elec;
cfg_layout.projection = 'polar';
cfg_layout.rotate     = 90;     % 若朝向不对，可改 90/180/270
layout = ft_prepare_layout(cfg_layout);

%% ===== 3) 定义 0–600ms 每50ms 一个时间窗（12个） =====
t0 = 0.0;   % seconds
t1 = 0.6;   % seconds
step = 0.05;
edges = t0:step:t1;            % 0,0.05,...,0.6 (13个边界)
nWin = numel(edges)-1;         % 12

% 每个窗 [edges(k), edges(k+1)]
win = [edges(1:end-1); edges(2:end)]';  % nWin x 2

%% ===== 4) 统一颜色范围 zlim（推荐：用0-600ms内最大绝对值） =====
tmask = (GA_all.time >= t0) & (GA_all.time <= t1);
ZLIM_ALL = max(abs(GA_all.avg(:, tmask)), [], 'all');
ZLIM_ALL = [-ZLIM_ALL, ZLIM_ALL];   % 对称范围

%% ===== 5) 画 12 张 topo 到一个 figure（3×4） =====
fh1 = figure('Color','w','Position',[100 100 1100 650]);
sgtitle(sprintf('Grand Average ERP — Topographies (0–600 ms, 50 ms bins)'), ...
    'FontWeight','bold','FontSize',12);

% （可选）如果你还想高亮某些电极，填这里（channel label cellstr）
% highlight_ch = clusters.Attend_N2b;   % 例如你之前那种
highlight_ch = [];  % 默认不高亮

for k = 1:nWin
    ax = subplot(3,4,k);
    pos = get(ax, 'Position');
    newPos = [pos(1) pos(2) pos(3)+0.01 pos(4)+0.01];
    axDraw = subplot('Position', newPos);

    cfg = [];
    cfg.layout      = layout;
    cfg.comment     = 'no';
    cfg.electrodes  = 'off';
    cfg.xlim        = win(k,:);        % 秒
    cfg.avgovertime = 'yes';
    cfg.zlim        = ZLIM_ALL;        % 统一颜色范围
    cfg.colorbar    = 'no';
    cfg.figure      = 'gcf';
    cfg.axes        = axDraw;
   

    % 可选：只高亮指定电极为黑点
    if ~isempty(highlight_ch)
        cfg.highlight        = 'on';
        cfg.highlightchannel = highlight_ch;
        cfg.highlightsymbol  = '.';
        cfg.highlightsize    = 5;
        cfg.highlightcolor   = [0 0 0];
    end

    ft_topoplotER(cfg, GA_all);

    title(sprintf('%.0f–%.0f ms', win(k,1)*1000, win(k,2)*1000), ...
        'FontWeight','bold','FontSize',10);
    axis(axDraw, 'square');
end

colormap(fh1, 'jet');

% 公共 colorbar（位置按你风格放右侧）
cb1 = colorbar('Position', [0.93 0.25 0.012 0.55]);
cb1.FontSize = 8;
ylabel(cb1, 'Amplitude (\muV)', 'FontWeight','bold', 'FontSize',8);
