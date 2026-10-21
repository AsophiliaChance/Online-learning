"""Figure 2 (revised): bar charts of model-predicted production accuracy and
posterior alpha power, with individual observations and within-group
significance brackets.

Modified from plot_model_trajectories_and_changes.py. Inputs are the same
MATLAB model exports used by that script, plus the within-group contrast
tables (Post-Pre for production, B3-B1 for alpha). The script does not refit
any model.

Run using the physio environment:
S:\\Program\\anaconda\\envs\\physio\\python.exe plot_figure2_bars.py
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

COLORS = {'Online': '#0072B2', 'Offline': '#D55E00'}
MARKERS = {'Online': 'o', 'Offline': 's'}
GROUPS = ['Online', 'Offline']


def fmt_p(p):
    return 'p < .001' if p < .001 else 'p = ' + f'{p:.3f}'.removeprefix('0')


def bar_panel(ax, means, time_col, levels, labels, mean, lower, upper,
              individual, ind_time, ind_value, scale=1.0, step=.62, gap=.55, width=.62):
    """Bars = model-predicted values with 95% CI; points = individual data.
    Online on the left, Offline on the right. Returns x positions per group."""
    rng = np.random.default_rng(2026)
    n = len(levels)
    xpos = {}
    shades = np.linspace(.35, .95, n)          # light -> dark across time
    for g_i, group in enumerate(GROUPS):
        x0 = g_i * (n * step + gap)
        xs = x0 + np.arange(n) * step
        xpos[group] = xs
        rows = means.loc[means.Group.eq(group)].set_index(time_col).loc[levels]
        y = rows[mean].to_numpy() * scale
        err = np.vstack([y - rows[lower].to_numpy() * scale,
                         rows[upper].to_numpy() * scale - y])
        for k in range(n):
            ax.bar(xs[k], y[k], width=width, color=COLORS[group], alpha=shades[k],
                   edgecolor=COLORS[group], linewidth=1.2, zorder=2)
        ax.errorbar(xs, y, yerr=err, fmt='none', ecolor='#202020',
                    elinewidth=2.8, capsize=5, capthick=2.8, zorder=5)
        pts = individual.loc[individual.Group.eq(group)]
        for k in range(n):
            v = pts.loc[pts[ind_time].eq(k + 1), ind_value].to_numpy() * scale
            ax.scatter(xs[k] + rng.uniform(-.17, .17, len(v)), v, s=96,
                       marker=MARKERS[group], facecolor='white',
                       edgecolor=COLORS[group], linewidth=1.0, alpha=.9, zorder=4)
    ticks = np.concatenate([xpos[g] for g in GROUPS])
    ax.set_xticks(ticks, labels * 2)
    for group in GROUPS:
        ax.text(xpos[group].mean(), -.17, group, transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=44, fontweight='bold',
                color=COLORS[group])
    ax.set_xlim(-.45, ticks.max() + .45)
    ax.spines[['top', 'right']].set_visible(False)
    ax.axhline(0, color='#202020', linewidth=1)
    return xpos


def bracket(ax, x1, x2, y, h, text):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color='#202020',
            linewidth=1.4, zorder=6)
    ax.text((x1 + x2) / 2, y + h * 1.4, text, ha='center', va='bottom',
            fontsize=40, zorder=6)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(r'D:\online_learning'))
    parser.add_argument('--project', type=Path,
                        default=Path(r'F:\Documents\ChatGPT\online'))
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    root, project = args.root, args.project
    output = args.output_dir or project / 'figure2_bars'
    output.mkdir(parents=True, exist_ok=True)

    # ---- model predictions (same exports as the original Figure 2)
    production = pd.read_csv(root / 'production/production_simple_effects_revised/'
                             'production_adjusted_means_group_session.csv')
    alpha = pd.read_csv(root / 'training_psychophysiology/alpha_VOFF_simple_effects/'
                        'alpha_adjusted_means_Group_Block.csv')
    alpha.Group = alpha.Group.replace({'Group2': 'Online', 'Group1': 'Offline'})
    # ---- individual observations
    raw_p = pd.read_csv(root / 'production/production_simple_effects_revised/'
                        'production_simple_effects_long_table.csv')
    raw_p = raw_p.groupby(['subj', 'group', 'session_num'], as_index=False).agg(
        correct=('correct', 'sum'), trials=('nTrials', 'sum'))
    raw_p['Accuracy'] = raw_p.correct / raw_p.trials
    raw_p = raw_p.rename(columns={'group': 'Group', 'session_num': 'SessionNum'})
    raw_a = pd.read_csv(root / 'training_psychophysiology/alpha_VOFF_LMM_TypeIII_randomEffects/'
                        'alpha_raw_subject_points.csv')
    # ---- within-group contrasts (Holm across the two groups)
    wp = pd.read_csv(project / 'production_main_supplementary_results/results_20261018_212433/'
                     'production_within_group_Post_minus_Pre.csv').set_index('Group')
    wa = pd.read_csv(project / 'alpha_analysis_effect_sizes/results_20261018_192905/'
                     'alpha_within_group_B3_B1_ratios.csv').set_index('Group')

    assert len(production) == 8 and len(alpha) == 6
    # Alpha bars and the within-group contrast must describe the same model.
    a = alpha.pivot(index='Block', columns='Group', values='AdjustedMean')
    for g in GROUPS:
        np.testing.assert_allclose(a.loc['B3', g] - a.loc['B1', g],
                                   wa.loc[g, 'Estimate'], atol=1e-6)

    plt.rcParams.update({'font.family': ['Times New Roman', 'Liberation Serif'],
                         'font.size': 40, 'axes.titlesize': 44, 'axes.labelsize': 40,
                         'xtick.labelsize': 36, 'ytick.labelsize': 38,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig = plt.figure(figsize=(26, 11))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=.32)
    sub = grid[0, 1].subgridspec(1, 2, width_ratios=[.74, .26], wspace=.02)
    axA = fig.add_subplot(grid[0, 0])
    axB = fig.add_subplot(sub[0, 0])
    axT = fig.add_subplot(sub[0, 1])

    # ---- A: production accuracy
    xs = bar_panel(axA, production, 'Session', ['test1', 'test2', 'test3', 'test4'],
                   ['Pre', 'Int.\n1', 'Int.\n2', 'Post'],
                   'AdjustedProbability', 'CI_Lower_probability', 'CI_Upper_probability',
                   raw_p, 'SessionNum', 'Accuracy', scale=100)
    axA.set_ylabel('Accuracy (%)')
    axA.set_ylim(0, 118)
    axA.set_yticks(np.arange(0, 101, 20))
    for g in GROUPS:
        p = wp.loc[g, 'p_Holm2']
        if p < .05:
            bracket(axA, xs[g][0], xs[g][3], 103, 3, f'* {fmt_p(p)}')

    # ---- B: posterior alpha power
    xs = bar_panel(axB, alpha, 'Block', ['B1', 'B2', 'B3'],
                   ['B1', 'B2', 'B3'],
                   'AdjustedMean', 'CI_Lower', 'CI_Upper',
                   raw_a, 'BlockNum', 'AlphaPower')
    axB.set_ylabel('Power (log10)')
    top = raw_a.AlphaPower.max()
    lo = raw_a.AlphaPower.min()
    span = top - lo
    axB.set_ylim(lo - .08 * span, top + .30 * span)
    for g in GROUPS:
        p = wa.loc[g, 'p_Holm2']
        if p < .05:
            bracket(axB, xs[g][0], xs[g][2], top + .06 * span, .04 * span, f'** {fmt_p(p)}')

    # ---- topography next to panel B (same EEGLAB map as before)
    topo_path = root / ('spectral_theta_alpha/alpha_topography_condition_subject_average/'
                        'alpha_topography_only_panel_c_600dpi.tif')
    with Image.open(topo_path) as image:
        topo = np.asarray(image.convert('RGB'))
    from scipy.ndimage import label, find_objects
    chroma = np.max(topo, axis=2).astype(float) - np.min(topo, axis=2).astype(float)
    components, _ = label(chroma > 35)
    largest = np.argmax(np.bincount(components.ravel())[1:]) + 1
    sy, sx = find_objects(components)[largest - 1]
    topo = topo[max(0, sy.start - 10):sy.stop + 10, max(0, sx.start - 10):sx.stop + 10]
    axT.imshow(topo, interpolation='none')
    axT.set_axis_off()
    axT.set_anchor('W')

    for ax, letter, title in [(axA, 'A', 'Production accuracy'),
                              (axB, 'B', 'Posterior α power')]:
        ax.set_title(title, pad=22, loc='left', fontweight='bold')
        ax.text(-.15 if ax is axA else -.2, 1.06, letter, transform=ax.transAxes,
                fontweight='bold', fontsize=48)
    fig.subplots_adjust(left=.085, right=.99, bottom=.20, top=.85)

    stem = output / 'Figure2_bars'
    fig.savefig(stem.with_suffix('.png'), dpi=400, facecolor='white')
    fig.savefig(stem.with_suffix('.pdf'), facecolor='white')
    fig.savefig(stem.with_suffix('.tif'), dpi=600, facecolor='white',
                pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
    print(stem.with_suffix('.png'))


if __name__ == '__main__':
    main()
