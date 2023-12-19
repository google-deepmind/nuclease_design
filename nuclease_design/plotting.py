# Copyright 2023 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Library for figure plotting."""

from typing import Sequence

from matplotlib import gridspec
import matplotlib.pyplot as plt
from nuclease_design import amino_acids
from nuclease_design import constants
from nuclease_design import utils
import numpy as np
import pandas as pd
import seaborn as sns


WT_COLOR = sns.color_palette()[-3]  # grey
ML_COLOR = sns.color_palette()[1]  # orange
DE_COLOR = sns.color_palette()[0]  # blue
HR_COLOR = sns.color_palette()[2]  # green
ZERO_SHOT_COLOR = ML_COLOR
ML3_COLOR = ML_COLOR
ML2_COLOR = ML_COLOR
EPPCR_COLOR = DE_COLOR
DE3_COLOR = DE_COLOR

DEFAULT_MIN_NUM_OBSERVATIONS_PER_BIN = 20


def _set_fontsize(ax, fontsize) -> None:
  """Sets fontsize for `ax`."""
  for item in (
      [ax.title, ax.xaxis.label, ax.yaxis.label]
      + ax.get_xticklabels()
      + ax.get_yticklabels()
  ):
    item.set_fontsize(fontsize)

  legend = ax.get_legend()
  if not legend:
    return
  plt.setp(legend.get_texts(), fontsize=fontsize)
  plt.setp(legend.get_title(), fontsize=fontsize)


def plot_hit_rate_per_num_mutations(
    df: pd.DataFrame,
    reference_name: str,
    *,
    random_seed: np.random.RandomState,
    hue_feature: str = 'sublibrary_name',
    hue_order: Sequence[str] | None = None,
    palette: sns.color_palette(),
    group_cols: Sequence[str] = ('sublibrary_name', 'num_mutations'),
    num_bootstraps: int = utils.DEFAULT_NUM_BOOTSTRAPS,
    expected_false_discovery_rate: float = utils.EXPECTED_FDR,
    min_num_observations_per_bin: int = DEFAULT_MIN_NUM_OBSERVATIONS_PER_BIN,
    ax: plt.Axes | None = None,
    fontsize: int = 20,
):
  """Plots hit rates.

  Args:
    df: DataFrame with pvalue columns.
    reference_name: a string indicating the reference sequence for the "hit"
      comparison. (e.g. 'stop', 'wt', 'a73r')
    random_seed: controls randomness.
    hue_feature: the feature to use to color the hit rate plot, e.g. 'sequence
      type'.
    hue_order: the order to select for the palette.
    palette: the color palette to use for.
    group_cols: the features to group together for FDR correction (done
      independently on each group).
    num_bootstraps: the number of times to resample the dataset (to produce
      confidence bands).
    expected_false_discovery_rate: the FDR rate at which to control.
    min_num_observations_per_bin: The minimum number of samples at a given
      radius. Groupings with fewer samples are filtered from the plot.
    ax: Optional Axes for plotting.
    fontsize: Font size.

  Returns:
    A matplotlib Axes object
  """
  random_state = np.random.RandomState(random_seed)
  df = df.copy()
  df['pvalue'] = df.apply(
      utils.get_pvalue,
      reference_name=reference_name,
      axis=1,
  )
  # filter to minimum size
  df = (
      df.groupby(group_cols)
      .filter(lambda df: len(df) > min_num_observations_per_bin)
      .reset_index()
  )
  hit_rate_df = (
      df.groupby(group_cols)
      .apply(
          utils.get_bootstrapped_hitrate_df,
          pval_col='pvalue',
          expected_false_discovery_rate=expected_false_discovery_rate,
          random_state=random_state,
          num_bootstraps=num_bootstraps,
      )
      .reset_index()
  )
  ax = sns.lineplot(
      data=hit_rate_df,
      x='num_mutations',
      y='hit_rate',
      hue=hue_feature,
      hue_order=hue_order,
      errorbar='sd',
      err_style='band',
      palette=palette,
      ax=ax,
  )
  mean_hit_rate_df = (
      hit_rate_df.groupby(group_cols).agg('mean').reset_index(drop=False)
  )
  sns.scatterplot(
      data=mean_hit_rate_df,
      x='num_mutations',
      y='hit_rate',
      hue=hue_feature,
      hue_order=hue_order,
      palette=palette,
      s=50,
      legend=False,
  )
  ax.set_xticks(range(hit_rate_df.num_mutations.max() + 1))
  reference_to_title_str = {
      'wt': 'WT',
      'neg_control': '0',
      'a73r': 'A73R',
  }
  ax.set_title(
      f'Hit Rate for Activity > {reference_to_title_str[reference_name]}'
  )
  ax.set_ylim(bottom=0.0)
  ax.set_xlim(left=1)
  _set_fontsize(ax, fontsize)
  ax.legend(title=hue_feature, fontsize=fontsize, title_fontsize=fontsize)
  return ax


def plot_hit_rate_per_num_mutations_with_histogram(
    df: pd.DataFrame,
    reference_name: str,
    *,
    random_seed: np.random.RandomState,
    hue_feature: str = 'sublibrary_name',
    hue_order: Sequence[str] | None = None,
    palette: sns.color_palette(),
    group_cols: Sequence[str] = ('sublibrary_name', 'num_mutations'),
    num_bootstraps: int = utils.DEFAULT_NUM_BOOTSTRAPS,
    expected_false_discovery_rate: float = utils.EXPECTED_FDR,
    min_num_observations_per_bin: int = DEFAULT_MIN_NUM_OBSERVATIONS_PER_BIN,
    fontsize: int = 20) -> plt.Axes:
  """Plots hit rates including histogram for num_mutations below.

  Args:
    df: DataFrame with pvalue columns.
    reference_name: a string indicating the reference sequence for the "hit"
      comparison. (e.g. 'stop', 'wt', 'a73r')
    random_seed: controls randomness.
    hue_feature: the feature to use to color the hit rate plot, e.g. 'sequence
      type'.
    hue_order: the order to select for the palette.
    palette: the color palette to use for.
    group_cols: the features to group together for FDR correction (done
      independently on each group).
    num_bootstraps: the number of times to resample the dataset (to produce
      confidence bands).
    expected_false_discovery_rate: the FDR rate at which to control.
    min_num_observations_per_bin: The minimum number of samples at a given
      radius. Groupings with fewer samples are filtered from the plot.
    fontsize: Font size.

  Returns:
    A matplotlib Axes object
  """

  num_histplots = df[hue_feature].nunique()
  with sns.axes_style('ticks', {'axes.grid': False}):
    height_ratios = [
        15,
    ] + [
        2,
    ] * num_histplots
    gs = gridspec.GridSpec(
        nrows=1 + num_histplots, ncols=1, height_ratios=height_ratios
    )
    ax0 = plt.subplot(gs[0])
    plot_hit_rate_per_num_mutations(
        df,
        reference_name,
        random_seed=random_seed,
        hue_feature=hue_feature,
        hue_order=hue_order,
        palette=palette,
        group_cols=group_cols,
        num_bootstraps=num_bootstraps,
        expected_false_discovery_rate=expected_false_discovery_rate,
        min_num_observations_per_bin=min_num_observations_per_bin,
        ax=ax0,
        fontsize=fontsize,
    )
    ax0.legend(title=hue_feature, fontsize=fontsize, title_fontsize=fontsize)
    _set_fontsize(ax0, fontsize)
    for i, grouping in enumerate(hue_order):
      ax = plt.subplot(gs[i + 1], sharex=ax0)
      sns.histplot(
          data=df[df[hue_feature] == grouping],
          x='num_mutations',
          ax=ax,
          discrete=True,
          hue=hue_feature,
          element='step',
          hue_order=hue_order,
          palette=palette,
          stat='count',
      )
      ax.set_xlabel('')
      ax.set_ylabel('')
      ax.set(yticklabels=[])
      ax.tick_params('y', labelbottom=False)
      ax.tick_params('x', labelbottom=False)
      ax.legend_ = None
      ax.grid(False)
      _set_fontsize(ax, fontsize)

  sns.despine()
  plt.tight_layout(pad=1)
  remove_underscores_from_axis_labels(ax0)

  return ax0


def plot_overall_hit_rate(
    df: pd.DataFrame,
    reference_name: str,
    random_seed: int,
    expected_false_discovery_rate: float = utils.EXPECTED_FDR,
    num_bootstraps: int = utils.DEFAULT_NUM_BOOTSTRAPS) -> None:
  """Plot barplot of overall hit rates.

  Args:
    df: DataFrame with pvalue columns.
    reference_name: a string indicating the reference sequence for the "hit"
      comparison. (e.g. 'neg_control', 'wt', 'a73r')
    random_seed: controls randomness.
    expected_false_discovery_rate: the FDR rate at which to control.
    num_bootstraps: the number of times to resample the dataset (to produce
      confidence bands).
  """
  random_state = np.random.RandomState(random_seed)

  df = df.copy()
  df['pvalue'] = df.apply(
      utils.get_pvalue,
      reference_name=reference_name,
      axis=1,
  )

  hit_rate_df = (
      df.groupby('sublibrary_name')
      .apply(
          utils.get_bootstrapped_hitrate_df,
          pval_col='pvalue',
          expected_false_discovery_rate=expected_false_discovery_rate,
          random_state=random_state,
          num_bootstraps=num_bootstraps,
      )
      .reset_index()
  )
  sns.barplot(data=hit_rate_df, x='sublibrary_name', y='hit_rate')
  plt.title(
      f'hit rate for activity > {reference_name} (efdr ='
      f' {expected_false_discovery_rate})'
  )
  plt.xticks(rotation=45)
  plt.ylim(top=1.0)


def remove_underscores_from_axis_labels(ax: plt.Axes) -> None:
  """Removes underscores from axis labels."""
  ax.set_xlabel(ax.get_xlabel().replace('_', ' '))
  ax.set_ylabel(ax.get_ylabel().replace('_', ' '))


def plot_purified_protein_activity(
    df,
    library_order=('WT', 'epPCR', 'ML2', 'DE3', 'ML3'),
    custom_palette=(WT_COLOR, EPPCR_COLOR, ML2_COLOR, DE3_COLOR, ML3_COLOR),
) -> plt.Axes:
  """Plot purified protein activity.

  Args:
    df: A `DataFrame` with columns "library", "fold_change_activity", where
    "fold_change_activity" is normalized to the wild type.
    library_order: The left-to-right order for the bar chart.
    custom_palette: The colors for the systems in `library_order`.

  Returns:
    plt.Axes object
  """
  with sns.axes_style('ticks'):
    ax = sns.barplot(
        data=df,
        x='library',
        order=library_order,
        y='fold_change_activity',
        dodge=False,
        errorbar='sd',
        palette=custom_palette,
    )

  plt.xlabel('')
  plt.ylabel('Nuclease Activity Fold-Improvement')
  plt.ylim(bottom=0)
  plt.yticks([1, 5, 10, 15, 20, 25])
  plt.title('Activity of Best Variant')
  sns.despine()
  return ax


def make_diversity_overlay_plot(
    df,
    post_sort_population,
    hue_feature,
    hue_order,
    palette=sns.color_palette(),
    xticks_max=26,
):
  """Makes multi-element plot showing library diversity."""
  opacity_map = {'post-sort': 1.0, 'pre-sort': 0.2}
  interval_opacity_map = {'post-sort': 0.2, 'pre-sort': 0.1}

  with sns.axes_style('ticks'):
    df = df.sort_values([hue_feature, 'population']).copy()
    populations_to_plot = ['pre-sort', post_sort_population]
    plotdf = df[df['population'].isin(populations_to_plot)]

    ax = sns.lineplot(
        data=plotdf,
        x='max_intra_cluster_hamming_distance',
        y='num_clusters',
        hue=hue_feature,
        hue_order=hue_order,
        style='population',
        errorbar=('ci', 95),
        palette=palette,
    )

    ax.set(yscale='log')
    remove_underscores_from_axis_labels(ax)
    plt.ylabel('Number of Clusters')
    plt.xlabel('Cluster Diameter')
    sns.move_legend(ax, 'upper left', bbox_to_anchor=(1, 1))
    plt.ylim(bottom=1.0)
    plt.xlim(left=0.0, right=xticks_max - 1)
    ax.spines[['right', 'top']].set_visible(False)

    opacity_labels = [
        'post-sort',
        'pre-sort',
        'post-sort',
        'pre-sort',
    ]  # TODO(neilthomas) make this not manual

    for i, (line, poly) in enumerate(zip(ax.lines, ax.collections)):
      opacity_label = opacity_labels[i]
      line.set_alpha(opacity_map[opacity_label])
      poly.set_alpha(interval_opacity_map[opacity_label])

    plt.xticks(list(range(0, xticks_max, 5)))
  sns.despine()


def _add_mutations_to_count_matrix(counts, mutations, val, aa_to_index):
  for _, one_indexed_position, new_value in mutations:
    counts[one_indexed_position, aa_to_index[new_value]] = val
  return counts


def make_mutation_heatmap(
    mutations,
    ax,
    use_x_ticks,
    use_y_ticks,
    inclusive_one_indexed_position_start,
    exclusive_one_indexed_position_end,
    color,
    parent_sequence,
):
  """Makes a single mutation heatmap."""
  aa_sorted = sorted(
      amino_acids.AA, key=amino_acids.AA_TO_ISOELECTRIC_POINT.get, reverse=True
  )
  aa_to_index = {aa: i for i, aa in enumerate(aa_sorted)}
  positions = list(
      range(
          inclusive_one_indexed_position_start,
          exclusive_one_indexed_position_end,
      )
  )
  counts = np.zeros(shape=(max(positions) + 1, len(amino_acids.AA)))
  counts = _add_mutations_to_count_matrix(counts, mutations, 1, aa_to_index)

  df = pd.DataFrame(counts.T, index=aa_to_index.keys())

  light_color = sns.light_palette(color)[0]
  cmap = sns.color_palette([light_color, color], as_cmap=True)

  sns.heatmap(
      data=df,
      ax=ax,
      xticklabels=5 if use_x_ticks else False,
      yticklabels=use_y_ticks,
      linewidths=0.1,
      cbar=False,
      cmap=cmap,
  )
  ax.set(xlabel='', ylabel='')
  if parent_sequence is not None:
    parent_plot_indexes = []
    for one_indexed_position in range(
        inclusive_one_indexed_position_start, exclusive_one_indexed_position_end
    ):
      parent_aa = parent_sequence[one_indexed_position - 1]
      parent_plot_indexes.append((one_indexed_position, aa_to_index[parent_aa]))

    xs = [index[0] + 0.5 for index in parent_plot_indexes]
    ys = [index[1] + 0.5 for index in parent_plot_indexes]
    ax.plot(xs, ys, marker='.', color=WT_COLOR, linestyle='None')


def _select_hits_at_threshold(df, threshold):
  if threshold == 'pre-sort':
    return df
  return utils.select_hit_rows(df, threshold, utils.EXPECTED_FDR)


def make_mutation_heatmap_grid(
    df_a,
    df_b,
    a_name,
    b_name,
    parent_sequence,
    fiducial='wt',
    palette=sns.color_palette(),
) -> None:
  """Makes a grid of mutation heatmaps."""
  _, axs = plt.subplots(
      nrows=2, ncols=2, figsize=(65, 15), sharex=True, dpi=1000
  )
  for i, (lib_name, df) in enumerate(
      [(a_name, utils.filter_to_g4_positions(df_a)), (b_name, df_b)]
  ):
    for j, (pop_name, threshold) in enumerate(
        [('full', 'pre-sort'), (f'activity > {fiducial}', fiducial)]
    ):
      ax = axs[i][j]
      df_to_plot = _select_hits_at_threshold(df, threshold)
      make_mutation_heatmap(
          df_to_plot['mutations'].explode().dropna(),
          ax=ax,
          use_y_ticks=(j == 0),
          use_x_ticks=(i == 1),
          inclusive_one_indexed_position_start=constants.G4_VARIABLE_REGION_START,
          exclusive_one_indexed_position_end=constants.G4_VARIABLE_REGION_END,
          color=palette[[a_name, b_name].index(lib_name)],
          parent_sequence=parent_sequence,
      )
      ax.yaxis.set_tick_params(labelsize=30)
      ax.xaxis.set_tick_params(labelsize=50)
      n_mut = df_to_plot['mutations'].explode().nunique()
      lib_name_for_title = lib_name.replace('_', '-')
      if pop_name == 'full':
        title_prefix = f'full {lib_name_for_title} library'
      else:
        title_prefix = f'{lib_name_for_title} hits ( > {threshold})'
      ax.set_title(
          f'{title_prefix}: {len(df_to_plot)} variants, {n_mut} distinct'
          ' mutations',
          size=60,
      )
  plt.xlim(
      left=constants.G4_VARIABLE_REGION_START,
      right=constants.G4_VARIABLE_REGION_END,
  )
  plt.tight_layout()


def plot_subsampling_hit_rate(
    df, reference_name, score_specs, ax=None, fontsize=20, title=None
) -> None:
  """Plots hit rate when sub-sampling by zero-shot score."""
  reference_description = (
      'functional' if reference_name == 'neg_control' else '> WT'
  )
  df = df.copy()
  pval_col = utils.get_pvalue_column_name('g1', reference_name)
  df['is_functional'] = utils.select_hits(df[pval_col], utils.EXPECTED_FDR)

  rows = []
  for model_name, score_column_for_sorting in score_specs:
    sorted_df = df.sort_values(score_column_for_sorting, ascending=False)
    total_num_hits = sorted_df['is_functional'].sum()
    base_rate = sorted_df['is_functional'].mean()
    for n in np.linspace(50, len(sorted_df), 100):
      hits = sorted_df.head(int(n))['is_functional']
      hit_rate = hits.mean()
      fraction_of_total_hits_recovered = hits.sum() / total_num_hits
      enrichment = hit_rate / base_rate
      rows.append(
          dict(
              library_subsampling_fraction=(n / len(sorted_df)),
              library_size=n,
              hit_rate=hit_rate,
              enrichment=enrichment,
              fraction_of_total_hits_recovered=fraction_of_total_hits_recovered,
              model=model_name,
          )
      )
  ranking_df = pd.DataFrame(rows)
  ranking_df['library_subsampling_percentage'] = (
      1 - ranking_df['library_subsampling_fraction']
  ) * 100
  ranking_df['hit_recovery_percentage'] = (
      ranking_df['fraction_of_total_hits_recovered'] * 100
  )
  with sns.axes_style('ticks'):
    sns.scatterplot(
        data=ranking_df,
        x='library_subsampling_percentage',
        y='hit_recovery_percentage',
        hue='model' if len(score_specs) > 1 else None,
        ax=ax,
    )
    ax = ax if ax else plt.gca()

    lims = np.linspace(0, 100, 100)
    ax.plot(lims[::-1], lims, 'k-', alpha=0.25, zorder=0, label='random')

    ax.set_xlabel('% reduction in library size')
    ax.set_ylabel(f'% of {reference_description}\nvariants recovered')
    if title:
      ax.set_title(title)
    if len(score_specs) > 1:
      ax.legend(loc='lower left')
    _set_fontsize(ax, fontsize=fontsize)


def plot_zero_shot_histograms(
    eppcr_df, score_column, ax=None, title=None, use_legend=True, fontsize=15
) -> None:
  """Plots zero-shot histograms."""
  pval_col = utils.get_pvalue_column_name('g1', 'neg_control')
  eppcr_df['activity'] = pd.Series(
      utils.select_hits(eppcr_df[pval_col], utils.EXPECTED_FDR)
  ).apply(lambda label: 'functional' if label else 'non-functional')
  color_palette = {'non-functional': 'red', 'functional': 'green'}
  with sns.axes_style('ticks'):
    sns.histplot(
        data=eppcr_df,
        x=score_column,
        hue='activity',
        common_norm=False,
        multiple='dodge',
        stat='probability',
        palette=color_palette,
        bins=50,
        hue_order=['non-functional', 'functional'],
        ax=ax,
        legend=use_legend,
    )
    if title:
      ax.set_title(title)
    ax = ax if ax else plt.gca()
    ax.set_ylabel('proportion')
    ax.set_xlabel('zero-shot model score (WT score = 0)')
    _set_fontsize(ax, fontsize)
