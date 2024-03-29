# Copyright 2024 DeepMind Technologies Limited
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

"""Data loading functions."""

import multiprocessing.pool
from typing import Sequence, Tuple

import numpy as np
import pandas as pd

from nuclease_design import data_utils
from nuclease_design import utils

MutationTuple = data_utils.MutationTuple


def load_count_df(filepath: str, data_dir: str) -> pd.DataFrame:
  """Loads count data from an unprocessed CSV file.

  Args:
    filepath: The path to the input CSV to load. The input CSV should have the
    following columns: 'nuc_mut_string' (str), 'mut_string' (str), 'num_AA_mut'
      (int), 'counts' (int)
    data_dir: The directory to load from.

  Returns:
    A DataFrame with the columns: 'nuc_mutations' (MutationTuple),
    'mutations' (MutationTuple), and 'counts' (int).

  Raises:
    ValueError: If any row in the input CSV has a discrepancy between the number
      of mutations in the mutation string and the number in the
      num_AA_mut column.

    ValueError: If any expected column is not in the input.
  """
  with utils.open_file(filepath, 'r', data_dir) as f:
    df = pd.read_csv(f)

  expected_input_cols = ['nuc_mut_string', 'mut_string', 'num_AA_mut', 'counts']
  if not set(expected_input_cols).issubset(df.keys()):
    raise ValueError('Expected input columns missing in {}'.format(df.keys()))

  df['mut_string'] = df['mut_string'].fillna(value='')
  df['nuc_mut_string'] = df['nuc_mut_string'].fillna(value='')
  df['nuc_mutations'] = df.nuc_mut_string.apply(data_utils.parse_mutants)
  df['mutations'] = df.mut_string.apply(data_utils.parse_mutants)

  if (df.num_AA_mut != df['mutations'].apply(len)).any():
    raise ValueError(
        'Parsed mutation string does have the expected number of mutations.')

  output_cols = ['nuc_mutations', 'mutations', 'counts']
  return df[output_cols]


def group_by_dna_seq(df: pd.DataFrame) -> pd.DataFrame:
  """Aggregates identical nucleic acid sequences.

  Args:
    df: A DataFrame with the columns: "nuc_mutations" (MutationTuple),
      "mutations" (MutationTuple), and "counts" (int).

  Returns:
    A DataFrame with the total counts of each dna variant. The returned
    DataFrame has a "mutations" (MutationTuple) and "counts" (int) column.
  """
  agg_df = df.groupby('nuc_mutations').agg({
      'counts': 'sum',
      'mutations': 'first'
  })
  return agg_df.reset_index()


def _outer_join_dfs(list_of_dfs: Sequence[pd.DataFrame]) -> pd.DataFrame:
  """Returns a DataFrame containing all variants from `list_of_dfs`."""
  df = pd.DataFrame(columns=['nuc_mutations', 'mutations'])

  for df_to_merge in list_of_dfs:
    df = df.merge(df_to_merge, on=('nuc_mutations', 'mutations'), how='outer')

  return df


def load_ngs_counts_from_paths(
    names_and_paths: Sequence[Tuple[str, str]],
    data_dir: str,
) -> pd.DataFrame:
  """Loads unprocessed data from multiple experiments into a single DataFrame.

  Args:
    names_and_paths: A list of (name, path), where `path` points to a CSV file
      to load, and `name` is the desired name of the column in the output
      corresponding to counts from that experiment.
    data_dir: Directory to load data from.
  Returns:
    A DataFrame with one row per nucleic acid variant, and one column of counts
    per specified experiment.
  """
  mut_cols = ['nuc_mutations', 'mutations']

  def load(name_and_path):
    name, path = name_and_path
    df = load_count_df(path, data_dir)
    df = group_by_dna_seq(df)
    df = df.rename(columns=dict(counts=name))
    # Drop any unused columns.
    cols = mut_cols + [name]
    return df[cols]

  with multiprocessing.pool.ThreadPool(len(names_and_paths)) as pool:
    dfs = list(pool.map(load, names_and_paths))
    pool.close()
    pool.join()

  df = _outer_join_dfs(dfs)
  names, _ = zip(*names_and_paths)
  count_cols = list(names)

  # If a variant didn't appear in a column, it is considered to have 0 counts.
  df[count_cols] = df[count_cols].fillna(value=0)
  return df


def get_abundance(counts: Sequence[float]) -> np.ndarray:
  """Returns the abundance given a sequence of counts."""
  counts = np.array(counts)
  total_count = counts.sum()
  return counts / total_count


def get_enrichment_factor(
    input_abundance: Sequence[float], post_sort_abundance: Sequence[float]
) -> np.ndarray:
  """Returns the enrichment factor."""
  return np.array(post_sort_abundance) / np.array(input_abundance)


def _get_ef_column_name(count_col_name: str) -> str:
  return 'ef_' + count_col_name.removeprefix('read_count_')


def get_enrichment_factor_df(
    outer_join_df: pd.DataFrame,
    input_count_cols: Sequence[str],
    post_sort_count_cols: Sequence[str],
    count_threshold: int = 10,
    group_aa_seqs: bool = True,
) -> pd.DataFrame:
  """Computes enrichment factors.

  In an experiment (or "sort"), variants are sorted according to their ability
  to reach some activity threshold. The enrichment factor for a variant in a
  post-sort column is computed by taking the ratio of
  (abundance in the post-sort column) : (abundance in the input columns).
  The "abundance" of a variant is the proportion of all counts ascribed to the
  variant. Before computing abundances, drop all nucleotide variants that have
  counts < `count_threshold`.

  Args:
    outer_join_df: A DataFrame with counts for each experiment being analyzed,
      as well as a 'mutations' and 'nuc_mutations' column.
    input_count_cols: Columns to use for aggregating input abundances.
    post_sort_count_cols: Columns to use for computing post-sort enrichment
      factors.
    count_threshold: A nucleic acid variant must have a count >= count threshold
      in at least one input_count_col in order to be kept.
    group_aa_seqs: If True, report enrichment factors for amino acid sequences,
      else, for nucleic acid sequences.

  Returns:
     A DataFrame with one row per variant, and columns corresponding to
     enrichment factors computed on the post-sort columns. Enrichment factor
     columns have "ef_" prepended to the post-sort column name.
  """
  df = outer_join_df.copy()

  max_input_read_count = df[input_count_cols].max(axis=1)
  df = df[max_input_read_count >= count_threshold]

  # consider AA or DNA variants
  if group_aa_seqs:
    df = df.drop(columns=['nuc_mutations'])
    df = df.groupby('mutations').agg('sum').reset_index()
    df = df.sort_values(by='mutations', key=lambda col: [len(s) for s in col])
  else:
    if df['nuc_mutations'].nunique() != df.shape[0]:
      raise ValueError('Expected nuc_mutations column to be unique')
    df = df.sort_values(
        by='nuc_mutations', key=lambda col: [len(s) for s in col]
    )

  # sum all input counts
  input_counts = df[input_count_cols].sum(axis=1)
  input_abundance = get_abundance(input_counts)
  if input_abundance.min() <= 0:
    raise ValueError(
        f'Invalid minimum input abundance: {input_abundance.min()}'
    )

  # compute enrichment factors
  for post_sort_count_col in post_sort_count_cols:
    ef_col_name = _get_ef_column_name(post_sort_count_col)
    df[ef_col_name] = get_enrichment_factor(
        input_abundance, get_abundance(df[post_sort_count_col])
    )
  return df.reset_index(drop=True)


def get_synonym_df(outer_join_df: pd.DataFrame,
                   mutations: Sequence[MutationTuple],
                   input_count_cols: Sequence[str],
                   post_sort_count_cols: Sequence[str],
                   count_threshold: int = 10) -> pd.DataFrame:
  """Gets enrichment factors for synonymous variants."""
  df = get_enrichment_factor_df(
      outer_join_df,
      input_count_cols,
      post_sort_count_cols,
      count_threshold=count_threshold,
      group_aa_seqs=False)
  return df[df['mutations'] == mutations]


def get_stop_codon_df(
    outer_join_df: pd.DataFrame,
    input_count_cols: Sequence[str],
    post_sort_count_cols: Sequence[str],
    count_threshold: int = 10,
) -> pd.DataFrame:
  """Gets enrichment factors for variants to stop codons."""
  df = get_enrichment_factor_df(
      outer_join_df,
      input_count_cols,
      post_sort_count_cols,
      count_threshold=count_threshold,
      group_aa_seqs=True,
  )
  introduces_stop_codon = df['mutations'].apply(
      data_utils.introduces_stop_codon
  )
  return df[introduces_stop_codon].copy()
