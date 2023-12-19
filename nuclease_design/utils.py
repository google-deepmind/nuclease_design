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

"""Utilities for nuclease publication."""

import functools
import io
import json
from os import path
from typing import Sequence

from nuclease_design import amino_acids
from nuclease_design import constants
from nuclease_design import data_utils
from nuclease_design import metrics
import numpy as np
import pandas as pd
import requests
from statsmodels.stats import multitest


AA_STR_TO_INT = dict(((aa, i) for i, aa in enumerate(amino_acids.AA)))
DEFAULT_NUM_BOOTSTRAPS = 15

# False discovery rate
EXPECTED_FDR = 0.1


def _read_from_gcs(file_path: str):
  response = requests.get(file_path)
  return io.BytesIO(response.content)


def open_file(filename: str, mode: str, data_dir: str):
  """Layer of indirection to allow for using native python `open`."""
  full_filename = path.join(data_dir, filename)

  if data_dir == constants.GCS_DATA_DIR:
    if mode != 'r':
      return ValueError('mode must be "r" when reading from GCS.')
    return _read_from_gcs(full_filename)
  return open(full_filename, mode)


def concat_unique_values(items):
  return tuple(sorted(set(items)))


def _load_mutations_to_sublibrary_df(data_dir):
  with open_file(constants.MUTATIONS_TO_SUBLIBRARY_PATH, 'r', data_dir) as f:
    return decode_df(pd.read_csv(f))


def get_mutations_to_sublibraries(
    generation, data_dir: str = constants.DATA_DIR
) -> dict[str, Sequence[str]]:
  df = _load_mutations_to_sublibrary_df(data_dir)
  df = df[df['generation'] == generation].copy()
  grouped_df = df.groupby('mutations').agg({
      'sublibrary_name': concat_unique_values,
  })
  return grouped_df.to_dict()['sublibrary_name']


def encode_df(df: pd.DataFrame) -> pd.DataFrame:
  """Encodes `df` columns containing nested data structures."""

  df = df.copy()
  if 'mutations' in df.columns:
    df['mutations'] = df['mutations'].apply(json.dumps)
  if 'sublibrary_names' in df.columns:
    df['sublibrary_names'] = df['sublibrary_names'].apply(json.dumps)
  if 'generations' in df.columns:
    df['generations'] = df['generations'].apply(json.dumps)
  return df


def decode_df(df: pd.DataFrame) -> pd.DataFrame:
  """Decodes `df` columns containing nested data structures."""
  df = df.copy()
  if 'mutations' in df.columns:
    df['mutations'] = (
        df['mutations']
        .apply(json.loads)
        .apply(lambda muts: tuple(tuple(mut) for mut in muts))
    )
  if 'sublibrary_names' in df.columns:
    df['sublibrary_names'] = (
        df['sublibrary_names'].apply(json.loads).apply(tuple)
    )
  if 'generations' in df.columns:
    df['generations'] = df['generations'].apply(json.loads).apply(tuple)
  return df


def load_fiducial_data(
    fiducial: str, generation: str, data_dir: str = constants.DATA_DIR
) -> pd.DataFrame:
  with open_file(
      constants.get_fiducial_filename(fiducial, generation), 'r', data_dir
  ) as f:
    return decode_df(pd.read_csv(f))


def load_data(
    generation: str,
    include_sequence=True,
    parent=constants.FULL_REFERENCE_SEQ,
    data_dir: str = constants.DATA_DIR,
) -> pd.DataFrame:
  """Loads processed enrichment factor data for `generation`.

  Args:
    generation: Which generation to load data for.
    include_sequence: Whether to include full-length protein sequences in the
      output `DataFrame`. Otherwise, the sequence is represented implicitly
      by the 'mutations' column.
    parent: The parent sequence used for constructing the full-length sequences
      from 'mutations'.
    data_dir: Directory to load data from.

  Returns:
    A `DataFrame` containing the enrichment factor data for `generation`.
  """
  with open_file(
      constants.get_processed_data_file(generation), 'r', data_dir
  ) as f:
    df = decode_df(pd.read_csv(f))
  df['num_mutations'] = df['mutations'].apply(len)
  if include_sequence:
    df['sequence'] = df.mutations.apply(
        data_utils.apply_mutations_to_parent, parent=parent
    )
  return df


def load_all_data(data_dir: str = constants.DATA_DIR) -> pd.DataFrame:
  return pd.concat(
      [
          load_data(generation, data_dir=data_dir).assign(generation=generation)
          for generation in ['g1', 'g2', 'g3', 'g4']
      ],
      ignore_index=True,
  )


def load_landscape(data_dir: str = constants.DATA_DIR) -> pd.DataFrame:
  with open_file(constants.LANDSCAPE_PATH, 'r', data_dir) as f:
    return decode_df(pd.read_csv(f))


def get_activity_label(reference_name: str) -> str:
  if reference_name == 'stop':
    return 'neg_control'
  else:
    return reference_name


@functools.cache
def _encode_sequence(sequence: str) -> list[int]:
  ints = []
  for aa in sequence:
    ints.append(AA_STR_TO_INT[aa])
  return ints


def num_sequence_clusters(
    sequences: Sequence[str],
    max_intra_cluster_hamming_distance: int,
) -> int:
  """Returns the number of sequence clusters.

  See the docstring for metrics.num_clusters for
  more details on clustering. This method returns clusters built by "farthest
  point" inter-cluster distance, i.e. d(u, v) =  max_ij(d(u[i], v[j]).

  This is consistent with UniRef90 clustering. From
  https://www.uniprot.org/help/uniref:
    each cluster is composed of sequences that have at least 90% sequence
    identity to and 80% overlap with the longest sequence (a.k.a. seed sequence)
    of the cluster.

  Args:
    sequences: a sequence of strings.
    max_intra_cluster_hamming_distance: the distance cutoff for clustering. If
      two sequences are in the same cluster, they have distance <=
      `max_intra_cluster_hamming_distance`.
  """
  int_encoded_sequences = np.vstack([_encode_sequence(s) for s in sequences])
  pdist = metrics.pairwise_hamming_distance(int_encoded_sequences)
  return metrics.num_clusters(pdist, max_intra_cluster_hamming_distance)


def select_hits(
    independent_pvals: Sequence[float], expected_false_discovery_rate: float
) -> list[float]:
  return multitest.multipletests(
      pvals=independent_pvals,
      alpha=expected_false_discovery_rate,
      method='fdr_bh',
  )[0]


def select_hit_rows(
    df: pd.DataFrame, reference_name: str, expected_false_discovery_rate: float
) -> pd.DataFrame:
  assert df['generation'].nunique() == 1
  generation = df['generation'].iloc[0]
  pval_col_name = get_pvalue_column_name(generation, reference_name)
  is_hit = select_hits(df[pval_col_name], expected_false_discovery_rate)
  return df[is_hit].copy()


def get_hit_rate(df, pval_col, expected_false_discovery_rate):
  return np.mean(select_hits(df[pval_col], expected_false_discovery_rate))


def get_hit_rate_stats(df, pval_col, expected_false_discovery_rate):
  is_hit = select_hits(df[pval_col], expected_false_discovery_rate)
  return dict(
      hit_rate=np.mean(is_hit),
      num_hits=np.sum(is_hit),
      library_size=len(is_hit),
  )


def get_pvalue_column_name(generation: str, ref: str) -> str:
  gate_name = constants.GENERATION_AND_REF_NAME_TO_GATE_NAME[generation][ref]
  return f'pvalue_gen_{generation}_ref_{ref}_{gate_name}_right'


def get_enrichment_factor_column_name(generation: str, ref: str) -> str:
  gate_name = constants.GENERATION_AND_REF_NAME_TO_GATE_NAME[generation][ref]
  return f'ef_{gate_name}_{generation}'


def get_pvalue(
    row,
    reference_name: str,
):
  pval_col_name = get_pvalue_column_name(row['generation'], reference_name)
  return row[pval_col_name]


def _bootstrap(df, random_state):
  return df.sample(frac=1, random_state=random_state, replace=True)


def get_bootstrapped_hitrate_df(
    df,
    pval_col,
    expected_false_discovery_rate,
    random_state,
    num_bootstraps=DEFAULT_NUM_BOOTSTRAPS,
):
  """Returns a long df from `num_bootstraps` calls to `get_hit_rate_stats`."""
  rows = []
  for _ in range(num_bootstraps):
    rows.append(
        get_hit_rate_stats(
            _bootstrap(df, random_state),
            pval_col=pval_col,
            expected_false_discovery_rate=expected_false_discovery_rate,
        )
    )
  return pd.DataFrame(rows)


def get_bootstrapped_hit_rate_stats(
    df: pd.DataFrame,
    reference_name: str,
    random_seed: int,
    expected_false_discovery_rate: float = EXPECTED_FDR,
) -> pd.DataFrame:
  """Returns bootstrapped hit rate mean/std for each sublibrary_name.

  Args:
    df: DataFrame with pvalue columns.
    reference_name: a string indicating the reference sequence for the "hit"
      comparison. (e.g. 'neg_control', 'wt', 'a73r')
    random_seed: controls randomness.
    expected_false_discovery_rate: the FDR rate at which to control.
  """
  df = df.copy()
  random_state = np.random.RandomState(random_seed)
  df['pvalue'] = df.apply(
      get_pvalue,
      reference_name=reference_name,
      axis=1,
  )

  hit_rate_df = df.groupby('sublibrary_name').apply(
      get_bootstrapped_hitrate_df,
      pval_col='pvalue',
      expected_false_discovery_rate=expected_false_discovery_rate,
      random_state=random_state,
  )
  agg_dict = {
      key: ('mean', 'std') for key in ['num_hits', 'hit_rate', 'library_size']
  }
  return hit_rate_df.groupby('sublibrary_name').agg(agg_dict)


def filter_to_g4_positions(df: pd.DataFrame) -> pd.DataFrame:
  """Selects rows from `df` where all mutations are in the G4 design region."""
  return df[
      df['mutations'].apply(
          lambda muts: all(
              pos >= constants.G4_VARIABLE_REGION_START
              and pos < constants.G4_VARIABLE_REGION_END
              for _, pos, _ in muts
          )
      )
  ]


def _get_sublibrary_df(df: pd.DataFrame, sublibrary_name: str) -> pd.DataFrame:
  includes_sublibrary = df['sublibrary_names'].apply(
      lambda t: sublibrary_name in t)
  return (df[includes_sublibrary]
          .assign(sublibrary_name=sublibrary_name)
          .drop(columns='sublibrary_names'))


def expand_sublibraries(df: pd.DataFrame) -> pd.DataFrame:
  """Expands each `df` row for each value in 'sublibrary_names'.

  Example:
  >>> df = pd.DataFrame([dict(sequence=ABC, sublibrary_names=('lib1', 'lib2')])
  >>> expand_sublibraries(df).to_csv()
  sequence,sublibrary_name
  ABC,lib1
  ABC,lib2

  Args:
    df: A `DataFrame` containing a 'sublibrary_names' column.

  Returns:
    A `DataFrame` containing a 'sublibrary_name' column, where
    each row of `df` is replicated for each value in `row.sublibrary_names`.
  """
  sublibrary_names = df['sublibrary_names'].explode().unique()
  return pd.concat(
      [
          _get_sublibrary_df(df, sublibrary_name)
          for sublibrary_name in sublibrary_names
      ],
      ignore_index=True,
  )
