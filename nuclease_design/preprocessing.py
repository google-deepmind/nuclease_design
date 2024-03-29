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

"""Library for data preprocessing."""

import functools
from typing import Sequence, Tuple

import pandas as pd
from statsmodels.stats import weightstats

from nuclease_design import constants
from nuclease_design import data_utils
from nuclease_design import preprocessing_utils
from nuclease_design import utils


def get_pvalue(
    enrichment_factor: float, null_distribution_efs: Sequence[float]
) -> float:
  """Returns the p-value from a pooled t-test."""
  return weightstats.ttest_ind(
      [
          enrichment_factor,
      ],
      null_distribution_efs,
      alternative='larger',
      # We assume the variance of alternative distribution is the same
      # as the variance of the null. We only have one sample from the
      # alternative distribution so this is forced.
      usevar='pooled',
  )[1]


def _drop_introduces_stop_codons(df: pd.DataFrame) -> pd.DataFrame:
  introduces_stop_codon = df['mutations'].apply(
      data_utils.introduces_stop_codon
  )
  print(
      '{} variants introduced stop codons'.format(introduces_stop_codon.sum())
  )
  return df[~introduces_stop_codon]


def _drop_mutates_stop_codon(df: pd.DataFrame) -> pd.DataFrame:
  mutates_stop_codon = df['mutations'].apply(data_utils.mutates_stop_codon)
  print(
      '{} variants mutated the WT stop codon'.format(mutates_stop_codon.sum())
  )
  return df[~mutates_stop_codon]


def _get_fiducial_df_fn_for_fiducial(fiducial_name):
  if fiducial_name == 'neg_control':
    return preprocessing_utils.get_stop_codon_df
  else:
    fiducial_mutation_tuple = constants.FIDUCIAL_NAME_TO_MUTATION_TUPLE[
        fiducial_name
    ]
    return functools.partial(
        preprocessing_utils.get_synonym_df, mutations=fiducial_mutation_tuple
    )


def preprocess_generation(
    names_and_paths: Sequence[tuple[str, str]],
    input_names: Sequence[str],
    post_sort_names: Sequence[str],
    count_threshold: int,
    fiducials: Sequence[str],
    generation: str,
    data_dir: str,
) -> Tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
  """Returns a df of EFs and a dict mapping fiducials to EFs.

  Args:
    names_and_paths: A list of (name, path), where path points to a CSV file to
      load, and name is the desired name of the column in the output
      corresponding to counts from that experiment.
    input_names: Columns to use for aggregating input abundances.
    post_sort_names: Columns to use for computing post-sort enrichment factors.
    count_threshold: A nucleic acid variant must have a count >= count threshold
      in at least one input_count_col in order to be kept.
    fiducials: A sequence of fiducials to use for comparison. Must be a subset
      of ['neg_control', 'wt', 'a73r', 'double', 'quad']
    generation: The string indicating the generation.
    data_dir: The directory to load from.

  Returns: A tuple whose first element is a dataframe mapping from variants to
    enrichment factors and fiducial-comparison p-values, and whose second
      element
    is a dictionary {fiducial_name: dataframe} where the dataframes including
      DNA-level enrichment factors
    for DNA replicates of the corresponding fiducial. For 'neg_control', the
      dataframe includes
    DNA variants that introduced premature stop codons. For all other fiducials,
      the dataframe
    includes synonymous DNA variants of the fiducial variant.
  """
  count_df = preprocessing_utils.load_ngs_counts_from_paths(
      names_and_paths, data_dir
  )
  ef_df = preprocessing_utils.get_enrichment_factor_df(
      count_df, input_names, post_sort_names, count_threshold=count_threshold
  )
  # remove stop codon mutations from variant DataFrame
  ef_df = _drop_introduces_stop_codons(ef_df)
  ef_df = _drop_mutates_stop_codon(ef_df)

  # get fiducial enrichment factors
  fiducial_ef_dfs = {}
  for fiducial_name in fiducials:
    get_fiducial_df_fn = _get_fiducial_df_fn_for_fiducial(fiducial_name)
    fiducial_ef_dfs[fiducial_name] = get_fiducial_df_fn(
        outer_join_df=count_df,
        input_count_cols=input_names,
        post_sort_count_cols=post_sort_names,
        count_threshold=count_threshold,
    )

  # compute pvalues
  for fiducial_name in fiducials:
    ef_col = utils.get_enrichment_factor_column_name(generation, fiducial_name)
    pvalue_col = utils.get_pvalue_column_name(generation, fiducial_name)
    pvalues = ef_df[ef_col].apply(
        get_pvalue, null_distribution_efs=fiducial_ef_dfs[fiducial_name][ef_col]
    )
    ef_df[pvalue_col] = pvalues
    is_hit_col = f'activity_greater_than_{fiducial_name}'
    ef_df[is_hit_col] = utils.select_hits(pvalues, utils.EXPECTED_FDR)
  pvalue_df = ef_df.copy()
  pvalue_df['num_mutations'] = pvalue_df['mutations'].apply(len)

  # add sublibrary names
  mutations_to_sublibraries = utils.get_mutations_to_sublibraries(generation)
  pvalue_df['sublibrary_names'] = pvalue_df['mutations'].apply(
      mutations_to_sublibraries.get
  )
  return pvalue_df, fiducial_ef_dfs
