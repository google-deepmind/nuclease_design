# Copyright 2025 Google LLC
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

"""Utilities for selecting batches."""

import collections
import itertools
from typing import Any, cast, Sequence

import more_itertools
import numpy as np
import pandas as pd

from zero_shot_sampling import alignment


def _flatten(list_of_lists):
  return list(itertools.chain(*list_of_lists))


def _get_weighted_interleaved_iterator(
    iterators, weights, random_state, n_copies=1000
):
  """Randomly interleaves `iterators` with sampling weights `weights`."""
  # To approximate weighted random sampling from the iterators, we make an
  # expanded list of references to the iterators, where the number of copies
  # of iterators[i] is proportional to weights[i], and then shuffle their order.
  weights = np.asarray(weights)
  counts = n_copies * weights / weights.sum()

  expanded_iterators = _flatten(
      [[iterator] * int(count) for iterator, count in zip(iterators, counts)]
  )
  random_state.shuffle(expanded_iterators)
  return more_itertools.roundrobin(*expanded_iterators)


def _get_constrained_iterator(items, filter_fn, update_fn):
  for item in items:
    if filter_fn(item):
      update_fn(item)
      yield item


def select_top_from_groups_with_mutation_usage_constraint(
    df: pd.DataFrame,
    num_to_select: int,
    group_col: str,
    mutations_col: str,
    group_proportions: dict[Any, float],
    sort_column: str,
    sort_ascending: bool,
    random_state: np.random.RandomState,
    max_mutation_usage: int,
    mutation_allowlist: Sequence[tuple[int, str]] = (),
    mutation_denylist: Sequence[tuple[int, str]] = (),
    max_mutation_usage_for_denylist: int | None = None,
) -> pd.DataFrame:
  """Selects the top-scoring rows per group while capping mutation usage.

  This function selects a subset of rows from `df` with the highest total score
  that satisfy the following
  constraints:
  1. The number of times that any given mutation appears in the output is at
     most `max_mutation_usage`.
  2. The proportions of rows in each group, given by `group_col` induced by
  `group_proportions` are
     respected.
  3. `num_to_select` rows are selected in total.

  Note that it may be impossible to satisfy all three constraints, or our greedy
  algorithm may not be able to find a solution. In this case, fewer than
  `num_to_select` rows may be returned.

  Args:
    df: A `DataFrame` with columns `group_col`, `mutations_col`, and
      `sort_column`.
    num_to_select: The number of rows to select.
    group_col: The name of the column in `df` that groups the rows.
    mutations_col: The name of the column in `df` that contains the mutations,
      where each mutation is a list of `alignment.Edit`.
    group_proportions: A dictionary mapping the values of `group_col` to
      their respective proportions in the output. They do not need to sum to 1,
      as the normalization will be done inside this function.
    sort_column: The name of the column in `df` that should be used to sort the
      rows.
    sort_ascending: Whether to sort the rows in ascending order of
      `sort_column`.
    random_state: A `RandomState` for reproducibility.
    max_mutation_usage: The maximum number of times that any given mutation can
      appear in the output.
    mutation_allowlist: A sequence of `alignment.Edit` tuples that are not
      subject to the mutation usage constraint.
    mutation_denylist: A sequence of `alignment.Edit` tuples that are subject to
      a different mutation usage constraint given by
      `max_mutation_usage_for_denylist`.
    max_mutation_usage_for_denylist: See above.

  Returns:
    A `DataFrame` with the selected rows.
  """
  df = df.sample(frac=1.0, random_state=random_state).sort_values(
      sort_column, ascending=sort_ascending
  )

  df = cast(pd.DataFrame, df)
  mutation_usage_counter = collections.Counter()

  def filter_fn(row):

    usages_of_restricted_mutations = [
        mutation_usage_counter[mutation]
        for mutation in row[mutations_col]
        if mutation not in mutation_allowlist
    ]
    if not usages_of_restricted_mutations:
      return True
    if max(usages_of_restricted_mutations) >= max_mutation_usage:
      return False
    if mutation_denylist:
      usages_of_denylist_mutations = [
          mutation_usage_counter[mutation]
          for mutation in row[mutations_col]
          if mutation in mutation_denylist
      ]
      if (
          usages_of_denylist_mutations
          and max(usages_of_denylist_mutations)
          > max_mutation_usage_for_denylist
      ):
        return False
    return True

  def update_fn(row):
    mutation_usage_counter.update(row[mutations_col])

  iterators = []
  weights = []
  for key, group_df in df.groupby(group_col):
    if key in group_proportions:
      weights.append(group_proportions[key])
      iterators.append(
          _get_constrained_iterator(
              group_df.to_dict(orient='records'), filter_fn, update_fn
          )
      )

  if not iterators:
    raise ValueError(
        'No groups were found in the input dataframe that are present in'
        ' group_proportions. Please check your input dataframe and'
        ' group_proportions.'
    )
  iterator = _get_weighted_interleaved_iterator(
      iterators, weights, random_state
  )

  selected_rows = itertools.islice(iterator, num_to_select)

  to_return = pd.DataFrame(selected_rows)

  actual_mutation_counts = to_return[mutations_col].explode().value_counts()
  actual_mutation_counts = actual_mutation_counts[
      ~actual_mutation_counts.index.isin(mutation_allowlist)
  ]
  assert actual_mutation_counts.max() <= max_mutation_usage, (
      'The maximum number of mutations in the selected library is'
      f' {actual_mutation_counts.max()}, which is greater than the maximum'
      f' mutation usage allowed: {max_mutation_usage}.'
  )

  assert to_return[group_col].isin(list(group_proportions)).all()
  return to_return


def _sort_edits_by_position(edits):
  return tuple(sorted(edits, key=lambda x: (x[1], x[0], x[2])))


def _remove_edits(edits, substitutions_to_remove):
  return [edit for edit in edits if edit not in substitutions_to_remove]


def _edit_positions_are_unique(edits: tuple[alignment.Edit, ...]) -> bool:
  positions = [position for _, position, _ in edits]
  return len(set(positions)) == len(positions)


def expand_edits_by_removing_heavy_hitter_substitutions(
    edit_lists: Sequence[tuple[alignment.Edit, ...]],
    top_n: int = 10,
    max_num_substitutions_to_remove: int = 3,
) -> list[tuple[alignment.Edit, ...]]:
  """Expands a list of edits by removing heavy hitter substitutions.

  Args:
    edit_lists: A list of lists of alignment.Edit tuples.
    top_n: The `top_n` most frequent substitutions will be considered
      heavy-hitters.
    max_num_substitutions_to_remove: The maximum number of heavy hitter
      substitutions to remove from a given element of `edit_lists`.

  Returns:
    A list of lists of alignment.Edit tuples that contains the input
    `edit_lists` along with expanded versions of the edits where
    heavy-hitter substitutions have been removed using the above approach.
  """
  for edits in edit_lists:
    alignment.validate_edit_positions(edits)
  edit_lists = [_sort_edits_by_position(edits) for edits in edit_lists]
  all_edits = pd.Series(edit_lists).explode().dropna()
  most_frequent_substitutions = (
      all_edits[all_edits.apply(alignment.is_substitution)]
      .value_counts()
      .head(top_n)
      .index
  )

  edit_list_to_expand = set(edit_lists)
  to_return = list(edit_lists)
  for substitutions_to_remove in itertools.chain(*([
      itertools.combinations(most_frequent_substitutions, n_comb)
      for n_comb in range(1, max_num_substitutions_to_remove + 1)
  ])):
    to_return.extend([
        _remove_edits(edits, substitutions_to_remove)
        for edits in edit_list_to_expand
    ])

  to_return = [
      edits for edits in to_return if _edit_positions_are_unique(edits)
  ]
  to_return = [_sort_edits_by_position(edits) for edits in to_return if edits]
  return list(set(to_return))
