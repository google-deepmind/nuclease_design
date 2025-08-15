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

import itertools

from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import pandas as pd

from zero_shot_sampling import alignment
from zero_shot_sampling import batch_selection_utils


class BatchSelectionUtilsTest(parameterized.TestCase):

  def assertProportionsClose(self, actual_values, expected_proportions):
    expected_proportions = pd.Series(expected_proportions)
    expected_normalized_proportions = (
        expected_proportions / expected_proportions.sum()
    )

    actual_normalized_proportions = actual_values.value_counts(normalize=True)
    self.assertSetEqual(
        set(actual_normalized_proportions.index),
        set(expected_normalized_proportions.index),
    )
    pd.testing.assert_series_equal(
        actual_normalized_proportions,
        expected_normalized_proportions,
        check_names=False,
        check_like=True,
        atol=0.05,
    )

  @parameterized.parameters(
      dict(
          group_proportions={'a': 1, 'b': 1, 'c': 1},
          max_mutation_usage=1000,
          expected_score_proportions={2.0: 1},
      ),
      dict(
          group_proportions={'a': 7, 'b': 4},
          max_mutation_usage=1000,
          expected_score_proportions={2.0: 1},
      ),
      dict(
          group_proportions={'a': 1, 'b': 1},
          max_mutation_usage=250,
          expected_score_proportions={2.0: 1, 1.0: 1},
      ),
      dict(
          group_proportions={'a': 1, 'b': 1},
          max_mutation_usage=10,
          expected_score_proportions={2.0: 1, 1.0: 1},
          expected_output_size=40,
      ),
  )
  def test_select_top_from_groups_with_mutation_usage_constraint_group_proportions(
      self,
      group_proportions,
      max_mutation_usage,
      expected_score_proportions,
      expected_output_size=1000,
  ):

    df = pd.DataFrame({
        'sequence': ['s1', 's2', 's3', 's4', 's5', 's6'],
        'group': ['a', 'a', 'b', 'b', 'c', 'c'],
        'mutations': [
            (('X', 0, 'Y'),),
            (('Z', 1, 'W'),),
            (('X', 0, 'S'),),
            (('Z', 1, 'T'),),
            (('X', 0, 'R'),),
            (('Z', 1, 'Q'),),
        ],
        'score': [1.0, 2.0, 2.0, 1.0, 1.0, 2.0],
    })

    df = pd.concat([df] * 10000).sample(
        frac=1, random_state=np.random.RandomState(0)
    )

    actual_output = batch_selection_utils.select_top_from_groups_with_mutation_usage_constraint(
        df,
        num_to_select=1000,
        group_col='group',
        mutations_col='mutations',
        group_proportions=group_proportions,
        sort_column='score',
        sort_ascending=False,
        random_state=np.random.RandomState(0),
        max_mutation_usage=max_mutation_usage,
    )

    self.assertLen(actual_output, expected_output_size)
    self.assertProportionsClose(actual_output['group'], group_proportions)
    self.assertProportionsClose(
        actual_output['score'], expected_score_proportions
    )

  def test_expand_edits_by_removing_heavy_hitter_substitutions(self):
    heavy_hitter_substitutions = [('A', 0, 'B'), ('C', 1, 'D')]
    other_edits = [('A', 2, '-'), ('C', 3, 'E')]

    all_edits = [tuple([edit]) for edit in heavy_hitter_substitutions] * 5
    all_edits.extend([
        tuple([other_edit, heavy_hitter_edit])
        for other_edit, heavy_hitter_edit in itertools.product(
            other_edits, heavy_hitter_substitutions
        )
    ])

    expected_output = all_edits + [tuple([edit]) for edit in other_edits]
    expected_output = [
        batch_selection_utils._sort_edits_by_position(edits)
        for edits in expected_output
    ]

    actual_output = batch_selection_utils.expand_edits_by_removing_heavy_hitter_substitutions(
        all_edits,
        top_n=2,
        max_num_substitutions_to_remove=1,
    )

    for edits in actual_output:
      alignment.validate_edit_positions(edits)

    self.assertSetEqual(set(actual_output), set(expected_output))


if __name__ == '__main__':
  absltest.main()
