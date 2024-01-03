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

"""Tests for data_utils."""

from absl.testing import absltest
from absl.testing import parameterized

from nuclease_design import data_utils


class DataUtilsTest(parameterized.TestCase):

  @parameterized.parameters(
      ('A45C, A56T', (('A', 45, 'C'), ('A', 56, 'T'))),
      ('A45C,A56T', (('A', 45, 'C'), ('A', 56, 'T'))),
      ('A45C,    A56T', (('A', 45, 'C'), ('A', 56, 'T'))),
      ('M99V, P88R, A132L', (('P', 88, 'R'), ('M', 99, 'V'), ('A', 132, 'L'))),
      ('*426A', (('*', 426, 'A'),)),
      ('A9*', (('A', 9, '*'),)),
      ('', ()))
  def test_parse_mutants(self, input_string, expected_output):
    actual_output = data_utils.parse_mutants(input_string)
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(
      ('nan'),
      ('AC'),
      ('a47k'),
      ('A1C,'),
      ('23AK'),
      ('lemonadelemonade'))
  def test_parse_mutants_raises(self, input_string):
    with self.assertRaisesRegex(ValueError, 'Malformed input'):
      data_utils.parse_mutants(input_string)

  @parameterized.parameters(
      ((('A', 45, 'C'), ('A', 56, '*')), True),
      ((('A', 45, 'C'), ('A', 56, 'T')), False),
      ((('*', 45, 'C'), ('A', 56, 'T')), False),
      ((('C', 1, '*'),), True),
      ((), False)
  )
  def test_introduces_stop_codon(self, mutations, expected_output):
    actual_output = data_utils.introduces_stop_codon(mutations)
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(
      ((('*', 45, 'C'),), True),
      ((('A', 45, 'C'), ('A', 56, 'T')), False),
      ((('*', 45, 'C'), ('A', 56, 'T')), True),
      ((('C', 1, '*'),), False),
      ((), False)
  )
  def test_mutates_stop_codon(self, mutations, expected_output):
    actual_output = data_utils.mutates_stop_codon(mutations)
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(((('A', 1, 'B'),), 'AA', 'BA'),
                            ((('B', 10, 'A'),), 'AAAAAAAAAB', 'AAAAAAAAAA'))
  def test_apply_mutations_to_parent(self, mutations, parent, expected_output):
    self.assertEqual(
        data_utils.apply_mutations_to_parent(mutations, parent),
        expected_output)

  @parameterized.parameters(((('A', 1, 'B'),), 'CA'),
                            ((('B', 10, 'A'),), 'AAAAAAAAAA'))
  def test_apply_bad_mutations_to_parent(self, mutations, parent):
    with self.assertRaisesRegex(ValueError, 'Invalid Mutation'):
      data_utils.apply_mutations_to_parent(mutations, parent)


if __name__ == '__main__':
  absltest.main()
