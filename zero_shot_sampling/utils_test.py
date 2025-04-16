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

from unittest import mock
from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import pandas as pd
from zero_shot_sampling import utils


class UtilsTest(parameterized.TestCase):

  @parameterized.parameters((1, 1), (2, 2), (None, 4), (5, 4))
  def test_read_fasta_as_df(self, num_records, expected_output_length):
    ids = ['u0', 'u1', 'u2', 'u3']
    descriptions = ['u0 d0', 'u1 d1', 'u2 d2', 'u3 d3']
    sequences = ['AB', 'Cd', 'E.', '-G']
    fasta_contents = '>u0 d0\nAB\n>u1 d1\nCd\n>u2 d2\nE.\n>u3 d3\n-G\n'
    filename = self.create_tempfile(content=fasta_contents).full_path
    expected_df = pd.DataFrame(
        dict(
            id=ids[:expected_output_length],
            description=descriptions[:expected_output_length],
            sequence=sequences[:expected_output_length]))

    actual_df = utils.read_fasta_as_df(filename, num_records=num_records)

    pd.testing.assert_frame_equal(actual_df, expected_df)

  @parameterized.parameters((5, 5), (5, 1), (5, 2), (5, 6), (5, 12),
                            (5, 6, 3), (5, 12, 5))
  def test_batch_apply(self, batch_size, num_inputs, input_dim=2):
    def fn(inputs):
      return np.power(inputs + 1, 2)

    inputs = np.random.randint(low=0, high=8, size=(num_inputs, input_dim))
    unbatched_output = fn(inputs)
    batch_fn = mock.Mock(wraps=fn)
    batched_output = utils.batch_apply(batch_fn, inputs, batch_size)

    np.testing.assert_array_equal(unbatched_output, batched_output)

    self.assertEqual(batch_fn.call_count, np.ceil(num_inputs / batch_size))
    for call in batch_fn.call_args_list:
      self.assertLen(call[0][0], batch_size)

  @parameterized.parameters((5, (10, 2, 1)), (3, (11, 3, 4)))
  def test_batch_apply_multidim(self, batch_size, data_shape):
    input_data = np.ones(data_shape)

    def fn(inputs):
      return np.sum(inputs.shape[1:]) * inputs

    unbatched_output = np.array(fn(input_data))
    batched_output = utils.batch_apply(fn, input_data, batch_size)
    np.testing.assert_array_equal(unbatched_output, batched_output)


if __name__ == '__main__':
  absltest.main()
