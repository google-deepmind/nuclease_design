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

"""Tests for metrics.py."""

from absl.testing import absltest
from absl.testing import parameterized
from nuclease_design import metrics
import numpy as np


class MetricsTest(parameterized.TestCase):

  @parameterized.named_parameters(
      dict(
          testcase_name='2_clusters',
          sequences=[[0, 0, 0, 0], [1, 2, 0, 0], [1, 2, 0, 0]],
          max_intra_cluster_hamming_distance=1,
          expected_num_clusters=2,
      ),
      dict(
          testcase_name='3_clusters',
          sequences=[[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0]],
          max_intra_cluster_hamming_distance=1,
          expected_num_clusters=3,
      ),
      dict(
          testcase_name='1_cluster',
          sequences=[[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0]],
          max_intra_cluster_hamming_distance=2,
          expected_num_clusters=1,
      ),
  )
  def test_num_clusters(
      self, sequences, max_intra_cluster_hamming_distance, expected_num_clusters
  ):
    sequences = np.vstack(sequences)
    pdist = metrics.pairwise_hamming_distance(sequences)
    num_clusters = metrics.num_clusters(
        pdist,
        max_intra_cluster_hamming_distance=max_intra_cluster_hamming_distance,
    )
    self.assertEqual(num_clusters, expected_num_clusters)

if __name__ == '__main__':
  absltest.main()
