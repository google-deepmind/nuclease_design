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

"""Metrics for evaluating design performance."""

import numpy as np
from scipy.cluster import hierarchy as sp_hierarchy
from scipy.spatial import distance as spd


def num_clusters(pdist: np.ndarray, max_intra_cluster_hamming_distance: int):
  """Returns the number of clusters given a pairwise distance matrix.

  This function uses hierarchical clustering based on complete linkage. First,
    we build a tree (dendrogram) between sequences by connecting two nodes
    (initially, one node per sequence) if the distance separating them is the
    smallest of all pairwise distances. After being connected, two nodes form
    a new node whose distance to any other node `n` is the max over all
    distances of subnodes to `n`. Nodes are joined iteratively until all
    sequences are connected.
    See the docstring of scipy.cluster.hierarchy.linkage for more details.
  Clusters are formed by the largest nodes that contain sequences whose
    maximum pairwise distance is no larger than `distance_threshold`. See
    the docstring of scipy.cluster.hierarchy.fcluster for more details.

  Args:
    pdist: np.array of shape [n_sequences, n_sequences] containing pairwise
      distances between sequences.
    max_intra_cluster_hamming_distance: The maximum hamming distance of two
      sequences in a cluster (inclusive).

  Returns:
    The number of clusters.
  """
  # Get the upper-triangular part of pdist.
  pdist = spd.squareform(pdist, checks=False)

  lk = sp_hierarchy.linkage(pdist, method='complete')
  clustering = sp_hierarchy.fcluster(
      lk, t=max_intra_cluster_hamming_distance, criterion='distance'
  )
  return np.max(clustering)


def pairwise_hamming_distance(array: np.ndarray) -> np.ndarray:
  """Returns a pairwise Hamming distance matrix for the input array.

  Args:
    array: NxL integer encoded array.

  Returns:
    NxN square matrix of pairwise Hamming distances.
  """
  l = array.shape[1]
  return l * spd.squareform(spd.pdist(array, metric='hamming'))
