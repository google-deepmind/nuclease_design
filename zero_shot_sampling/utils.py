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

"""Utilities for working with data processing and modeling."""

import io
import itertools
import os

from Bio import SeqIO
import numpy as np
import pandas as pd
import requests



_GCS_PREFIX = 'https://storage.googleapis.com/'


def _read_from_gcs(file_path: str):
  response = requests.get(file_path)
  return io.BytesIO(response.content)


def open_file(filename: str, mode: str) -> ...:
  """Layer of indirection for loading data from multiple locations."""

  if filename.startswith(_GCS_PREFIX):
    if mode != 'r':
      raise ValueError('mode must be "r" when reading from GCS.')
    return _read_from_gcs(filename)
  return open(filename, mode)  # pylint: disable=unreachable


def make_directory(directory: str) -> None:
  """Creates a directory if it does not already exist."""
  os.makedirs(directory, exist_ok=True)  # pylint: disable=unreachable


def read_fasta_as_df(
    filename: str, num_records: int | None = None
) -> pd.DataFrame:
  """Reads records from a FASTA file.

  Args:
    filename: File name of the input FASTA file.
    num_records: The number of records to read. If None, all records are read.

  Returns:
    A `pd.DataFrame` with columns 'id', 'description', and 'sequence'.
    Given a fasta record
        >id...
        sequence
    each row is a tuple `(id, description, sequence)`, where `description`
    represents all characters after '>', i.e. 'id...'
  """
  lines = []
  with open(filename, 'r') as f:
    for record in itertools.islice(SeqIO.parse(f, 'fasta'), num_records):
      lines.append(
          dict(
              id=record.id,
              description=record.description,
              sequence=str(record.seq),
          )
      )

  return pd.DataFrame(lines)


def _batchify(inputs, batch_size):
  """Reshapes and pads inputs to include an additional batch dimension.

  The inputs can be of arbitrary length. They length does not need to be a
  multiple of batch_size, in which case padding will be added.

  Args:
    inputs: An np.ndarray or iterable of np.ndarray of shape [input_size, ...].
    batch_size:
      The size of the batches to group the data into.

  Returns:
    batched_inputs: np.ndarray of size [num_batches, batch_size, ...],
    where num_batches is ceil(input_size / batch_size).
    pad_size: Number of examples in the final batch that are padding. We use
      copies of inputs[0] as padding.
  """

  inputs = np.asarray(inputs)

  pad_size = -len(inputs) % batch_size
  if pad_size:
    padding = np.tile(inputs[:1], [pad_size] + (inputs.ndim - 1) * [1])
    padded_inputs = np.concatenate([inputs, padding], axis=0)
  else:
    padded_inputs = inputs
  batched_shape = (-1, batch_size) + padded_inputs.shape[1:]
  batched_inputs = np.reshape(padded_inputs, batched_shape)
  return batched_inputs, pad_size


def batch_apply(fn, inputs, batch_size):
  """Applies fn() to inputs in batches of size batch_size.

  The inputs can be of arbitrary length. They length does not need to be a
  multiple of batch_size. Padding will be added (and then removed) such that
  fn() is always called on inputs of size exactly batch_size. fn() is assumed
  to operate independently across the batch dimension of its inputs (e.g.
  computing predictions of a model on inputs) instead of performing an operation
  where the batch elements interact (e.g. performing a gradient step of a model
  on a batch of inputs).

  Args:
    fn: The function to map across the inputs.
    inputs: An np.ndarray or iterable of np.ndarray. fn() is mapped
      along the first dimension.
    batch_size:
      The size of the batches to evaluate fn() on.

  Returns:
    np.ndarray where outputs[i] = fn(inputs[i])
  """

  batched_inputs, pad_size = _batchify(inputs, batch_size)
  results = np.concatenate([fn(batch) for batch in batched_inputs])
  if pad_size:
    results = results[:-pad_size]
  return results
