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

"""Tools for sequence alignment."""

import dataclasses
import json
import re
from typing import Callable, Mapping, Sequence

import numpy as np
from numpy import typing as npt
import pandas as pd

from zero_shot_sampling import amino_acids
from zero_shot_sampling import domains
from zero_shot_sampling import utils

# An MSA contains amino acid characters and gap characters, which are used to
# denote that a particular sequence (row) has no token in a given column
# (sequence position) of the alignment. Sometimes, such as when training
# profile HMMs, columns are treated separately as being either 'match' or
# 'insert' columns. Depending on the type of column, gaps are denoted with
# different characters.

Edit = tuple[str, int, str]

# Deletions that occur in match states.
_MATCH_GAP_CHARACTER = amino_acids.GAP
# Deletions that occur in insert states.
_INSERT_GAP_CHARACTER = '.'

_GAP_CHARACTERS = (_MATCH_GAP_CHARACTER, _INSERT_GAP_CHARACTER)
_GAP_CHARACTER_REGEX = re.compile('[%s]' % ''.join(_GAP_CHARACTERS))


_INSERTION = 'insertion'
_DELETION = 'deletion'
_SUBSTITUTION = 'substitution'
_EDIT_TYPES = (_INSERTION, _DELETION, _SUBSTITUTION)


def _get_edit_type(new_value: str, parent_value: str) -> str:
  if new_value == parent_value:
    raise ValueError('New value and parent value must be different.')
  if is_gap(new_value) and not is_gap(parent_value):
    return _DELETION
  elif not is_gap(new_value) and is_gap(parent_value):
    return _INSERTION
  elif not is_gap(new_value) and not is_gap(parent_value):
    return _SUBSTITUTION
  else:
    raise ValueError(
        f'New value and parent value are both gap characters: {new_value} and'
        f' {parent_value}.'
    )


def get_edit_counts(
    aligned_sequence: str, aligned_parent_sequence: str
) -> Mapping[str, int]:
  """Returns a mapping from edit type to number of edits of that type."""
  if len(aligned_sequence) != len(aligned_parent_sequence):
    raise ValueError(
        'Aligned sequence and aligned parent sequence must have the same'
        ' length.'
    )
  to_return = {}
  edit_types = [
      _get_edit_type(new_val, parent_val)
      for new_val, parent_val in zip(aligned_sequence, aligned_parent_sequence)
      if new_val != parent_val
  ]

  for edit_type in _EDIT_TYPES:
    to_return[edit_type] = sum(
        observed_edit_type == edit_type for observed_edit_type in edit_types
    )
  return to_return


def is_substitution(edit: Edit) -> bool:
  value, _, parent_value = edit
  return not is_gap(value) and not is_gap(parent_value)


def get_edits(seq: str, parent_seq: str) -> tuple[Edit, ...]:
  return tuple([
      (parent_val, pos, val)
      for pos, (val, parent_val) in enumerate(zip(seq, parent_seq))
      if val != parent_val
  ])


def validate_edit_positions(edits: tuple[Edit, ...]):
  positions = [position for _, position, _ in edits]
  if len(set(positions)) != len(positions):
    raise ValueError(
        'Edits must have unique positions. Found duplicate positions:'
        f' {set(positions)} for {edits}.'
    )


def apply_edits(edits: tuple[Edit, ...], parent_seq: str) -> str:
  validate_edit_positions(edits)
  parent_seq = list(parent_seq)
  for edit in edits:
    old_value, position, new_value = edit
    if parent_seq[position] != old_value:
      raise ValueError(
          f'{edit} is not a valid edit for {parent_seq}, which has value'
          f' {parent_seq[position]} at position {position}.'
      )
    parent_seq[position] = new_value
  return ''.join(parent_seq)


def get_unaligned_seq(aligned_seq: str) -> str:
  """Removes gap characters from sequence and converts to upper case."""
  return re.sub(_GAP_CHARACTER_REGEX, '', aligned_seq).upper()


def is_gap(token: str) -> bool:
  return token in _GAP_CHARACTERS


@dataclasses.dataclass(frozen=True)
class MSA:
  parent: str
  homologs: tuple[str, ...]


def load_msa_from_fasta(fasta_file: str) -> MSA:
  """Loads an MSA from a fasta file where the first entry is the parent."""
  msa_df = utils.read_fasta_as_df(fasta_file)
  parent = msa_df.sequence.iloc[0]
  homologs = tuple(msa_df.sequence.iloc[1:])
  return MSA(parent, homologs)


def _convert_match_states(aligned_seq, match_positions):
  chars = []
  for i, char in enumerate(aligned_seq):
    if i in match_positions:
      char = '-' if is_gap(char) else char.upper()
    else:
      char = '.' if is_gap(char) else char.lower()
    chars.append(char)
  return ''.join(chars)


def convert_parent_residues_to_match_states(msa: MSA) -> MSA:
  """Sets all non-gap characters in the parent to match states.

  Each column of an MSA can either be a 'match state' or an 'insert state'.
  Match states are denoted by upper-case letters and '-' and insert states are
  denoted by lower-case letters and '.'. This function transforms the MSA such
  that every column corresponding to anon-gap character in the parent sequence
  is a match state and all other columns are insert states.

  Args:
    msa: An MSA.

  Returns:
    An MSA with the parent and homologs transformed in the above manner.
  """
  non_gap_parent_positions = [
      i for i, char in enumerate(msa.parent) if not is_gap(char)
  ]
  match_positions = non_gap_parent_positions
  parent = _convert_match_states(msa.parent, match_positions)
  homologs = [
      _convert_match_states(seq, match_positions) for seq in msa.homologs
  ]
  return MSA(parent, tuple(homologs))


def get_aligned_start_and_end_indices(
    aligned_parent: str,
    start_index_in_unaligned_parent: int,
    end_index_in_unaligned_parent: int,
) -> tuple[int, int]:
  """Returns the indices of the design region in the aligned parent sequence.

  Example:
    aligned_parent: 'K-A--CS'
    start_index_in_unaligned_parent: 1
    end_index_in_unaligned_parent: 3
    output: (2, 6)

  Explanation of example:
    The design region in the unaligned parent sequence is the substring
    'AC' which corresponds to the positions [2, 3, 4, 5] in the aligned parent
    sequence, which corresponds to start and end indices of 2 and 6.

  Args:
    aligned_parent: The aligned parent sequence.
    start_index_in_unaligned_parent: The index in the un-aligned parent sequence
      where the design region starts (inclusive).
    end_index_in_unaligned_parent: The index in the un-aligned parent sequence
      where the design region ends (exclusive).

  Returns:
    start_index_in_aligned_parent: The index in the aligned parent sequence
      where the design region starts (inclusive).
    end_index_in_aligned_parent: The index in the aligned parent sequence
      where the design region ends (exclusive).
  """

  aligned_index_to_unaligned_index = (
      np.cumsum([not is_gap(char) for char in aligned_parent]) - 1
  )

  in_design_region = np.logical_and(
      aligned_index_to_unaligned_index >= start_index_in_unaligned_parent,
      aligned_index_to_unaligned_index < end_index_in_unaligned_parent,
  )
  array_indices_in_design_region = in_design_region.nonzero()[0]
  start_index_in_unaligned_parent = array_indices_in_design_region[0]
  end_index_in_unaligned_parent = array_indices_in_design_region[-1] + 1
  return int(start_index_in_unaligned_parent), int(
      end_index_in_unaligned_parent
  )


def get_design_region_feasibility_fn(
    aligned_parent_seq: str,
    max_num_insertions: int,
    max_num_deletions: int,
    max_num_substitutions: int,
    design_region_start: int,
    design_region_end: int,
) -> Callable[[str], bool]:
  """Returns a function that checks if a sequence is meets the given constraints.

  To determine whether a sequence `seq` is feasible, we first
  select the sub-sequence of `seq` corresponding to the design region. Then, we
  compute the number of insertions, deletions, and substitutions in the design
  region relative to the design region of the aligned parent sequence. We
  require that these counts are less than or equal to the corresponding
  arguments to this function.

  Args:
    aligned_parent_seq: The aligned parent sequence.
    max_num_insertions: The maximum number of insertions allowed in the design
      region (inclusive).
    max_num_deletions: The maximum number of deletions allowed in the design
      region (inclusive).
    max_num_substitutions: The maximum number of substitutions allowed in the
      design region (inclusive).
    design_region_start: The index in the un-aligned parent sequence where the
      design region starts (inclusive).
    design_region_end: The index in the un-aligned parent sequence where the
      design region ends (exclusive).

  Returns:
    A function that takes an aligned sequence and returns a boolean indicating
    whether the sequence meets the above constraints.
  """
  aligned_design_region_start, aligned_design_region_end = (
      get_aligned_start_and_end_indices(
          aligned_parent_seq, design_region_start, design_region_end
      )
  )

  def _select_design_region_seq(aligned_sequence):
    return aligned_sequence[
        aligned_design_region_start:aligned_design_region_end
    ]

  aligned_design_region_parent_seq = _select_design_region_seq(
      aligned_parent_seq
  )

  def is_feasible(
      aligned_sequence,
  ):

    design_region_sequence = _select_design_region_seq(aligned_sequence)
    edit_counts = get_edit_counts(
        design_region_sequence, aligned_design_region_parent_seq
    )
    return (
        (edit_counts['insertion'] <= max_num_insertions)
        and (edit_counts['deletion'] <= max_num_deletions)
        and (edit_counts['substitution'] <= max_num_substitutions)
    )

  return is_feasible


def is_match(token: str) -> bool:
  return token == _MATCH_GAP_CHARACTER or token.isupper()


def get_spliced_sequence_df(
    aligned_seqs: Sequence[str],
    aligned_parent_seq: str,
    design_region_start: int,
    design_region_end: int) -> pd.DataFrame:
  """Splices the design region of `aligned_seqs` into `aligned_parent_seq`.

  For each sequence `seq` in `aligned_seqs`, we do the following:

  1) Determine the indices in `aligned_parent_seq` corresponding to the design
    region (this is non-trivial because `design_region_start` and
    `design_region_end` are w.r.t the un-aligned version of `parent_seq`).
  2) Determine the sub-sequence of `seq` corresponding to the design region
  3) Splice this sub-sequence into `aligned_parent_seq` at the appropriate
    location.
  4) Compute the edits in the design region w.r.t. `aligned_parent_seq`.

  Args:
    aligned_seqs: A list of aligned sequences.
    aligned_parent_seq: The aligned parent sequence.
    design_region_start: The start of the design region in the unaligned parent
      sequence.
    design_region_end: The end of the design region in the unaligned parent
      sequence.

  Returns:
    A `DataFrame` with the following columns:
    - aligned_sequence: The aligned sequence with the design region spliced into
      `aligned_parent_seq`.
    - edits: A tuple of (parent_value, position, value) corresponding to each
      edit in the design region.
    - substitutions: A list of (parent_value, position, value) corresponding to
      each substitution in the design region.
    - num_edits: The number of edits in the design region.
    - sequence: The unaligned version of `aligned_sequence`.
  """

  aligned_design_region_start, aligned_design_region_end = (
      get_aligned_start_and_end_indices(
          aligned_parent_seq, design_region_start, design_region_end
      )
  )
  rows = []
  for seq in aligned_seqs:
    spliced_aligned_seq = (
        aligned_parent_seq[:aligned_design_region_start]
        + seq[aligned_design_region_start:aligned_design_region_end]
        + aligned_parent_seq[aligned_design_region_end:]
    )
    edits = get_edits(spliced_aligned_seq, aligned_parent_seq)
    rows.append(
        dict(
            aligned_sequence=spliced_aligned_seq,
            edits=edits,
            substitutions=[
                edit for edit in edits if is_substitution(edit)
            ],
            num_edits=len(edits),
            sequence=get_unaligned_seq(spliced_aligned_seq),
        )
    )
  return pd.DataFrame(rows)


def _get_match_indices(aligned_seq: str) -> list[int]:
  return [i for i, token in enumerate(aligned_seq) if is_match(token)]


def _scatter_non_gap_indices(aligned_seq: str,
                             fill_value: int = -1) -> np.ndarray:
  to_return = np.full(len(aligned_seq), fill_value)
  idx = np.logical_not(np.isin(list(aligned_seq), _GAP_CHARACTERS))
  to_return[idx] = np.arange(idx.sum())
  return to_return


def _seqs_as_array(seqs: npt.ArrayLike) -> np.ndarray:
  if isinstance(seqs[0], str):
    return np.asarray([list(s) for s in seqs])
  return np.asarray(seqs)


class MatchStateSelector:
  """Helper for selecting sequence positions corresponding to match states."""

  def __init__(self, aligned_parent: str):
    """Constructs an instance of this class.

    Args:
      aligned_parent: A string corresponding to row in a multiple sequence
        alignment. Upper-case letters and '-' denote match states and lower-case
        letters and '.' denote insert states.
    """
    self._match_indices_alignment = _get_match_indices(aligned_parent)
    if not self._match_indices_alignment:
      raise ValueError('Input aligned_parent has 0 match states.')

    scattered_indices = _scatter_non_gap_indices(aligned_parent)
    self._match_indices_unaligned_parent = scattered_indices[
        self._match_indices_alignment]

  @property
  def num_match_positions(self) -> int:
    """Returns number of match states in `aligned_parent`."""
    return len(self._match_indices_alignment)

  def select_from_aligned_seqs(self, seqs: npt.ArrayLike) -> np.ndarray:
    """Selects match state positions from aligned `seqs`.

    Args:
      seqs: List of N sequences that are aligned to the `aligned_parent` used
        to construct the class. They can either be string or int encoded.

    Returns:
      np.ndarray of either strings or ints with shape
      `[N, self.num_match_positions]`.
    """
    seqs_array = _seqs_as_array(seqs)
    return seqs_array[:, self._match_indices_alignment]

  def select_from_unaligned_seqs(self, seqs: npt.ArrayLike) -> np.ndarray:
    """Selects match state positions from un-aligned `seqs`.

    Args:
      seqs: List of N string or int encoded sequences with length equal to the
        number of non-gap characters in the `aligned_parent` used to create
        the class. It is assumed that the positions in each seq in `seqs`
        corresponding to match states are identical to the positions of match
        states in the `aligned_parent` used to construct the class. This is
        appropriate when each element of `seqs` contains a small number of
        mutations from `aligned_parent`.

    Returns:
      List of N lists of either strings or ints, where each inner list
      is of length `self.num_match_positions`.
    """
    seqs_array = _seqs_as_array(seqs)
    return seqs_array[:, self._match_indices_unaligned_parent]


class MatchStateEncoder:
  """Wrapper for encoding unaligned seqs using a model built on match states."""

  def __init__(self, encode_fn, match_state_selector, batch_size=256):
    """Creates an instance of this class.

    Args:
      encode_fn: A function from a batch of int or string encoded sequences to
        an np.ndarray. The function will be called on the output of
        `match_state_selector.select_from_unaligned_seqs`.
      match_state_selector: A `MatchStateSelector`.
      batch_size: Batch size used when calling `encode_fn`
    """
    self._encode_fn = encode_fn
    self._match_state_selector = match_state_selector
    self._batch_size = batch_size

  def __call__(self, seqs: npt.ArrayLike) -> np.ndarray:
    """Selects the match states from `x` and applies `encode_fn` in batches."""
    match_state_seqs = self._match_state_selector.select_from_unaligned_seqs(
        seqs)
    return utils.batch_apply(
        self._encode_fn, match_state_seqs, batch_size=self._batch_size)


@dataclasses.dataclass(frozen=True)
class HomologConfig:
  """Container for objects used when modeling aligned homologous sequences.

  It is common for statistical models of homologs to only be defined at
  match state positions in the MSA. This container contains objects used to
  construct such models and process data for model fitting.

  Attributes:
    match_state_selector: MatchStateSelector used to select match state
      positions from aligned sequences.
    match_state_domain: Domain for the match state positions in the MSA.
    aligned_parent: String representation of the full parent sequence, including
      at non-match positions.
  """
  match_state_selector: MatchStateSelector
  match_state_domain: domains.DiscreteDomain
  aligned_parent: str

  @property
  def match_state_parent(self) -> str:
    """Returns the parent sequence only at match states."""
    return ''.join(self.match_state_selector.select_from_aligned_seqs(
        [self.aligned_parent]
    )[0])


def get_homolog_config_from_aligned_parent(
    aligned_parent: str, tokens: Sequence[str] = amino_acids.AA_GAP
) -> HomologConfig:
  """Builds a HomologConfig using an aligned parent sequence.

  Args:
    aligned_parent: A string representing an aligned sequence. Positions with
      capital letters or '-' denote match states and lower-case letters or '.'
      denote insert states.
    tokens: The list of tokens for the Vocabulary used to encode the sequences.

  Returns:
    A HomologConfig
  """

  match_state_selector = MatchStateSelector(aligned_parent)

  msa_vocab = domains.Vocabulary(tokens=tokens)
  match_state_domain = domains.DiscreteDomain.from_shared_vocab(
      match_state_selector.num_match_positions, msa_vocab
  )

  return HomologConfig(
      match_state_selector=match_state_selector,
      match_state_domain=match_state_domain,
      aligned_parent=aligned_parent,
  )


def get_homolog_config_from_a2m_file(
    a2m_file: str, tokens: Sequence[str] = amino_acids.AA_GAP
) -> HomologConfig:
  """Builds a HomologConfig using the first record in an a2m fasta file."""
  homolog_df = utils.read_fasta_as_df(a2m_file, num_records=1)
  aligned_parent = homolog_df.sequence.iloc[0]
  return get_homolog_config_from_aligned_parent(aligned_parent, tokens=tokens)


def save_homolog_config(homolog_config: HomologConfig, filename: str) -> None:
  """Saves a HomologConfig to a file."""
  with utils.open_file(filename, 'w') as f:
    json.dump(
        dict(aligned_parent=homolog_config.aligned_parent), f
    )


def load_homolog_config(filename: str) -> HomologConfig:
  """Loads a HomologConfig from a file."""
  with utils.open_file(filename, 'r') as f:
    aligned_parent = json.load(f)['aligned_parent']
  return get_homolog_config_from_aligned_parent(aligned_parent)
