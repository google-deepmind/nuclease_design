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

"""Utility functions for preprocessing data."""

import re
from typing import Iterable, Tuple


# Mutation strings look like 'A45C'
_MUTATION_REGEX = re.compile(r'^[A-Z*]\d+[A-Z*]$')
_EMPTY_MUTANT_STRING = ''
_STOP_CODON = '*'

# By convention, the wildtype has no mutations and is thus an empty tuple
WILDTYPE_MUTATION_TUPLE = ()

# A MutationTuple is assumed to have the mutation position 1-indexed.
# e.g. ('A', 1, 'B')
# Using one-indexed mutations introduces vulnerability for bugs, but we found
# that it also reduced the chance of mis-communications.
MutationTuple = Tuple[str, int, str]


def _parse_mutant(mutant_string: str) -> MutationTuple:
  """Parses a mutation string into a tuple of (initial_value, index, final_value)."""
  if not _MUTATION_REGEX.match(mutant_string):
    raise ValueError(f'Malformed input: {mutant_string}')
  return (mutant_string[0], int(mutant_string[1:-1]), mutant_string[-1])


def parse_mutants(mutants_string: str) -> Iterable[MutationTuple]:
  """Parses a string encoding of multiple mutations.

  Args:
    mutants_string : A string, with comma delimited mutations, e.g., 'A45C,
      A56T'. The input '' is interpreted as the empty set of mutations.
      In keeping with comp-bio conventions, positions are one-indexed.

  Raises:
    ValueError: If the input string is malformed.

  Returns:
    A tuple of parsed mutation tuples, sorted by their position.
    e.g., (('A', 45, 'C'), ('A', 56, 'T')).
    The empty set of mutations is returned for an empty string input
    >>> parse_mutants('')
    ()
  """
  mutants_string = str(mutants_string)
  if mutants_string == _EMPTY_MUTANT_STRING:
    return WILDTYPE_MUTATION_TUPLE

  mutations = mutants_string.split(',')
  mutations = [_parse_mutant(s.replace(' ', '')) for s in mutations]
  # sort by mutation position
  mutations = tuple(sorted(mutations, key=lambda seq: seq[1]))
  return mutations


def introduces_stop_codon(mutations: Iterable[MutationTuple]) -> bool:
  """Returns True if any mutation in `mutations` introduces a stop codon.

  Args:
    mutations: An iterable of mutation tuples, e.g., (('A', 45, 'C'), ('A', 56,
      'T')).

  Returns:
    True if any mutation introduces a stop codon, False otherwise.

  """
  return any(m[2] == _STOP_CODON and m[0] != _STOP_CODON for m in mutations)


def mutates_stop_codon(mutations: Iterable[MutationTuple]) -> bool:
  """Returns True if any mutation in `mutations` mutates a stop codon.

  Example usage:
  >>> mutation = (('A', 45, 'C'), ('*', 56, 'T'))
  >>> mutates_stop_codon(mutation)
  True

  Args:
    mutations: An iterable of mutation tuples, e.g., (('A', 45, 'C'), ('*', 56,
      'T')).

  Returns:
    True if any mutation mutates a stop codon, False otherwise.
  """
  return any(m[0] == _STOP_CODON and m[2] != _STOP_CODON for m in mutations)


def apply_mutations_to_parent(mutations: Iterable[MutationTuple],
                              parent: str) -> str:
  """Mutates a parent sequence.

  Args:
    mutations: A list of MutationTuples.
    parent: The sequence to apply mutations to.

  Returns:
    The mutated parent sequence.

  Raises:
    AssertionError: If the expected starting string is not found in the
    parent sequence.
  """
  new_seq = list(parent)
  for mutation in mutations:
    old, position, new = mutation
    if position <= 0:
      raise ValueError('Invalid Mutation: position {}'.format(position))
    # Mutation tuples are 1 indexed.
    idx = position - 1
    if parent[idx] != old:
      raise ValueError(
          'Invalid Mutation: Mutation start {} does not match parent {}'.format(
              old, parent[idx]))
    new_seq[idx] = new
  return ''.join(new_seq)
