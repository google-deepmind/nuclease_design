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

"""Specifications for discrete domains."""

from typing import Sequence

import numpy as np


_DEFAULT_VOCAB_NAME = 'vocab'


class Vocabulary:
  """Basic vocabulary used to represent output tokens for domains."""

  def __init__(self, tokens: Sequence[str] | int, name: str | None = None):
    """A token vocabulary.

    Args:
      tokens: An list of tokens to put in the vocab. If an int, will be
        interpreted as the number of tokens and '0', ..., 'tokens-1' will be
        used as tokens.
      name: An optional name for the Vocabulary.
    """
    if isinstance(tokens, (int, np.integer)):
      tokens = range(tokens)
    tokens = [str(token) for token in tokens]
    self._set_tokens(tokens)
    self._name = name if name else ''

  def _set_tokens(self, tokens):
    """Set self._tokens and related indices from list of token strings."""
    self._tokens = tokens
    self._token_ids = list(range(len(tokens)))
    self._id_to_token = dict(
        zip(self._token_ids, self._tokens))
    self._token_to_id = dict(
        zip(self._tokens, self._token_ids))

  def __len__(self):
    return len(self._tokens)

  @property
  def tokens(self) -> list[str]:
    """Return the tokens of the vocabulary."""
    return list(self._tokens)

  @property
  def name(self) -> str:
    """Return the name of the vocabulary."""
    return self._name

  @property
  def token_ids(self) -> list[int]:
    """Return the tokens ids of the vocabulary."""
    return list(self._token_ids)

  def is_valid(self, value: int) -> bool:
    """Tests if a value is a valid token id and returns a bool."""
    return value in self._token_ids

  def are_valid(self, values: Sequence[int]) -> np.ndarray:
    """Tests if values are valid token ids and returns an array of bools."""
    return np.array([self.is_valid(value) for value in values])

  def encode_token(self, token: str) -> int:
    """Maps a single character to an int."""
    return self._token_to_id[token]

  def encode(self, tokens: Sequence[str]) -> Sequence[int]:
    """Maps an iterable of string tokens to a list of integer token ids."""
    if isinstance(tokens, bytes):
      tokens = tokens.decode('utf-8')
    return [self.encode_token(token) for token in tokens]

  def decode_token(self, int_token: int) -> str:
    """Maps a single int to a character."""
    return self._id_to_token[int_token]

  def decode(
      self, values: Sequence[int]
  ) -> str:
    """Maps an iterable of integer token ids to strings."""
    tokens = []
    for value in values:
      value = int(value)  # Requires if value is a scalar tensor.
      tokens.append(self.decode_token(value))
    return ''.join(tokens)


class DiscreteDomain:
  """A domain for fixed-length sequences of categorical variables.

  Each position in the sequence has a separate Vocabulary.
  These can have with varying cardinalities.

  Example usage:
    vocabs = [Vocabulary(2, name='a'), Vocabulary(['x', 'y', 'z'], name='b')]
    domain = domains.DiscreteDomain(vocabs)

    domain.decode(domain.sample_uniformly(2))
    >> [['0', 'z'], ['1', 'y']]

    domain.encode([['0', 'y'], ['0', 'z']])
    >> [[0, 1], ['0, 2]]

    domain.vocab_sizes
    >> [2, 3]
  """

  def __init__(self, vocabs):
    """Creates an instance of this class.

    Args:
      vocabs: Sequence of Vocabularies for each position in the sequence. The
        Vocabularies can have differing cardinalities. All Vocabularies must
        have a non-None `name` field.
    """
    self._vocab_names = [str(vocab.name) for vocab in vocabs]
    if len(set(self._vocab_names)) != len(self._vocab_names):
      raise ValueError(
          f'Input Vocabulary names must be unique: {self._vocab_names}'
      )
    self._vocabs = vocabs
    self._vocab_sizes = [len(vocab) for vocab in vocabs]
    self._length = len(vocabs)

  @classmethod
  def from_shared_vocab(
      cls, length: int, vocab: Vocabulary
  ) -> 'DiscreteDomain':
    """Returns a `DiscreteDomain` of length `length` with a shared vocab.

    Args:
      length: The length of the domain.
      vocab: The vocabulary to use for all positions.

    Returns:
      A `DiscreteDomain` based on these vocabularies.
    """
    vocabs = []
    for i in range(length):
      name = '%s_%d' % (_DEFAULT_VOCAB_NAME, i)
      vocabs.append(Vocabulary(tokens=vocab.tokens, name=name))
    return cls(vocabs)

  @classmethod
  def from_vocab_sizes(cls, vocab_sizes: Sequence[int]) -> 'DiscreteDomain':
    names = ['p%d' % i for i in range(len(vocab_sizes))]
    return cls([Vocabulary(vocab_size, name=name) for vocab_size, name in
                zip(vocab_sizes, names)])

  def __str__(self):
    return '%s: length=%d vocab_sizes=%s' % (
        self.__class__.__name__, self.length, self.vocab_sizes)

  @property
  def length(self) -> int:
    return self._length

  def __len__(self):
    return self._length

  @property
  def vocab_sizes(self) -> Sequence[int]:
    """Returns the number of possible values at each sequence position."""
    return self._vocab_sizes

  @property
  def vocab_names(self) -> Sequence[str]:
    return self._vocab_names

  @property
  def vocabs(self) -> Sequence[Vocabulary]:
    return self._vocabs

  def _encode_sequence(self, sequence):
    """Maps a list of string tokens to a list of lists of integer token ids."""
    return [
        vocab.encode_token(token)
        for vocab, token in zip(self._vocabs, sequence)
    ]

  def encode(self, sequences: Sequence[str]) -> Sequence[Sequence[int]]:
    return [self._encode_sequence(sequence) for sequence in sequences]

  def _decode_sequence(self, sequence: Sequence[int], sep: str):
    """Maps list of lists of integer token ids to list of strings."""
    values = [
        vocab.tokens[token] for vocab, token in zip(self._vocabs, sequence)
    ]
    return sep.join([str(value) for value in values])

  def decode(
      self, sequences: Sequence[Sequence[int]], sep: str = ''
  ) -> Sequence[str]:
    return [
        self._decode_sequence(sequence, sep) for sequence in sequences
    ]

  def is_valid(self, sequence: Sequence[int]) -> bool:
    return len(sequence) == self.length and all(
        vocab.is_valid(token) for vocab, token in zip(self._vocabs, sequence))

  def are_valid(self, sequences: Sequence[Sequence[int]]) -> np.ndarray:
    """Tests if the given sequences are valid for this domain."""
    return np.array([self.is_valid(sequence) for sequence in sequences])

  def sample_uniformly(
      self, num_samples: int, random_state: np.random.RandomState | None = None
  ) -> np.ndarray:
    """Returns `num_sequences` sequences from the domain."""
    random_state = random_state or np.random.RandomState()
    return random_state.randint(
        0, self._vocab_sizes, (num_samples, self._length), dtype=np.int32)
