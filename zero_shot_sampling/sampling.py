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

"""Utilities for sampling sequences from a VAE model."""

from typing import Callable, Sequence

import numpy as np
import tensorflow as tf

from zero_shot_sampling import alignment
from zero_shot_sampling import vae


def _get_tf_seed(random_state: np.random.RandomState) -> np.ndarray:
  """Uses `random_state` to samples two integer seeds."""
  return random_state.randint(np.iinfo('int32').max, size=2)


def _sample_int_variant_batch_from_vae(
    num_samples: int, model: vae.VAE, parent_seq: Sequence[int],
    random_state: np.random.RandomState) -> np.ndarray:
  """Samples from `P(x|z)`, where `z ~ P(z|parent_seq)`."""
  parent_seq = np.asarray(parent_seq)
  batch_input = tf.expand_dims(parent_seq, 0)
  batch_latents = model.sample_latents_from_posterior(
      batch_input, num_samples, seed=_get_tf_seed(random_state))
  latents = tf.reshape(batch_latents, [num_samples, -1])
  return model.sample_structures_from_latents(
      latents, seed=_get_tf_seed(random_state)).numpy()


def sample_variants_from_vae(
    num_samples: int,
    vae_model: vae.VAE,
    homolog_config: alignment.HomologConfig,
    random_state: np.random.RandomState,
    sampling_batch_size: int = 512,
) -> list[str]:
  """Samples variants of an aligned parent sequence from a VAE.

  Instead of sampling from the raw VAE distribution, which may produce samples
  with very high diversity, we sample from
  P(x|z)`, where `z ~ P(z|parent_seq)`, where `parent_seq` is the subsequence of
  the aligned parent sequence defining `homolog_config` at its
  match state positions.

  Args:
    num_samples: The number of samples to return.
    vae_model: The VAE model to use for sampling.
    homolog_config: A `alignment.HomologConfig` object.
    random_state: A `np.random.RandomState` for controlling randomness.
    sampling_batch_size: The batch size to use for sampling from the VAE.

  Returns:
    A list of `num_samples` samples from the VAE.
  """
  encoded_match_state_parent = homolog_config.match_state_domain.encode(
      [homolog_config.match_state_parent]
  )[0]

  samples = []
  while len(samples) < num_samples:
    int_samples = _sample_int_variant_batch_from_vae(
        sampling_batch_size,
        vae_model,
        encoded_match_state_parent,
        random_state=random_state,
    )
    samples.extend(homolog_config.match_state_domain.decode(int_samples))

  return samples[:num_samples]


def rejection_sample(
    sample_fn: Callable[[np.random.RandomState], str],
    feasibility_fn: Callable[[str], bool],
    num_samples: int,
    random_state: np.random.RandomState,
    max_num_proposals: int = int(1e6),
    verbose: bool = True,
) -> list[str]:
  """Generic rejection sampler.

  Args:
    sample_fn: A function that takes a `np.random.RandomState` and returns a
      sampled string.
    feasibility_fn: A function that takes a string and returns a boolean.
      Samples for which `feasibility_fn(sample) == True` are considered
      accepted.
    num_samples: The number of samples to return.
    random_state: A `np.random.RandomState` for controlling randomness.
    max_num_proposals: The maximum number of proposals to consider. If this
      budget is exceeded, an error will be raised.
    verbose: Whether to print diagnostic information about the acceptance rate
      of the sampling.

  Returns:
    A list of accepted samples.

  Raises:
    RuntimeError: If the maximum number of proposals is exceeded without
      accepting enough samples.
  """
  accepted_samples = []
  num_proposals = 0
  while len(accepted_samples) < num_samples:
    if num_proposals > max_num_proposals:
      raise RuntimeError(
          f'Only {len(accepted_samples)} samples were accepted out of a budget'
          f' of max_num_proposals = {max_num_proposals}'
      )
    samples = sample_fn(random_state)
    num_proposals += len(samples)

    for sample in samples:
      if feasibility_fn(sample):
        accepted_samples.append(sample)

  if verbose:
    print(
        f'Accepted {len(accepted_samples)} samples out of'
        f' {num_proposals} total.'
    )
  return accepted_samples
