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

"""Utilities for fitting VAE models on MSA data."""

import functools
import json
from os import path
from typing import Any, Mapping, Sequence

import numpy as np
import tensorflow as tf

from zero_shot_sampling import alignment
from zero_shot_sampling import domains
from zero_shot_sampling import utils
from zero_shot_sampling import vae


def train_vae(
    msa: alignment.MSA,
    homolog_config: alignment.HomologConfig,
    training_config: vae.VAETrainingConfig,
    architecture_hparams: Mapping[str, Any],
) -> tuple[vae.VAE, tf.keras.callbacks.History]:
  """Trains a VAE on the data from `MSA`.

  Args:
    msa: An alignment.MSA. The model will be trained on the parent and all
      homologs.
    homolog_config: A HomologConfig object used to process the MSA data.
    training_config: A VAETrainingConfig object containing the training
      hyper-parameters.
    architecture_hparams: A dictionary containing the non-default values for
      architecture hparams. See the `vae.VAE` constructor for more details.

  Returns:
    model: A `vae.VAE`
    history: A tf.keras.callbacks.History object containing the training
      history.
  """
  model = vae.VAE(
      domain=homolog_config.match_state_domain,
      optimizer=training_config.optimizer_fn(),
      **architecture_hparams
  )

  train_seqs = [msa.parent] + list(msa.homologs)

  match_state_train_seqs = (
      homolog_config.match_state_selector.select_from_aligned_seqs(train_seqs)
  )

  int_train_seqs = np.array(
      homolog_config.match_state_domain.encode(match_state_train_seqs)
  )

  history = model.fit(
      int_train_seqs,
      epochs=training_config.epochs,
      batch_size=training_config.batch_size,
      verbose=1,
  )

  return model, history


def save(
    model: vae.VAE,
    model_path: str,
    homolog_config: alignment.HomologConfig,
    architecture_hparams: Mapping[str, Any],
) -> None:
  """Saves a VAE model and its metadata to disk."""
  utils.make_directory(model_path)

  model.save_weights(path.join(model_path, 'model.weights.h5'))

  with utils.open_file(path.join(model_path, 'architecture'), 'w') as f:
    json.dump(architecture_hparams, f)

  alignment.save_homolog_config(
      homolog_config, path.join(model_path, 'homolog_config')
  )


def load(
    model_path: str,
) -> tuple[vae.VAE, alignment.HomologConfig]:
  """Saves a VAE model and its corresponding HomologConfig from disk."""
  homolog_config = alignment.load_homolog_config(
      path.join(model_path, 'homolog_config')
  )

  with utils.open_file(path.join(model_path, 'architecture'), 'r') as f:
    architecture_hparams = json.load(f)

  model = vae.VAE(
      domain=homolog_config.match_state_domain,
      **architecture_hparams,
  )

  # Build the model with a dummy input before loading weights
  placeholder_input = np.zeros((1, len(homolog_config.match_state_parent)), dtype=np.int32)
  _ = model.model(placeholder_input)

  model.load_weights(path.join(model_path, 'model.weights.h5'))

  return model, homolog_config


def get_likelihood_scores(
    seqs: Sequence[str],
    model: vae.VAE,
    homolog_config: alignment.HomologConfig,
    num_posterior_samples: int = 1000,
    batch_size: int = 16,
) -> list[float]:
  """Computes the likelihood scores of `seqs` under `model`.

  Args:
    seqs: A sequence of strings representing aligned sequences.
    model: A `vae.VAE`.
    homolog_config: A `HomologConfig` that the model was trained with respect
      to.
    num_posterior_samples: The number of posterior samples to use when computing
      the log-likelihood. See `vae.VAE.log_probability` for more details.
    batch_size: The batch size used when calling `model.log_probability`.

  Returns:
    A list of floats containing the likelihood scores of `seqs` under `model`.
  """
  for seq in seqs:
    if len(seq) != len(homolog_config.match_state_parent):
      raise ValueError(
          'All sequences must have the same length as'
          ' homolog_config.match_state_parent Required:'
          f' {len(homolog_config.match_state_parent)}, got: {len(seq)}'
      )
  encoder = alignment.MatchStateEncoder(
      functools.partial(
          model.log_probability, num_posterior_samples=num_posterior_samples
      ),
      homolog_config.match_state_selector,
      batch_size=batch_size,
  )

  unaligned_domain = domains.DiscreteDomain.from_shared_vocab(
      length=len(homolog_config.aligned_parent),
      vocab=homolog_config.match_state_domain.vocabs[0],
  )

  return encoder(unaligned_domain.encode(seqs)).tolist()
