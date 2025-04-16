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

"""Variational Auto-Encoder."""

import dataclasses
from typing import Callable, Sequence

import numpy as np
from numpy import typing as npt
import tensorflow as tf
import tensorflow_probability as tfp

from zero_shot_sampling import domains
from zero_shot_sampling import layers

_MIN_FLOAT = np.finfo(np.float32).min

# Castable to a Tensor.
TensorLike = npt.ArrayLike | tf.Tensor

# Seed used to control randomness in Tensorflow.
# Guide here: https://www.tensorflow.org/api_docs/python/tf/random/set_seed
# The seed should either be an int or a length-2 TensorOrArray.
TFSeed = TensorLike


def identity_loss(y_true: TensorLike, y_pred: TensorLike) -> TensorLike:
  """Identity loss function."""
  del y_true
  return y_pred


def summed_categorical_crossentropy(
    y_true: TensorLike, y_pred: TensorLike
) -> TensorLike:
  """Summed categorical cross-entropy loss function."""
  return tf.reduce_sum(
      tf.keras.losses.categorical_crossentropy(y_true, y_pred), axis=-1)


def mean_categorical_crossentropy(
    y_true: TensorLike, y_pred: TensorLike
) -> TensorLike:
  """Mean categorical cross-entropy loss function."""
  return tf.reduce_mean(
      tf.keras.losses.categorical_crossentropy(y_true, y_pred), axis=-1)


def importance_log_ratio(
    z_mean: tf.Tensor, z_log_std: tf.Tensor, latent_vectors: tf.Tensor
) -> tf.Tensor:
  """Computes the log of the prior/posterior importance ratio.

  Args:
    z_mean: A numpy array or tensor of shape (batch_size, num_latents)
      containing the mean of the posterior.
    z_log_std: A numpy array or tensor of shape (batch_size, num_latents)
      containing the log_std of the posterior.
    latent_vectors: A numpy array or tensor of shape (batch_size, num_latents)
      containing the latent vectors.

  Returns:
    A tensor of shape (batch_size,) containing the logarithm of the
      prior/posterior ratio.
  """
  prior = tfp.distributions.Normal(loc=tf.zeros(z_mean.shape),
                                   scale=tf.ones(z_log_std.shape))
  prior_log_prob = prior.log_prob(latent_vectors)
  posterior = tfp.distributions.Normal(loc=z_mean,
                                       scale=tf.exp(z_log_std))
  posterior_log_prob = posterior.log_prob(latent_vectors)
  ratio = prior_log_prob - posterior_log_prob
  # This is a sum over the dimensions in the latent space.
  # We want log (p(z)/q(z))
  # = log (\prod_i p(z_i) / \prod_i q(z_i))
  # = log (\prod_i p(z_i)) - log (\prod_i q(z_i))
  # = sum_i log p(z_i) - sum_i log q(z_i) = sum_i (log p(z_i) - log q(z_i)).
  ratio = tf.reduce_sum(ratio, axis=1)
  return ratio


def build_decoder_output_model(
    domain: domains.DiscreteDomain,
) -> tf.keras.Model:
  """Builds a Keras model with the output layers of the VAE decoder."""
  padded_vocab_size = max(domain.vocab_sizes)
  model = tf.keras.Sequential([
      tf.keras.layers.Dense(domain.length * padded_vocab_size),
      tf.keras.layers.Reshape((domain.length, padded_vocab_size),),
      tf.keras.layers.Dense(padded_vocab_size)
  ])
  if len(set(domain.vocab_sizes)) > 1:
    mask = layers.get_mask_for_discrete_domain(domain)
    model.add(layers.get_masking_layer(mask, mask_value=_MIN_FLOAT,
                                       name='mask_logits'))

  model.add(tf.keras.layers.Softmax(axis=-1))

  return model


class VAE:
  """Variational Auto-Encoder."""

  def __init__(
      self,
      domain: domains.DiscreteDomain,
      fc_layer_params: Sequence[int] = (64,),
      num_latents: int = 32,
      activation_fn: str = 'elu',
      optimizer: tf.keras.optimizers.Optimizer | None = None,
      loss_fn: ... = summed_categorical_crossentropy,
      tf_compile: bool = True,
      name: str = 'Vae',
  ):
    """Creates an instance of this class.

    Args:
      domain: An instance of a `domains.DiscreteDomain`.
      fc_layer_params: A sequence of ints specifying the number of units in each
        layer of the fully_connected_encoder. The decoder uses the same sizes.
      num_latents: The number of latent dimensions.
      activation_fn: A keras-compatible string name for an activation function.
      optimizer: Either a `tf.keras.optimizers.Optimizer` or None, in which
        case `tf.keras.optimizers.Adam(learning_rate=0.001)` is used.
      loss_fn: The reconstruction loss function, e.g.
        `summed_categorical_crossentropy` or `mean_categorical_crossentropy`.
      tf_compile: Whether the fit() and log_probability() functions should be
        compiled using tf.function.
      name: The name of the Keras model.
    """
    self._domain = domain

    # Encoder
    encoder_input = tf.keras.layers.Input((domain.length,), dtype=tf.int32)
    one_hot_encoder_input = tf.one_hot(encoder_input, max(domain.vocab_sizes))
    x = layers.EncodingLayer(
        conv_layer_params=None,
        fc_layer_params=fc_layer_params,
        activation_fn=activation_fn,
        dropout_rate=0,
    )(one_hot_encoder_input)

    mu = tf.keras.layers.Dense(num_latents)(x)
    log_std = tf.keras.layers.Dense(num_latents)(x)

    encoder = tf.keras.Model(inputs=encoder_input, outputs=[mu, log_std],
                             name='%s/encoder' % name)

    # Decoder
    decoder_input = tf.keras.layers.Input((num_latents,))
    x = layers.EncodingLayer(
        conv_layer_params=None,
        fc_layer_params=fc_layer_params[::-1],
        activation_fn=activation_fn,
        dropout_rate=0,
    )(decoder_input)
    x = build_decoder_output_model(domain)(x)
    decoder = tf.keras.Model(inputs=decoder_input, outputs=x,
                             name='%s/decoder' % name)

    # VAE
    eps = tf.random.normal(shape=tf.shape(mu))
    std = tf.exp(log_std)
    z = eps * std + mu
    x_pred = decoder(z)

    kl = -0.5 * tf.reduce_sum(1 + log_std - mu**2 - tf.exp(log_std), axis=1)

    vae = tf.keras.models.Model(inputs=[encoder_input],
                                outputs=[x_pred, kl],
                                name='%s/base' % name)

    loss = [loss_fn, identity_loss]
    optimizer = optimizer or tf.keras.optimizers.Adam(learning_rate=0.001)
    vae.compile(optimizer=optimizer, loss=loss, run_eagerly=(not tf_compile))
    if tf_compile:
      self._log_probability_as_tensor = tf.function(
          self._log_probability_as_tensor)

    self._encoder = encoder
    self._decoder = decoder
    self._model = vae
    self._num_latents = num_latents
    self._loss_fn = loss_fn

  @property
  def domain(self) -> domains.DiscreteDomain:
    return self._domain

  @property
  def model(self) -> tf.keras.Model:
    """Returns the underlying tf.Keras model."""
    return self._model

  def _get_targets(self, structures: TensorLike) -> tf.Tensor:
    return tf.one_hot(structures, max(self._domain.vocab_sizes))

  def sample_from_prior(
      self,
      num_samples: int,
      mean: float = 0.,
      stddev: float = 1.,
      seed: TFSeed | None = None) -> tf.Tensor:
    """Samples from the prior distribution over latent variables.

    Args:
      num_samples: An int indicating the number of samples to return.
      mean: The mean of the normal distribution.
      stddev: A standard deviation of the normal distribution. Note that we used
        stddev equal to 1 for training, but we might want to draw samples closer
        to/further from the mean.
      seed: A `TFSeed` for controlling randomness.

    Returns:
      A numpy array of shape (num_samples, num_latents).
    """

    return tfp.distributions.Normal(
        loc=mean, scale=stddev).sample(
            sample_shape=(num_samples, self._num_latents), seed=seed)

  def sample_latents_from_posterior(
      self,
      structures: TensorLike,
      num_samples_per_structure: int,
      seed: TFSeed | None = None) -> tf.Tensor:
    """Samples latent variables from posterior given observed structures.

    Args:
      structures: A numpy array or tensor of shape (batch_size,
        sequence_length).
      num_samples_per_structure: Number of posterior samples per structure.
      seed: A `TFSeed` for controlling randomness.

    Returns:
      A (batch_size, num_samples_per_structure, num_latents) np array.
    """
    mu, log_std = self._encoder(structures)
    std = tf.exp(log_std)

    # Tile because each posterior sample for a given structure uses the same
    # mean and stddev.
    def _tile(x):
      return tf.tile(tf.expand_dims(x, 1), [1, num_samples_per_structure, 1])

    mu = _tile(mu)
    std = _tile(std)
    return tfp.distributions.Normal(loc=mu, scale=std).sample(seed=seed)

  def sample_structures_from_latents(
      self,
      latents: TensorLike,
      greedy: bool = False,
      seed: TFSeed | None = None) -> tf.Tensor:
    """Samples structures given latents.

    Args:
      latents: A numpy array or tensor of shape (batch_size, num_latents).
      greedy: Whether to discretize the decoder output deterministically by
        using argmax instead of sampling from a multinomial distribution.
      seed: A `TFSeed` for controlling randomness.

    Returns:
      (batch_size, sequence_dim) int array of samples.
    """
    probs = self._decoder(latents)
    tf.debugging.assert_all_finite(probs,
                                   'Probabilities include NaNs! %s' % probs)
    if greedy:
      return tf.argmax(probs, -1)
    else:
      return tfp.distributions.Categorical(probs=probs).sample(seed=seed)

  def sample(
      self,
      num_samples: int,
      greedy: bool = False,
      latents_seed: TFSeed | None = None,
      structures_seed: TFSeed | None = None,
  ) -> tf.Tensor:
    """Draws `num_samples` from the model.

    Args:
      num_samples: Number of samples to draw.
      greedy: Whether to discretize the decoder output deterministically
        by using argmax instead of sampling from a multinomial
        distribution.
      latents_seed: A `TFSeed` for controlling randomness when sampling
        latent factors from the prior.
      structures_seed: A `TFSeed` for controlling randomness when
        sampling structures conditional on the latent factors.

    Returns:
      An int Tensor containing sampled structures.
    """
    latents = self.sample_from_prior(num_samples, seed=latents_seed)
    return self.sample_structures_from_latents(
        latents, greedy, seed=structures_seed)

  def fit(
      self,
      structures: TensorLike,
      sample_weight: TensorLike | None = None,
      **kwargs,
  ) -> tf.keras.callbacks.History:
    inputs = np.asarray(structures)
    targets = self._get_targets(structures)
    if sample_weight is not None:
      # Duplicate weights since the model has two outputs (structures, kl).
      sample_weight = [sample_weight] * 2
    return self.model.fit(
        inputs,
        [targets, tf.zeros((len(structures),))],
        sample_weight=sample_weight,
        **kwargs,
    )

  def _conditional_log_probability(
      self, structures: TensorLike, latents: TensorLike
  ) -> tf.Tensor:
    """Computes the log-probability of sequences conditioned on latent vectors.

    Args:
      structures: A numpy array or tensor of shape
        (batch_size, sequence_length).
      latents: A numpy array or tensor of shape (batch_size, num_latents)
        containing the latent vectors.

    Returns:
      A tensor of shape (batch_size,) containing the conditional
        log-probabilities per sequence.
    """
    probs = self._decoder(latents)
    return -self._loss_fn(
        self._get_targets(structures), probs)

  def _log_probability_as_tensor(
      self,
      structures: TensorLike,
      num_posterior_samples: int,
      seed: TFSeed | None = None,
  ) -> tf.Tensor:
    """Computes log probabilities and can be compiled with `tf.function`."""

    structures = tf.convert_to_tensor(structures)
    batch_size = structures.shape[0]

    def _tile(array):
      return tf.tile(array, (num_posterior_samples, 1))

    z_mean, z_log_std = self._encoder(structures)
    z_mean_tiled = _tile(z_mean)
    z_log_std_tiled = _tile(z_log_std)
    z_tiled = tfp.distributions.Normal(
        z_mean_tiled, tf.exp(z_log_std_tiled)).sample(seed=seed)

    samples_tiled = _tile(structures)
    cond_log_prob_tiled = self._conditional_log_probability(
        samples_tiled, z_tiled)
    log_ratio_tiled = importance_log_ratio(z_mean_tiled, z_log_std_tiled,
                                           z_tiled)
    log_prob_tiled = cond_log_prob_tiled + log_ratio_tiled

    # Reshape the flattened tensor into shape (num_posterior_samples,
    # batch_size), consistent with the description in the _posterior_parameters
    # method and aggregate over the posterior samples, returning an array of
    # shape (batch_size,): one log-probability value per input sequence.
    log_prob_tiled = tf.reshape(log_prob_tiled,
                                (num_posterior_samples, batch_size))
    log_prob = tf.math.reduce_logsumexp(log_prob_tiled, axis=0)
    log_prob -= np.log(num_posterior_samples)

    return log_prob

  def log_probability(
      self,
      structures: TensorLike,
      num_posterior_samples: int = 1000,
      seed: TFSeed | None = None,
  ) -> np.ndarray:
    r"""Computes the log probability via importance sampling.

    Uses the following importance sampling-based approximation:
      log \sum_i p(x, z_i)/q(z_i),
        where z_i ~ q(z_i), i = 1 ... num_posterior_samples.

    Args:
      structures: A numpy array or tensor of shape (batch_size, sequence_length)
        or a list of (sequence_length,) arrays.
      num_posterior_samples: An int indicating the number of structures from the
        posterior to use for computing the importance sampling estimate.
      seed: A `TFSeed` for controlling randomness.

    Returns:
      A numpy array of shape (batch_size,) containing log probabilities.
    """
    return self._log_probability_as_tensor(
        structures, num_posterior_samples, seed
    ).numpy()

  def save_weights(self, path: str, **kwargs) -> None:
    """Saves the weights of the model to `path`.

    Args:
      path: The path to save the weights to.
      **kwargs: See tf.keras.Model.save_weights().
    """
    self.model.save_weights(path, **kwargs)

  def load_weights(self, path: str, **kwargs) -> None:
    """Loads the weights of the model from `path`.

    Args:
      path: The path to load the weights from.
      **kwargs: See tf.keras.Model.load_weights().
    """
    self.model.load_weights(path, **kwargs)


@dataclasses.dataclass(frozen=True)
class VAETrainingConfig:
  epochs: int
  batch_size: int
  optimizer_fn: Callable[[], tf.keras.optimizers.Optimizer]
