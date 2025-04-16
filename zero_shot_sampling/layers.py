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

"""Utilities for building neural network layers."""

from collections.abc import Sequence

import numpy as np
import tensorflow as tf

from zero_shot_sampling import domains


def _conv_stack(
    conv_layer_params: Sequence[tuple[int, int, int]],
    kernel_initializer: tf.keras.initializers.Initializer,
    activation: str,
) -> tf.keras.Model:
  """Builds a stack of conv layers.

  Args:
    conv_layer_params: A list of tuples (filters, kernel_size, strides). Each
      entry results in a convolutional layer with the specified parameters. See
      the Tensorflow documentation of convolutional layers for more details.
    kernel_initializer: Initializer for creating the conv layers.
    activation: The activation function between layers.

  Returns:
    A sequential Keras model.
  """

  model = tf.keras.models.Sequential()
  for conv_layer_param in conv_layer_params:
    model.add(
        tf.keras.layers.Conv1D(
            filters=conv_layer_param[0],
            kernel_size=conv_layer_param[1],
            strides=conv_layer_param[2],
            activation=activation,
            padding='same',
            kernel_initializer=kernel_initializer,
        )
    )
  return model


def _dense_stack(
    units: Sequence[int],
    activation: str,
    kernel_initializer: tf.keras.initializers.Initializer,
    dropout: float | None,
) -> tf.keras.Model:
  """Builds a stack of dense/dropout layers.

  Args:
    units: A sequence of ints corresponding to the number of units of each
      dense layer in the stack.
    activation: The activation function of dense layers.
    kernel_initializer: Initializer for creating the dense layers.
    dropout: If not None, will add dropout layers with this rate after each
      dense layer.

  Returns:
    A sequential Keras model.
  """

  model = tf.keras.models.Sequential()
  for layer_units in units:
    model.add(
        tf.keras.layers.Dense(
            layer_units,
            activation=activation,
            kernel_initializer=kernel_initializer,
        )
    )
    if dropout is not None:
      model.add(tf.keras.layers.Dropout(rate=dropout))
  return model


class EncodingLayer(tf.keras.layers.Layer):
  """A layer for encoding sequence data.

  Applies a stack of convolutional and fully-connected layers.
  After the convolutional layers, per-position features are converted to a
  single feature vector for the entire sequence by flattening. Therefore, this
  model is not strictly convolutional, since the dense layers have
  position-specific parameters.
  """

  def __init__(self,
               conv_layer_params: Sequence[tuple[int, int, int]] | None,
               fc_layer_params: Sequence[int] | None,
               activation_fn: str,
               dropout_rate: float | None):
    """Creates an instance of this class.

    Args:
      conv_layer_params: A list of tuples (filters, kernel_size, strides). Each
        entry results in a convolutional layer with the specified parameters.
      fc_layer_params: List of parameters for fully-connected layers, where each
        item is the number of units in the layer.
      activation_fn: A keras-compatible string specifying an activation
        function.
      dropout_rate: A dropout layer with this dropout rate will be added after
        each dense layer.
    """
    super().__init__()
    self._components = []

    kernel_initializer = tf.compat.v1.variance_scaling_initializer(
        scale=2.0, mode='fan_in', distribution='truncated_normal'
    )

    if conv_layer_params:
      self._components.append(
          _conv_stack(
              conv_layer_params=conv_layer_params,
              activation=activation_fn,
              kernel_initializer=kernel_initializer,
          )
      )

    self._components.append(tf.keras.layers.Flatten())

    if fc_layer_params:
      self._components.append(
          _dense_stack(
              units=fc_layer_params,
              dropout=dropout_rate,
              activation=activation_fn,
              kernel_initializer=kernel_initializer,
          )
      )

  def compute_output_shape(self, input_shape: ...) -> ...:
    output_shape = input_shape
    for component in self._components:
      output_shape = component.compute_output_shape(output_shape)
    return output_shape

  def call(self, x, **kwargs):
    for layer in self._components:
      x = layer(x, **kwargs)
    return x


def get_mask_for_discrete_domain(domain: domains.DiscreteDomain) -> np.ndarray:
  """Returns a binary mask with True indicating valid token indices.

    `mask[i, j] == True` at position `i` and token index `j` corresponds to
    `j < domain.vocab_sizes[i]`, i.e. `j` is a valid token
    index at position `i`.

  Args:
    domain: A DiscreteDomain.

  Returns:
    A [domain.length, max(domain.vocab_sizes)] bool array.
  """
  mask = np.zeros(
      shape=(domain.length, max(domain.vocab_sizes)), dtype=bool)
  for i in range(domain.length):
    mask[i, :domain.vocab_sizes[i]] = True
  return mask


def get_masking_layer(
    mask: np.ndarray, mask_value: float, name: str | None = None
) -> tf.keras.layers.Layer:
  """Returns a keras layer that maps masked values to `mask_value`.

  The layer inputs a Tensor of size [batch_size, N, M] and outputs a Tensor
  of the same size that is a copy of the input, except that,
  `output[b, i, j] = max_value` for positions i, j where `mask` is False.

  Args:
    mask: An [N, M] bool array.
    mask_value: Value to place in positions that are False in `mask`.
    name: An optional name for the returned layer.

  Returns:
    A Keras layer that performs the masking.
  """

  mask = tf.constant(np.expand_dims(mask, 0), dtype=tf.float32)

  def apply_mask(logits):
    return  mask * logits + (1. - mask) * mask_value

  return tf.keras.layers.Lambda(apply_mask, name=name)
