# Copyright 2024 DeepMind Technologies Limited
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

"""Utilities for loading and using G4 CNN models."""

from os import path
from typing import Callable, Sequence

import numpy as np
from scipy import special
import tensorflow as tf

from nuclease_design import amino_acids
from nuclease_design import constants
from nuclease_design import utils


_ENSEMBLE_SIZE = 5


# The 4-class classification model has labels with the following interpretation.
# Let WT(alpha) be the alpha-percentile of the distribution of enrichment
# factors for WT synonyms. For a given observed enrichment factor EF, we have:
# <WT: EF < WT(alpha)
# WT: WT(alpha) <= EF < WT(1 - alpha)
# >WT: WT(1 - alpha) <= EF < A73R(alpha)
# >=A73R: EF >= A73R(alpha)

# Note that these are different from constants.LANDSCAPE_ACTIVITY_LEVELS.
# This is because at the time of modeling we had a different labeling approach
# that is more permissive (e.g., we don't do a false-discovery rate correction).

OUTPUT_CLASSES = ('<WT', 'WT', '>WT', '>=A73R')
_NUM_OUTPUT_CLASSES = len(OUTPUT_CLASSES)


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
               conv_layer_params: Sequence[tuple[int, int, int]],
               fc_layer_params: Sequence[int],
               dropout_layer_params: float | None):
    """Creates an instance of this class.

    Args:
      conv_layer_params: A list of tuples (filters, kernel_size, strides). Each
        entry results in a convolutional layer with the specified parameters.
      fc_layer_params: List of parameters for fully-connected layers, where each
        item is the number of units in the layer.
      dropout_layer_params: A dropout layer with this dropout rate will be added
        after each dense layer.
    """
    super().__init__()
    self._components = []

    kernel_initializer = tf.compat.v1.variance_scaling_initializer(
        scale=2.0, mode='fan_in', distribution='truncated_normal'
    )

    self._components.append(
        _conv_stack(
            conv_layer_params=conv_layer_params,
            activation='relu',
            kernel_initializer=kernel_initializer,
        )
    )

    self._components.append(tf.keras.layers.Flatten())

    self._components.append(
        _dense_stack(
            units=fc_layer_params,
            dropout=dropout_layer_params,
            activation='relu',
            kernel_initializer=kernel_initializer,
        )
    )

  def compute_output_shape(self, input_shape):
    output_shape = input_shape
    for component in self._components:
      output_shape = component.compute_output_shape(output_shape)
    return output_shape

  def call(self, x, **kwargs):
    for layer in self._components:
      x = layer(x, **kwargs)
    return x


def get_model(
    sequence_length: int,
    vocab_size: int,
    num_output_classes: int,
    conv_layer_params: Sequence[tuple[int, int, int]] = (
        (32, 5, 1),
        (32, 5, 1),
        (32, 5, 1),
    ),
    fc_layer_params: Sequence[int] = (64,),
    dropout_layer_params: float | None = 0.05,
) -> tf.keras.Model:
  """Builds a multi-class classification model with one `EncodingLayer`.

  Inputs are assumed to be int-encoded sequences. Outputs are predicted
  multi-class probabilities.

  Args:
    sequence_length: The name of the input sequences.
    vocab_size: The number of possible values that each position in the sequence
      can take on.
    num_output_classes: The number of output classes.
    conv_layer_params: A list of tuples (filters, kernel_size, strides). Each
      entry results in a convolutional layer with the specified parameters.
    fc_layer_params: List of parameters for fully-connected layers, where each
      item is the number of units in the layer.
    dropout_layer_params: A dropout layer with this dropout rate will be added
      after each dense layer.

  Returns:
    A `tf.keras.Model`.
  """
  inputs = tf.keras.layers.Input(sequence_length, dtype=tf.int32)
  encoded_inputs = tf.one_hot(inputs, vocab_size)
  embeddings = EncodingLayer(
      conv_layer_params=conv_layer_params,
      fc_layer_params=fc_layer_params,
      dropout_layer_params=dropout_layer_params,
  )(encoded_inputs)
  outputs = tf.keras.layers.Dense(
      units=num_output_classes, activation='softmax'
  )(embeddings)
  return tf.keras.Model(inputs=inputs, outputs=outputs)


def int_encode_sequences(sequences: list[str]) -> np.ndarray:
  residue_to_index = {
      residue: index for index, residue in enumerate(amino_acids.AA)
  }

  def _encode(sequence):
    return [residue_to_index[residue] for residue in sequence]

  return np.vstack([_encode(s) for s in sequences])


def get_cnn_filename(model_index: int) -> str:
  return f'cnn_models/model_{model_index}.h5'


def copy_all_models_to_local_dir(local_data_dir: str) -> None:
  for index in range(_ENSEMBLE_SIZE):
    utils.copy_file(get_cnn_filename(index), constants.DATA_DIR, local_data_dir)


def load_cnn_model(model_index: int, data_dir: str) -> tf.keras.Model:
  model = get_model(
      sequence_length=constants.VARIABLE_REGION_LENGTH,
      vocab_size=len(amino_acids.AA),
      num_output_classes=_NUM_OUTPUT_CLASSES,
  )
  try:
    model.load_weights(path.join(data_dir, get_cnn_filename(model_index)))
    return model
  except FileNotFoundError as exc:
    raise ValueError(
        'CNN models can only be loaded from a local file. '
        'You should use utils.copy_all_models_to_local_dir first.') from exc


def load_cnn_ensemble(data_dir: str) -> list[tf.keras.Model]:
  return [
      load_cnn_model(model_index=model_index, data_dir=data_dir)
      for model_index in range(_ENSEMBLE_SIZE)
  ]


def trim_sequence_to_variable_region(sequence: str) -> str:
  return sequence[constants.VARIABLE_REGION_START_INDEX:]


def get_predict_fn(
    model: tf.keras.Model,
) -> Callable[[Sequence[str]], np.ndarray]:
  def predict_fn(sequences: Sequence[str]) -> np.ndarray:
    sequences = [trim_sequence_to_variable_region(s) for s in sequences]
    int_sequences = int_encode_sequences(sequences)
    return model.predict(int_sequences)

  return predict_fn


def get_ensemble_predict_fn(
    models: Sequence[tf.keras.Model],
) -> Callable[[Sequence[str]], np.ndarray]:
  predict_fns = [get_predict_fn(m) for m in models]

  def predict(sequences: Sequence[str]) -> np.ndarray:
    predictions = [predict_fn(sequences) for predict_fn in predict_fns]
    predictions = np.stack(predictions, axis=-1)
    log_probs = np.log(predictions + 1e-6)
    return special.softmax(np.mean(log_probs, axis=-1))

  return predict
