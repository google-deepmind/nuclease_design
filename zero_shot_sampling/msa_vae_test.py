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

from absl.testing import absltest
import numpy as np
import tensorflow as tf

from zero_shot_sampling import alignment
from zero_shot_sampling import msa_vae
from zero_shot_sampling import vae


class MsaVaeTest(absltest.TestCase):

  def setUp(self):
    super().setUp()
    self._msa = alignment.MSA(parent='A-C', homologs=('ATC', 'TT-'))
    self._homolog_config = alignment.get_homolog_config_from_aligned_parent(
        self._msa.parent
    )
    self._training_config = vae.VAETrainingConfig(
        epochs=1,
        batch_size=2,
        optimizer_fn=lambda: tf.keras.optimizers.Adam(learning_rate=0.01),
    )
    self._architecture_hparams = dict(
        fc_layer_params=(2,),
        num_latents=2,
    )

  def assertNestedAlmostEqual(self, nested_array_1, nested_array_2):
    for array_1, array_2 in zip(
        tf.nest.flatten(nested_array_1), tf.nest.flatten(nested_array_2)
    ):
      np.testing.assert_allclose(array_1, array_2)

  def test_train_vae(self):
    model, history = msa_vae.train_vae(
        self._msa,
        self._homolog_config,
        self._training_config,
        self._architecture_hparams)

    self.assertIsInstance(model, vae.VAE)
    self.assertLen(history.history['loss'], 1)

  def test_save_and_load(self):
    model_path = self.create_tempdir().full_path

    model, _ = msa_vae.train_vae(
        self._msa,
        self._homolog_config,
        self._training_config,
        self._architecture_hparams,
    )

    msa_vae.save(
        model=model,
        model_path=model_path,
        homolog_config=self._homolog_config,
        architecture_hparams=self._architecture_hparams,
    )
    loaded_model, loaded_homolog_config = msa_vae.load(model_path)
    self.assertNestedAlmostEqual(
        model.model.get_weights(), loaded_model.model.get_weights()
    )
    self.assertEqual(
        self._homolog_config.aligned_parent,
        loaded_homolog_config.aligned_parent,
    )

  def test_get_likelihood_scores(self):
    model, _ = msa_vae.train_vae(
        self._msa,
        self._homolog_config,
        self._training_config,
        self._architecture_hparams,
    )
    # Note that additional test coverage for the values of the likelihood scores
    # is provided in vae_test.py.
    likelihood_scores = msa_vae.get_likelihood_scores(
        seqs=['A-C', 'ATC', 'TT-'],
        model=model,
        homolog_config=self._homolog_config,
        num_posterior_samples=2, batch_size=3
    )
    self.assertLen(likelihood_scores, 3)

    with self.assertRaisesRegex(ValueError, 'must have the same length'):
      _ = msa_vae.get_likelihood_scores(
          seqs=['A-CTT'],
          model=model,
          homolog_config=self._homolog_config,
      )


if __name__ == '__main__':
  absltest.main()
