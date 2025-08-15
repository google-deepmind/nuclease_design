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

from unittest import mock

from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import tensorflow as tf

from zero_shot_sampling import domains
from zero_shot_sampling import vae


def _build_vae(
    domain,
    tf_compile=False,
    num_latents=4,
    fc_layer_params=(5,),
    learning_rate=0.01,
    **kwargs,
):
  # Use tf_compile=False as the default because it avoids the time overhead
  # of compilation.
  return vae.VAE(
      domain=domain,
      tf_compile=tf_compile,
      num_latents=num_latents,
      fc_layer_params=fc_layer_params,
      optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate),
      **kwargs,
  )


class VaeTest(parameterized.TestCase):

  @parameterized.parameters(([2, 3, 2], 1), ([3, 3, 3], 13))
  def test_encoder_outputs(self, vocab_sizes, num_latents):
    num_samples = 10
    domain = domains.DiscreteDomain.from_vocab_sizes(vocab_sizes)
    model = _build_vae(domain=domain, num_latents=num_latents)
    x = domain.sample_uniformly(num_samples)
    y = model._encoder(x)

    self.assertLen(y, 2)
    mu, log_std = y
    np.testing.assert_array_equal(
        (num_samples, num_latents), mu.shape.as_list()
    )
    np.testing.assert_array_equal(
        (num_samples, num_latents), log_std.shape.as_list()
    )

  @parameterized.parameters(([2, 2, 2], 1), ([3, 3, 3], 13))
  def test_decoder_outputs(self, vocab_sizes, num_latents):
    num_samples = 10
    domain = domains.DiscreteDomain.from_vocab_sizes(vocab_sizes)
    model = _build_vae(domain=domain, num_latents=num_latents)
    eps = tf.random.uniform((num_samples, num_latents))
    y = model._decoder(eps)
    np.testing.assert_array_equal(
        (num_samples, domain.length, max(domain.vocab_sizes)), y.shape.as_list()
    )

  @parameterized.parameters(([2, 3, 2], 1), ([3, 3, 3], 13))
  def test_vae_outputs(self, vocab_sizes, num_latents):
    num_samples = 10
    domain = domains.DiscreteDomain.from_vocab_sizes(vocab_sizes)
    model = _build_vae(domain=domain, num_latents=num_latents)
    x_pred, kl = model.model(domain.sample_uniformly(num_samples))
    np.testing.assert_array_equal(
        (num_samples, domain.length, max(domain.vocab_sizes)),
        x_pred.shape.as_list(),
    )
    np.testing.assert_array_equal((num_samples,), kl.shape.as_list())

  @parameterized.parameters((range(2, 6),), ([3, 3, 3, 3],))
  def test_fit(self, vocab_sizes):
    num_samples = 4
    domain = domains.DiscreteDomain.from_vocab_sizes(vocab_sizes)
    model = _build_vae(
        domain=domain,
        fc_layer_params=(2,),
        num_latents=8,
        learning_rate=0.001
    )

    structures = domain.sample_uniformly(num_samples)
    model.fit(structures, epochs=5, batch_size=num_samples)
    latents = model.sample_latents_from_posterior(
        structures, num_samples_per_structure=1
    )[:, 0, :]

    samples = model.sample_structures_from_latents(latents)
    np.testing.assert_array_equal(
        (num_samples, domain.length),
        samples.shape.as_list(),
    )

  def test_fit_from_list_inputs(self):
    domain = domains.DiscreteDomain.from_vocab_sizes([3, 4, 5])
    model = _build_vae(domain=domain)
    model.fit(list(domain.sample_uniformly(5)), epochs=1)

  @parameterized.parameters((True,), (False,))
  def test_sample_structures_from_latents(self, greedy):
    num_latents = 3
    num_samples = 5
    domain = domains.DiscreteDomain.from_vocab_sizes([3, 3, 5, 2, 1],)
    latents = np.random.normal(size=(num_samples, num_latents))
    model = _build_vae(domain=domain, num_latents=num_latents)

    structures = model.sample_structures_from_latents(latents, greedy)
    np.testing.assert_array_equal(
        (num_samples, domain.length), structures.shape
    )
    self.assertTrue(all(domain.are_valid(structures)))

  @parameterized.parameters(((0, 1),), ((1, 2),))
  def test_sample_structures_from_latents_deterministic(self, seed):
    num_latents = 3
    domain = domains.DiscreteDomain.from_vocab_sizes([3, 3, 5, 2, 1],)
    latents = np.random.uniform(size=(25, num_latents))
    model = _build_vae(domain=domain, num_latents=num_latents)

    def _sample(seed):
      return model.sample_structures_from_latents(latents, seed=seed)

    np.testing.assert_array_equal(_sample(seed), _sample(seed))
    different_seed = np.array(seed) + 1
    self.assertTrue(np.any(_sample(seed) != _sample(different_seed)))

  @parameterized.parameters(((0, 1),), ((1, 2),))
  def test_sample_from_prior_deterministic(self, seed):
    domain = domains.DiscreteDomain.from_vocab_sizes([3, 3, 5, 2, 1],)
    model = _build_vae(domain=domain)

    def _sample(seed):
      return model.sample_from_prior(num_samples=25, seed=seed)

    np.testing.assert_array_equal(_sample(seed), _sample(seed))
    different_seed = np.array(seed) + 1
    self.assertTrue(np.any(_sample(seed) != _sample(different_seed)))

  @parameterized.parameters(
      ([3, 3, 5, 2, 1], 5, False),
      ([3, 3, 5, 2, 1], 5, True),
  )
  def test_sample(self, vocab_sizes, num_samples, greedy):
    domain = domains.DiscreteDomain.from_vocab_sizes(vocab_sizes)
    model = _build_vae(domain=domain)
    structures = model.sample(num_samples, greedy=greedy)
    np.testing.assert_array_equal(
        (num_samples, domain.length), structures.shape
    )
    self.assertTrue(all(domain.are_valid(structures)))

  @parameterized.parameters(
      ((0, 1), (2, 3),),
      ((1, 2), (3, 4)))
  def test_sample_deterministic(self, latents_seed, structures_seed):
    domain = domains.DiscreteDomain.from_vocab_sizes([3, 3, 5, 2, 1])
    model = _build_vae(domain=domain)

    def _sample(latents_seed, structures_seed):
      return model.sample(
          25, latents_seed=latents_seed, structures_seed=structures_seed
      )

    np.testing.assert_array_equal(
        _sample(latents_seed, structures_seed),
        _sample(latents_seed, structures_seed),
    )

    different_latents_seed = np.array(latents_seed) + 1
    self.assertTrue(
        np.any(
            _sample(latents_seed, structures_seed)
            != _sample(different_latents_seed, structures_seed)
        )
    )

  @parameterized.parameters((2, 2, 5), (3, 1, 4))
  def test_sample_latents_from_posterior(self, batch_size, num_latents,
                                         num_samples):
    domain = domains.DiscreteDomain.from_vocab_sizes([3, 3, 5, 2, 1])
    model = _build_vae(domain=domain, num_latents=num_latents)
    inputs = domain.sample_uniformly(batch_size)
    latents = model.sample_latents_from_posterior(inputs, num_samples)
    np.testing.assert_array_equal(
        latents.shape, (batch_size, num_samples, num_latents)
    )

  @parameterized.parameters(((0, 1),), ((1, 2),))
  def test_sample_latents_from_posterior_deterministic(self, seed):
    domain = domains.DiscreteDomain.from_vocab_sizes([3, 3, 5, 2, 1])
    inputs = domain.sample_uniformly(3)
    model = _build_vae(domain=domain)

    def _sample(seed):
      return model.sample_latents_from_posterior(inputs, 4, seed=seed)

    np.testing.assert_array_equal(_sample(seed), _sample(seed))
    different_seed = np.array(seed) + 1
    self.assertTrue(np.any(_sample(seed) != _sample(different_seed)))

  def assertNestedAlmostEqual(self, nested_array_1, nested_array_2):
    for array_1, array_2 in zip(
        tf.nest.flatten(nested_array_1), tf.nest.flatten(nested_array_2)
    ):
      np.testing.assert_allclose(array_1, array_2)

  def test_save_load_weights(self):
    domain = domains.DiscreteDomain.from_vocab_sizes(range(2, 6))
    x = domain.sample_uniformly(10)

    def _get_vae():
      return vae.VAE(domain=domain, fc_layer_params=(32, 16), num_latents=8)

    vae1 = _get_vae()
    vae1.fit(x)
    path = self.create_tempfile().full_path
    vae1.model.save_weights(path)

    vae2 = _get_vae()
    vae2.model.load_weights(path)

    self.assertNestedAlmostEqual(
        vae1.model.get_weights(), vae2.model.get_weights()
    )
    np.testing.assert_allclose(
        vae1.log_probability(x, seed=(0, 0)),
        vae2.log_probability(x, seed=(0, 0)),
    )

  def test_importance_log_ratio_posterior_equal_prior(self):
    batch_size = 5
    latent_size = 3
    z_mean = np.zeros((batch_size, latent_size), dtype=np.float32)
    z_log_std = np.zeros((batch_size, latent_size), dtype=np.float32)
    z_sample = np.random.normal(size=(batch_size, latent_size))
    importance_ratio = np.exp(
        vae.importance_log_ratio(z_mean, z_log_std, z_sample))
    np.testing.assert_allclose(importance_ratio, [1.0] * batch_size)

  def test_importance_log_ratio(self):
    batch_size = 5
    latent_size = 3
    z_mean = np.random.normal(size=(batch_size, latent_size))
    z_mean = z_mean.astype(dtype=np.float32)
    z_log_std = np.zeros((batch_size, latent_size), dtype=np.float32)
    # Given the same variance, the probability for the posterior mean is higher
    # under the posterior than under the prior.
    importance_ratio = np.exp(
        vae.importance_log_ratio(z_mean, z_log_std, z_mean))
    np.testing.assert_array_less(importance_ratio, 1.0)

  @parameterized.parameters(
      (range(2, 6),),
      ([3, 3, 3, 3],),
      ([3, 3, 3, 3], True)
      )
  def test_log_probability(self, vocab_sizes, tf_compile=False):
    num_samples = 10
    domain = domains.DiscreteDomain.from_vocab_sizes(vocab_sizes)
    model = _build_vae(domain=domain, tf_compile=tf_compile)
    pos_samples = tf.ones((num_samples, domain.length), dtype=tf.int32)
    neg_samples = tf.zeros((num_samples, domain.length), dtype=tf.int32)
    model.fit(pos_samples, epochs=100, verbose=False)

    log_probs_pos = model.log_probability(pos_samples)
    log_probs_neg = model.log_probability(neg_samples)
    np.testing.assert_array_equal((num_samples,), log_probs_neg.shape)
    np.testing.assert_array_equal((num_samples,), log_probs_pos.shape)
    self.assertTrue(np.all(log_probs_neg < log_probs_pos))

  def test_log_probability_from_list_inputs(self):
    domain = domains.DiscreteDomain.from_vocab_sizes(range(2, 6))
    model = _build_vae(domain=domain)
    pos_samples = tf.ones((4, domain.length), dtype=tf.int32)
    model.fit(pos_samples, epochs=2, verbose=False)
    model.log_probability(list(domain.sample_uniformly(5)))

  def test_log_probability_of_invalid_inputs(self):
    domain = domains.DiscreteDomain.from_vocab_sizes(range(2, 6))
    model = _build_vae(domain=domain)
    for i in [2, 3, 4]:
      log_probability = model.log_probability(np.full((1, domain.length), i))
      np.testing.assert_array_less(log_probability, -5)

  @parameterized.parameters(((0, 1),), ((1, 2),))
  def test_log_probability_deterministic(self, seed):
    domain = domains.DiscreteDomain.from_vocab_sizes(range(2, 6))
    model = vae.VAE(domain)
    structures = domain.sample_uniformly(10)

    def _get_log_probs():
      return model.log_probability(structures, seed=seed)

    np.testing.assert_allclose(_get_log_probs(), _get_log_probs())

  @parameterized.parameters(
      (vae.summed_categorical_crossentropy,),
      (vae.mean_categorical_crossentropy,),
  )
  def test_loss_fn(self, loss_fn):
    mock_loss_fn = mock.Mock(wraps=loss_fn)
    mock_loss_fn.__name__ = 'foo'
    domain = domains.DiscreteDomain.from_vocab_sizes(range(2, 6))
    model = _build_vae(domain=domain, loss_fn=mock_loss_fn)
    model.fit(domain.sample_uniformly(10))
    mock_loss_fn.assert_called()


if __name__ == '__main__':
  absltest.main()
