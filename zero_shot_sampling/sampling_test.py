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
from absl.testing import parameterized
import numpy as np

from zero_shot_sampling import alignment
from zero_shot_sampling import domains
from zero_shot_sampling import sampling
from zero_shot_sampling import vae


class TestSampling(parameterized.TestCase):

  def _get_domain(self):
    vocab = domains.Vocabulary(['A', 'B', 'C', '-'])
    return domains.DiscreteDomain.from_shared_vocab(length=4, vocab=vocab)

  def _get_vae(self, domain):
    return vae.VAE(domain, fc_layer_params=(5,), num_latents=3)

  def test_sample_int_variant_batch_from_vae_valid_samples(self):
    domain = self._get_domain()
    model = self._get_vae(domain)
    parent = domain.sample_uniformly(1)[0]
    samples = sampling._sample_int_variant_batch_from_vae(
        25, model, parent, np.random.RandomState(0))
    self.assertLen(samples, 25)
    self.assertTrue(np.all(domain.are_valid(samples)))

  @parameterized.parameters((0,), (1,))
  def test_sample_int_variant_batch_from_vae_deterministic(self, random_seed):
    domain = self._get_domain()
    model = self._get_vae(domain)
    parent = domain.sample_uniformly(1)[0]

    def _sample(seed):
      random_state = np.random.RandomState(seed)
      return sampling._sample_int_variant_batch_from_vae(
          25, model, parent, random_state)

    np.testing.assert_array_equal(_sample(random_seed), _sample(random_seed))
    # Assert that there is at least one difference.
    self.assertTrue(np.any(_sample(random_seed) != _sample(random_seed + 1)))

  def test_sample_variants_from_vae_valid_samples(self):

    homolog_config = alignment.get_homolog_config_from_aligned_parent('AaC')
    model = self._get_vae(homolog_config.match_state_domain)
    samples = sampling.sample_variants_from_vae(
        25, model, homolog_config, np.random.RandomState(0))
    int_samples = homolog_config.match_state_domain.encode(samples)
    self.assertTrue(
        np.all(homolog_config.match_state_domain.are_valid(int_samples))
    )

  @parameterized.parameters(
      (lambda sample: sample == 'A', 10),
      (lambda sample: sample in ['B', 'C'], 25),
  )
  def test_rejection_sample_feasibility(
      self,
      feasibility_fn,
      num_samples,
  ):
    def sample_fn(random_state):
      return random_state.choice(['A', 'B', 'C'])

    samples = sampling.rejection_sample(
        sample_fn=sample_fn,
        feasibility_fn=feasibility_fn,
        num_samples=num_samples,
        random_state=np.random.RandomState(0),
    )
    self.assertLen(samples, num_samples)
    self.assertTrue(all(feasibility_fn(s) for s in samples))

  def test_rejection_sample_raises_when_budget_exceeded(
      self,
  ):
    def sample_fn(random_state):
      return random_state.choice(['A', 'B', 'C'])

    feasibility_fn = lambda sample: sample == 'D'
    with self.assertRaisesRegex(RuntimeError, 'Only 0 samples were accepted'):
      _ = sampling.rejection_sample(
          sample_fn=sample_fn,
          feasibility_fn=feasibility_fn,
          num_samples=10,
          random_state=np.random.RandomState(0),
      )


if __name__ == '__main__':
  absltest.main()
