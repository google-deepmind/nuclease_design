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
from zero_shot_sampling import alignment
from zero_shot_sampling import amino_acids
from zero_shot_sampling import domains


class AlignmentTest(parameterized.TestCase):

  @parameterized.parameters(
      (
          'A-C',
          'A-C',
          dict(
              insertion=0,
              deletion=0,
              substitution=0,
          ),
      ),
      (
          'A-C',
          'ACC',
          dict(
              insertion=0,
              deletion=1,
              substitution=0,
          ),
      ),
      (
          'ACC',
          'A--',
          dict(
              insertion=2,
              deletion=0,
              substitution=0,
          ),
      ),
      (
          'B-C',
          'AC-',
          dict(
              insertion=1,
              deletion=1,
              substitution=1,
          ),
      ),
  )
  def test_get_edit_counts(
      self,
      aligned_sequence: str,
      aligned_parent_sequence: str,
      expected_output: dict[str, int],
  ):
    actual_output = alignment.get_edit_counts(
        aligned_sequence, aligned_parent_sequence
    )
    self.assertDictEqual(actual_output, expected_output)

  def test_get_edit_counts_mismatched_lengths_raises(self):
    with self.assertRaisesRegex(ValueError, 'must have the same length'):
      alignment.get_edit_counts('A-C', 'AC')

  @parameterized.parameters(
      ('A-C', 'A-C', ()),
      ('A-C', 'ACC', (('C', 1, '-'),)),
      ('A--', 'AB-', (('B', 1, '-'),)),
      ('A-C', 'CC-', (('C', 0, 'A'), ('C', 1, '-'), ('-', 2, 'C'))),
  )
  def test_get_edits(self, query_seq, parent_seq, expected_output):
    actual_output = alignment.get_edits(query_seq, parent_seq)
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(
      ((('A', 0, 'B'),), 'A-C', 'B-C'),
      ((('A', 0, 'B'), ('-', 1, 'D')), 'A-C', 'BDC')
      )
  def test_apply_edits(self, edits, parent_seq, expected_output):
    actual_output = alignment.apply_edits(edits, parent_seq)
    self.assertEqual(actual_output, expected_output)

  def test_apply_edits_duplicate_positions_raises(self):
    with self.assertRaisesRegex(ValueError, 'duplicate positions'):
      alignment.apply_edits(
          (('A', 0, 'B'), ('A', 0, 'C')),
          'A-C',
      )

  @parameterized.parameters(
      dict(
          fasta_text='>i1\nA-C\n>i2\nACC',
          expected_parent='A-C',
          expected_homologs=('ACC',),
      ),
      dict(
          fasta_text='>i1\nA-C.\n>i2\nACCt\n>i3\nA-Cc',
          expected_parent='A-C.',
          expected_homologs=('ACCt', 'A-Cc'),
      ),
  )
  def test_load_msa_from_fasta(
      self, fasta_text, expected_parent, expected_homologs
  ):
    fasta_file = self.create_tempfile(content=fasta_text).full_path
    msa = alignment.load_msa_from_fasta(fasta_file=fasta_file)

    expected_msa = alignment.MSA(
        parent=expected_parent, homologs=expected_homologs
    )
    self.assertEqual(msa, expected_msa)

  @parameterized.parameters(
      dict(
          msa=alignment.MSA(parent='A-C', homologs=('ACC',)),
          expected_output=alignment.MSA(parent='A.C', homologs=('AcC',)),
      ),
      dict(
          msa=alignment.MSA(parent='a-C', homologs=('ACC',)),
          expected_output=alignment.MSA(parent='A.C', homologs=('AcC',)),
      ),
  )
  def test_convert_parent_residues_to_match_states(self, msa, expected_output):
    actual_output = alignment.convert_parent_residues_to_match_states(msa=msa)
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(('A', False),
                            ('M', False),
                            ('-', True),
                            ('.', True),
                            ('a', False))
  def test_is_gap(self, token, expected_output):
    actual_output = alignment.is_gap(token)
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(('A', True),
                            ('M', True),
                            ('-', True),
                            ('.', False),
                            ('a', False))
  def test_is_match(self, token, expected_output):
    actual_output = alignment.is_match(token)
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(('A-C', [0, 1, 2]),
                            ('aA-C', [1, 2, 3]),
                            ('aA-C.A', [1, 2, 3, 5]),
                            ('a.c', []))
  def test_get_match_indices(self, aligned_seq, expected_output):
    actual_output = alignment._get_match_indices(aligned_seq)
    np.testing.assert_array_equal(actual_output, expected_output)

  @parameterized.parameters(('A-C', 'AC'),
                            ('aA-C', 'AAC'),
                            ('aA-C.A', 'AACA'),
                            ('a.c', 'AC'))
  def test_get_unaligned_seq(self, aligned_seq, expected_output):
    actual_output = alignment.get_unaligned_seq(aligned_seq)
    np.testing.assert_array_equal(actual_output, expected_output)

  @parameterized.parameters(
      ('A-C', [0, -1, 1]),
      ('aA-C', [0, 1, -1, 2]),
      ('.aA-C.A', [-1, 0, 1, -1, 2, -1, 3]),
      ('.aA-C.A', [0, 0, 1, 0, 2, 0, 3], 0),
  )
  def test_scatter_non_gap_indices(
      self, aligned_seq, expected_output, fill_value=-1
  ):
    actual_output = alignment._scatter_non_gap_indices(aligned_seq, fill_value)
    np.testing.assert_array_equal(actual_output, expected_output)

  @parameterized.parameters(
      (
          'ACC',
          0,
          3,
          (0, 3),
      ),
      (
          'A-C',
          0,
          2,
          (0, 3),
      ),
      (
          'A-C',
          0,
          1,
          (0, 2),
      ),
      ('K-A--CS', 1, 3, (2, 6)),
      ('ACC-K', 0, 4, (0, 5))
  )
  def test_get_aligned_start_and_end_indices(
      self,
      aligned_parent,
      start_index_in_unaligned_parent,
      end_index_in_unaligned_parent,
      expected_output,
  ):
    actual_output = alignment.get_aligned_start_and_end_indices(
        aligned_parent,
        start_index_in_unaligned_parent,
        end_index_in_unaligned_parent,
    )
    self.assertEqual(actual_output, expected_output)

  @parameterized.parameters(
      (
          dict(
              aligned_parent_seq='ACC-K',
              max_num_insertions=1,
              max_num_deletions=0,
              max_num_substitutions=2,
              design_region_start=0,
              design_region_end=4,
          ),
          ['A-CCK', 'ATC-K', 'A--TK', 'AKSTK', 'AC-CK'],
          [False, True, False, True, False],
      ),
      (
          dict(
              aligned_parent_seq='ACC-K',
              max_num_insertions=1,
              max_num_deletions=0,
              max_num_substitutions=1,
              design_region_start=0,
              design_region_end=4,
          ),
          ['A-CCK', 'ATC-K', 'A--TK', 'AKSTK', 'AC-CK'],
          [False, True, False, False, False],
      ),
      (
          dict(
              aligned_parent_seq='ACC-K',
              max_num_insertions=1,
              max_num_deletions=0,
              max_num_substitutions=1,
              design_region_start=2,
              design_region_end=4,
          ),
          ['A-CCK', 'ATC-K', 'A--TK', 'AKSTK', 'AC-CK'],
          [True, True, False, True, False],
      )
  )
  def test_get_design_region_feasibility_fn(
      self, feasibility_fn_kwargs, input_seqs, expected_output
  ):
    feasibility_fn = alignment.get_design_region_feasibility_fn(
        **feasibility_fn_kwargs
    )
    actual_output = [feasibility_fn(seq) for seq in input_seqs]
    self.assertEqual(tuple(actual_output), tuple(expected_output))

  @parameterized.parameters(
      (
          'ACC-K',
          ['A-CCK', 'ATC-K', 'A--TK', 'AKSTK', 'AC-CK'],
          0,
          3,
          ['ACCK', 'ATCK', 'ATK', 'AKSTK', 'ACCK'],
      ),
      (
          'ACC-K',
          ['A-CCK', 'ATC-K', 'A--TK', 'A-C-T'],
          2,
          4,
          ['ACCCK', 'ACCK', 'ACTK', 'ACCT'],
      ),
      (
          'ACC-K',
          ['A-CCK', 'ATC-K', 'A--TS', 'A-C-T'],
          2,
          3,
          ['ACCCK', 'ACCK', 'ACTK', 'ACCK'],
      ),
  )
  def test_get_spliced_sequence_df(
      self,
      aligned_parent_seq,
      aligned_seqs,
      design_region_start,
      design_region_end,
      expected_output
  ):
    actual_output_df = alignment.get_spliced_sequence_df(
        aligned_seqs,
        aligned_parent_seq,
        design_region_start,
        design_region_end,
    )
    np.testing.assert_array_equal(actual_output_df['sequence'], expected_output)


class MatchStateSelectorTest(parameterized.TestCase):

  @parameterized.parameters(('A-C', 3),
                            ('aA-C', 3),
                            ('aA-C.A', 4))
  def test_num_match_positions(self, aligned_parent, expected_output):
    match_state_selector = alignment.MatchStateSelector(aligned_parent)
    np.testing.assert_array_equal(
        match_state_selector.num_match_positions, expected_output
    )

  @parameterized.parameters(
      ('AaC', ['-BC', '---', 'ABC'], ['-C', '--', 'AC']),
      ('aA-C.A', ['.A-BcC', 'mB-B.C'], ['A-BC', 'B-BC']),
  )
  def test_select_from_aligned_seqs_aa_inputs(self, aligned_parent, input_seqs,
                                              expected_output_strings):
    match_state_selector = alignment.MatchStateSelector(aligned_parent)
    actual_outputs = match_state_selector.select_from_aligned_seqs(input_seqs)
    expected_outputs = [list(s) for s in expected_output_strings]
    np.testing.assert_array_equal(actual_outputs, expected_outputs)

  @parameterized.parameters(
      ('AaC', [[0, 1, 2], [3, 3, 3], [2, 1, 0]], [[0, 2], [3, 3], [2, 0]]))
  def test_select_from_aligned_seqs_int_inputs(self, aligned_parent, input_seqs,
                                               expected_outputs):
    match_state_selector = alignment.MatchStateSelector(aligned_parent)
    actual_outputs = match_state_selector.select_from_aligned_seqs(input_seqs)
    np.testing.assert_array_equal(actual_outputs, expected_outputs)

  @parameterized.parameters(('AaC', ['ABC', 'ACD'], ['AC', 'AD']),
                            ('aAAC', ['ABCD', 'DCBA'], ['BCD', 'CBA']),
                            ('aA-.A', ['ABCDE', 'EDCBA'], ['BCE', 'DCA']))
  def test_select_from_unaligned_seqs_aa_inputs(self, aligned_parent,
                                                input_seqs,
                                                expected_output_strings):
    match_state_selector = alignment.MatchStateSelector(aligned_parent)
    actual_outputs = match_state_selector.select_from_aligned_seqs(input_seqs)
    expected_outputs = [list(s) for s in expected_output_strings]
    np.testing.assert_array_equal(actual_outputs, expected_outputs)

  @parameterized.parameters(
      ('AaC', [[0, 1, 2], [2, 1, 0]], [[0, 2], [2, 0]]),)
  def test_select_from_unaligned_seqs_int_inputs(self, aligned_parent,
                                                 input_seqs, expected_outputs):
    match_state_selector = alignment.MatchStateSelector(aligned_parent)
    actual_outputs = match_state_selector.select_from_unaligned_seqs(input_seqs)
    np.testing.assert_array_equal(actual_outputs, expected_outputs)

  def test_no_match_states_raises(self):
    with self.assertRaisesRegex(ValueError, 'match states'):
      alignment.MatchStateSelector('a.c')


class MatchStateEncoderTest(parameterized.TestCase):

  def _test_call(self, seqs, encode_fn, aligned_reference, batch_size):
    match_state_selector = alignment.MatchStateSelector(aligned_reference)
    mock_encode_fn = mock.Mock(wraps=encode_fn)
    encoder = alignment.MatchStateEncoder(
        mock_encode_fn, match_state_selector, batch_size=batch_size)
    encoder(seqs)
    self.assertEqual(mock_encode_fn.call_count,
                     np.ceil(len(seqs) / batch_size))
    for args, _ in mock_encode_fn.call_args_list:
      np.testing.assert_array_equal(args[0].shape, (batch_size, 2))

  @parameterized.parameters((10, 3,),
                            (12, 5))
  def test_call_int_inputs(self, num_inputs, batch_size):
    encode_fn = lambda seqs: np.sum(seqs, axis=1)
    seqs = np.random.uniform(size=(num_inputs, 3))
    self._test_call(seqs, encode_fn, 'AaC', batch_size)

  @parameterized.parameters((10, 3,),
                            (12, 5))
  def test_call_string_inputs(self, num_inputs, batch_size):
    encode_fn = lambda seqs: np.float32([seq[0] == 'A' for seq in seqs])
    domain = domains.DiscreteDomain.from_shared_vocab(
        length=3, vocab=domains.Vocabulary(tokens=['A', 'B', 'C'])
    )
    seqs = domain.decode(domain.sample_uniformly(num_inputs))
    self._test_call(seqs, encode_fn, 'AaC', batch_size)


class HomologConfigTest(parameterized.TestCase):

  def _validate_homolog_config(self, homolog_config,
                               expected_match_state_parent, expected_tokens):
    self.assertLen(expected_match_state_parent,
                   homolog_config.match_state_selector.num_match_positions)
    self.assertLen(homolog_config.match_state_domain,
                   len(expected_match_state_parent))
    for vocab in homolog_config.match_state_domain.vocabs:
      np.testing.assert_array_equal(vocab.tokens, expected_tokens)
    self.assertEqual(
        homolog_config.match_state_parent, expected_match_state_parent)

  @parameterized.parameters(('AaC', 'AC'), ('Aa-', 'A-'),
                            ('aAAC.', 'AAC', ['C', 'A', 'B']))
  def test_get_homolog_config_from_aligned_parent(
      self,
      input_seq,
      expected_match_state_parent,
      tokens=amino_acids.AA_GAP):
    homolog_config = alignment.get_homolog_config_from_aligned_parent(
        input_seq, tokens=tokens)
    self._validate_homolog_config(homolog_config,
                                  expected_match_state_parent, tokens)

  @parameterized.parameters(('AaC', 'AC'), ('Aa-', 'A-'),
                            ('aAAC.', 'AAC', ['C', 'A', 'B']))
  def test_get_homolog_config_from_a2m_file(self,
                                            input_seq,
                                            expected_match_state_parent,
                                            tokens=amino_acids.AA_GAP):
    fasta_str = '>header\n%s' % input_seq
    temp_file = self.create_tempfile(content=fasta_str).full_path
    homolog_config = alignment.get_homolog_config_from_a2m_file(
        temp_file, tokens)
    self._validate_homolog_config(homolog_config,
                                  expected_match_state_parent, tokens)

  def test_save_and_load_homolog_config(self):
    full_parent = 'AaC-Kt'
    homolog_config = alignment.get_homolog_config_from_aligned_parent(
        full_parent
    )
    temp_file = self.create_tempfile().full_path
    alignment.save_homolog_config(homolog_config, temp_file)
    loaded_homolog_config = alignment.load_homolog_config(temp_file)
    self.assertEqual(homolog_config.aligned_parent,
                     loaded_homolog_config.aligned_parent)
    self.assertEqual(homolog_config.aligned_parent, full_parent)


if __name__ == '__main__':
  absltest.main()
