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

"""Tests for preprocessing.py."""

from absl.testing import absltest
from absl.testing import parameterized
import pandas as pd

from nuclease_design import preprocessing_utils


_MOCK_DATA_DIR = ""


def _get_mock_csvs():
  input_1_csv_string = ("""nuc_mut_string,mut_string,counts,num_AA_mut
      A1C,,30,0
      A2C,,10,0
      A5C,"A1T, V2W",30,2
      A3C,"A1*",10,1
      A4C,"*2A",10,1
      A6C,"A5B",5,1""")

  input_2_csv_string = ("""nuc_mut_string,mut_string,counts,num_AA_mut
      A1C,,10,0
      A11C,,5,0
      A12C,"A1V",5,1""")

  sort_1_csv_string = ("""nuc_mut_string,mut_string,counts,num_AA_mut
      A1C,,1,0
      A1C,,1,0""")

  sort_2_csv_string = ("""nuc_mut_string,mut_string,counts,num_AA_mut
      A1C,,10,0
      A2C,,10,0
      A5C,"A1T, V2W",10,2
      A3C,"A1*",10,1
      A4C,"*2A",10,1""")

  name_to_csv_str = {
      "input1": input_1_csv_string,
      "input2": input_2_csv_string,
      "sort1": sort_1_csv_string,
      "sort2": sort_2_csv_string,
  }
  return name_to_csv_str


class EnrichmentFactorTest(parameterized.TestCase):

  def _write_mock_csvs(self):
    name_to_csv_str = _get_mock_csvs()
    name_to_filepath = {}
    for name, csv_string in name_to_csv_str.items():
      filepath = self.create_tempfile().full_path
      with open(filepath, "w") as f:
        f.write(csv_string)
      name_to_filepath[name] = filepath
    return name_to_filepath

  @parameterized.named_parameters(
      dict(
          testcase_name="keep stops",
          count_threshold=8,
          group_aa_seqs=True,
          expected_output=pd.DataFrame(
              dict(
                  mutations=[
                      (),
                      (("*", 2, "A"),),
                      (("A", 1, "*"),),
                      (("A", 1, "T"), ("V", 2, "W")),
                  ],
                  input1=[40.0, 10.0, 10.0, 30.0],
                  input2=[10.0, 0.0, 0.0, 0.0],
                  sort1=[2.0, 0.0, 0.0, 0.0],
                  sort2=[20.0, 10.0, 10.0, 10.0],
                  # input_abundance=[1/2, 1/10, 1/10, 3/10],
                  ef_sort1=[2.0, 0.0, 0.0, 0.0],
                  ef_sort2=[4 / 5, 2.0, 2.0, 2 / 3],
              )
          ),
      ),
      dict(
          testcase_name="lower threshold",
          count_threshold=6,
          group_aa_seqs=True,
          expected_output=pd.DataFrame(
              dict(
                  mutations=[
                      (),
                      (("*", 2, "A"),),
                      (("A", 1, "*"),),
                      (("A", 1, "T"), ("V", 2, "W")),
                  ],
                  input1=[40.0, 10.0, 10.0, 30.0],
                  input2=[10.0, 0.0, 0.0, 0.0],
                  sort1=[2.0, 0.0, 0.0, 0.0],
                  sort2=[20.0, 10.0, 10.0, 10.0],
                  # input_abundance=[1/2, 1/10, 1/10, 3/10],
                  ef_sort1=[2.0, 0.0, 0.0, 0.0],
                  ef_sort2=[4 / 5, 2.0, 2.0, 2 / 3],
              )
          ),
      ),
      dict(
          testcase_name="dna variants",
          count_threshold=10,
          group_aa_seqs=False,
          expected_output=pd.DataFrame(
              dict(
                  nuc_mutations=[
                      (("A", 1, "C"),),
                      (("A", 2, "C"),),
                      (("A", 3, "C"),),
                      (("A", 4, "C"),),
                      (("A", 5, "C"),),
                  ],
                  mutations=[
                      (),
                      (),
                      (("A", 1, "*"),),
                      (("*", 2, "A"),),
                      (("A", 1, "T"), ("V", 2, "W")),
                  ],
                  input1=[30.0, 10.0, 10.0, 10.0, 30.0],
                  input2=[10.0, 0.0, 0.0, 0.0, 0.0],
                  sort1=[2.0, 0.0, 0.0, 0.0, 0.0],
                  sort2=[10.0, 10.0, 10.0, 10.0, 10.0],
                  # input_abundance=[4/10, 1/10, 1/10, 1/10, 3/10],
                  ef_sort1=[10 / 4, 0.0, 0.0, 0.0, 0.0],
                  ef_sort2=[1 / 2, 2.0, 2.0, 2.0, 2 / 3],
              )
          ),
      ),
  )
  def test_compute_ef(self, expected_output, count_threshold, group_aa_seqs):
    name_to_filename = self._write_mock_csvs()

    outer_join_df = preprocessing_utils.load_ngs_counts_from_paths(
        tuple(name_to_filename.items()), data_dir=_MOCK_DATA_DIR)
    processed_df = preprocessing_utils.get_enrichment_factor_df(
        outer_join_df,
        ["input1", "input2"],
        ["sort1", "sort2"],
        count_threshold=count_threshold,
        group_aa_seqs=group_aa_seqs,
    )

    pd.testing.assert_frame_equal(
        processed_df, expected_output, check_like=True)

  def test_load_ngs_counts_from_paths(self):
    name_to_filename = self._write_mock_csvs()
    outer_join_df = preprocessing_utils.load_ngs_counts_from_paths(
        tuple(name_to_filename.items()), data_dir=_MOCK_DATA_DIR)

    expected_keys = ("mutations", "nuc_mutations", "input1", "input2", "sort1",
                     "sort2")
    output_keys = outer_join_df.keys()
    for key in expected_keys:
      self.assertIn(key, set(output_keys))

  def test_get_stop_codon_df(self):
    name_to_filename = self._write_mock_csvs()
    outer_join_df = preprocessing_utils.load_ngs_counts_from_paths(
        tuple(name_to_filename.items()), data_dir=_MOCK_DATA_DIR)

    df = preprocessing_utils.get_stop_codon_df(
        outer_join_df,
        ["input1", "input2"],
        ["sort1", "sort2"],
        count_threshold=10,
    )
    expected_mutations = ((("A", 1, "*"),),)
    self.assertCountEqual(df.mutations, expected_mutations)

  @parameterized.named_parameters(
      dict(
          testcase_name="wt",
          mutations=(),
          expected_nuc_mutations=((("A", 1, "C"),), (("A", 2, "C"),))),
      dict(
          testcase_name="double",
          mutations=(("A", 1, "T"), ("V", 2, "W")),
          expected_nuc_mutations=((("A", 5, "C"),),))
  )
  def test_get_synonym_df(self, mutations, expected_nuc_mutations):
    name_to_filename = self._write_mock_csvs()
    outer_join_df = preprocessing_utils.load_ngs_counts_from_paths(
        tuple(name_to_filename.items()), data_dir=_MOCK_DATA_DIR)
    df = preprocessing_utils.get_synonym_df(
        outer_join_df,
        mutations,
        ["input1", "input2"],
        ["sort1", "sort2"],
        count_threshold=10,
    )
    self.assertCountEqual(df.nuc_mutations, expected_nuc_mutations)


if __name__ == "__main__":
  absltest.main()
