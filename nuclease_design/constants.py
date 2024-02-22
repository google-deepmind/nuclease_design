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

"""Constants for data."""

from os import path

FULL_REFERENCE_SEQ = 'MIKKWAVHLLFSALVLLGLSGGAAYSPQHAEGAARYDDVLYFPASRYPETGAHISDAIKAGHADVCTIERSGADKRRQESLKGIPTKPGFDRDEWPMAMCEEGGKGASVRYVSSSDNRGAGSWVGNRLNGYADGTRILFIVQ'

GCS_DATA_DIR = 'https://storage.googleapis.com/nuclease_design'
DATA_DIR = GCS_DATA_DIR


G4_PATHS = (
    (
        'read_count_0_input_g4',
        'OB268_reprocessed.csv',
    ),
    (
        'read_count_1_70_g4',
        'OB269_reprocessed.csv',
    ),
    (
        'read_count_2_90_g4',
        'OB270_reprocessed.csv',
    ),
    (
        'read_count_3_98_g4',
        'OB271_reprocessed.csv',
    ),
    (
        'read_count_4_99.5_g4',
        'OB272_reprocessed.csv',
    ),
)


G3_PATHS = (
    (
        'read_count_0_input_g3',
        'ob222_processed_updated.csv',
    ),
    (
        'read_count_1_59_g3',
        'ob231_processed.csv',
    ),
    (
        'read_count_1_80_g3',
        'ob230_processed.csv',
    ),
    (
        'read_count_2_95_g3',
        'ob241_processed.csv',
    ),
    (
        'read_count_2_99_g3',
        'ob239_processed.csv',
    ),
)


G2_PATHS = (
    (
        'read_count_1_input_deep_g2',
        'ob188_processed.csv',
    ),
    (
        'read_count_1_input_g2',
        'ob188_ns1_processed.csv',
    ),
    (
        'read_count_1_86_g2',
        'ob193_processed.csv',
    ),
    (
        'read_count_1_97.5_g2',
        'ob194_processed.csv',
    ),
    (
        'read_count_2_input_g2',
        'ob188_ns2_processed.csv',
    ),
    (
        'read_count_2_93_g2',
        'ob196_processed.csv',
    ),
)

G1_PATHS = (
    (
        'read_count_1_high_g1',
        'Sort_1_high_gate_output_ob142.csv',
    ),
    (
        'read_count_1_input_g1',
        'Sort_1_input_ob133_3_1.csv',
    ),
    (
        'read_count_1_low_g1',
        'Sort_1_low_gate_output_ob141.csv',
    ),
    (
        'read_count_2_high_g1',
        'Sort_2_high_gate_output_ob145.csv',
    ),
    (
        'read_count_2_input_g1',
        'Sort_2_input_ob133_3_2.csv',
    ),
    (
        'read_count_2_low_g1',
        'Sort_2_low_gate_ob146.csv',
    ),
    (
        'read_count_3_high_g1',
        'Sort_3_high_gate_ob149.csv',
    ),
    (
        'read_count_3_input_g1',
        'Sort_3_input_ob133_3_3.csv',
    ),
    (
        'read_count_3_low_g1',
        'Sort_3_low_gate_ob150.csv',
    ),
    (
        'read_count_1_input_reseq_g1',
        'Sort_1_input_ob133_3_1_reseq.csv',
    ),
    (
        'read_count_2_high_reseq_g1',
        'Sort_2_high_gate_ob145_reseq.csv',
    ),
    (
        'read_count_2_input_reseq_g1',
        'Sort_2_input_ob133_3_2_reseq.csv',
    ),
)


def _get_full_path_tuple(generation, paths):
  return tuple(
      (name, path.join('raw_count_data', generation, base_path))
      for name, base_path in paths
  )


G4_INFO = {
    'names_and_paths': _get_full_path_tuple('G4', G4_PATHS),
    'count_threshold': 50,
    'input_names': ['read_count_0_input_g4'],
    'post_sort_names': [
        'read_count_1_70_g4',
        'read_count_2_90_g4',
        'read_count_3_98_g4',
        'read_count_4_99.5_g4',
    ],
    'fiducials': [
        'neg_control',
        'wt',
        'a73r',
        'a73r_d74s',
        'a63p_a73r_d74h_i84y',
    ],
    'generation': 'g4',
}

G3_INFO = {
    'names_and_paths': _get_full_path_tuple('G3', G3_PATHS),
    'count_threshold': 50,
    'input_names': ['read_count_0_input_g3'],
    'post_sort_names': [
        'read_count_1_59_g3',
        'read_count_1_80_g3',
        'read_count_2_95_g3',
        'read_count_2_99_g3',
    ],
    'fiducials': ['neg_control', 'wt', 'a73r'],
    'generation': 'g3',
}

G2_INFO = {
    'names_and_paths': _get_full_path_tuple('G2', G2_PATHS),
    'count_threshold': 100,
    'input_names': [
        'read_count_1_input_deep_g2',
        'read_count_1_input_g2',
        'read_count_2_input_g2',
    ],
    'post_sort_names': [
        'read_count_1_86_g2',
        'read_count_2_93_g2',
        'read_count_1_97.5_g2',
    ],
    'fiducials': ['neg_control', 'wt'],
    'generation': 'g2',
}

G1_INFO = {
    'names_and_paths': _get_full_path_tuple('G1', G1_PATHS),
    'count_threshold': 20,
    'input_names': [
        'read_count_1_input_g1',
        'read_count_2_input_g1',
        'read_count_3_input_g1',
        'read_count_1_input_reseq_g1',
        'read_count_2_input_reseq_g1',
    ],
    'post_sort_names': [
        'read_count_1_high_g1',
        'read_count_2_high_g1',
        'read_count_3_high_g1',
        'read_count_2_high_reseq_g1',
    ] + ['read_count_1_low_g1', 'read_count_2_low_g1', 'read_count_3_low_g1'],
    'fiducials': ['neg_control', 'wt'],
    'generation': 'g1',
}

GENERATION_AND_REF_NAME_TO_GATE_NAME = {
    'g1': {
        'neg_control': '2_high_reseq',
        'wt': '2_high_reseq',
    },
    'g2': {
        'neg_control': '1_86',
        'wt': '2_93',
    },
    'g3': {
        'neg_control': '1_59',
        'wt': '1_80',
        'a73r': '1_80',
    },
    'g4': {
        'neg_control': '1_70',
        'wt': '2_90',
        'a73r': '4_99.5',
        'a73r_d74s': '4_99.5',
        'a63p_a73r_d74h_i84y': '4_99.5',
    },
}

FIDUCIAL_NAME_TO_MUTATION_TUPLE = {
    'wt': (),
    'a73r': (('A', 73, 'R'),),
    'a73r_d74s': (('A', 73, 'R'), ('D', 74, 'S')),
    'a63p_a73r_d74h_i84y': (
        ('A', 63, 'P'),
        ('A', 73, 'R'),
        ('D', 74, 'H'),
        ('I', 84, 'Y'),
    ),
}


ALIGNMENT_PATH = path.join('alignments', 'F1BV52_BACLI_e3_n8_mFirst.a2m')

# N-terminal one indexed 33 110 (inclusive)
G4_N_DESIGN_SEQUENCE_START = 32  # (zero-indexed, inclusive)
G4_N_DESIGN_SEQUENCE_END = 110  # (zero-indexed, exclusive)

# C-terminal one indexed 65 142 (inclusive)
G4_C_DESIGN_SEQUENCE_START = 64  # (zero-indexed, inclusive)
G4_C_DESIGN_SEQUENCE_END = 142  # (zero-indexed, exclusive)

_DESIGN_REGION_INDICES = [
    G4_C_DESIGN_SEQUENCE_START,
    G4_C_DESIGN_SEQUENCE_END,
    G4_N_DESIGN_SEQUENCE_START,
    G4_N_DESIGN_SEQUENCE_END,
]

G4_VARIABLE_REGION_START = min(_DESIGN_REGION_INDICES)
G4_VARIABLE_REGION_END = max(_DESIGN_REGION_INDICES) + 1


def get_processed_data_file(generation: str) -> str:
  return path.join('processed_data', f'{generation}.csv')


def get_fiducial_filename(
    fiducial: str, generation: str) -> str:
  return path.join(
      'processed_fiducial_data',
      f'{generation}_{fiducial}.csv',
  )


NUM_CLUSTERS_G4_PATH = path.join('analysis', 'num_clusters_g4_df.csv')

NUM_CLUSTERS_ZERO_SHOT_PATH = path.join(
    'analysis', 'num_clusters_zero_shot_df.csv'
)

ZERO_SHOT_SCORE_CSV = path.join('analysis', 'zero_shot_scores.csv')


MUTATIONS_TO_SUBLIBRARY_PATH = path.join(
    'library_designs', 'mutation_sublibrary_df.csv'
)


PLATE_GENOTYPES_PATH = path.join('plate_data', 'genotypes.csv')

TIME_SERIES_ACTIVITY_XML_PATH = path.join(
    'plate_data',
    '2021-09-07-purified-protein-comparison-activity.xml',
)

LANDSCAPE_PATH = path.join('processed_data', 'landscape.csv')


# Secondary structure coordinates taken from NucB crystal structure paper
# https://academic.oup.com/view-large/figure/107884442/gkx1170fig3.jpg
def get_secondary_structure(position: int) -> str | None:
  """Maps a one-indexed int position to its secondary structure element."""
  if position in range(38, 42 + 1):
    return 'beta1'
  elif position in range(47, 60 + 1):
    return 'alpha1'
  elif position in range(64, 68 + 1):
    return 'beta2'
  elif position in range(72, 82 + 1):
    return 'alpha2'
  elif position in range(89, 95 + 1):
    return 'beta3'
  elif position in range(108, 113 + 1):
    return 'beta4'
  elif position in range(114, 129 + 1):
    return 'alpha3'
  elif position in range(135, 140 + 1):
    return 'beta4'
  else:
    return None


# Inclusive indexes
INDEX_TO_CORE_BASE_LABEL = {
    (38, 42): 'core',
    (64, 48): 'core',
    (72, 88): 'base',
    (89, 113): 'core',
    (114, 129): 'base',
    (135, 140): 'core',
}


VARIABLE_REGION_START_INDEX = 24
VARIABLE_REGION_LENGTH = len(FULL_REFERENCE_SEQ) - VARIABLE_REGION_START_INDEX

LANDSCAPE_ACTIVITY_LEVELS = (
    'non-functional',
    'activity_greater_than_0',
    'activity_greater_than_WT',
    'activity_greater_than_A73R',
)


NON_CAMPAIGN_SUBLIBRARIES = (
    'g2_g1_plate_assay_variants',
    'g2_wt_synonyms',
    'g2_stratified_sample',
    'g3_wt_synonyms',
    'g3_g1_stratified_sample',
    'g3_g2_plate_assay_variants',
    'g3_a73r_synonyms',
    'g4_mbo_seeds',
    'g4_wt_synonyms',
    'g4_stratified_sample',
    'g4_g3_plate_assay_variants',
    'g4_a73r_synonyms',
    'g4_g3_hit_constituents',
    'g4_homolog_graft',
    'g4_double_synonyms',
    'g4_quad_synonyms',
    'prosar+screen_g2_redux',
)
