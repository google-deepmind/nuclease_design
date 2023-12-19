# Copyright 2023 DeepMind Technologies Limited
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

"""Defines amino acids."""


# Single letter, three letter, and full amino acid names.
AA_NAMES = (
    ('A', 'Ala', 'alanine'),
    ('R', 'Arg', 'arginine'),
    ('N', 'Asn', 'asparagine'),
    ('D', 'Asp', 'aspartic acid'),
    ('C', 'Cys', 'cysteine'),
    ('E', 'Glu', 'glutamic acid'),
    ('Q', 'Gln', 'glutamine'),
    ('G', 'Gly', 'glycine'),
    ('H', 'His', 'histidine'),
    ('I', 'Ile', 'isoleucine'),
    ('L', 'Leu', 'leucine'),
    ('K', 'Lys', 'lysine'),
    ('M', 'Met', 'methionine'),
    ('F', 'Phe', 'phenylalanine'),
    ('P', 'Pro', 'proline'),
    ('S', 'Ser', 'serine'),
    ('T', 'Thr', 'threonine'),
    ('W', 'Trp', 'tryptophan'),
    ('Y', 'Tyr', 'tyrosine'),
    ('V', 'Val', 'valine'),
    # Extended AAs
    ('B', 'Asx', 'asparagine or aspartic acid'),
    ('Z', 'Glx', 'glutamine or glutamic acid'),
    ('X', 'Xaa', 'Any'),
    ('J', 'Xle', 'Leucine or isoleucine'),
)

# Indices of standard amino acids in `AA_NAMES`.
STANDARD_INDICES = tuple(range(20))
# Indices of extended amino acids in `AA_NAMES`.
EXTENDED_INDICES = tuple(range(20, 24))

# Single letter codes of standard amino acids.
STANDARD_AA = tuple(AA_NAMES[i][0] for i in STANDARD_INDICES)
AA = STANDARD_AA


# Isoelectric points (higher pI ~ more positively charged)
# https://www.vanderbilt.edu/AnS/Chemistry/Rizzo/stuff/AA/AminoAcids.html
AA_TO_ISOELECTRIC_POINT = {
    'A': 6.11,   # Alanine
    'R': 10.76,  # Arginine
    'N': 5.43,   # Asparagine
    'D': 2.98,   # Aspartic Acid
    'C': 5.15,   # Cysteine
    'E': 3.08,   # Glutamic Acid
    'Q': 5.65,   # Glutamine
    'G': 6.06,   # Glycine
    'H': 7.64,   # Histidine
    'I': 6.04,   # Isoleucine
    'L': 6.04,   # Leucine
    'K': 9.47,   # Lysine
    'M': 5.71,   # Methionine
    'F': 5.76,   # Phenylalanine
    'P': 6.30,   # Proline
    'S': 5.70,   # Serine
    'T': 5.60,   # Threonine
    'W': 5.88,   # Tryptophan
    'Y': 5.63,   # Tyrosine
    'V': 6.02    # Valine
}
