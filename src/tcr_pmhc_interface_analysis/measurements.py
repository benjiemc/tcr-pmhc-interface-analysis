import numpy as np
import pandas as pd
from python_pdb.chemistry import MOLECULAR_WEIGHTS

from tcr_pmhc_interface_analysis.chi_atoms import CHI_ATOMS
from tcr_pmhc_interface_analysis.utils import get_coords


def get_distances(vec1, vec2):
    xx = np.square(vec1[:, 0] - vec2[:, 0])
    yy = np.square(vec1[:, 1] - vec2[:, 1])
    zz = np.square(vec1[:, 2] - vec2[:, 2])

    return np.sqrt(xx + yy + zz)


def get_distance(vec1, vec2):
    xx = np.square(vec1[0] - vec2[0])
    yy = np.square(vec1[1] - vec2[1])
    zz = np.square(vec1[2] - vec2[2])

    return np.sqrt(xx + yy + zz)


def compute_residue_com(residue: pd.DataFrame) -> np.array:
    '''Compute the centre of mass for a residue.'''
    weights = residue['element'].map(MOLECULAR_WEIGHTS)
    coords = get_coords(residue)

    return np.average(coords, weights=weights, axis=0)


def measure_chi_angle(residue_df: pd.DataFrame, number: int = 1) -> float:
    res_name = residue_df.iloc[0]['residue_name']
    residue_chi_atoms = CHI_ATOMS[res_name][number]

    atom_positions = [residue_df.query("atom_name == @atom")[['pos_x', 'pos_y', 'pos_z']].iloc[0].to_numpy()
                      for atom in residue_chi_atoms]

    b1 = atom_positions[1] - atom_positions[0]
    b2 = atom_positions[2] - atom_positions[1]
    b3 = atom_positions[3] - atom_positions[1]

    b1_norm = b1 / np.linalg.norm(b1)
    b2_norm = b2 / np.linalg.norm(b2)
    b3_norm = b3 / np.linalg.norm(b3)

    n1 = np.cross(b1_norm, b2_norm)
    n2 = np.cross(b2_norm, b3_norm)

    m1 = np.cross(n1, b2_norm)

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return np.arctan2(x, y)
