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


def calculate_angle(v1, v2):
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)

    cos_theta = np.dot(v1, v2) / (n1 * n2)

    cos_theta = min(cos_theta, 1)
    cos_theta = max(cos_theta, -1)

    return np.arccos(cos_theta)


def calculate_dihedral_angle(a, b, c, d):
    ba = a - b
    bc = c - b
    cd = d - c

    u = np.cross(ba, bc)
    v = np.cross(cd, bc)

    w = np.cross(u, v)

    angle = calculate_angle(u, v)

    if calculate_angle(bc, w) > 0.001:
        angle = -angle

    return angle


def calculate_phi_psi_angles(residue: pd.DataFrame,
                             prev_residue: pd.DataFrame,
                             next_residue: pd.DataFrame) -> tuple[float, float]:
    '''Calculate the dihedral (phi and psi) angles for a given residue.'''
    c_prev_pos = prev_residue.query("atom_name == 'C'")[['pos_x', 'pos_y', 'pos_z']].iloc[0]

    n_pos = residue.query("atom_name == 'N'")[['pos_x', 'pos_y', 'pos_z']].iloc[0]
    ca_pos = residue.query("atom_name == 'CA'")[['pos_x', 'pos_y', 'pos_z']].iloc[0]
    c_pos = residue.query("atom_name == 'C'")[['pos_x', 'pos_y', 'pos_z']].iloc[0]

    n_next_pos = next_residue.query("atom_name == 'N'")[['pos_x', 'pos_y', 'pos_z']].iloc[0]

    phi_angle = calculate_dihedral_angle(c_prev_pos, n_pos, ca_pos, c_pos)
    psi_angle = calculate_dihedral_angle(n_pos, ca_pos, c_pos, n_next_pos)

    return (phi_angle, psi_angle)
