import argparse
import itertools
import os

import pandas as pd
from python_pdb.aligners import align_pandas_structure
from python_pdb.comparisons import rmsd
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_structure_tools.measurements import compute_residue_com, get_distance, get_distances, measure_chi_angle
from tcr_pmhc_structure_tools.processing import annotate_tcr_df
from tcr_pmhc_structure_tools.utils import get_coords

parser = argparse.ArgumentParser()

parser.add_argument('input', help='path to data directory')
parser.add_argument('--output', '-o', help='path to output file')


def main():
    args = parser.parse_args()

    complexes = [complex_id for complex_id in os.listdir(args.input)
                 if os.path.isdir(os.path.join(complex_id, args.input))]
    num_complexes = len(complexes)

    # Info
    complex_ids = []
    apo_ranks = []
    holo_ranks = []
    cdrs = []
    chain_types = []
    residue_names = []
    residue_seq_ids = []
    residue_insert_codes = []

    # Measurements
    ca_distances = []
    rmsds = []
    chi_angle_changes = []
    com_distances = []

    for num, complex_id in enumerate(complexes, 1):
        print(complex_id, '-', num, 'of', num_complexes, flush=True)

        complex_path = os.path.join(args.input, complex_id)

        for rank_dir in os.listdir(complex_path):
            print(rank_dir)

            _, holo_rank, _, apo_rank = rank_dir.split('_')

            rank_path = os.path.join(complex_path, rank_dir)
            tcr_pmhc_path = os.path.join(rank_path, 'tcr_pmhc.pdb')
            tcr_path = os.path.join(rank_path, 'tcr.pdb')

            structures = []
            for path in tcr_path, tcr_pmhc_path:
                with open(path, 'r') as fh:
                    structure_df = parse_pdb_to_pandas(fh.read())

                structure_df = annotate_tcr_df(structure_df, alpha_chain_id='D', beta_chain_id='E')
                structure_df['backbone'] = structure_df['atom_name'].map(
                    lambda atom_name: (atom_name == 'N' or atom_name == 'CA' or atom_name == 'C' or atom_name == 'O')
                )

                structures.append(structure_df)

            tcr_df, tcr_pmhc_df = structures

            # Query desired residues/atoms
            for chain_type, cdr_number in itertools.product(['alpha_chain', 'beta_chain'], [1, 2, 3]):
                apo_cdr = tcr_df.query("cdr == @cdr_number and chain_type == @chain_type")
                holo_cdr = tcr_pmhc_df.query("cdr == @cdr_number and chain_type == @chain_type")

                apo_backbone_cdr_coords = get_coords(apo_cdr.query('backbone'))
                holo_backbone_cdr_coords = get_coords(holo_cdr.query('backbone'))

                holo_cdr = align_pandas_structure(holo_backbone_cdr_coords, apo_backbone_cdr_coords, holo_cdr)

                # Compute CA distance
                distance = get_distances(get_coords(apo_cdr.query("atom_name == 'CA'")),
                                         get_coords(holo_cdr.query("atom_name == 'CA'")))

                cdr_rmsds = []
                cdr_angle_changes = []
                cdr_com_distances = []

                cdr_residue_names = []
                cdr_seq_ids = []
                cdr_insert_codes = []

                group1 = apo_cdr.groupby(['residue_name', 'residue_seq_id', 'residue_insert_code'], dropna=False)
                group2 = holo_cdr.groupby(['residue_name', 'residue_seq_id', 'residue_insert_code'], dropna=False)

                for ((res_name, seq_id, insert_code), res1), (_, res2) in zip(group1, group2):
                    # Compute Residue RMSD
                    try:
                        cdr_rmsds.append(rmsd(get_coords(res1), get_coords(res2)))
                    except ValueError:
                        print('Mismatched number of atoms in residue:')
                        print(apo_rank, holo_rank, res_name, seq_id, insert_code)
                        cdr_rmsds.append(None)

                    # Compute chi angle changes
                    if res_name == 'GLY' or res_name == 'ALA':
                        cdr_angle_changes.append(None)
                    else:
                        try:
                            cdr_angle_changes.append(measure_chi_angle(res1) - measure_chi_angle(res2))
                        except IndexError:
                            print('Missing atoms needed to calculate chi angle:')
                            print(apo_rank, holo_rank, res_name, seq_id, insert_code)
                            cdr_angle_changes.append(None)

                    # Compute C.O.M changes
                    cdr_com_distances.append(get_distance(compute_residue_com(res1), compute_residue_com(res2)))

                    cdr_residue_names.append(res_name)
                    cdr_seq_ids.append(seq_id)
                    cdr_insert_codes.append(insert_code)

                # append values to lists
                num_residues = len(cdr_residue_names)

                residue_names += cdr_residue_names
                residue_seq_ids += cdr_seq_ids
                residue_insert_codes += cdr_insert_codes

                ca_distances += list(distance)
                rmsds += cdr_rmsds
                chi_angle_changes += cdr_angle_changes
                com_distances += cdr_com_distances

                complex_ids += [complex_id] * num_residues
                apo_ranks += [apo_rank] * num_residues
                holo_ranks += [holo_rank] * num_residues

                cdrs += [cdr_number] * num_residues
                chain_types += [chain_type] * num_residues

    pd.DataFrame({
        'complex_id': complex_ids,
        'apo_rank': apo_ranks,
        'holo_rank': holo_ranks,
        'cdr': cdrs,
        'chain_type': chain_types,
        'residue_name': residue_names,
        'residue_seq_id': residue_seq_ids,
        'residue_insert_code': residue_insert_codes,
        'ca_distance': ca_distances,
        'rmsd': rmsds,
        'chi_angle_change': chi_angle_changes,
        'com_distance': com_distances,
    }).to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
