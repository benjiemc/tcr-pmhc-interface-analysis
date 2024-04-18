import argparse
import os

import numpy as np
import pandas as pd
from dtaidistance.dtw_ndim import distance_fast
from python_pdb.aligners import align_pandas_structure
from python_pdb.parsers import parse_pdb_to_pandas
from tcr_structure_tools.cdr_numbering import assign_cdr_number
from tcr_structure_tools.utils import get_coords

parser = argparse.ArgumentParser(description='Create voxel representation of TCR-CDR loops')

parser.add_argument('--stcrdab-path',
                    default='/project/koohylab/shared/tcr_data/raw_DONOTMODIFY/structure/STCRDab_all_2022-11-10')
parser.add_argument('--output-dir', required=True)
parser.add_argument('--number-of-anchors', type=int, default=5, help='number of anchors to include in alignment')


def find_anchors(cdr_df: pd.DataFrame,
                 structure_df: pd.DataFrame,
                 num_anchors: int = 1) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''Get the anchors of a cdr loop. DOES NOT SUPPORT MULTIPLE MODELS.'''
    cdr_start_info = tuple(cdr_df.iloc[0][['chain_id', 'residue_seq_id', 'residue_insert_code']].tolist())
    cdr_end_info = tuple(cdr_df.iloc[-1][['chain_id', 'residue_seq_id', 'residue_insert_code']].tolist())

    residues = list(structure_df.groupby(['chain_id', 'residue_seq_id', 'residue_insert_code'], dropna=False))

    for i, (residue_info, residue) in enumerate(residues):
        if pd.isna(residue_info[-1]):
            residue_info = (residue_info[0], residue_info[1], None)

        if cdr_start_info == residue_info:
            start_anchor = pd.concat([res for _, res in residues[i - num_anchors: i]])

        if cdr_end_info == residue_info:
            end_anchor = pd.concat([res for _, res in residues[(i + 1): (i + 1) + num_anchors]])
            break

    return start_anchor, end_anchor


def main():
    args = parser.parse_args()
    stcrdab_summary = pd.read_csv(os.path.join(args.stcrdab_path, 'db_summary.dat'), delimiter='\t')
    selected_stcrdab = stcrdab_summary.copy()

    # Resolution better than 3.50 Å
    selected_stcrdab['resolution'] = pd.to_numeric(selected_stcrdab['resolution'], errors='coerce')
    selected_stcrdab = selected_stcrdab.query("resolution < 3.50")

    # alpha-beta TCRs
    selected_stcrdab = selected_stcrdab.query("TCRtype == 'abTCR'")

    # General clean: drop columns that don't contain anything useful
    selected_stcrdab = selected_stcrdab.loc[:, selected_stcrdab.nunique() > 1]
    selected_stcrdab = selected_stcrdab.dropna(axis=1, how='all')

    # Reset Index
    selected_stcrdab = selected_stcrdab.reset_index(drop=True)

    structure_names = []
    cdrs_with_anchors = {
        'alpha_chain': {1: [], 2: [], 3: []},
        'beta_chain': {1: [], 2: [], 3: []},
    }

    for _, row in selected_stcrdab.iterrows():
        structure_name = f'{row.pdb}_{row.Achain}{row.Bchain}'
        structure_names.append(structure_name)

        print('Collecting loops and anchors from', structure_name, flush=True)

        with open(os.path.join(args.stcrdab_path, 'imgt', row.pdb + '.pdb'), 'r') as fh:
            structure_df = parse_pdb_to_pandas(fh.read())

        tcr_df = structure_df.query(
            "record_type == 'ATOM' and (chain_id == @row.Achain or chain_id == @row.Bchain)"
        ).copy()

        tcr_df['chain_type'] = tcr_df['chain_id'].map({
            row.Achain: 'alpha_chain',
            row.Bchain: 'beta_chain',
        })
        tcr_df['cdr'] = tcr_df['residue_seq_id'].map(assign_cdr_number)
        tcr_backbone_df = tcr_df.query("atom_name == 'N' or atom_name == 'CA' or atom_name == 'C' or atom_name == 'O'")

        for chain in ('alpha_chain', 'beta_chain'):
            for cdr in 1, 2, 3:
                cdr_backbone_df = tcr_backbone_df.query('chain_type == @chain and cdr == @cdr').copy()
                start_anchor, end_anchor = find_anchors(cdr_backbone_df, tcr_backbone_df, args.number_of_anchors)

                cdrs_with_anchors[chain][cdr].append(pd.concat([start_anchor, cdr_backbone_df, end_anchor]))

    with open(os.path.join(args.output_dir, 'structure_names.txt'), 'w') as fh:
        fh.write('\n'.join(structure_names))
        fh.write('\n')

    for chain in ('alpha_chain', 'beta_chain'):
        for cdr in 1, 2, 3:
            print('Working on', chain, cdr, flush=True)

            cdr_loops = cdrs_with_anchors[chain][cdr]
            num_loops = len(cdr_loops)
            distance_matrix = np.zeros((num_loops, num_loops))

            for i in range(num_loops):
                for j in range(i + 1, num_loops):
                    loop_with_anchor_1 = cdr_loops[i]
                    loop_with_anchor_2 = cdr_loops[j]

                    anchor_coords_1 = get_coords(loop_with_anchor_1.query('cdr.isnull()'))
                    anchor_coords_2 = get_coords(loop_with_anchor_2.query('cdr.isnull()'))

                    # Align
                    loop_with_anchor_2 = align_pandas_structure(anchor_coords_2, anchor_coords_1, loop_with_anchor_2)

                    # Calculate DTW
                    loop_coords_1 = get_coords(loop_with_anchor_1.query('cdr.notnull()'))
                    loop_coords_2 = get_coords(loop_with_anchor_2.query('cdr.notnull()'))

                    distance = distance_fast(loop_coords_1.astype(np.double), loop_coords_2.astype(np.double))

                    # Add to distance matrix
                    distance_matrix[i, j] = distance

            # Write distance matrix
            name = f"cdr{cdr}_{chain.split('_')[0]}_distance_matrix.txt"
            np.savetxt(os.path.join(args.output_dir, name), distance_matrix)


if __name__ == '__main__':
    main()
