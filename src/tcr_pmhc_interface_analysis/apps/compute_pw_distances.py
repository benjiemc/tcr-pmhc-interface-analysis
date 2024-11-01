'''Compute the pairwise DTW distance between all loops in the STCRDab.'''
import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
from dtaidistance.dtw_ndim import distance_fast
from python_pdb.aligners import align_pandas_structure
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_interface_analysis.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_interface_analysis.processing import annotate_tcr_pmhc_df, find_anchors
from tcr_pmhc_interface_analysis.utils import get_coords

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('stcrdab', help='Path to STCRDab')
parser.add_argument('--output', '-o', help='Path to output location')
parser.add_argument('--resolution-cutoff', type=float, default=3.50,
                    help='Maximum resolution allowed from the structures (Default: 3.50)')
parser.add_argument('--number-of-anchors', type=int, default=5,
                    help='number of anchors to include in alignment (Default: 5)')
parser.add_argument('--compress-output', action='store_true', help='compress the output matrices using gzip')

add_logging_arguments(parser)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    stcrdab_summary = pd.read_csv(os.path.join(args.stcrdab, 'db_summary.dat'), delimiter='\t')

    stcrdab_summary['resolution'] = pd.to_numeric(stcrdab_summary['resolution'], errors='coerce')
    stcrdab_summary = stcrdab_summary.query("resolution <= @args.resolution_cutoff")

    stcrdab_summary = stcrdab_summary.query("TCRtype == 'abTCR'")

    # Reset Index
    stcrdab_summary = stcrdab_summary.reset_index(drop=True)

    structure_names = []
    cdrs_with_anchors = {
        'alpha_chain': {1: [], 2: [], 3: []},
        'beta_chain': {1: [], 2: [], 3: []},
    }

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    for _, row in stcrdab_summary.iterrows():
        structure_name = f'{row.pdb}_{row.Achain}{row.Bchain}'
        structure_names.append(structure_name)

        logger.info('Collecting loops and anchors from %s', structure_name)

        with open(os.path.join(args.stcrdab, 'imgt', row.pdb + '.pdb'), 'r') as fh:
            structure_df = parse_pdb_to_pandas(fh.read())

        structure_df = annotate_tcr_pmhc_df(structure_df, alpha_chain_id=row.Achain, beta_chain_id=row.Bchain)

        tcr_df = structure_df.query("chain_type == 'alpha_chain' or chain_type == 'beta_chain'").copy()
        tcr_backbone_df = tcr_df.query("atom_name == 'N' or atom_name == 'CA' or atom_name == 'C' or atom_name == 'O'")

        for chain in ('alpha_chain', 'beta_chain'):
            for cdr in 1, 2, 3:
                cdr_backbone_df = tcr_backbone_df.query('chain_type == @chain and cdr == @cdr').copy()
                start_anchor, end_anchor = find_anchors(cdr_backbone_df, tcr_backbone_df, args.number_of_anchors)

                cdrs_with_anchors[chain][cdr].append(pd.concat([start_anchor, cdr_backbone_df, end_anchor]))

    with open(os.path.join(args.output, 'structure_names.txt'), 'w') as fh:
        fh.write('\n'.join(structure_names))
        fh.write('\n')

    for chain in ('alpha_chain', 'beta_chain'):
        for cdr in 1, 2, 3:
            logger.info('Working on %s %d', chain, cdr)

            cdr_loops = cdrs_with_anchors[chain][cdr]
            num_loops = len(cdr_loops)
            distance_matrix = np.zeros((num_loops, num_loops))

            for i in range(num_loops):
                for j in range(i + 1, num_loops):
                    loop_with_anchor_1 = cdr_loops[i]
                    loop_with_anchor_2 = cdr_loops[j]

                    anchor_coords_1 = get_coords(loop_with_anchor_1.query('cdr.isnull()'))
                    anchor_coords_2 = get_coords(loop_with_anchor_2.query('cdr.isnull()'))

                    loop_with_anchor_2 = align_pandas_structure(anchor_coords_2, anchor_coords_1, loop_with_anchor_2)

                    loop_coords_1 = get_coords(loop_with_anchor_1.query('cdr.notnull()'))
                    loop_coords_2 = get_coords(loop_with_anchor_2.query('cdr.notnull()'))

                    distance = distance_fast(loop_coords_1.astype(np.double), loop_coords_2.astype(np.double))

                    distance_matrix[i, j] = distance

            distance_matrix = np.maximum(distance_matrix, distance_matrix.transpose())

            logger.info('Writing distance matrix')
            name = f"cdr{cdr}_{chain.split('_')[0]}_distance_matrix.txt"

            if args.compress_output:
                name += '.gz'

            np.savetxt(os.path.join(args.output, name), distance_matrix)


if __name__ == '__main__':
    main()
