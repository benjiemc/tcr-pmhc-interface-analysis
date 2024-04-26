import argparse
import itertools
import os

import pandas as pd
from python_pdb.aligners import align_pandas_structure
from python_pdb.comparisons import rmsd
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_structure_tools.processing import annotate_tcr_df
from tcr_pmhc_structure_tools.utils import get_coords

parser = argparse.ArgumentParser()

parser.add_argument('input', help='path to data directory')
parser.add_argument('--output', '-o', help='path to output file')
parser.add_argument('--align-loops', action='store_true',
                    help='perform an alignment on the loops before computing RMSD.')


def main():
    args = parser.parse_args()

    complexes = [complex_id for complex_id in os.listdir(args.input)
                 if os.path.isdir(os.path.join(complex_id, args.input))]
    num_complexes = len(complexes)

    complex_ids = []
    apo_ranks = []
    holo_ranks = []
    chain_types = []
    cdrs = []
    rmsds = []

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

            for chain_type, cdr_num in itertools.product(('alpha_chain', 'beta_chain'), (1, 2, 3)):
                tcr_cdr_backbone = tcr_df.query('cdr == @cdr_num and chain_type == @chain_type and backbone')
                tcr_pmhc_cdr_backbone = tcr_pmhc_df.query('cdr == @cdr_num and chain_type == @chain_type and backbone')

                tcr_cdr_backbone_coords = get_coords(tcr_cdr_backbone)
                tcr_pmhc_cdr_backbone_coords = get_coords(tcr_pmhc_cdr_backbone)

                if args.align_loops:
                    tcr_pmhc_cdr_backbone = align_pandas_structure(tcr_pmhc_cdr_backbone_coords,
                                                                   tcr_cdr_backbone_coords,
                                                                   tcr_pmhc_cdr_backbone)
                    tcr_pmhc_cdr_backbone_coords = get_coords(tcr_pmhc_cdr_backbone)

                complex_ids.append(complex_id)
                apo_ranks.append(apo_rank)
                holo_ranks.append(holo_rank)

                chain_types.append(chain_type)
                cdrs.append(cdr_num)

                rmsds.append(rmsd(tcr_cdr_backbone_coords, tcr_pmhc_cdr_backbone_coords))

    pd.DataFrame({
        'complex_id': complex_ids,
        'apo_rank': apo_ranks,
        'holo_rank': holo_ranks,
        'chain_type': chain_types,
        'cdr': cdrs,
        'rmsd': rmsds,
    }).to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
