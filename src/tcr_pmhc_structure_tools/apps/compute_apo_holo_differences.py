import argparse
import glob
import itertools
import logging
import os

import pandas as pd
from python_pdb.aligners import align_pandas_structure
from python_pdb.comparisons import rmsd
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_structure_tools.apps._log import setup_logger
from tcr_pmhc_structure_tools.processing import annotate_tcr_df
from tcr_pmhc_structure_tools.utils import get_coords

logger = logging.getLogger()

parser = argparse.ArgumentParser()

parser.add_argument('input', help='path to data directory')
parser.add_argument('--output', '-o', help='path to output file')
parser.add_argument('--align-loops', action='store_true',
                    help='perform an alignment on the loops before computing RMSD.')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error'], default='warning',
                    help="Level to log messages at (Default: 'warning')")


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    summary_path, = glob.glob(os.path.join(args.input, '*summary.csv'))
    summary_df = pd.read_csv(summary_path)

    complexes = [complex_id for complex_id in os.listdir(args.input)
                 if os.path.isdir(os.path.join(args.input, complex_id))]
    num_complexes = len(complexes)

    complex_ids = []
    structure_x_names = []
    structure_y_names = []
    chain_types = []
    cdrs = []
    rmsds = []

    for num, complex_id in enumerate(complexes, 1):
        logger.info('%s - %d of %d', complex_id, num, num_complexes)

        complex_path = os.path.join(args.input, complex_id)
        complex_pdb_files = [file_ for file_ in os.listdir(complex_path) if file_.endswith('.pdb')]
        complex_summary = summary_df[summary_df['file_name'].isin(complex_pdb_files)]

        complex_tcrs = complex_summary.query("structure_type == 'tcr' or structure_type == 'tcr_pmhc'")
        comparisons = complex_tcrs.merge(complex_tcrs, how='cross')
        comparisons = comparisons.query('file_name_x != file_name_y')

        for _, comparison in comparisons.iterrows():
            logger.debug('Computing changes between %s and %s', comparison['file_name_x'], comparison['file_name_y'])

            structures = []
            for suffix in '_x', '_y':
                with open(os.path.join(complex_path, comparison['file_name' + suffix]), 'r') as fh:
                    structure_df = parse_pdb_to_pandas(fh.read())

                structure_df = annotate_tcr_df(structure_df,
                                               alpha_chain_id=comparison['alpha_chain' + suffix],
                                               beta_chain_id=comparison['beta_chain' + suffix])
                structure_df['backbone'] = structure_df['atom_name'].map(
                    lambda atom_name: (atom_name == 'N' or atom_name == 'CA' or atom_name == 'C' or atom_name == 'O')
                )

                structures.append(structure_df)

            structure_x, structure_y = structures

            for chain_type, cdr_num in itertools.product(('alpha_chain', 'beta_chain'), (1, 2, 3)):
                tcr_cdr_backbone_x = structure_x.query('cdr == @cdr_num and chain_type == @chain_type and backbone')
                tcr_cdr_backbone_y = structure_y.query('cdr == @cdr_num and chain_type == @chain_type and backbone')

                tcr_cdr_backbone_coords_x = get_coords(tcr_cdr_backbone_x)
                tcr_cdr_backbone_coords_y = get_coords(tcr_cdr_backbone_y)

                if args.align_loops:
                    tcr_cdr_backbone_x = align_pandas_structure(tcr_cdr_backbone_coords_x,
                                                                tcr_cdr_backbone_coords_y,
                                                                tcr_cdr_backbone_x)
                    tcr_cdr_backbone_coords_x = get_coords(tcr_cdr_backbone_x)

                complex_ids.append(complex_id)
                structure_x_names.append(comparison['file_name_x'])
                structure_y_names.append(comparison['file_name_y'])

                chain_types.append(chain_type)
                cdrs.append(cdr_num)

                rmsds.append(rmsd(tcr_cdr_backbone_coords_x, tcr_cdr_backbone_coords_y))

    pd.DataFrame({
        'complex_id': complex_ids,
        'structure_x_name': structure_x_names,
        'structure_y_name': structure_y_names,
        'chain_type': chain_types,
        'cdr': cdrs,
        'rmsd': rmsds,
    }).to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
