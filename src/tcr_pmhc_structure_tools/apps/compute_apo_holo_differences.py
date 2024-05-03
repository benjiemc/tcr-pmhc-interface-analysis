import argparse
import glob
import itertools
import logging
import os
import re

import pandas as pd
from python_pdb.aligners import align_pandas_structure
from python_pdb.comparisons import rmsd
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_structure_tools.apps._log import setup_logger
from tcr_pmhc_structure_tools.measurements import compute_residue_com, get_distance, measure_chi_angle
from tcr_pmhc_structure_tools.processing import annotate_tcr_df
from tcr_pmhc_structure_tools.utils import get_coords

logger = logging.getLogger()

parser = argparse.ArgumentParser()

parser.add_argument('input', help='path to data directory')
parser.add_argument('--output', '-o', help='path to output file')
parser.add_argument('--align-loops', action='store_true',
                    help='perform an alignment on the loops before computing RMSD.')
parser.add_argument('--per-residue', action='store_true',
                    help='Perform measurements on each residue individually as oppose to the loop as a whole.')
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

    info = {
        'complex_id': [],
        'structure_x_name': [],
        'structure_y_name': [],
        'chain_type': [],
        'cdr': [],
    }
    if args.per_residue:
        info['residue_name'] = []
        info['residue_seq_id'] = []
        info['residue_insert_code'] = []

    measurements = {}
    measurements['rmsd'] = []

    if args.per_residue:
        measurements['ca_distance'] = []
        measurements['chi_angle_change'] = []
        measurements['com_distance'] = []

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
                logger.debug('Computing differences in CDR-%s%d', 'A' if chain_type == 'alpha_chain' else 'B', cdr_num)

                tcr_cdr_x = structure_x.query('cdr == @cdr_num and chain_type == @chain_type')
                tcr_cdr_y = structure_y.query('cdr == @cdr_num and chain_type == @chain_type')

                tcr_cdr_backbone_x = tcr_cdr_x.query('backbone')
                tcr_cdr_backbone_y = tcr_cdr_y.query('backbone')

                tcr_cdr_backbone_coords_x = get_coords(tcr_cdr_backbone_x)
                tcr_cdr_backbone_coords_y = get_coords(tcr_cdr_backbone_y)

                if args.align_loops:
                    tcr_cdr_x = align_pandas_structure(tcr_cdr_backbone_coords_x,
                                                       tcr_cdr_backbone_coords_y,
                                                       tcr_cdr_x)
                    tcr_cdr_backbone_x = tcr_cdr_x.query('backbone')
                    tcr_cdr_backbone_coords_x = get_coords(tcr_cdr_backbone_x)

                if args.per_residue:
                    cdr_ca_distances = []
                    cdr_rmsds = []
                    cdr_angle_changes = []
                    cdr_com_distances = []

                    cdr_residue_names = []
                    cdr_seq_ids = []
                    cdr_insert_codes = []

                    common_columns = ['residue_name', 'residue_seq_id', 'residue_insert_code', 'atom_name']
                    cdr_comparison = pd.merge(tcr_cdr_x, tcr_cdr_y, how='inner', on=common_columns)
                    residues = cdr_comparison.groupby(['residue_seq_id', 'residue_insert_code', 'residue_name'],
                                                      dropna=False)

                    for (seq_id, insert_code, res_name), residue in residues:
                        res_x_columns = residue.filter(regex=r'_x$').columns.tolist()
                        res_y_columns = residue.filter(regex=r'_y$').columns.tolist()

                        res_x = residue[res_x_columns + common_columns]
                        res_y = residue[res_y_columns + common_columns]

                        res_x = res_x.rename(columns=lambda col_name: re.sub('_x$', '', col_name))
                        res_y = res_y.rename(columns=lambda col_name: re.sub('_y$', '', col_name))

                        # CA Distance
                        cdr_ca_distances.append(get_distance(get_coords(res_x.query("atom_name == 'CA'").iloc[0]),
                                                             get_coords(res_y.query("atom_name == 'CA'").iloc[0])))
                        # Compute Residue RMSD
                        try:
                            cdr_rmsds.append(rmsd(get_coords(res_x), get_coords(res_y)))
                        except ValueError:
                            logger.warning('Mismatched number of atoms in residue: %s %d%s',
                                           res_name,
                                           seq_id,
                                           insert_code if insert_code else '')
                            cdr_rmsds.append(None)

                        # Compute chi angle changes
                        if res_name == 'GLY' or res_name == 'ALA':
                            cdr_angle_changes.append(None)
                        else:
                            try:
                                cdr_angle_changes.append(measure_chi_angle(res_x) - measure_chi_angle(res_y))
                            except IndexError:
                                logger.warning('Missing atoms needed to calculate chi angle: %s %d%s',
                                               res_name,
                                               seq_id,
                                               insert_code if pd.notnull(insert_code) else '')
                                cdr_angle_changes.append(None)

                        # Compute C.O.M changes
                        cdr_com_distances.append(get_distance(compute_residue_com(res_x), compute_residue_com(res_y)))

                        cdr_residue_names.append(res_name)
                        cdr_seq_ids.append(seq_id)
                        cdr_insert_codes.append(insert_code)

                    num_residues = len(cdr_residue_names)

                    info['residue_name'] += cdr_residue_names
                    info['residue_seq_id'] += cdr_seq_ids
                    info['residue_insert_code'] += cdr_insert_codes

                    measurements['ca_distance'] += cdr_ca_distances
                    measurements['rmsd'] += cdr_rmsds
                    measurements['chi_angle_change'] += cdr_angle_changes
                    measurements['com_distance'] += cdr_com_distances

                    info['complex_id'] += [complex_id] * num_residues
                    info['structure_x_name'] += [comparison['file_name_x']] * num_residues
                    info['structure_y_name'] += [comparison['file_name_y']] * num_residues

                    info['cdr'] += [cdr_num] * num_residues
                    info['chain_type'] += [chain_type] * num_residues

                else:
                    info['complex_id'].append(complex_id)
                    info['structure_x_name'].append(comparison['file_name_x'])
                    info['structure_y_name'].append(comparison['file_name_y'])

                    info['chain_type'].append(chain_type)
                    info['cdr'].append(cdr_num)

                    measurements['rmsd'].append(rmsd(tcr_cdr_backbone_coords_x, tcr_cdr_backbone_coords_y))

    logger.info('Outputing results...')
    pd.DataFrame(info | measurements).to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
