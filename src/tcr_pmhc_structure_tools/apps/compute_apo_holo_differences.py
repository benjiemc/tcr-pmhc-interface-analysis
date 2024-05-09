import argparse
import glob
import logging
import os
import re

import numpy as np
import pandas as pd
from python_pdb.aligners import align_pandas_structure
from python_pdb.comparisons import rmsd
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_structure_tools.apps._log import setup_logger
from tcr_pmhc_structure_tools.measurements import compute_residue_com, get_distance, measure_chi_angle
from tcr_pmhc_structure_tools.processing import annotate_tcr_pmhc_df
from tcr_pmhc_structure_tools.utils import get_coords

logger = logging.getLogger()

parser = argparse.ArgumentParser()

parser.add_argument('input', help='path to data directory')
parser.add_argument('--output', '-o', help='path to output file')
parser.add_argument('--select-entities', choices=['tcr', 'pmhc'])
parser.add_argument('--align-entities', action='store_true',
                    help='perform an alignment on the selected entities before computing RMSD.')
parser.add_argument('--per-residue', action='store_true',
                    help='Perform measurements on each residue individually as oppose to the entity as a whole.')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error'], default='warning',
                    help="Level to log messages at (Default: 'warning')")


def split_merge(merged_df: pd.DataFrame,
                common_columns: list[str],
                suffixes: tuple[str, str] = ('_x', '_y')) -> pd.DataFrame:
    '''Split a dataframe back after a merge operation based on the columns used to merge and the suffixes.'''
    dfs = []
    for suffix in suffixes:
        columns = merged_df.filter(regex=f'{suffix}$').columns.tolist()
        df = merged_df[columns + common_columns]

        dfs.append(df.rename(columns=lambda col_name: re.sub(f'{suffix}$', '', col_name)))

    return tuple(dfs)


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
        'entity': [],
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

        comparison_structures = complex_summary.query("structure_type == @args.select_entities or state == 'holo'")
        comparisons = pd.merge(comparison_structures, comparison_structures, how='cross')
        comparisons['comparison'] = comparisons.apply(lambda row: '-'.join(sorted([row.file_name_x, row.file_name_y])),
                                                      axis='columns')
        comparisons = comparisons.drop_duplicates('comparison')
        comparisons = comparisons.drop('comparison', axis='columns')
        comparisons = comparisons.query('file_name_x != file_name_y')

        for _, comparison in comparisons.iterrows():
            logger.debug('Computing changes between %s and %s', comparison['file_name_x'], comparison['file_name_y'])

            structures = []
            for suffix in '_x', '_y':
                with open(os.path.join(complex_path, comparison['file_name' + suffix]), 'r') as fh:
                    structure_df = parse_pdb_to_pandas(fh.read())

                chains = comparison.filter(like='chain').filter(regex=f'{suffix}$').replace({np.nan: None}).tolist()
                structure_df = annotate_tcr_pmhc_df(structure_df, *chains)

                structure_df['backbone'] = structure_df['atom_name'].map(
                    lambda atom_name: (atom_name == 'N' or atom_name == 'CA' or atom_name == 'C' or atom_name == 'O')
                )

                structures.append(structure_df)

            structure_x, structure_y = structures

            if args.select_entities == 'tcr':
                entity_columns = ['chain_type', 'cdr']

            elif args.select_entities == 'pmhc':
                entity_columns = ['chain_type']

            structure_common_columns = entity_columns + ['residue_name',
                                                         'residue_seq_id',
                                                         'residue_insert_code',
                                                         'atom_name']

            structure_comparison = pd.merge(structure_x, structure_y, how='inner', on=structure_common_columns)

            # Necessary to avoid pandas warning
            if len(entity_columns) == 1:
                entity_columns = entity_columns[0]

            for entity_name, selected_entity in structure_comparison.groupby(entity_columns):
                logger.debug('Computing differences in %s', entity_name)
                entity_x, entity_y = split_merge(selected_entity, structure_common_columns)

                entity_backbone_x = entity_x.query('backbone')
                entity_backbone_y = entity_y.query('backbone')

                entity_backbone_coords_x = get_coords(entity_backbone_x)
                entity_backbone_coords_y = get_coords(entity_backbone_y)

                if args.align_entities:
                    entity_x = align_pandas_structure(entity_backbone_coords_x,
                                                      entity_backbone_coords_y,
                                                      entity_x)
                    entity_backbone_x = entity_x.query('backbone')
                    entity_backbone_coords_x = get_coords(entity_backbone_x)

                if args.per_residue:
                    entity_ca_distances = []
                    entity_rmsds = []
                    entity_angle_changes = []
                    entity_com_distances = []

                    entity_residue_names = []
                    entity_seq_ids = []
                    entity_insert_codes = []

                    residue_common_columns = ['residue_name', 'residue_seq_id', 'residue_insert_code', 'atom_name']
                    entity_comparison = pd.merge(entity_x, entity_y, how='inner', on=residue_common_columns)
                    residues = entity_comparison.groupby(['residue_seq_id', 'residue_insert_code', 'residue_name'],
                                                         dropna=False)

                    for (seq_id, insert_code, res_name), residue in residues:
                        res_x, res_y = split_merge(residue, residue_common_columns)

                        # CA Distance
                        entity_ca_distances.append(get_distance(get_coords(res_x.query("atom_name == 'CA'").iloc[0]),
                                                                get_coords(res_y.query("atom_name == 'CA'").iloc[0])))
                        # Compute Residue RMSD
                        try:
                            entity_rmsds.append(rmsd(get_coords(res_x), get_coords(res_y)))
                        except ValueError:
                            logger.warning('Mismatched number of atoms in residue: %s %d%s',
                                           res_name,
                                           seq_id,
                                           insert_code if insert_code else '')
                            entity_rmsds.append(None)

                        # Compute chi angle changes
                        if res_name == 'GLY' or res_name == 'ALA':
                            entity_angle_changes.append(None)
                        else:
                            try:
                                entity_angle_changes.append(measure_chi_angle(res_x) - measure_chi_angle(res_y))
                            except IndexError:
                                logger.warning('Missing atoms needed to calculate chi angle: %s %d%s',
                                               res_name,
                                               seq_id,
                                               insert_code if pd.notnull(insert_code) else '')
                                entity_angle_changes.append(None)

                        # Compute C.O.M changes
                        entity_com_distances.append(get_distance(compute_residue_com(res_x),
                                                                 compute_residue_com(res_y)))

                        entity_residue_names.append(res_name)
                        entity_seq_ids.append(seq_id)
                        entity_insert_codes.append(insert_code)

                    num_residues = len(entity_residue_names)

                    info['residue_name'] += entity_residue_names
                    info['residue_seq_id'] += entity_seq_ids
                    info['residue_insert_code'] += entity_insert_codes

                    measurements['ca_distance'] += entity_ca_distances
                    measurements['rmsd'] += entity_rmsds
                    measurements['chi_angle_change'] += entity_angle_changes
                    measurements['com_distance'] += entity_com_distances

                    info['complex_id'] += [complex_id] * num_residues
                    info['structure_x_name'] += [comparison['file_name_x']] * num_residues
                    info['structure_y_name'] += [comparison['file_name_y']] * num_residues

                    info['entity'] += [entity_name] * num_residues

                else:
                    info['complex_id'].append(complex_id)
                    info['structure_x_name'].append(comparison['file_name_x'])
                    info['structure_y_name'].append(comparison['file_name_y'])

                    info['entity'].append(entity_name)

                    measurements['rmsd'].append(rmsd(entity_backbone_coords_x, entity_backbone_coords_y))

    if args.select_entities == 'tcr':
        info['chain_type'] = [chain_type for chain_type, _ in info['entity']]
        info['cdr'] = [int(cdr) for _, cdr in info['entity']]
        info.pop('entity')

    elif args.select_entities == 'pmhc':
        info['chain_type'] = info['entity']
        info.pop('entity')

    logger.info('Outputing results...')

    output_columns = ['complex_id', 'structure_x_name', 'structure_y_name']

    if args.select_entities == 'tcr':
        output_columns += ['chain_type', 'cdr']

    elif args.select_entities == 'pmhc':
        output_columns += ['chain_type']

    if args.per_residue:
        output_columns += ['residue_name', 'residue_seq_id', 'residue_insert_code']

    output_columns.append('rmsd')
    if args.per_residue:
        output_columns += ['ca_distance', 'chi_angle_change', 'com_distance']

    pd.DataFrame(info | measurements).loc[:, output_columns].to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
