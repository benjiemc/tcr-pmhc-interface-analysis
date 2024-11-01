'''Cluster TCR CDR Loop structures and assign them as canonical or pseudo clusters.

Adding cannonical conformation information to clusters
------------------------------------------------------

Based on the definition from Comparative Analysis of the CDR Loops of Antigen Receptors
(https://doi.org/10.3389/fimmu.2019.02454), canonical clusters are any cluster where there are more than two unique
sequences within the density clusters. The other clusters will be refered to as pseudo-clusters, as these may just be
the effect of the same loop finding the same conformation.

'''
import argparse
import logging
import os
import sys

import hdbscan
import numpy as np
import pandas as pd
from python_pdb.formats.residue import THREE_TO_ONE_CODE
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_interface_analysis.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_interface_analysis.processing import annotate_tcr_pmhc_df

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('structure_names', help='path to the structure names file')
parser.add_argument('distance_matrices', nargs='+', help='paths to the distance matrices')
parser.add_argument('--output', '-o', help='output path')
parser.add_argument('--assign-cluster-types', action='store_true',
                    help='assign cluster types (requires --stcrdab-path input)')
parser.add_argument('--stcrdab-path', required=False, help='path to the STCRDab')

add_logging_arguments(parser)


def get_sequence(df):
    residue_names = df.drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code'])['residue_name']
    return ''.join(residue_names.map(THREE_TO_ONE_CODE).tolist())


def get_cdr_sequences(names: list[str], stcrdab_path: str) -> pd.DataFrame:
    '''Get CDR Sequences for a list of structure names (format: <pdb_id>_<alpha_chain_id><beta_chain_id>).'''
    def load_cdrs(pdb_id, alpha_chain_id, beta_chain_id):
        with open(os.path.join(stcrdab_path, 'imgt', pdb_id + '.pdb'), 'r') as fh:
            structure_df = parse_pdb_to_pandas(fh.read())

        structure_df = annotate_tcr_pmhc_df(structure_df, alpha_chain_id, beta_chain_id)
        tcr_df = structure_df.query('chain_type.notnull()')

        cdr_sequences = tcr_df.groupby(['chain_type', 'cdr']).apply(get_sequence)
        cdr_sequences.name = 'sequence'

        index = cdr_sequences.index.to_flat_index()
        index = index.map(lambda index: f"cdr_{'a' if index[0] == 'alpha_chain' else 'b'}{int(index[1])}_sequence")
        cdr_sequences.index = index

        return cdr_sequences

    structures = pd.DataFrame({'name': names})

    structures[['pdb_id', 'chains']] = structures['name'].str.split('_').apply(pd.Series)
    structures[['alpha_chain_id', 'beta_chain_id']] = structures['chains'].apply(list).apply(pd.Series)

    structures = pd.concat([structures,
                            structures.apply(lambda row: load_cdrs(row.pdb_id, row.alpha_chain_id, row.beta_chain_id),
                                             axis=1)], axis=1)
    structures = structures.melt(id_vars=['name'],
                                 value_vars=['cdr_a1_sequence', 'cdr_a2_sequence', 'cdr_a3_sequence',
                                             'cdr_b1_sequence', 'cdr_b2_sequence', 'cdr_b3_sequence'],
                                 value_name='sequence')

    structures[['chain_type', 'cdr']] = structures['variable'].map(
        lambda name: tuple(name.split('_')[1])
    ).apply(pd.Series)

    structures['chain_type'] = structures['chain_type'].map(
        lambda letter: 'alpha_chain' if letter == 'a' else 'beta_chain'
    )

    return structures


def assign_cluster_types(df: pd.DataFrame, min_uniq: int = 2) -> pd.Series:
    '''Assign clusters as canonical or pseudo'''
    cluster_types = df.query("cluster != 'noise'").groupby(
        ['chain_type', 'cdr', 'cluster'],
    )['sequence'].agg(lambda sequences: 'canonical' if sequences.nunique() > min_uniq else 'pseudo')
    cluster_types.name = 'cluster_type'

    return cluster_types


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Loading structure names from %s', args.structure_names)
    with open(args.structure_names, 'r') as fh:
        structure_names = [line.strip() for line in fh.readlines()]

    logger.info('Computing Cluster')
    df = pd.DataFrame()
    for path in args.distance_matrices:
        name = os.path.basename(path).split('.')[0].replace('_distance_matrix', '')

        cdr, chain = name.split('_')
        cdr = cdr.replace('cdr', '')

        logger.info('Loading CDR%s%s distance matrix', cdr, chain)

        cdr_distance_matrix = np.loadtxt(path)

        logger.info('Clustering loops')
        cdr_clusters = hdbscan.HDBSCAN(min_cluster_size=5, metric='precomputed').fit_predict(cdr_distance_matrix)

        cdr_df = pd.DataFrame({
            'name': structure_names,
            'cluster': cdr_clusters,
        })
        cdr_df['chain_type'] = chain + '_chain'
        cdr_df['cdr'] = cdr

        df = pd.concat([df, cdr_df])

    df['cluster'] = df['cluster'].apply(str)
    df['cluster'] = df['cluster'].replace('-1', 'noise')

    if args.assign_cluster_types:
        logger.info('Assigning cluster types')
        structures = get_cdr_sequences(df['name'].unique(), args.stcrdab_path)

        df = df.merge(structures[['name', 'chain_type', 'cdr', 'sequence']],
                      how='left',
                      on=['name', 'chain_type', 'cdr'])

        cluster_types = assign_cluster_types(df)

        df = df.merge(cluster_types.reset_index(), how='left', on=['chain_type', 'cdr', 'cluster'])

    logger.info('Outputting clusters to %s', args.output)
    df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
