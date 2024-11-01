'''Select apo and holo structures from the STCRDab and Histo.fyi.'''
import argparse
import glob
import logging
import os
import sys

import numpy as np
import pandas as pd
from python_pdb.entities import Structure
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_interface_analysis.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_interface_analysis.histo_fyi_utils import (PMHC_CLASS_I_URL, TCR_PMHC_CLASS_I_URL, fetch_structure,
                                                         retrieve_data_from_api)
from tcr_pmhc_interface_analysis.missing_residues import (get_raw_structures_with_missing_residues,
                                                          screen_pmhcs_for_missing_residues,
                                                          screen_tcrs_for_missing_residues)
from tcr_pmhc_interface_analysis.stcrdab_utils import (get_ab_tcr_mhc_class_Is_from_stcrdab, get_ab_tcrs_from_stcrdab,
                                                       get_stcrdab_sequences)

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('stcrdab', help='Path to STCRDab')
parser.add_argument('--resolution-cutoff', type=float, default=3.50,
                    help='Maximum resolution allowed from the structures (Default: 3.50)')
parser.add_argument('--add-mhcs', action='store_true', help='Add MHC apo forms to all holo structures')
parser.add_argument('--drop-duplicate-ids', action='store_true',
                    help='Only keep one copy of the structure from a pdb id')
parser.add_argument('--output', '-o', help='Path to output location')

add_logging_arguments(parser)


def select_apo_holo(group: pd.DataFrame) -> bool:
    '''Select groups with both *apo* and *holo* confomations'''
    return 'apo' in group['state'].unique().tolist() and 'holo' in group['state'].unique().tolist()


def screen_quality(df: pd.DataFrame,
                   structure_type: str,
                   resolution_cutoff: float = 3.50,
                   stcrdab_path: str | None = None) -> pd.DataFrame:
    '''Screen structures for quality.

    Resolution must be below threshold and structures can not be missing residues in key domains. For TCRs, this is the
     variable domain, for pMHCs, this is the peptide.

    Args:
        df: dataframe containing structure information
        structure_type: either 'tcr', 'pmhc', or 'tcr-pmhc'
        resolution_cuttoff: maximum allowed resolution of structures in ångstroms (Default: 3.50)

    Returns:
        Dataframe with structures that meet these criteria

    '''
    selected_df = df.copy()

    # Resolution Screen
    selected_df = selected_df[selected_df['resolution'] <= resolution_cutoff]

    # Missing Residues Screen
    raw_structures = get_raw_structures_with_missing_residues(df['pdb_id'].unique().tolist(), stcrdab_path)

    match structure_type:
        case 'tcr':
            valid_structures = screen_tcrs_for_missing_residues(selected_df, raw_structures)

        case 'pmhc':
            valid_structures = screen_pmhcs_for_missing_residues(selected_df, raw_structures)

        case 'tcr-pmhc':
            valid_structures_tcr = screen_tcrs_for_missing_residues(selected_df, raw_structures)
            valid_structures_pmhc = screen_pmhcs_for_missing_residues(selected_df, raw_structures)
            valid_structures = [valid_tcr and valid_pmhc
                                for valid_tcr, valid_pmhc in zip(valid_structures_tcr, valid_structures_pmhc)]

    selected_df = selected_df[valid_structures]

    return selected_df


def main():
    args = parser.parse_args()
    setup_logger(logger, level=args.log_level)

    logger.info('Loading STCRDab Summary')
    stcrdab_summary = pd.read_csv(os.path.join(args.stcrdab, 'db_summary.dat'), delimiter='\t')

    stcrdab_summary['resolution'] = pd.to_numeric(stcrdab_summary['resolution'], errors='coerce')
    stcrdab_summary['file_path_imgt'] = stcrdab_summary['pdb'].map(
        lambda pdb_id: os.path.join(args.stcrdab, 'imgt', pdb_id + '.pdb')
    )
    stcrdab_summary['file_path_raw'] = stcrdab_summary['pdb'].map(
        lambda pdb_id: os.path.join(args.stcrdab, 'raw', pdb_id + '.pdb')
    )
    stcrdab_summary = stcrdab_summary.rename({'pdb': 'pdb_id', 'Achain': 'alpha_chain', 'Bchain': 'beta_chain'},
                                             axis='columns')

    stcrdab_summary['chains'] = stcrdab_summary[['alpha_chain',
                                                 'beta_chain',
                                                 'antigen_chain',
                                                 'mhc_chain1',
                                                 'mhc_chain2']].apply(lambda chains: '-'.join(chains.dropna()), axis=1)

    logger.info('Retrieving unbound abTCRs from STCRDab')
    stcrdab_tcrs = get_ab_tcrs_from_stcrdab(stcrdab_summary)
    stcrdab_tcrs = get_stcrdab_sequences(stcrdab_tcrs, structure_type='tcr')
    stcrdab_tcrs['state'] = 'apo'
    stcrdab_tcrs['structure_type'] = 'tcr'

    if args.add_mhcs:
        logger.info('Retrieving unbound pMHCs from histo.fyi')
        histo_pmhcs = retrieve_data_from_api(PMHC_CLASS_I_URL)
        histo_pmhcs['state'] = 'apo'
        histo_pmhcs['structure_type'] = 'pmhc'

    logger.info('Retrieving TCR-pMHC')
    stcrdab_tcr_pmhcs = get_ab_tcr_mhc_class_Is_from_stcrdab(stcrdab_summary)
    stcrdab_tcr_pmhcs = get_stcrdab_sequences(stcrdab_tcr_pmhcs,  structure_type='tcr-pmhc')
    histo_tcr_pmhcs = retrieve_data_from_api(TCR_PMHC_CLASS_I_URL)

    merged_tcr_pmhcs = stcrdab_tcr_pmhcs.merge(histo_tcr_pmhcs,
                                               how='outer',
                                               left_on=['pdb_id', 'antigen_chain', 'mhc_chain1'],
                                               right_on=['pdb_id', 'antigen_chain', 'mhc_chain1'],
                                               suffixes=['_stcrdab', '_histo'])

    merged_tcr_pmhcs['mhc_chain2'] = merged_tcr_pmhcs['mhc_chain2_stcrdab'].combine_first(
        merged_tcr_pmhcs['mhc_chain2_histo'],
    )
    merged_tcr_pmhcs = merged_tcr_pmhcs.drop(['mhc_chain2_stcrdab', 'mhc_chain2_histo'], axis='columns')

    merged_tcr_pmhcs['resolution_stcrdab'] = merged_tcr_pmhcs['resolution_stcrdab'].fillna(np.inf)
    merged_tcr_pmhcs['resolution_histo'] = merged_tcr_pmhcs['resolution_histo'].fillna(np.inf)
    merged_tcr_pmhcs['resolution'] = merged_tcr_pmhcs.apply(
        lambda row: np.min([row.resolution_stcrdab, row.resolution_histo]),
        axis=1,
    )
    merged_tcr_pmhcs = merged_tcr_pmhcs.drop(['resolution_stcrdab', 'resolution_histo'], axis='columns')

    merged_tcr_pmhcs['state'] = 'holo'
    merged_tcr_pmhcs['structure_type'] = 'tcr_pmhc'

    logger.info('Checking quality of structures')

    logger.info('Screening TCRs')
    stcrdab_tcrs = screen_quality(stcrdab_tcrs, 'tcr', args.resolution_cutoff, args.stcrdab)

    if args.add_mhcs:
        logger.info('Screening pMHCs')
        histo_pmhcs = screen_quality(histo_pmhcs, 'pmhc', args.resolution_cutoff)

    logger.info('Screening TCR-pMHCs')
    merged_tcr_pmhcs = screen_quality(merged_tcr_pmhcs, 'tcr-pmhc', args.resolution_cutoff, args.stcrdab)

    if args.drop_duplicate_ids:
        logger.info('Removing duplicate PDB IDs')
        stcrdab_tcrs = stcrdab_tcrs.drop_duplicates('pdb_id')
        merged_tcr_pmhcs = merged_tcr_pmhcs.drop_duplicates('pdb_id')

        if args.add_mhcs:
            histo_pmhcs = histo_pmhcs.drop_duplicates('pdb_id')

    logger.info('Finding examples of apo and holo structures')
    apo_holo_tcrs = pd.concat([stcrdab_tcrs, merged_tcr_pmhcs])
    apo_holo_tcrs = apo_holo_tcrs.groupby('cdr_sequences_collated').filter(select_apo_holo)

    if args.add_mhcs:
        apo_holo_pmhcs = pd.concat([histo_pmhcs, merged_tcr_pmhcs])
        apo_holo_pmhcs = apo_holo_pmhcs.groupby(['mhc_slug', 'peptide_sequence']).filter(select_apo_holo)

        apo_holo = pd.concat([apo_holo_tcrs, apo_holo_pmhcs])

    else:
        apo_holo = apo_holo_tcrs

    apo_holo = apo_holo.drop_duplicates()

    # Clean DataFrame
    apo_holo['chains'] = apo_holo['chains'].combine_first(apo_holo['chains_stcrdab'])
    apo_holo['chains'] = apo_holo['chains'].combine_first(apo_holo['chains_histo'])
    apo_holo = apo_holo.drop(['chains_stcrdab', 'chains_histo'], axis='columns')

    apo_holo['file_name'] = (apo_holo['pdb_id'] + '_' + apo_holo['chains'] + '_' + apo_holo['structure_type'] + '.pdb')

    logger.info('Exporting')
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    logger.info('Writing summary file')
    apo_holo[['file_name',
              'pdb_id',
              'structure_type',
              'state',
              'alpha_chain',
              'beta_chain',
              'antigen_chain',
              'mhc_chain1',
              'mhc_chain2',
              'cdr_sequences_collated',
              'peptide_sequence',
              'mhc_slug']].to_csv(os.path.join(args.output, 'apo_holo_summary.csv'), index=False)

    logger.info('Collecting structures')

    stcrdab_pdb_ids = [path.rsplit('/', 1)[-1].split('.')[0]
                       for path in glob.glob(os.path.join(args.stcrdab, 'imgt', '*.pdb'))]
    num_structures = len(apo_holo)

    for num, (_, row) in enumerate(apo_holo.iterrows(), 1):
        logger.debug('Exporting PDB ID %s - %d of %d', row.pdb_id, num, num_structures)

        match row.structure_type:
            case 'tcr':
                with open(os.path.join(args.stcrdab, 'imgt', f'{row.pdb_id}.pdb'), 'r') as fh:
                    structure_df = parse_pdb_to_pandas(fh.read())

                tcr_df = structure_df.query('chain_id == @row.alpha_chain or chain_id == @row.beta_chain')
                output_text = str(Structure.from_pandas(tcr_df))

            case 'pmhc':
                output_text = fetch_structure(row.pdb_id, row.assembly_number)

            case 'tcr_pmhc':
                if row.pdb_id not in stcrdab_pdb_ids:
                    output_text = fetch_structure(row.pdb_id, row.assembly_number)

                else:
                    with open(os.path.join(args.stcrdab, 'imgt', f'{row.pdb_id}.pdb'), 'r') as fh:
                        structure_df = parse_pdb_to_pandas(fh.read())

                    tcr_pmhc_df = structure_df.query(('chain_id == @row.alpha_chain '
                                                      'or chain_id == @row.beta_chain '
                                                      'or chain_id == @row.antigen_chain '
                                                      'or chain_id == @row.mhc_chain1 '
                                                      'or chain_id == @row.mhc_chain2'))

                    output_text = str(Structure.from_pandas(tcr_pmhc_df))

        with open(os.path.join(args.output, row.file_name), 'w') as fh:
            fh.write(output_text)


if __name__ == '__main__':
    main()
