import argparse
import os
import glob
import logging
import shutil

import pandas as pd
from pymol import cmd

from tcr_pmhc_structure_tools.apps._log import setup_logger
from tcr_pmhc_structure_tools.imgt_numbering import IMGT_CDR, IMGT_VARIABLE_DOMAIN

logger = logging.getLogger()

parser = argparse.ArgumentParser()

parser.add_argument('structures', help='')
parser.add_argument('--output', '-o', help='')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error'], default='warning',
                    help="Level to log messages at (Default: 'warning')")


def get_floor_selection() -> str:
    '''Create pymol selection for the floor of the MHC binding groove.'''
    mhc_floor_residues = [
        (1, 14),  # A1
        (18, 28),  # B1
        (31, 38),  # C1
        (1001, 1014),  # A2
        (1018, 1028),  # B2
        (1031, 1038),  # C2
        (1042, 1049),  # D2
    ]
    return ' or '.join([f'resi {high}-{low}' for high, low in mhc_floor_residues])


def get_framework_selection() -> str:
    '''Create pymol selection for the framework region of a TCR based on IMGT numbering.'''
    indices = [index for index in IMGT_VARIABLE_DOMAIN if index not in IMGT_CDR]
    return ' or '.join([f'resi {index}' for index in indices])


def align_tcr(mobile_path: str, alpha_chain_id: str, beta_chain_id: str, name: str = 'tcr') -> None:
    '''Align a group of TCRs to a reference structure.'''
    cmd.load(mobile_path, name)

    framework_selection = (f'{name} and '
                           f'(chain {alpha_chain_id} or chain {beta_chain_id}) and '
                           f'({get_framework_selection()})')
    cmd.select(f'{name}-Fw', framework_selection)

    cmd.align(f'{name}-Fw', 'target-Fw')


def align_pmhc(mobile_path: str, mhc_chain_id: str, name: str = 'pmhc') -> None:
    '''Align a (TCR-)pMHC to a reference structure already loaded in the environment.'''
    cmd.load(mobile_path, name)

    floor_selection = (f'{name} and '
                       f'chain {mhc_chain_id} and '
                       f'({get_floor_selection()})')
    cmd.select(f'{name}-floor', floor_selection)

    cmd.align(f'{name}-floor', 'target-floor')


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    summary_path, = glob.glob(os.path.join(args.structures, '*summary.csv'))
    summary_df = pd.read_csv(summary_path)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    holo_structures = summary_df.query("state == 'holo' and structure_type == 'tcr_pmhc'")
    num_holo_structures = len(holo_structures)

    for num, (_, holo_structure) in enumerate(holo_structures.iterrows(), 1):
        holo_name = holo_structure.file_name.replace('.pdb', '')
        logger.info('Aligning %s - %d of %d', holo_name, num, num_holo_structures)

        apo_tcrs = summary_df.query(('cdr_sequences_collated == @holo_structure.cdr_sequences_collated '
                                     "and state == 'apo' and structure_type == 'tcr'"))

        apo_pmhcs = summary_df.query(('peptide_sequence == @holo_structure.peptide_sequence '
                                      'and mhc_slug == @holo_structure.mhc_slug '
                                      "and state == 'apo' and structure_type == 'pmhc'"))

        cmd.load(os.path.join(args.structures, holo_structure.file_name), 'target')
        framework_selection = ('target and '
                               f'(chain {holo_structure.alpha_chain} or chain {holo_structure.beta_chain}) and '
                               f'({get_framework_selection()})')
        cmd.select('target-Fw', framework_selection)

        floor_selection = ('target and '
                           f'chain {holo_structure.mhc_chain1} and '
                           f'({get_floor_selection()})')
        cmd.select('target-floor', floor_selection)

        group_path = os.path.join(args.output, holo_name)
        if not os.path.exists(group_path):
            os.mkdir(group_path)

        logger.info('Aliging %d TCR(s)', len(apo_tcrs))
        for _, apo_structure in apo_tcrs.iterrows():
            apo_name = apo_structure.file_name.replace('.pdb', '')
            logger.debug('Aligning %s', apo_name)
            tcr_path = os.path.join(args.structures, apo_structure.file_name)
            align_tcr(tcr_path, apo_structure.alpha_chain, apo_structure.beta_chain, apo_name)
            cmd.save(os.path.join(group_path, apo_structure.file_name), apo_name)

        logger.info('Aliging %d pMHC(s)', len(apo_pmhcs))
        for _, apo_structure in apo_pmhcs.iterrows():
            apo_name = apo_structure.file_name.replace('.pdb', '')
            logger.debug('Aligning %s', apo_name)
            pmhc_path = os.path.join(args.structures, apo_structure.file_name)
            align_pmhc(pmhc_path, apo_structure.mhc_chain1, apo_name)
            cmd.save(os.path.join(group_path, apo_structure.file_name), apo_name)

        cmd.save(os.path.join(group_path, holo_structure.file_name), 'target')
        cmd.reinitialize()

    logger.info('Copying summary file')
    shutil.copy2(summary_path, args.output)


if __name__ == '__main__':
    main()
