import argparse
import glob
import logging
import os

import pandas as pd
from pymol import cmd

from tcr_pmhc_interface_analysis.apps._log import setup_logger
from tcr_pmhc_interface_analysis.imgt_numbering import IMGT_CDR, IMGT_VARIABLE_DOMAIN

logger = logging.getLogger()

parser = argparse.ArgumentParser()

parser.add_argument('structures', help='path to the structures to align')
parser.add_argument('--output', '-o', help='path to output the aligned files')
parser.add_argument('--only-holo', action='store_true',
                    help=('only align the holo structures based on either TCR CDR sequences'
                          ' or mhc allele and peptide sequence.'))
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

    if args.only_holo:
        for cdr_sequence, group in holo_structures.groupby('cdr_sequences_collated'):
            if len(group) == 1:
                logger.warning('Skipping CDR %s as there is only 1 holo form', cdr_sequence)
                continue

            logger.info('Aligning TCRs with %s CDRs', cdr_sequence)

            reference_info = group.iloc[0]
            logger.debug('Reference structure is %s', reference_info.file_name)

            cmd.load(os.path.join(args.structures, reference_info.file_name), 'target')
            framework_selection = ('target and '
                                   f'(chain {reference_info.alpha_chain} or chain {reference_info.beta_chain}) and '
                                   f'({get_framework_selection()})')
            cmd.select('target-Fw', framework_selection)

            group_path = os.path.join(args.output, cdr_sequence)
            if not os.path.exists(group_path):
                os.mkdir(group_path)

            for _, structure in group.iterrows():
                name = structure.file_name.replace('.pdb', '')
                logger.debug('Aligning %s', name)
                tcr_path = os.path.join(args.structures, structure.file_name)
                align_tcr(tcr_path, structure.alpha_chain, structure.beta_chain, name)
                cmd.save(os.path.join(group_path, structure.file_name), name)

            cmd.reinitialize()

        for (mhc_slug, peptide_sequence), group in holo_structures.groupby(['mhc_slug', 'peptide_sequence']):
            if len(group) == 1:
                logger.warning('Skipping MHC (%s) with %s as there is only 1 holo form', mhc_slug, peptide_sequence)
                continue

            logger.info('Aligning MHCs (%s) with %s', mhc_slug, peptide_sequence)

            reference_info = group.iloc[0]
            logger.debug('Reference structure is %s', reference_info.file_name)

            cmd.load(os.path.join(args.structures, reference_info.file_name), 'target')

            floor_selection = ('target and '
                               f'chain {reference_info.mhc_chain1} and '
                               f'({get_floor_selection()})')
            cmd.select('target-floor', floor_selection)

            group_path = os.path.join(args.output, mhc_slug + '_' + peptide_sequence)
            if not os.path.exists(group_path):
                os.mkdir(group_path)

            for _, structure in group.iterrows():
                name = structure.file_name.replace('.pdb', '')
                logger.debug('Aligning %s', name)
                pmhc_path = os.path.join(args.structures, structure.file_name)
                align_pmhc(pmhc_path, structure.mhc_chain1, name)
                cmd.save(os.path.join(group_path, structure.file_name), name)

            cmd.reinitialize()

    else:
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
    file_names_in_output = [file_name.split('/')[-1]
                            for file_name in glob.glob(os.path.join(args.output, '**/*.pdb'), recursive=True)]
    output_summary_df = summary_df[summary_df['file_name'].isin(file_names_in_output)]
    output_summary_df.to_csv(os.path.join(args.output, summary_path.split('/')[-1]), index=False)


if __name__ == '__main__':
    main()
