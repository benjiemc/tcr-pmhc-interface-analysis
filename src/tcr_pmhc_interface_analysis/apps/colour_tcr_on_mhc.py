import argparse
import logging
import re
import sys

import pandas as pd

from tcr_pmhc_interface_analysis.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

CDR_COLOURS = {
    'CDR-A1': 'violetpurple',
    'CDR-A2': 'lightmagenta',
    'CDR-A3': 'marine',
    'CDR-B1': 'slate',
    'CDR-B2': 'pink',
    'CDR-B3': 'cyan',
}


def pymol_session_file(arg_value, pat=re.compile(r'\.pse$')):
    if not pat.search(arg_value):
        raise argparse.ArgumentTypeError("invalid output file name, the pymol session must end in '.pse'")
    return arg_value


parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('mhc_path', help='path to the template MHC PDB file to colour with the contact information')
parser.add_argument('--output', '-o', required=True, type=pymol_session_file, help='name of output pymol session')
parser.add_argument('--contacts-path', required=True, help='path to the CSV file containg contacts information')
parser.add_argument('--mhc-chain-id', required=True, help='chain id of the main MHC chain')
parser.add_argument('--antigen-chain-id', required=True, help='chain id of the antigen')
parser.add_argument('--num-contacts-cutoff', type=int, default=100,
                    help=('number of contacts needed to highlight the MHC residue with contact information '
                          '(Default: 100)'))
parser.add_argument('--percentage-contacts-cutoff', type=float, default=None,
                    help=('A percentage cutoff instead of the raw count numbers, will override `--num-contacts-cutoff`'
                          ' (Default: None)'))
parser.add_argument('--dominant-peptide-contacts', required=False, nargs='+',
                    help='list of dominat CDRs for each peptide contact postion (CDR-A3, CDR-A3, ..., CDR-B3)')

add_logging_arguments(parser)


def select_dominant(group: pd.DataFrame) -> str:
    '''Select dominant CDR from counts.'''
    return group.sort_values('count', ascending=False).iloc[0]['cdr_name']


def aggregate_residues(resis: pd.Series) -> str:
    '''Create pymol selection for residue ids'''
    resis = [f'resi {resi}' for resi in resis.tolist()]
    return ' or '.join(resis)


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    try:
        from pymol import cmd

    except ImportError:
        logger.error('PyMOL not found. Please install: https://www.pymol.org/')
        sys.exit(1)

    cmd.load(args.mhc_path)

    cmd.hide('everything', f'chain {args.antigen_chain_id}')
    cmd.color('grey80', 'all')
    cmd.show('surface', f'chain {args.mhc_chain_id}' if args.dominant_peptide_contacts is not None else 'all')
    cmd.show('sticks', f'chain {args.mhc_chain_id} and not (name c,n)')
    cmd.set('transparency', 0.25, f'chain {args.mhc_chain_id}')

    contacts_count = pd.read_csv(args.contacts_path)

    if args.percentage_contacts_cutoff is not None:
        logger.info('Using percentage cutoff of %f', args.percentage_contacts_cutoff)
        percentages = (contacts_count.groupby(['resi_mhc'])['count'].transform('sum')
                       / contacts_count['count'].sum()
                       * 100)
        contacts_count = contacts_count[percentages > args.percentage_contacts_cutoff]

    else:
        logger.info('Using count cutoff of %d', args.num_contacts_cutoff)
        contacts_count = contacts_count[contacts_count['count'] >= args.num_contacts_cutoff]

    tcr_mhc_contacts = contacts_count.groupby('resi_mhc').apply(select_dominant)
    tcr_mhc_contacts.name = 'cdr_name'
    tcr_mhc_contacts = tcr_mhc_contacts.reset_index()

    tcr_mhc_contacts = tcr_mhc_contacts.groupby('cdr_name')['resi_mhc'].apply(aggregate_residues)

    for cdr_name, mhc_residues in tcr_mhc_contacts.items():
        cmd.select(cdr_name, f'chain {args.mhc_chain_id} and ({mhc_residues})')
        cmd.deselect()
        cmd.color(CDR_COLOURS[cdr_name], cdr_name)

    if args.dominant_peptide_contacts is not None:
        cmd.show('spheres', f'chain {args.antigen_chain_id}')
        for resi, cdr_name in enumerate(args.dominant_peptide_contacts, 1):
            cmd.color(CDR_COLOURS[cdr_name], f'chain {args.antigen_chain_id} and resi {resi}')

    cmd.save(args.output)


if __name__ == '__main__':
    main()
