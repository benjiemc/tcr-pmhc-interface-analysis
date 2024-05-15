import argparse
import csv
import logging
import re

from pymol import cmd

from tcr_pmhc_interface_analysis.apps._log import setup_logger

logger = logging.getLogger()

CDR_COLOURS = {
    'CDR1alpha': 'slate',
    'CDR2alpha': 'lightpink',
    'CDR3alpha': 'marine',
    'CDR1beta': 'yellow',
    'CDR2beta': 'purple',
    'CDR3beta': 'cyan',
}


def pymol_session_file(arg_value, pat=re.compile(r'\.pse$')):
    if not pat.search(arg_value):
        raise argparse.ArgumentTypeError("invalid output file name, the pymol session must end in '.pse'")
    return arg_value


parser = argparse.ArgumentParser()

parser.add_argument('mhc_path', help='path to the template MHC PDB file to colour with the contact information')
parser.add_argument('--output', '-o', required=True, type=pymol_session_file, help='name of output pymol session')
parser.add_argument('--contacts-path', required=True, help='path to the CSV file containg contacts information')
parser.add_argument('--num-contacts-cutoff', type=int, default=100,
                    help='number of contacts needed to highlight the MHC residue with contact information')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error'], default='warning',
                    help="Level to log messages at (Default: 'warning')")


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    cmd.load(args.mhc_path)

    cmd.hide('everything', 'chain C')
    cmd.color('grey80', 'all')
    cmd.show('surface')
    cmd.show('sticks', 'all and not (name c,n)')
    cmd.set('transparency', 0.25, 'all')

    with open(args.contacts_path, 'r') as fh:
        reader = csv.reader(fh)
        next(reader)  # skip header

        cdr_mhc_contacts = {}

        for _, cdr_name, mhc_residue, count in reader:
            if int(count) >= args.num_contacts_cutoff:
                if cdr_name not in cdr_mhc_contacts:
                    cdr_mhc_contacts[cdr_name] = f'resi {mhc_residue}'

                else:
                    cdr_mhc_contacts[cdr_name] += f' or resi {mhc_residue}'

    for cdr_name, mhc_residues in cdr_mhc_contacts.items():
        cmd.select(cdr_name, f'chain A and ({mhc_residues})')
        cmd.deselect()
        cmd.color(CDR_COLOURS[cdr_name], cdr_name)

    cmd.save(args.output)


if __name__ == '__main__':
    main()
