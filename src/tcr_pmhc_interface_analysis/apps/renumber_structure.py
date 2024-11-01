'''Renumber TCR structure following either IMGT or Aho numbering.

Requirements:
    - ANARCI: https://github.com/oxpig/ANARCI

'''
import argparse
import logging
import sys

from python_pdb.formats.residue import THREE_TO_ONE_CODE
from python_pdb.parsers import parse_pdb, stringify_structure

from tcr_pmhc_interface_analysis.apps._log import add_logging_arguments, setup_logger
from tcr_pmhc_interface_analysis.utils import get_header

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('structure', help='path to the pdb structure file')
parser.add_argument('--output', '-o', help='name of output structure file')

add_logging_arguments(parser)


def main():
    '''Entry point for script.'''
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    try:
        import anarci

    except ImportError:
        logger.error('Anarci not found. Please follow install instructions here: https://github.com/oxpig/ANARCI')
        sys.exit(1)

    with open(args.structure, 'r') as fh:
        pdb_contents = fh.read()

    structure = parse_pdb(pdb_contents, silent=True)
    header = get_header(pdb_contents)

    for model in structure:
        sequences = {chain.name: [THREE_TO_ONE_CODE[res.name] for res in chain if res.name in THREE_TO_ONE_CODE]
                     for chain in model}
        sequences = {chain: ''.join(sequence) for chain, sequence in sequences.items()}

        chain_map = []

        for chain_id, sequence in sequences.items():
            numbering, chain_type = anarci.number(sequence)

            if not numbering:
                logger.info('Chain ID %s not identified by ANARCI', chain_id)
                continue

            if len(numbering) == 2:
                numbering = numbering[0]
                chain_type = chain_type[0]

                logger.warning(('Multiple possible chain annotations found for chain id %s. '
                                'Defaulting to first (%sCHAIN)'), chain_id, chain_type)

            numbering = [(seq_id, insert_code) for (seq_id, insert_code), res_name in numbering if res_name != '-']

            residues = [res for res in model[chain_id].get_residues()]
            num_residues_not_numbered = len(residues) - len(numbering)

            next_seq_id = numbering[-1][0] + 1

            for _ in range(num_residues_not_numbered):
                numbering.append((next_seq_id, ' '))
                next_seq_id += 1

            for (seq_id, insert_code), res in zip(numbering, residues):
                res.seq_id = seq_id
                res.insert_code = insert_code if insert_code.strip() != '' else None

            chain_map.append((chain_type, chain_id))

    with open(args.output, 'w') as fh:
        fh.write(('REMARK     Renumbered using IMGT numbering provided by ANARCI '
                  '(DOI: 10.1093/bioinformatics/btv552)'))
        fh.write('\n')
        fh.write(f"REMARK     {' '.join([f'{chain_type}CHAIN={chain_id}' for chain_type, chain_id in chain_map])}")
        fh.write('\n')

        if len(header) > 0:
            fh.write(header)
            fh.write('\n')

        fh.write(stringify_structure(structure))


if __name__ == '__main__':
    main()
