'''Download all of the STCRDab structures to a specified output path.'''
import argparse
import logging
import os
import sys

import requests

from tcr_pmhc_interface_analysis.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('output', help='path to the downloaded data directory eg. some/path/stcrdab')

add_logging_arguments(parser)

STCRDAB_BASE_URL = 'https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred'
STCRDAB_SUMMARY_FILE_URL = f'{STCRDAB_BASE_URL}/summary/all'
STCRDAB_IMGT_STRUCTURE_BASE_URL = f'{STCRDAB_BASE_URL}/pdb/%s'
STCRDAB_RAW_STRUCTURE_BASE_URL = f'{STCRDAB_IMGT_STRUCTURE_BASE_URL}?raw=true'


def main():
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    logger.info('Creating output dir')
    os.mkdir(args.output)

    logger.info('Downloading summary file')
    suumary_req = requests.get(STCRDAB_SUMMARY_FILE_URL)
    summary_file_contents = suumary_req.text

    with open(os.path.join(args.output, 'db_summary.dat'), 'w') as fh:
        fh.write(summary_file_contents)

    logger.info('Downloading PDB files')
    os.mkdir(os.path.join(args.output, 'imgt'))
    os.mkdir(os.path.join(args.output, 'raw'))

    pdb_ids = sorted(list(set([line.split('\t')[0] for line in summary_file_contents.strip().split('\n')[1:]])))

    for pdb_id in pdb_ids:
        logger.info('Downloading PDB ID: %s', pdb_id)

        for download_type, url in ('imgt', STCRDAB_IMGT_STRUCTURE_BASE_URL), ('raw', STCRDAB_RAW_STRUCTURE_BASE_URL):
            pdb_req = requests.get(url % pdb_id)

            with open(os.path.join(args.output, download_type, pdb_id + '.pdb'), 'w') as fh:
                fh.write(pdb_req.text)


if __name__ == '__main__':
    main()
