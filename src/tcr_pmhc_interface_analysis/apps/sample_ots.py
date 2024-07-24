'''Sample OTS sequences and select columns.'''
import argparse
import glob
import logging
import os

import pandas as pd

from tcr_pmhc_interface_analysis.apps._log import setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(description=__doc__)

DEFAULT_COLUMNS = ['cdr1_aa_alpha', 'cdr2_aa_alpha', 'cdr3_aa_alpha',
                   'cdr1_aa_beta', 'cdr2_aa_beta', 'cdr3_aa_beta',
                   'v_call_alpha', 'v_call_beta',
                   'j_call_alpha', 'j_call_beta']

parser.add_argument('path', help='path to OTS files')
parser.add_argument('--sample-size', '-n', type=int, default=1000, help='number of samples (default: 1000)')
parser.add_argument('--seed', default=None, type=int, help='Seed for sampling (Default: None)')
parser.add_argument('--redundant', action='store_true', help='do not remove redundant TCRs based on input columns')
parser.add_argument('--columns',
                    nargs='+',
                    default=DEFAULT_COLUMNS,
                    help=f"relevant columns to include (Default: {', '.join(DEFAULT_COLUMNS)})")
parser.add_argument('--output', '-o', required=True, help='path to output csv')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error'], default='warning',
                    help="Level to log messages at (Default: 'warning')")


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    ots = []
    for file_ in glob.glob(os.path.join(args.path, '*.csv*')):
        logger.info('Collecting %s', os.path.basename(file_))
        ots.append(pd.read_csv(file_, skiprows=1))

    ots = pd.concat(ots)

    if not args.redundant:
        logger.info('Removing redundant entries')
        ots = ots.drop_duplicates(args.columns)

    logger.info('Sampling %d entries', args.sample_size)
    ots_sample = ots.sample(args.sample_size, random_state=args.seed)

    logger.info('Writing output to %s', args.output)
    ots_sample[args.columns].to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
