'''Sample OTS sequences and select columns.'''
import argparse
import glob
import gzip
import json
import logging
import os
import sys

import numpy as np
import pandas as pd

from tcr_pmhc_interface_analysis.apps._log import add_logging_arguments, setup_logger

logger = logging.getLogger()

parser = argparse.ArgumentParser(prog=f'python -m {sys.modules[__name__].__spec__.name}',
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

DEFAULT_COLUMNS = ['cdr1_aa_alpha', 'cdr2_aa_alpha', 'cdr3_aa_alpha',
                   'cdr1_aa_beta', 'cdr2_aa_beta', 'cdr3_aa_beta',
                   'v_call_alpha', 'v_call_beta',
                   'j_call_alpha', 'j_call_beta']

parser.add_argument('path', help='path to OTS files')
parser.add_argument('--sample-size', type=int, default=1000, help='size of samples (Default: 1000)')
parser.add_argument('--num', '-n', type=int, default=1, help='number of samples to take (Default: 1)')
parser.add_argument('--seed', default=None, type=int, help='Seed for sampling (Default: None)')
parser.add_argument('--redundant', action='store_true', help='do not remove redundant TCRs based on input columns')
parser.add_argument('--columns',
                    nargs='+',
                    default=DEFAULT_COLUMNS,
                    help=f"relevant columns to include (Default: {', '.join(DEFAULT_COLUMNS)})")
parser.add_argument('--output', '-o', required=True, help='path to output csv')

add_logging_arguments(parser)


def main() -> None:
    args = parser.parse_args()
    setup_logger(logger, args.log_level)

    ots = []
    for file_ in sorted(glob.glob(os.path.join(args.path, '*.csv*')), key=os.path.basename):
        logger.info('Collecting %s', os.path.basename(file_))

        fh = gzip.open(file_, 'rt') if file_.endswith('.gz') else open(file_, 'r')
        meta = json.loads(fh.readline().replace('\"\"', '\"').strip('\"').rstrip().rstrip('\"'))
        fh.close()

        chunk_df = pd.read_csv(file_, skiprows=1)

        for key, value in meta.items():
            chunk_df[key] = value

        ots.append(chunk_df)

    ots = pd.concat(ots)

    if not args.redundant:
        logger.info('Removing redundant entries')
        ots = ots.drop_duplicates(args.columns)

    if args.seed:
        logger.info('Setting seed to %d', args.seed)
        np.random.seed(args.seed)

    samples = []
    for num in range(1, args.num + 1):
        logger.info('Taking sample %d of %d, sampling %d entries', num, args.num, args.sample_size)
        sample = ots.sample(args.sample_size)

        if args.num > 1:
            sample['sample_num'] = num

        samples.append(sample)

    ots_sample = pd.concat(samples)

    logger.info('Writing output to %s', args.output)
    ots_sample[args.columns + ['sample_num'] if args.num > 1 else args.columns].to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
