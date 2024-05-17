import logging
import re

import pandas as pd
from python_pdb.formats.residue import THREE_TO_ONE_CODE

logger = logging.getLogger(__name__)


def get_coords(df):
    return df[['pos_x', 'pos_y', 'pos_z']].to_numpy()


def get_sequence(df: pd.DataFrame) -> str:
    '''Get sequence from structre dataframe.'''
    df = df.copy()
    residues = df.drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code'])['residue_name']
    return ''.join(residues.map(lambda tlc: THREE_TO_ONE_CODE[tlc]).to_list())


def get_header(pdb_contents: str) -> str:
    '''Get the header lines from the contents of a pdb file.'''
    header = []
    for line in pdb_contents.split('\n'):
        record_type = line[0:6]
        if record_type in ("ATOM  ", "HETATM", "MODEL "):
            break

        header.append(line)

    return '\n'.join(header)


def mhc_code_to_slug(code: str) -> str:
    '''Convert MHC allele codes into slugs.

    >>> mhc_code_to_slug('HLA-A*02:01:59')
    'hla_a_02_01_59'

    '''
    slug = code.lower()
    slug = re.sub(r'[*:-]', '_', slug)

    return slug


def mhc_slug_to_code(slug: str) -> str:
    '''Convert mhc slugs into the allele codes.

    >>> mhc_slug_to_code('hla_a_02_01')
    'HLA-A*02:01'

    >>> mhc_slug_to_code('h2_kb')
    'H2-Kb'

    TODO: Make this better for other mouse alleles
    '''
    species = 'human' if slug.startswith('hla') else 'mouse' if slug.startswith('h2') else None

    if species == 'human':
        code = slug.upper()
        code = code.split('_')
        code = code[0] + '-' + code[1] + '*' + ':'.join(code[2:])

    elif species == 'mouse':
        code = slug.title()
        code = code.replace('_', '-')

    else:
        logger.error('Species not found, outputting slug')
        code = slug

    return code
