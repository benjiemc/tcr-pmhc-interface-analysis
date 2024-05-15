import pandas as pd
from python_pdb.formats.residue import THREE_TO_ONE_CODE


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
