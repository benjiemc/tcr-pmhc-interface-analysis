import numpy as np
import pandas as pd
from python_pdb.aligners import align_pandas_structure, align_sequences
from python_pdb.formats.residue import THREE_TO_ONE_CODE

from tcr_pmhc_structure_tools.imgt_numbering import IMGT_VARIABLE_DOMAIN  # noqa: F401


def align_tcrs(tcr_mobile_df: pd.DataFrame, tcr_target_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Align two dataframes containing TCR structures in a pandas dataframe, the structures are algined on the Fw region.
    '''
    mobile_coords = []
    target_coords = []

    for chain_type in ('alpha_chain', 'beta_chain'):
        fw_chain_sequences = []
        fw_chain_ca_coords = []

        for df in (tcr_mobile_df, tcr_target_df):
            fw_chain = df.query(('cdr.isnull() '
                                 'and residue_seq_id in @IMGT_VARIABLE_DOMAIN '
                                 'and chain_type == @chain_type'))

            fw_chain_ca_coords.append(fw_chain.query("atom_name == 'CA'")[['pos_x', 'pos_y', 'pos_z']].values)

            fw_chain_seq = fw_chain.drop_duplicates(
                ['residue_seq_id', 'residue_insert_code']
            )['residue_name'].map(lambda tlc: THREE_TO_ONE_CODE[tlc]).tolist()
            fw_chain_sequences.append(fw_chain_seq)

        fw_chain_seq_alignment, _ = align_sequences(*fw_chain_sequences)

        iter_mobile_ca_coords = iter(fw_chain_ca_coords[0])
        iter_target_ca_coords = iter(fw_chain_ca_coords[1])

        for res_id_mobile, res_id_target in fw_chain_seq_alignment:
            next_res_ca_coords_mobile = next(iter_mobile_ca_coords) if res_id_mobile != '-' else None
            next_res_ca_coords_target = next(iter_target_ca_coords) if res_id_target != '-' else None

            if next_res_ca_coords_mobile is not None and next_res_ca_coords_target is not None:
                mobile_coords.append(next_res_ca_coords_mobile)
                target_coords.append(next_res_ca_coords_target)

    mobile_coords = np.array(mobile_coords)
    target_coords = np.array(target_coords)

    return align_pandas_structure(mobile_coords, target_coords, tcr_mobile_df)
