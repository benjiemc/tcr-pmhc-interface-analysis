import pandas as pd

from tcr_pmhc_interface_analysis.imgt_numbering import IMGT_MHC_ABD, assign_cdr_number


def annotate_tcr_pmhc_df(structure_df: pd.DataFrame,
                         alpha_chain_id: str = None,
                         beta_chain_id: str = None,
                         antigen_chain_id: str = None,
                         mhc_chain1_id: str = None,
                         mhc_chain2_id: str = None) -> pd.DataFrame:
    '''Add chain_type and cdr columns to a structure dataframe.'''
    def assign_chain_type(chain_id: str) -> str | None:
        if chain_id == alpha_chain_id:
            return 'alpha_chain'

        if chain_id == beta_chain_id:
            return 'beta_chain'

        if chain_id == antigen_chain_id:
            return 'antigen_chain'

        if chain_id == mhc_chain1_id:
            return 'mhc_chain1'

        if chain_id == mhc_chain2_id:
            return 'mhc_chain2'

        return None

    structure_df = structure_df.copy()
    structure_df['chain_type'] = structure_df['chain_id'].map(assign_chain_type)

    structure_df['cdr'] = structure_df.apply(
        lambda row: (assign_cdr_number(row.residue_seq_id)
                     if row.chain_type == 'alpha_chain' or row.chain_type == 'beta_chain' else None),
        axis=1,
    )

    structure_df['mhc_abd'] = structure_df.apply(
        lambda row: row.residue_seq_id in IMGT_MHC_ABD and row.chain_type == 'mhc_chain1',
        axis=1,
    )

    return structure_df


def find_anchors(cdr_df: pd.DataFrame,
                 structure_df: pd.DataFrame,
                 num_anchors: int = 1) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''Get the anchors of a cdr loop. DOES NOT SUPPORT MULTIPLE MODELS.'''
    cdr_start_info = tuple(cdr_df.iloc[0][['chain_id', 'residue_seq_id', 'residue_insert_code']].tolist())
    cdr_end_info = tuple(cdr_df.iloc[-1][['chain_id', 'residue_seq_id', 'residue_insert_code']].tolist())

    residues = list(structure_df.groupby(['chain_id', 'residue_seq_id', 'residue_insert_code'], dropna=False))

    for i, (residue_info, residue) in enumerate(residues):
        if pd.isna(residue_info[-1]):
            residue_info = (residue_info[0], residue_info[1], None)

        if cdr_start_info == residue_info:
            start_anchor = pd.concat([res for _, res in residues[i - num_anchors: i]])

        if cdr_end_info == residue_info:
            end_anchor = pd.concat([res for _, res in residues[(i + 1): (i + 1) + num_anchors]])
            break

    return start_anchor, end_anchor
