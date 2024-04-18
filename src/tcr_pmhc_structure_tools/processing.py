import pandas as pd

from tcr_pmhc_structure_tools.imgt_numbering import assign_cdr_number


def annotate_tcr_df(structure_df: pd.DataFrame, alpha_chain_id: str, beta_chain_id: str) -> pd.DataFrame:
    '''Add chain_type and cdr columns to a structure dataframe.'''
    def assign_chain_type(chain_id: str) -> str | None:
        if chain_id == alpha_chain_id:
            return 'alpha_chain'

        if chain_id == beta_chain_id:
            return 'beta_chain'

        return None

    structure_df = structure_df.copy()
    structure_df['chain_type'] = structure_df['chain_id'].map(assign_chain_type)

    structure_df['cdr'] = structure_df.apply(
        lambda row: assign_cdr_number(row.residue_seq_id) if pd.notnull(row.chain_type) else None,
        axis=1,
    )
    
    return structure_df
