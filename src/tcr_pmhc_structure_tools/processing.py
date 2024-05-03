import pandas as pd

from tcr_pmhc_structure_tools.imgt_numbering import assign_cdr_number


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

    return structure_df
