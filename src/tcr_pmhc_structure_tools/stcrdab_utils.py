import pandas as pd
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_structure_tools.processing import annotate_tcr_df
from tcr_pmhc_structure_tools.utils import get_sequence


def _clean_data_frame(df: pd.DataFrame) -> pd.DataFrame:
    '''Drop all columns that contain NAs or non-unique values.'''
    df = df.copy()

    df = df.loc[:, df.nunique() > 1]
    df = df.dropna(axis=1, how='all')
    df = df.reset_index(drop=True)

    return df


def get_ab_tcrs_from_stcrdab(stcrdab_summary: pd.DataFrame) -> pd.DataFrame:
    '''Get unbound alpha-beta TCRs from STCRDab.'''
    selected_stcrdab = stcrdab_summary.copy()

    selected_stcrdab = selected_stcrdab.query("TCRtype == 'abTCR'")
    selected_stcrdab = selected_stcrdab.query('mhc_type.isnull() and antigen_type.isnull()')

    selected_stcrdab = _clean_data_frame(selected_stcrdab)

    return selected_stcrdab


def get_ab_tcr_mhc_class_Is_from_stcrdab(stcrdab_summary: pd.DataFrame) -> pd.DataFrame:
    '''Get bound abTCR-MHC Class Is from STCRDab.'''
    selected_stcrdab = stcrdab_summary.copy()

    selected_stcrdab = selected_stcrdab.query("mhc_type == 'MH1'")
    selected_stcrdab = selected_stcrdab.query("antigen_type == 'peptide'")

    selected_stcrdab = _clean_data_frame(selected_stcrdab)

    return selected_stcrdab


def get_stcrdab_sequences(stcrdab_summary: pd.DataFrame, structure_type: str) -> pd.DataFrame:
    '''Get the CDR sequences for STCRDab structures.

    Args:
        stcrdab_summary: dataframe with stcrdab structures and file paths
        structure_type: either 'tcr' or 'tcr-pmhc'

    Returns:
        dataframe with additional sequence columns

    '''
    stcrdab_summary = stcrdab_summary.copy()

    sequences = {
        'alpha_chain': {1: [], 2: [], 3: []},
        'beta_chain': {1: [], 2: [], 3: []},
    }

    if structure_type == 'tcr-pmhc':
        sequences['peptide'] = []
        sequences['mhc_chain_1'] = []
        sequences['mhc_chain_2'] = []

    for _, stcrdab_entry in stcrdab_summary.iterrows():
        with open(stcrdab_entry['file_path_imgt'], 'r') as fh:
            structure_df = parse_pdb_to_pandas(fh.read())

        structure_df = structure_df.query("record_type == 'ATOM'")
        structure_df = annotate_tcr_df(structure_df, stcrdab_entry['alpha_chain'], stcrdab_entry['beta_chain'])

        for chain_type in 'alpha_chain', 'beta_chain':
            for cdr in 1, 2, 3:
                cdr_df = structure_df.query('chain_type == @chain_type and cdr == @cdr')
                sequences[chain_type][cdr].append(get_sequence(cdr_df))

        if structure_type == 'tcr-pmhc':
            sequences['peptide'].append(get_sequence(structure_df.query('chain_id == @stcrdab_entry.antigen_chain')))
            sequences['mhc_chain_1'].append(get_sequence(structure_df.query('chain_id == @stcrdab_entry.mhc_chain1')))
            sequences['mhc_chain_2'].append(get_sequence(structure_df.query('chain_id == @stcrdab_entry.mhc_chain2')))

    stcrdab_summary['cdr_1_alpha_seq'] = sequences['alpha_chain'][1]
    stcrdab_summary['cdr_2_alpha_seq'] = sequences['alpha_chain'][2]
    stcrdab_summary['cdr_3_alpha_seq'] = sequences['alpha_chain'][3]

    stcrdab_summary['cdr_1_beta_seq'] = sequences['beta_chain'][1]
    stcrdab_summary['cdr_2_beta_seq'] = sequences['beta_chain'][2]
    stcrdab_summary['cdr_3_beta_seq'] = sequences['beta_chain'][3]

    stcrdab_summary['cdr_sequences_collated'] = stcrdab_summary[['cdr_1_alpha_seq',
                                                                 'cdr_2_alpha_seq',
                                                                 'cdr_3_alpha_seq',
                                                                 'cdr_1_beta_seq',
                                                                 'cdr_2_beta_seq',
                                                                 'cdr_3_beta_seq']].apply(
        lambda sequences: '-'.join(sequences.dropna()),
        axis=1,
    )

    if structure_type == 'tcr-pmhc':
        stcrdab_summary['peptide_seq'] = sequences['peptide']
        stcrdab_summary['mhc_chain_1_seq'] = sequences['mhc_chain_1']
        stcrdab_summary['mhc_chain_2_seq'] = sequences['mhc_chain_2']

    return stcrdab_summary
