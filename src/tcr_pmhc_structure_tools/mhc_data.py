import pandas as pd
import requests

HISTO_DATASETS_URL='https://api.histo.fyi/v1/sets'
TCR_PMHC_CLASS_I_URL = f'{HISTO_DATASETS_URL}/complex_types/class_i_with_peptide_and_alpha_beta_tcr'
PMHC_CLASS_I_URL = f'{HISTO_DATASETS_URL}/complex_types/class_i_with_peptide'


def _retrieve_data_from_api(url: str) -> pd.DataFrame:
    req = requests.get(url)
    
    pdb_ids = []
    peptide_sequences = []
    mhc_slugs = []
    chains = []
    resolutions = []
    
    for member in req.json()['set']['members']:
        pdb_ids.append(member['pdb_code'])
        peptide_sequences.append(member['peptide_sequence'])
        mhc_slugs.append(member['allele']['alpha']['slug'])
        chains.append([''.join(assembly['chains']) for assembly in member['assemblies'].values()])
        resolutions.append(member['resolution'])
    
    return pd.DataFrame({
        'pdb_id': pdb_ids,
        'peptide_sequence': peptide_sequences,
        'mhc_slug': mhc_slugs,
        'chains': chains,
    }).explode('chains')


def get_apo_mhc_ids(resolution_cutoff):
    '''Get the apo MHCs from a list of TCR-pMHC pdb ids.
    
    Data is collected from histo.fyi.
    
    Args:
        holo_pdb_ids: pdb ids of TCR-pMHC complexes optional
    
    Returns:
        dataframe with holo tcr-pmhc pdb_ids mapped to apo pdb_ids and chains
        
        Eg:
                    pdb_id_holo peptide_sequence     mhc_slug          pdb_id_apo  chains_apo
            0        7nme        QLPRLFPLL            hla_a_24_02       7nmd        ABC
            0        7nme        QLPRLFPLL            hla_a_24_02       7nmd        DEF
            1        7nmf        QLPRLFPLL            hla_a_24_02       7nmd        ABC
            1        7nmf        QLPRLFPLL            hla_a_24_02       7nmd        DEF
            2        7n5p       SSLCNFRAYV            h2_db             7n5q        ABC

    '''
    tcr_pmhc_data = _retrieve_data_from_api(TCR_PMHC_CLASS_I_URL)
    pmhc_data = _retrieve_data_from_api(PMHC_CLASS_I_URL)
    
    apo_holo_mhcs = tcr_pmhc_data.merge(pmhc_data,
                                        how='inner',
                                        on=['mhc_slug', 'peptide_sequence'],
                                        suffixes=('_holo', '_apo'))
    
    return apo_holo_mhcs[['pdb_id_holo', 'peptide_sequence', 'mhc_slug', 'pdb_id_apo', 'chains_apo']]