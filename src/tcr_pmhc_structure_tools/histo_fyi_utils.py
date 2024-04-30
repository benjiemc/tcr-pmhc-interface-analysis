import logging

import pandas as pd
import requests

logger = logging.getLogger(__name__)

HISTO_DATASETS_URL = 'https://api.histo.fyi/v1/sets'
TCR_PMHC_CLASS_I_URL = f'{HISTO_DATASETS_URL}/complex_types/class_i_with_peptide_and_alpha_beta_tcr'
PMHC_CLASS_I_URL = f'{HISTO_DATASETS_URL}/complex_types/class_i_with_peptide'
HISTO_STRUCTURE_BASE_URL = "https://coordinates.histo.fyi/structures/downloads/class_i/without_solvent"


def retrieve_data_from_api(url: str) -> pd.DataFrame:
    '''Get pMHC from API end point.

    Args:
        url: Either 'TCR_PMHC_CLASS_I_URL' or 'PMHC_CLASS_I_URL'.

    Returns:
         dataframe with (tcr-)pmhc pdb_ids, peptide sequences, mhc allel information, and complex ids.

        Eg:
                pdb_id     peptide_sequence mhc_slug            chains     resolution
            0   7s7f        LPASPAHQL        hla_b_07_02         ABCD       1.88
            1   7s8f        EPRSPSHSM        hla_b_07_02          ABC       1.80
            2   7s8e        EPRSPSHSM        hla_b_07_02          ABE       1.60
            3   7rzd        EPRSPSHSM        hla_b_07_02          ABC       1.82
            4   7s7d        EPRSPSHSM        hla_b_07_02          ABE       1.56
    '''
    req = requests.get(url)
    pages = req.json()['set']['pagination']['pages']

    pdb_ids = []
    peptide_sequences = []
    mhc_slugs = []
    chains = []
    antigen_chains = []
    mhc_chain1s = []
    mhc_chain2s = []
    assembly_numbers = []
    resolutions = []

    for page_number in pages:
        req = requests.get(f'{url}?page_number={page_number}')

        for member in req.json()['set']['members']:
            number_of_antigen_chains = len(member['assigned_chains']['peptide']['chains'])
            number_of_mhc_chain1s = len(member['assigned_chains']['class_i_alpha']['chains'])
            number_of_mhc_chain2s = len(member['assigned_chains']['beta2m']['chains'])
            number_of_assemblies = len(member['assemblies'])

            if len({number_of_antigen_chains, number_of_mhc_chain1s, number_of_mhc_chain2s, number_of_assemblies}) != 1:
                logger.warning('Skipping %s due to inconsistencies in chain annotations', member['pdb_code'])
                continue

            pdb_ids.append(member['pdb_code'])
            peptide_sequences.append(member['assigned_chains']['peptide']['sequence'])
            mhc_slugs.append(member['allele']['alpha']['slug'])
            chains.append(['-'.join(assembly['chains']) for _, assembly in member['assemblies'].items()])
            assembly_numbers.append([number for number, _ in member['assemblies'].items()])
            antigen_chains.append(member['assigned_chains']['peptide']['chains'])
            mhc_chain1s.append(member['assigned_chains']['class_i_alpha']['chains'])
            mhc_chain2s.append(member['assigned_chains']['beta2m']['chains'])
            resolutions.append(float(member['resolution']) if member['resolution'] else None)

    return pd.DataFrame({
        'pdb_id': pdb_ids,
        'peptide_sequence': peptide_sequences,
        'mhc_slug': mhc_slugs,
        'antigen_chain': antigen_chains,
        'mhc_chain1': mhc_chain1s,
        'mhc_chain2': mhc_chain2s,
        'chains': chains,
        'assembly_number': assembly_numbers,
        'resolution': resolutions,
    }).explode(['chains', 'antigen_chain', 'mhc_chain1', 'mhc_chain2', 'assembly_number']).reset_index()


def fetch_structure(pdb_id: str, assembly_number: int = 1, domain: str = 'all') -> str:
    '''Fetch a structure from the histo.fyi api

    Args:
        pdb_id: pdb id of the structure
        assembly_number: number if there are multiple of the same structure in the pdb file (Default: 1)
        domain: can be either 'abd' or 'peptide' or 'all' (Default: 'all')

    Returns:
        pdb file text

    Raises:
        ValueError: if a non-valid domain is given

    '''
    match domain:
        case 'peptide' | 'abd':
            req = requests.get(f'{HISTO_STRUCTURE_BASE_URL}/{pdb_id}_{assembly_number}_{domain}.pdb')
            return req.text

        case 'all':
            domains = [requests.get(f'{HISTO_STRUCTURE_BASE_URL}/{pdb_id}_{assembly_number}_{domain}.pdb').text
                       for domain in ('abd', 'peptide')]
            return '\n'.join(domains)

        case _:
            raise ValueError(f"Domain: {domain}, is not a valid selection. Use either 'peptide', 'abd' or 'all'.")
