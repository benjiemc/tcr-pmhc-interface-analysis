import glob
import logging
import os
import re

import numpy as np
import pandas as pd
import requests
from python_pdb.aligners import align_sequences
from python_pdb.parsers import parse_pdb_to_pandas

from tcr_pmhc_structure_tools.imgt_numbering import IMGT_CDR
from tcr_pmhc_structure_tools.utils import get_sequence, get_header

logger = logging.getLogger(__name__)


def get_missing_residues(header: str) -> list[dict]:
    '''
    Extract all missing residues from the header of a pdb file. Returns residue name, chain id, residue/sequence id,
    and insert id if available.
    '''
    lines = re.findall(r'^REMARK 465\s+\w?\s+(\w{3}) (\w)\s+(\d+)(\w)?', header, flags=re.MULTILINE)
    return [{'residue_name': res_name,
             'chain_id': chain,
             'residue_seq_id': int(seq_id),
             'residue_insert_code': insert_id} for res_name, chain, seq_id, insert_id in lines]


def get_missing_atoms(header: str) -> list[dict]:
    '''Extract residues with missing atoms from the header of a pdb file.'''
    lines = re.findall(r'^REMARK 470\s+\w?\s+(\w{3}) (\w)\s+(\d+)(\w)?', header, flags=re.MULTILINE)
    return [{'residue_name': res_name,
             'chain_id': chain,
             'residue_seq_id': int(seq_id),
             'residue_insert_code': insert_id} for res_name, chain, seq_id, insert_id in lines]


def screen_variable(chain: pd.DataFrame, raw_chain: pd.DataFrame) -> bool:
    '''Screen variable domains for missing residues'''
    raw_seq = get_sequence(raw_chain)
    seq = get_sequence(chain)
    alignment, _ = align_sequences(raw_seq, seq)

    raw_index = 0
    index = 0
    current_seq_id = 0

    for raw_res, res in alignment:
        if raw_res == '-':
            index += 1
            continue

        if res != '-':
            current_seq_id = chain.iloc[index]['residue_seq_id']

        if raw_chain.iloc[raw_index]['missing'] is True and current_seq_id in IMGT_CDR:
            return False

        if res == '-':
            raw_index += 1
            continue

        index += 1
        raw_index += 1

    return True


def get_missing_residues_and_atoms(header: str) -> pd.DataFrame:
    '''Create a dataframe of missing information for a structure.'''
    missing_residues = pd.DataFrame(get_missing_residues(header))
    missing_residues['missing'] = True
    missing_residues['missing'] = missing_residues['missing'].astype(bool)

    residues_missing_atoms = pd.DataFrame(get_missing_atoms(header))
    residues_missing_atoms['missing'] = True
    residues_missing_atoms['missing'] = residues_missing_atoms['missing'].astype(bool)

    return pd.concat([missing_residues, residues_missing_atoms]).reset_index(drop=True)


def add_missing_entities_to_structure(structure: pd.DataFrame, missing_entities: pd.DataFrame) -> pd.DataFrame:
    '''Create a new dataframe with line indicating missing entities.'''
    updated_structure = structure.copy()
    updated_structure['missing'] = False
    updated_structure['missing'] = updated_structure['missing'].astype(bool)

    for _, row in missing_entities.iterrows():
        chain_before = updated_structure.query("chain_id == @row.chain_id and residue_seq_id <= @row.residue_seq_id")
        chain_after = updated_structure.query("chain_id == @row.chain_id and residue_seq_id > @row.residue_seq_id")

        new_chain = pd.concat([chain_before, row.to_frame().T, chain_after])
        updated_structure = pd.concat([structure.query('chain_id != @row.chain_id'), new_chain]).reset_index(drop=True)

    return updated_structure


def screen_tcrs_for_missing_residues(df: pd.DataFrame, raw_structure_dfs: dict[pd.DataFrame]) -> np.ndarray:
    '''Disgard structures missing residues in TCR variable domains.'''
    valid_structures = []

    for _, entry in df.iterrows():
        logger.debug('Checking PDB ID %s', entry.pdb_id)

        # 1. check if there are any missing items
        raw_structure = raw_structure_dfs[entry.pdb_id]

        if len(raw_structure.query('missing == True')) == 0:
            logger.debug('No missing residues, adding to selection')
            valid_structures.append(True)
            continue

        if pd.isnull(entry['file_path_imgt']):
            logger.warning('Cannot check pdb id %s as there is no IMGT information available', entry.pdb_id)
            valid_structures.append(True)
            continue

        with open(entry['file_path_imgt'], 'r') as fh:
            structure = parse_pdb_to_pandas(fh.read())

        # 2. look if they are in TCR variable domains
        # 2a. alpha chain
        logger.debug('looking at alpha chain')
        alpha = (structure.query("record_type == 'ATOM' and chain_id == @entry.alpha_chain")
                          .drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())

        raw_alpha = (raw_structure.query("record_type == 'ATOM' and chain_id == @entry.alpha_chain")
                                  .drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())

        if not screen_variable(alpha, raw_alpha):
            logger.info('check failed, missing residue found in alpha chain of pdb id %s', entry.pdb_id)
            valid_structures.append(False)
            continue

        # 2b. beta chain
        logger.debug('looking at beta chain')
        beta = (structure.query("record_type == 'ATOM' and chain_id == @entry.beta_chain")
                         .drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())

        raw_beta = (raw_structure.query("record_type == 'ATOM' and chain_id == @entry.beta_chain")
                                 .drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())

        if not screen_variable(beta, raw_beta):
            logger.info('check failed, missing residue found in beta chain of pdb id %s', entry.pdb_id)
            valid_structures.append(False)
            continue

        logger.debug('all clear, adding to selection')
        valid_structures.append(True)

    return np.array(valid_structures, dtype=bool)


def screen_pmhcs_for_missing_residues(df: pd.DataFrame, raw_structures: dict[pd.DataFrame]) -> np.ndarray:
    '''Screen pmhcs for missing residues.'''
    valid_structures = []
    for _, entry in df.iterrows():
        raw_structure = raw_structures[entry.pdb_id]
        raw_peptide = raw_structure.query('chain_id == @entry.antigen_chain')

        if len(raw_peptide.query("missing == True")) > 0:
            logger.info('check failed, missing residue found in pdb id %s', entry.pdb_id)
            valid_structures.append(False)

        else:
            valid_structures.append(True)

    return np.array(valid_structures, dtype=bool)


def get_raw_structures_with_missing_residues(pdb_ids: list[str], stcrdab_path: str | None = None) -> dict[pd.DataFrame]:
    '''Get all raw structure files from the STCRDab or RCSB PDB.'''
    raw_structure_dfs = {}

    if stcrdab_path:
        stcrdab_raw_path = os.path.join(stcrdab_path, 'raw')
        local_pdb_ids = [path.rsplit('/', 1)[-1].replace('.pdb', '')
                         for path in glob.glob(os.path.join(stcrdab_raw_path, '*.pdb'))]

    else:
        local_pdb_ids = []

    for pdb_id in pdb_ids:
        if pdb_id in local_pdb_ids:
            with open(os.path.join(stcrdab_raw_path, f'{pdb_id}.pdb'), 'r') as fh:
                pdb_contents = fh.read()

        else:
            req = requests.get(f'https://files.rcsb.org/download/{pdb_id}.pdb')
            pdb_contents = req.text

        structure = parse_pdb_to_pandas(pdb_contents)
        header = get_header(pdb_contents)

        missing_entities = get_missing_residues_and_atoms(header)
        structure = add_missing_entities_to_structure(structure, missing_entities)

        raw_structure_dfs[pdb_id] = structure

    return raw_structure_dfs
