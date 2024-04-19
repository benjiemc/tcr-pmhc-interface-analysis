import argparse
import json
import os
import logging
import re

import pandas as pd
from python_pdb.aligners import align_sequences
from python_pdb.entities import Structure
from python_pdb.parsers import parse_pdb_to_pandas
from python_pdb.formats.residue import THREE_TO_ONE_CODE

from tcr_pmhc_structure_tools.apps._log import setup_logger
from tcr_pmhc_structure_tools.processing import annotate_tcr_df
from tcr_pmhc_structure_tools.imgt_numbering import IMGT_VARIABLE_DOMAIN

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser()
parser.add_argument('stcrdab', help='Path to STCRDab')
parser.add_argument('--resolution-cutoff', type=float, default=3.50,
                    help='Maximum resolution allowed from the structures')
parser.add_argument('--output', '-o', help='Path to output location')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error'], default='info',
                    help='Level to log at')


def get_header(pdb_contents: str) -> str:
	'''Get the header lines from the contents of a pdb file.'''
	header = []
	for line in pdb_contents.split('\n'):
		record_type = line[0:6]
		if record_type in ("ATOM  ", "HETATM", "MODEL "):
			break

		header.append(line)

	return '\n'.join(header)


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


def get_sequence(df: pd.DataFrame) -> str:
		return ''.join(df['residue_name'].map(lambda tlc: THREE_TO_ONE_CODE[tlc]).to_list())


def screen_variable(chain, raw_chain):
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
			
		if raw_chain.iloc[raw_index]['missing'] and current_seq_id in IMGT_VARIABLE_DOMAIN:
			return False

		if res == '-':
			raw_index += 1
			continue

		index +=1
		raw_index += 1
	
	return True


def select_tcrs(stcrdab_summary: pd.DataFrame, resolution_cutoff: float) -> pd.DataFrame:
	'''Select TCRs and TCR-pMHCS from the STCRDab'''
	selected_stcrdab = stcrdab_summary.copy()

	selected_stcrdab['resolution'] = pd.to_numeric(selected_stcrdab['resolution'], errors='coerce')
	selected_stcrdab = selected_stcrdab.query("resolution <= @resolution_cutoff")

	# alpha-beta TCRs
	selected_stcrdab = selected_stcrdab.query("TCRtype == 'abTCR'")

	# MHC class I bound or unbound
	selected_stcrdab = selected_stcrdab.query("(mhc_type == 'MH1') | (mhc_type.isnull() & antigen_type.isnull())")

	# peptide antigen
	selected_stcrdab = selected_stcrdab.query("antigen_type == 'peptide' or antigen_type.isnull()")

	# General clean: drop columns that don't contain anything useful
	selected_stcrdab = selected_stcrdab.loc[:, selected_stcrdab.nunique() > 1]
	selected_stcrdab = selected_stcrdab.dropna(axis=1, how='all')

	# Reset Index
	selected_stcrdab = selected_stcrdab.reset_index(drop=True)

	return selected_stcrdab


def get_missing_residues_and_atoms(header: str) -> pd.DataFrame:
	'''Create a dataframe of missing information for a structure.'''
	missing_residues = pd.DataFrame(get_missing_residues(header))
	missing_residues['missing'] = True
	
	residues_missing_atoms = pd.DataFrame(get_missing_atoms(header))
	residues_missing_atoms['missing'] = True
	
	return pd.concat([missing_residues, residues_missing_atoms]).reset_index(drop=True)


def add_missing_entities_to_structure(structure: pd.DataFrame, missing_entities: pd.DataFrame) -> pd.DataFame:
	'''Create a new dataframe with line indicating missing entities.'''
	updated_structure = structure.copy()
	updated_structure['missing'] = False

	for _, row in missing_entities.iterrows():
		chain_before = updated_structure.query("chain_id == @row.chain_id and residue_seq_id <= @row.residue_seq_id")
		chain_after = updated_structure.query("chain_id == @row.chain_id and residue_seq_id > @row.residue_seq_id")

		new_chain = pd.concat([chain_before, row.to_frame().T, chain_after])
		updated_structure = pd.concat([structure.query('chain_id != @row.chain_id'), new_chain]).reset_index(drop=True)

	return updated_structure


def screen_tcrs_for_missing_residues(stcrdab_summary: pd.DataFrame) -> pd.DataFrame:
	'''Disgard structures missing residues in important domains
 
	These domains include:
	  - TCR variable domains
	  - peptide

 	'''
	raw_structure_dfs = {}

	for _, pdb_id, path in stcrdab_summary[['pdb', 'file_path_raw']].drop_duplicates().itertuples():
		with open(path, 'r') as fh:
			pdb_contents = fh.read()
			
		structure = parse_pdb_to_pandas(pdb_contents)
		header = get_header(pdb_contents)
		
		missing_entities = get_missing_residues_and_atoms(header)
		structure = add_missing_entities_to_structure(structure, missing_entities)

		raw_structure_dfs[pdb_id] = structure
		
	selected_entries = []
	for _, entry in stcrdab_summary.iterrows():
		logger.debug('Checking PDB ID %s', entry.pdb)
		
		# 1. check if there are any missing items
		raw_structure = raw_structure_dfs[entry.pdb]
		
		if len(raw_structure.query('missing == True')) == 0:
			logger.debug('No missing residues, adding to selection')
			selected_entries.append(entry.to_frame().T)
			continue
		
		with open(entry['file_path_imgt'], 'r') as fh:
			structure = parse_pdb_to_pandas(fh.read())

		# 2. look if they are in TCR variable domains
		# 2a. alpha chain
		logger.debug('looking at alpha chain')
		alpha = (structure.query("record_type == 'ATOM' and chain_id == @entry.Achain")
						  .drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())
		
		raw_alpha = (raw_structure.query("record_type == 'ATOM' and chain_id == @entry.Achain")
								.drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())
		
		if not screen_variable(alpha, raw_alpha):
			logger.warning('check failed, missing residue found')
			continue
		
		# 2b. beta chain
		logger.debug('looking at beta chain')
		beta = (structure.query("record_type == 'ATOM' and chain_id == @entry.Bchain")
						.drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())
		
		raw_beta = (raw_structure.query("record_type == 'ATOM' and chain_id == @entry.Bchain")
								.drop_duplicates(['chain_id', 'residue_seq_id', 'residue_insert_code']).reset_index())
		
		if not screen_variable(beta, raw_beta):
			logger.warning('check failed, missing residue found')
			continue
		
		# 3. look if they are in the peptide
		logger.debug('looking at peptide')
		raw_peptide = raw_structure.query('chain_id == @entry.antigen_chain')
		
		if len(raw_peptide.query("missing == True")) > 0:
			logger.warning('check failed, missing residue found')
			continue
		
		logger.debug('all clear, adding to selection')
		selected_entries.append(entry.to_frame().T)

	return pd.concat(selected_entries)


def get_stcrdab_sequences(selected_stcrdab: pd.DataFrame) -> pd.DataFrame:
	'''Get the CDR sequences for STCRDab structures.'''
	selected_stcrdab = selected_stcrdab.copy()

	sequences = {
		'alpha_chain': {1: [], 2: [], 3: []},
		'beta_chain': {1: [], 2: [], 3: []},
		'peptide': [],
		'mhc_chain_1': [],
		'mhc_chain_2': [],
	}

	for _, stcrdab_entry in selected_stcrdab.iterrows():
		with open(stcrdab_entry['file_path_imgt'], 'r') as fh:
				structure_df = parse_pdb_to_pandas(fh.read())

		structure_df = structure_df.query("atom_type == 'ATOM'")
		structure_df = annotate_tcr_df(structure_df, stcrdab_entry['Achain'], stcrdab_entry['Bchain'])
  
		for chain_type in 'alpha_chain', 'beta_chain':
			for cdr in 1, 2, 3:
				cdr_df = structure_df.query('chain_type == @chain_type and cdr == @cdr')
				sequences[chain_type][cdr].append(get_sequence(cdr_df))
    
		if stcrdab_entry['state'] == 'holo':
			sequences['peptide'].append(get_sequence(structure_df.query('chain_id == @stcrdab_entry.antigen_chain')))
			sequences['mhc_chain_1'].append(get_sequence(structure_df.query('chain_id == @stcrdab_entry.mhc_chain1')))
			sequences['mhc_chain_2'].append(get_sequence(structure_df.query('chain_id == @stcrdab_entry.mhc_chain2')))

		else:   
			sequences['peptide'].append(None)
			sequences['mhc_chain_1'].append(None)
			sequences['mhc_chain_2'].append(None)

	selected_stcrdab['cdr_1_alpha_seq'] = sequences['alpha_chain'][1]
	selected_stcrdab['cdr_2_alpha_seq'] = sequences['alpha_chain'][2]
	selected_stcrdab['cdr_3_alpha_seq'] = sequences['alpha_chain'][3]

	selected_stcrdab['cdr_1_beta_seq'] = sequences['beta_chain'][1]
	selected_stcrdab['cdr_2_beta_seq'] = sequences['beta_chain'][2]
	selected_stcrdab['cdr_3_beta_seq'] = sequences['beta_chain'][3]

	selected_stcrdab['peptide_seq'] = sequences['peptide']

	selected_stcrdab['mhc_chain_1_seq'] = sequences['mhc_chain_1']
	selected_stcrdab['mhc_chain_2_seq'] = sequences['mhc_chain_2']
 
	return selected_stcrdab


def select_apo_holo(selected_stcrdab: pd.DataFrame) -> pd.DataFrame:
	'''Select groups with both *apo* and *holo* confomations'''
	apo_holo_dfs = []

	for sequences, tcr_group in selected_stcrdab.groupby(['cdr_1_alpha_seq',
                                                          'cdr_2_alpha_seq',
                                                          'cdr_3_alpha_seq',
                                                          'cdr_1_beta_seq',
                                                          'cdr_2_beta_seq',
                                                          'cdr_3_beta_seq']):
		# Screen out groups that don't have apo and holo forms
		if 'holo' not in tcr_group['state'].unique().tolist() or 'apo' not in tcr_group['state'].unique().tolist():
			continue
			
		tcr_group = tcr_group.copy()
		tcr_group['cdr_sequence_collated'] = '-'.join(sequences)
		
		apo_holo_dfs.append(tcr_group)

	return pd.concat(apo_holo_dfs).reset_index(drop=True)
	

def main():
	setup_logger(logger, level=args.log_level)
	args = parser.parse_args()

	logger.info('Loading STCRDab Summary')
	stcrdab_summary = pd.read_csv(os.path.join(args.stcrdab, 'db_summary.dat'), delimiter='\t')
	
	stcrdab_summary['file_path_imgt'] = stcrdab_summary['pdb'].map(
    	lambda pdb_id: os.path.join(args.stcrdab, 'imgt', pdb_id + '.pdb')
    )
	stcrdab_summary['file_path_raw'] = stcrdab_summary['pdb'].map(
    	lambda pdb_id: os.path.join(args.stcrdab, 'raw', pdb_id + '.pdb')
    )
	
	stcrdab_summary['chains'] = stcrdab_summary[['Achain', 'Bchain', 'antigen_chain', 'mhc_chain1']].apply(
    	lambda chains: '-'.join(chains.dropna()),
     	axis=1,
    )

	logger.info('Selecting high quality abTCRs')
	selected_stcrdab = select_tcrs(stcrdab_summary, args.resolution_cutoff)

	logger.info('Removing structures with missing residues')
	selected_stcrdab = screen_tcrs_for_missing_residues(selected_stcrdab)

	logger.info('Annotating the state of the TCR: unbound (*apo*) or bound (*holo*)')
	selected_stcrdab['state'] = selected_stcrdab.apply(
		lambda row: 'apo' if pd.isna(row.antigen_chain) and pd.isna(row.mhc_chain1) else 'holo',
		axis=1,
	)

	logger.info('Removing duplicate structures based on PDB ID')
	holo_pdb_ids = selected_stcrdab.query("state == 'holo'")['pdb'].unique().tolist()
	selected_stcrdab = selected_stcrdab.query("state == 'holo' or (state == 'apo' and pdb not in @holo_pdb_ids)")
	selected_stcrdab = selected_stcrdab.drop_duplicates('pdb')

	logger.info('Retrieving sequences from structures')
	selected_stcrdab = get_stcrdab_sequences(selected_stcrdab)

	logger.info('Finding apo and holo groups')
	apo_holo_tcrs = select_apo_holo(selected_stcrdab)

	logger.info('Exporting dataset')
	summary = {}
 
	if not os.path.exists(args.output):
		os.mkdir(args.output)

	for group_name, group_data in apo_holo_tcrs.groupby('cdr_sequence_collated'):
		output_path = os.path.join(args.output, group_name)
		
		if not os.path.exists(output_path):
			os.mkdir(output_path)
		
		summary[group_name] = []
		
		for _, entry in group_data.iterrows():
			# isolate structure
			with open(entry['file_path_imgt'], 'r') as fh:
				structure = parse_pdb_to_pandas(fh.read())
			
			structure = structure.query("record_type == 'ATOM'")
			structure = structure.query("chain_id in @entry.chains")
			
			with open(os.path.join(output_path, f"{entry.pdb}_{entry.chains}_{entry.state}.pdb"), 'w') as fh:
				fh.write(str(Structure.from_pandas(structure)))
			
			# add info to summary
			structure_summary = {'pdb_id': entry.pdb,
								'state': entry.state,
								'alpha_chain': entry.Achain,
								'beta_chain': entry.Bchain}
			
			if entry.state == 'holo':
				structure_summary['peptide_chain'] = entry['antigen_chain']
				structure_summary['mhc_chain'] = entry['mhc_chain1']
			
			summary[group_name].append(structure_summary)

	with open(os.path.join(args.output, 'summary.json'), 'w') as fh:
		json.dump(summary, fh, indent=1)


if __name__ == '__main__':
    main()