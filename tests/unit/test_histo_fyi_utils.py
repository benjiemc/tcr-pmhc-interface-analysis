import pytest

from tcr_pmhc_interface_analysis.histo_fyi_utils import fetch_structure, retrieve_data_from_api, PMHC_CLASS_I_URL


class TestFetchStructure:
    def test(self):
        structure_text = fetch_structure('3vfn', 1)
        assert len(structure_text) == 137309

    def test_peptide(self):
        structure_text = fetch_structure('3vfn', 1, domain='peptide')
        assert len(structure_text) == 8267

    def test_abd(self):
        structure_text = fetch_structure('3vfn', 1, domain='abd')
        assert len(structure_text) == 129041

    def test_value_error(self):
        with pytest.raises(ValueError):
            fetch_structure('3vfn', 1, domain='pinapple')


class TestRetrieveDataFromAPI:
    def test_pmhc(self):
        histo_pmhc_df = retrieve_data_from_api(PMHC_CLASS_I_URL)
        assert histo_pmhc_df.columns.tolist() == ['pdb_id',
                                                  'peptide_sequence',
                                                  'mhc_slug',
                                                  'antigen_chain',
                                                  'mhc_chain1',
                                                  'mhc_chain2',
                                                  'chains',
                                                  'assembly_number',
                                                  'resolution']
