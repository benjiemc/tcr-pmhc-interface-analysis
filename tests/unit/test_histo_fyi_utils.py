import pytest

from tcr_pmhc_structure_tools.histo_fyi_utils import fetch_structure


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
