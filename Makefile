.PHONY: all environment data analysis notebooks test lint

all: data analysis notebooks

environment:
	conda env create -f environment.yml
	pip install .

data: data/processed/apo-holo-tcr-pmhc-class-I

data/raw/stcrdab:
	python -m tcr_pmhc_structure_tools.apps.download_stcrdab $@

data/interim/apo-holo-tcr-pmhc-class-I: src/tcr_pmhc_structure_tools/apps/select_structures.py data/raw/stcrdab
	python -m tcr_pmhc_structure_tools.apps.select_structures -o $@ data/raw/stcrdab

data/interim/apo-holo-tcr-pmhc-class-I-imgt-numbered: src/tcr_pmhc_structure_tools/apps/renumber_structure.py data/interim/apo-holo-tcr-pmhc-class-I
	mkdir -p $@
	@for path in $(word 2,$^)/*.pdb; do \
		filename=$$(basename $$path); \
		echo Re-numbering $$filename; \
		python -m tcr_pmhc_structure_tools.apps.renumber_structure \
            --output $@/$$filename \
            $$path; \
	done
	cp $(word 2,$^)/apo_holo_summary.csv $@

data/processed/apo-holo-tcr-pmhc-class-I: src/tcr_pmhc_structure_tools/apps/align_tcr_pmhcs.py data/interim/apo-holo-tcr-pmhc-class-I-imgt-numbered
	mkdir -p $@
	python -m tcr_pmhc_structure_tools.apps.align_tcr_pmhcs -o $@ $(word 2,$^)

analysis: \
	data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_loop_align_results.csv \
	data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_fw_align_results.csv \
	data/processed/apo-holo-tcr-pmhc-class-I-comparisons/tcr_apo_per_res_holo_loop_align.csv \
	data/processed/apo-holo-tcr-pmhc-class-I-comparisons/pmhc_per_res_apo_holo.csv

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_loop_align_results.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities tcr --align-entities -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_fw_align_results.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities tcr -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/tcr_apo_per_res_holo_loop_align.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities tcr --align-entities --per-residue -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/pmhc_per_res_apo_holo.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities pmhc --per-residue -o $@ $^

notebooks: $(patsubst notebooks/%.ipynb,run_notebook_%,$(wildcard notebooks/*.ipynb))

run_notebook_%:
	@echo "Executing $*"
	@jupyter nbconvert --to notebook --execute --inplace notebooks/$*.ipynb

test:
	@pytest

lint:
	@echo "Running linting"
	@flake8
	@isort --check-only src