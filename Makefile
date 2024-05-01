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

analysis:
	@echo "Running analysis on data..."

notebooks: $(patsubst notebooks/%.ipynb,run_notebook_%,$(wildcard notebooks/*.ipynb))

run_notebook_%:
	@echo "Executing $*"
	@jupyter nbconvert --to notebook --execute --inplace notebooks/$*.ipynb

test:
	@pytest

lint:
	@echo "Running linting"
	@flake8