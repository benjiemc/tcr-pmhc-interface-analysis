.PHONY: all environment data analysis notebooks test lint

all: environment data analysis notebooks

environment:
	conda env create -f environment.yml
	pip install .

data:
	@echo "Downloading and processing data..."

data/raw/stcrdab:
	python -m tcr_pmhc_structure_tools.apps.download_stcrdab $@

data/interim/apo-holo-tcr-pmhc-class-I: src/tcr_pmhc_structure_tools/apps/select_structures.py data/raw/stcrdab
	python -m tcr_pmhc_structure_tools.apps.select_structures -o $@ data/raw/stcrdab

data/interim/apo-holo-tcr-pmhc-class-I-imgt-numbered: src/tcr_pmhc_structure_tools/apps/select_structures.py data/interim/apo-holo-tcr-pmhc-class-I
	mkdir $@
	@for path in $(word 2,$^)/*.pdb; do \
		filename=$$(basename $$path); \
		echo Re-numbering $$filename; \
		python -m tcr_pmhc_structure_tools.apps.renumber_structure \
            --output $@/$$filename \
            $$path; \
	done
	cp $(word 2,$^)/apo_holo_summary.csv $@

analysis:
	@echo "Running analysis on data..."

notebooks:
	@echo "Running notebooks..."

test:
	@echo "Testing package"

lint:
	flake8