.PHONY: all environment data analysis notebooks test lint

all: data analysis notebooks

environment:
	conda env create -f environment.yml
	pip install .

data: \
	data/processed/apo-holo-tcr-pmhc-class-I \
	data/processed/apo-holo-tcr-pmhc-class-I-holo-aligned \
	data/processed/structure-pw-distances

data/raw/stcrdab:
	python -m tcr_pmhc_structure_tools.apps.download_stcrdab $@

data/interim/apo-holo-tcr-pmhc-class-I: data/raw/stcrdab
	python -m tcr_pmhc_structure_tools.apps.select_structures -o $@ $^

data/interim/apo-holo-tcr-pmhc-class-I-imgt-numbered: data/interim/apo-holo-tcr-pmhc-class-I
	mkdir -p $@
	@for path in $(word 2,$^)/*.pdb; do \
		filename=$$(basename $$path); \
		echo Re-numbering $$filename; \
		python -m tcr_pmhc_structure_tools.apps.renumber_structure \
            --output $@/$$filename \
            $$path; \
	done
	head -n1 $^/apo_holo_summary.csv > $@/apo_holo_summary.csv
	ls $@/*.pdb | cut -d/ -f4 | xargs -I % grep % $^/apo_holo_summary.csv >> $@/apo_holo_summary.csv

data/processed/apo-holo-tcr-pmhc-class-I: data/interim/apo-holo-tcr-pmhc-class-I-imgt-numbered
	mkdir -p $@
	python -m tcr_pmhc_structure_tools.apps.align_tcr_pmhcs -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-holo-aligned: data/interim/apo-holo-tcr-pmhc-class-I-imgt-numbered/
	python -m tcr_pmhc_structure_tools.apps.align_tcr_pmhcs --only-holo -o $@ $^

data/processed/structure-pw-distances: data/raw/stcrdab
	python -m tcr_pmhc_structure_tools.apps.compute_pw_distances --compress-output -o $@ $^

analysis: $(wildcard data/processed/apo-holo-tcr-pmhc-class-I-comparisons/*.csv)

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_loop_align_results.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities tcr --align-entities -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_fw_align_results.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities tcr -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/tcr_per_res_apo_holo_loop_align.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities tcr --align-entities --per-residue -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/pmhc_per_res_apo_holo.csv: data/processed/apo-holo-tcr-pmhc-class-I
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities pmhc --per-residue -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_fw_align_holo.csv: data/processed/apo-holo-tcr-pmhc-class-I-holo-aligned
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --select-entities tcr -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/rmsd_cdr_loop_align_holo.csv: data/processed/apo-holo-tcr-pmhc-class-I-holo-aligned
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences --align-entities --select-entities tcr -o $@ $^

data/processed/apo-holo-tcr-pmhc-class-I-comparisons/pmhc_tcr_contact_apo_holo.csv: data/processed/apo-holo-tcr-pmhc-class-I data/processed/mhc_contacts.csv
	python -m tcr_pmhc_structure_tools.apps.compute_apo_holo_differences \
	-o $@ \
	--select-entities pmhc \
	--pmhc-tcr-contact-residues $(shell awk -F ',' '$$3 >= 100 { print $$2 }' $(word 2,$^) | tail -n +2 | sort | uniq | tr '\n' ' ') \
	$(word 1,$^)

data/processed/mhc_contacts.csv: run_notebook_Identify_contact_residues_on_MHC_Class_I_molecules

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