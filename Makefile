.PHONY: all environment data analysis notebooks test lint

all: environment data analysis notebooks

environment:
	conda env create -f environment.yml
	pip install .

data:
	@echo "Downloading and processing data..."
	python -m tcr_pmhc_structure_tools.apps.download_stcrdab data/raw/stcrdab
	python -m tcr_pmhc_structure_tools.apps.select_structures -o data/interim/apo-holo-tcr-pmhc-class-I data/raw/stcrdab

analysis:
	@echo "Running analysis on data..."

notebooks:
	@echo "Running notebooks..."

test:
	@echo "Testing package"

lint:
	flake8