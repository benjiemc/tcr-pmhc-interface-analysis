.PHONY: all environment data analysis notebooks test

all: environment data analysis notebooks

environment:
	conda env create -f environment.yml
	pip install .

data:
	@echo "Downloading and processing data..."
	python -m tcr_pmhc_structure_tools.apps.download_stcrdab data/raw/stcrdab

analysis:
	@echo "Running analysis on data..."

notebooks:
	@echo "Running notebooks..."

make test:
	@echo "Testing package"