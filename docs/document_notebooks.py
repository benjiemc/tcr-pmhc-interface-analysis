import glob
import os
import shutil
import sys

notebook_names = []

for notebook in glob.glob(os.path.join(sys.argv[1], '*.ipynb'), recursive=True):
    notebook_names.append(os.path.basename(notebook).rstrip('.ipynb'))
    shutil.copy(notebook, sys.argv[2])

notebooks_contents = ('Analysis Notebooks\n'
                      '==================\n'
                      '\n'
                      '.. toctree::\n'
                      '   :maxdepth: 1\n'
                      '\n') + '\n'.join(['   ' + name for name in sorted(notebook_names)])

with open(os.path.join(sys.argv[2], 'notebooks.rst'), 'w') as fh:
    fh.write(notebooks_contents)
