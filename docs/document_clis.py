import glob
import os
import subprocess
import sys

cli_apps = glob.glob('src/**/apps/*.py', recursive=True)
cli_apps = [app for app in cli_apps if not app.split('/')[-1].startswith('_')]

app_names = []

for app in cli_apps:
    app_name = app.strip('src/').replace('/', '.').rstrip('.py')

    app_doc = subprocess.run(['python', app, '--help'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL,
                             universal_newlines=True).stdout

    with open(os.path.join(sys.argv[1], app_name + '.rst'), 'w') as fh:
        fh.write(app_name)
        fh.write('\n')
        fh.write('=' * len(app_name))
        fh.write('\n')
        fh.write('.. code-block:: text\n')
        fh.write('\n')
        fh.write('\n'.join(['   ' + line for line in app_doc.split('\n')]))
        fh.write('\n')

    app_names.append(app_name)

apps_contents = ('Command Line Applications\n'
                 '=========================\n'
                 '\n'
                 '.. toctree::\n'
                 '   :maxdepth: 1\n'
                 '\n') + '\n'.join(['   ' + name for name in sorted(app_names)])

with open(os.path.join(sys.argv[1], 'apps.rst'), 'w') as fh:
    fh.write(apps_contents)
