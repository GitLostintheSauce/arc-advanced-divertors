#------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Run the OFT ARC V3A baseline example and save its equilibrium plot

Executes an unmodified copy of ARCV3A_baseline_ex.ipynb from the OFT fork in
work/arcv3a_baseline/ (the notebook downloads the public gEQDSK from Zenodo
record 19498373 there) and writes the flux-surface plot to figures/.

Usage (from the repo root, with .venv active):
    python scripts/run_arcv3a_baseline.py [--oft-root ~/Documents/code/OpenFUSIONToolkit]
'''
import argparse
import base64
import os
import shutil

import nbformat
from nbclient import NotebookClient

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE_DIR = os.path.join('src', 'examples', 'TokaMaker', 'ARC')
INPUT_FILES = ['ARCV3A_baseline_ex.ipynb', 'arc_v3a_mesh.h5']


def save_cell_png(nb, marker, out_path):
    '''! Save the first PNG output of the first code cell containing `marker`

    @param nb Executed notebook
    @param marker Source substring identifying the plotting cell
    @param out_path Destination PNG path
    @result True if an image was found and written
    '''
    for cell in nb.cells:
        if cell.cell_type != 'code' or marker not in cell.source:
            continue
        for output in cell.get('outputs', []):
            if 'image/png' in output.get('data', {}):
                with open(out_path, 'wb') as fid:
                    fid.write(base64.b64decode(output['data']['image/png']))
                return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--oft-root', default=os.path.expanduser('~/Documents/code/OpenFUSIONToolkit'),
                        help='Path to the OFT fork checkout')
    args = parser.parse_args()

    work_dir = os.path.join(REPO_ROOT, 'work', 'arcv3a_baseline')
    os.makedirs(work_dir, exist_ok=True)
    for filename in INPUT_FILES:
        shutil.copy(os.path.join(args.oft_root, EXAMPLE_DIR, filename), work_dir)

    nb = nbformat.read(os.path.join(work_dir, INPUT_FILES[0]), as_version=4)
    NotebookClient(nb, timeout=1200, resources={'metadata': {'path': work_dir}}).execute()
    nbformat.write(nb, os.path.join(work_dir, 'executed.ipynb'))

    out_path = os.path.join(REPO_ROOT, 'figures', 'arcv3a_baseline_equilibrium.png')
    if not save_cell_png(nb, 'plot_psi', out_path):
        raise RuntimeError('Equilibrium plot not found in executed notebook')
    print('Saved {0}'.format(out_path))


if __name__ == '__main__':
    main()
