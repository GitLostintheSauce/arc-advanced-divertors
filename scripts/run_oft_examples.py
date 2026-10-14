#------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Execute OFT TokaMaker example notebooks from the fork and dump their plots

Each notebook's example folder is copied to work/examples/<folder>/ and the
notebook is executed there unmodified. The executed notebook and every PNG
output (<notebook>_cellNN_M.png) are written alongside it.

Usage (from the repo root, with .venv active):
    python scripts/run_oft_examples.py [NOTEBOOK ...] [--oft-root PATH]

NOTEBOOK paths are relative to src/examples/TokaMaker in the fork. With none
given, the roadmap step 1.2 set is run.
'''
import argparse
import base64
import os
import shutil
import time

import nbformat
from nbclient import NotebookClient

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES_SUBDIR = os.path.join('src', 'examples', 'TokaMaker')
DEFAULT_NOTEBOOKS = [
    'fixed_boundary/fixed_boundary_ex1.ipynb',
    'ITER/ITER_baseline_ex.ipynb',
    'DIIID/DIIID_shapehelper_ex.ipynb',
    'ARC/ARCV3A_control_ex.ipynb',
]


def run_notebook(oft_root, rel_path, timeout):
    '''! Copy an example folder to work/ and execute one notebook in it

    @param oft_root Path to the OFT fork checkout
    @param rel_path Notebook path relative to src/examples/TokaMaker
    @param timeout Per-cell timeout [s]
    @result List of PNG paths written
    '''
    src_dir = os.path.join(oft_root, EXAMPLES_SUBDIR, os.path.dirname(rel_path))
    work_dir = os.path.join(REPO_ROOT, 'work', 'examples', os.path.dirname(rel_path))
    shutil.copytree(src_dir, work_dir, dirs_exist_ok=True)
    name = os.path.splitext(os.path.basename(rel_path))[0]

    nb = nbformat.read(os.path.join(work_dir, os.path.basename(rel_path)), as_version=4)
    try:
        # Some examples name a developer-local kernel (e.g. "oft_env"); always use this environment's
        NotebookClient(nb, timeout=timeout, kernel_name='python3',
                       resources={'metadata': {'path': work_dir}}).execute()
    finally:
        nbformat.write(nb, os.path.join(work_dir, name + '_executed.ipynb'))

    png_paths = []
    for icell, cell in enumerate(nb.cells):
        if cell.cell_type != 'code':
            continue
        iimg = 0
        for output in cell.get('outputs', []):
            if 'image/png' not in output.get('data', {}):
                continue
            png_path = os.path.join(work_dir, '{0}_cell{1:02d}_{2}.png'.format(name, icell, iimg))
            with open(png_path, 'wb') as fid:
                fid.write(base64.b64decode(output['data']['image/png']))
            png_paths.append(png_path)
            iimg += 1
    return png_paths


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('notebooks', nargs='*', default=DEFAULT_NOTEBOOKS)
    parser.add_argument('--oft-root', default=os.path.expanduser('~/Documents/code/OpenFUSIONToolkit'),
                        help='Path to the OFT fork checkout')
    parser.add_argument('--timeout', type=int, default=3600, help='Per-cell timeout [s]')
    args = parser.parse_args()

    failed = []
    for rel_path in args.notebooks:
        t0 = time.time()
        try:
            png_paths = run_notebook(args.oft_root, rel_path, args.timeout)
            print('OK    {0} ({1:.0f} s, {2} plots)'.format(rel_path, time.time() - t0, len(png_paths)), flush=True)
        except Exception as exc:
            failed.append(rel_path)
            print('FAIL  {0} ({1:.0f} s): {2}'.format(rel_path, time.time() - t0, str(exc).splitlines()[-1] if str(exc) else exc), flush=True)
    if failed:
        raise SystemExit('{0} notebook(s) failed'.format(len(failed)))


if __name__ == '__main__':
    main()
