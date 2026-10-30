#------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Check regenerated results against the expected values (used by CI)

  arc    Read work/arcv3a_baseline/executed.ipynb (from run_arcv3a_baseline.py) and check
         the equilibrium is diverted with Ip = 12 MA and |q95| close to 3.38.
  sparc  Compare every sparc/static_scan/results/*.json against the committed copy in
         git HEAD, number by number, within a relative tolerance.

Usage (from the repo root, with .venv active):
    python scripts/check_ci_results.py arc
    python scripts/check_ci_results.py sparc [--rtol 1e-4]
'''
import argparse
import glob
import json
import os
import re
import subprocess

import nbformat

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARC_NOTEBOOK = os.path.join(REPO_ROOT, 'work', 'arcv3a_baseline', 'executed.ipynb')
SPARC_RESULTS = os.path.join('sparc', 'static_scan', 'results')


def check_arc():
    '''! Check the ARC V3A baseline equilibrium statistics

    @result List of failure messages (empty if all checks pass)
    '''
    nb = nbformat.read(ARC_NOTEBOOK, as_version=4)
    text = ''
    for cell in nb.cells:
        if cell.cell_type == 'code' and 'print_info' in cell.source:
            text = ''.join(output.get('text', '') for output in cell.get('outputs', []))
            break
    print(text.split('  CS1U')[0])

    failures = []
    if not re.search(r'Topology\s*=\s*Diverted', text):
        failures.append('ARC baseline is not diverted')
    match = re.search(r'Toroidal Current \[A\]\s*=\s*(\S+)', text)
    if match is None or abs(float(match.group(1)) - 12.0E6) > 1.0E3:
        failures.append('ARC baseline Ip is not 12 MA')
    match = re.search(r'q_0, q_95\s*=\s*(\S+)\s+(\S+)', text)
    if match is None or abs(abs(float(match.group(2))) - 3.38) > 0.02:
        failures.append('ARC baseline |q95| is not 3.38 +/- 0.02')
    return failures


def compare(new, old, rtol, path, diffs):
    '''! Recursively compare two JSON values, recording numbers outside tolerance

    @param new Regenerated value
    @param old Committed value
    @param rtol Relative tolerance for numbers
    @param path Key path used in messages
    @param diffs List that mismatch messages are appended to
    @result Largest relative difference found
    '''
    if isinstance(old, bool) or isinstance(old, str) or old is None:
        if new != old:
            diffs.append('{0}: {1!r} != {2!r}'.format(path, new, old))
        return 0.0
    if isinstance(old, (int, float)):
        if not isinstance(new, (int, float)):
            diffs.append('{0}: {1!r} != {2!r}'.format(path, new, old))
            return 0.0
        rel = abs(new - old) / max(abs(old), 1.0E-12)
        if rel > rtol and abs(new - old) > 1.0E-9:
            diffs.append('{0}: {1!r} != {2!r} (rel {3:.2e})'.format(path, new, old, rel))
        return rel if abs(new - old) > 1.0E-9 else 0.0
    if isinstance(old, dict):
        if set(new) != set(old):
            diffs.append('{0}: keys differ'.format(path))
            return 0.0
        return max([compare(new[key], old[key], rtol, path + '.' + key, diffs) for key in old] + [0.0])
    if len(new) != len(old):
        diffs.append('{0}: length {1} != {2}'.format(path, len(new), len(old)))
        return 0.0
    return max([compare(a, b, rtol, '{0}[{1}]'.format(path, i), diffs) for i, (a, b) in enumerate(zip(new, old))] + [0.0])


def check_sparc(rtol):
    '''! Compare regenerated SPARC scan results with the committed ones

    @param rtol Relative tolerance for numbers
    @result List of failure messages (empty if all checks pass)
    '''
    failures = []
    paths = sorted(glob.glob(os.path.join(REPO_ROOT, SPARC_RESULTS, '*.json')))
    for path in paths:
        rel_path = os.path.relpath(path, REPO_ROOT)
        committed = subprocess.run(['git', 'show', 'HEAD:' + rel_path], cwd=REPO_ROOT,
                                   capture_output=True, text=True, check=True).stdout
        with open(path) as fid:
            new = json.load(fid)
        diffs = []
        max_rel = compare(new, json.loads(committed), rtol, os.path.basename(path), diffs)
        print('{0:6s} {1} (max rel diff {2:.2e})'.format('FAIL' if diffs else 'OK', rel_path, max_rel))
        failures += diffs
    if not paths:
        failures.append('no SPARC result files found')
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('check', choices=['arc', 'sparc'])
    parser.add_argument('--rtol', type=float, default=1.0E-4, help='Relative tolerance for the sparc check')
    args = parser.parse_args()

    failures = check_arc() if args.check == 'arc' else check_sparc(args.rtol)
    for message in failures:
        print('FAIL  ' + message)
    if failures:
        raise SystemExit('{0} check(s) failed'.format(len(failures)))
    print('All {0} checks passed'.format(args.check))


if __name__ == '__main__':
    main()
