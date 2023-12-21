# Licensed under GPL 2 - see LICENSE file
import os
import subprocess
import re

def check_no_warnings(out: subprocess.CompletedProcess) -> None:
    """Check that no warnings are included in the standard output
    """
    assert 'WARNING' not in out.stdout.decode('utf-8')
    assert 'Warning' not in out.stdout.decode('utf-8')
    assert len(out.stderr) == 0

def find_par_file(stem: str='marx') -> str:
    """Find the location of the marx.par file

    marx itself can parse the PFILE environment variable to find the
    location of the ``marx.par`` file, but when I call directly into the
    C interface, then it is easier to have the file location defined
    in Python to limit the C calls to the exact routines I want to test.

    Parameters
    ----------
    stem : str
        The name of the par file to find (default: marx)

    Returns
    -------
    parfile : str
        The full path to the par file
    """
    paths = re.split(r'[:;]', os.environ['PFILES'])
    for p in paths:
        if os.path.exists(os.path.join(p, stem + '.par')):
            return os.path.join(p, stem + '.par')
    raise FileNotFoundError(f"Could not find {stem}.par in {paths}")
