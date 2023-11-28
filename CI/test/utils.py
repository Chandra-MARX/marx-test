# Licensed under GPL 2 - see LICENSE file
import subprocess

def check_no_warnings(out: subprocess.CompletedProcess) -> None:
    """Check that no warnings are included in the standard output
    """
    assert 'WARNING' not in out.stdout.decode('utf-8')
    assert 'Warning' not in out.stdout.decode('utf-8')
    assert len(out.stderr) == 0
