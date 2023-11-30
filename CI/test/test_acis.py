# Licensed under GPL 2 - see LICENSE file
"""Test ACIS-specific functionality in MARX

Some of these test are done at the C level as unit tests
"""
from os import path
import os
from tempfile import TemporaryDirectory
import shutil
import sys
import pytest
from cffi import FFI
import numpy as np

@pytest.fixture(scope='module')
def acis_ccode():
    # figure out where marx is installed and assume that libraries and include
    # files are there, too
    marx_root = path.dirname(path.dirname(shutil.which('marx')))

    ffibuilder = FFI()

    # cdef() expects a single string declaring the C types, functions and
    # globals needed to use the shared object. It must be in valid C syntax.
    ffibuilder.cdef("""
    typedef ... Param_File_Type;
    Param_File_Type *marx_pf_parse_cmd_line (char *file, char *mode, int argc, char **argv);
    int marx_init_acis_s_rmf (Param_File_Type *p);
    int marx_map_energy_to_acis_pha (int ccd_id, int x, int y, double energy, short *phap);
    int marx_apply_acis_rmf (int ccd_id, float x, float y,
			 double energy, float *pip, short *phap);
    """)

    # set_source() gives the name of the python extension module to
    # produce, and some C source code as a string.  This C code needs
    # to make the declarated functions, types and globals available,
    # so it is often just the "#include".
    ffibuilder.set_source("_acis_cffi",
        """
        #include "marx.h"   // the C header of the library
        // Usually, this is set up when the parameter file is read
        // but since we just pick and choose a few functions, we need to
        // set it up manually
        int _Marx_Det_Extend_Flag = 0;
        """,
        libraries=['jdmath', 'jdfits', 'pfile', 'marx'],
        library_dirs=[path.join(marx_root, 'lib/')],
        include_dirs=[path.join(marx_root, 'include/')],
    )
    with TemporaryDirectory() as tmpdirname:
        ffibuilder.compile(tmpdir=tmpdirname)
        sys.path.append(tmpdirname)

        from _acis_cffi import lib, ffi
        yield lib, ffi


@pytest.mark.parametrize('ccd_id, x, y, energy, expected', 
                         [(4, 100, 200, 1.5, 372),
                          (7, 400, 234, 0.3, 43),
                          (7, 600, 534, 0.3, 43),
                          (7, 400, 234, 0.5, 85),
                          (6, 623, 1003, 6.8, 1455),
                         ])
def test_marx_map_energy_to_acis_pha_aciss(ccd_id, x, y, energy, expected, acis_ccode):
    """Test the C function marx_map_energy_to_acis_pha

    This test could be generalized to also text acis_i.
    The tested numbers are not from first principles, but set as a
    regressions test to match marx 5.5.0.
    """
    lib, ffi = acis_ccode
    # need to read a parameter file, because we need that as input
    # to marx_init_acis_s_rmf below
    marxpar = ffi.new("char[]",
                      (os.environ['PFILES'] + "/marx.par").encode('utf-8'))
    mode = ffi.new("char[]", b"rwL")
    # program name "marx" is the first command line argument
    p = ffi.new("char[]", b"marx")    # p is a 'char *'
    q = ffi.new("char **", p)
    pf = lib.marx_pf_parse_cmd_line(marxpar, mode, 1, q)
    assert 0 == lib.marx_init_acis_s_rmf(pf)

    phap = ffi.new("short *")
    lib.marx_map_energy_to_acis_pha(ccd_id, x, y, energy, phap)
    assert phap[0] == expected


def test_acis_fef_gaussian_rmf(acis_ccode):
    lib, ffi = acis_ccode
    # need to read a parameter file, because we need that as input
    # to marx_init_acis_s_rmf below
    marxpar = ffi.new("char[]",
                      (os.environ['PFILES'] + "/marx.par").encode('utf-8'))
    mode = ffi.new("char[]", b"rwL")
    # program name "marx" is the first command line argument
    p = ffi.new("char[]", b"marx")    # p is a 'char *'
    q = ffi.new("char **", p)
    pf = lib.marx_pf_parse_cmd_line(marxpar, mode, 1, q)
    assert 0 == lib.marx_init_acis_s_rmf(pf)

    pip = ffi.new("float *")
    phap = ffi.new("short *")

    pha = np.empty(1000)
    pi = np.empty(1000)
    for i in range(1000):
        lib.marx_apply_acis_rmf(7, 400, 234, 0.3, pip, phap)
        pha[i] = phap[0]
        pi[i] = pip[0]

    assert pha.mean() == pytest.approx(39.5, abs=.2)
    assert pha.std() == pytest.approx(9.2, rel=.1)
    assert pi.mean() == pytest.approx(0.23, abs=.05)
    assert pi.std() == pytest.approx(0.06, rel=0.2)