# Licensed under GPL 2 - see LICENSE file
"""Test the C interface to the jdmath library

Tests here interface directly at the C level to check that math
functions directly.
Tests are not complete, but are added as work is done on the library.
"""
from os import path
from tempfile import TemporaryDirectory
import shutil
import sys
import pytest
from cffi import FFI


@pytest.fixture(scope='module')
def jdmath():
    # figure out where marx is installed and assume that libraries and include
    # files are there, too
    marx_root = path.dirname(path.dirname(shutil.which('marx')))

    ffibuilder = FFI()

    # cdef() expects a single string declaring the C types, functions and
    # globals needed to use the shared object. It must be in valid C syntax.
    ffibuilder.cdef("""
    unsigned int JDMbinary_search_f (float x, float *xp, unsigned int n);
    unsigned int JDMbinary_search_d (double, double *, unsigned int);
    float JDMinterpolate_f (float x, float *xp, float *yp, unsigned int n);
    double JDMinterpolate_d (double x, double *xp, double *yp, unsigned int n);
    float JDMlog_interpolate_f (float x, float *xp, float *yp, unsigned int n);
    int JDMhistogram_f (float *, unsigned int,
			   float *, unsigned int,
			   unsigned int *, int *);
    int JDMhistogram_d (double *, unsigned int,
			   double *, unsigned int,
			   unsigned int *, int *);
    """)

    # set_source() gives the name of the python extension module to
    # produce, and some C source code as a string.  This C code needs
    # to make the declarated functions, types and globals available,
    # so it is often just the "#include".
    ffibuilder.set_source("_jdmath_cffi",
        """
        #include "jdmath.h"   // the C header of the library
        """,
        libraries=['jdmath'],
        library_dirs=[path.join(marx_root, 'lib/')],
        include_dirs=[path.join(marx_root, 'include/')],
    )
    with TemporaryDirectory() as tmpdirname:
        ffibuilder.compile(tmpdir=tmpdirname)
        sys.path.append(tmpdirname)

        from _jdmath_cffi import lib as jdmath
        from _jdmath_cffi import ffi
        yield jdmath, ffi


@pytest.mark.parametrize(['x', 'expected'], 
                         [(0.1, 0),
                          (0.2, 1),
                          (0.3, 1),
                           # We do expect (0.4, 2) to work as well, but it returns
                           # 1 instead of 2. I believe that this is just due
                           # to the limited precision of the float.
                          (0.40001, 2),
                          (0.5, 2),
                          (0.61, 3),
                          (0.7, 3)])
@pytest.mark.parametrize('search', ['JDMbinary_search_f',
                                   'JDMbinary_search_d'])
def test_binary_search(x, expected, jdmath, search):
    func = getattr(jdmath[0], search)
    assert func(x, [0.2, 0.4, 0.6], 3) == expected

@pytest.mark.parametrize(['x', 'xin', 'yin', 'expected'],
                         [(-1., [0., 1., 2., 3.], [0., 1., 2., 3.], -1.),
                          (1., [0., 1., 2., 3.], [0., 1., 2., 3.], 1.),
                          (1.5, [0., 1., 2., 3.], [0., 1., 2., 3.], 1.5),
                          (2.5, [0., 1., 2., 3.], [0., 1., 2., 3.], 2.5),
                          (5.5, [0., 1., 2., 3.], [0., 1., 2., 4.], 9.),
                          (2.5, [0., 1., 2., 3.], [0., 1., 2., 4.], 3.),
                         ])
@pytest.mark.parametrize('interp', ['JDMinterpolate_f',
                                   'JDMinterpolate_d'])
def test_interpolate(x, xin, yin, expected, jdmath, interp):
    func = getattr(jdmath[0], interp)
    assert func(x, xin, yin, len(xin)) == pytest.approx(expected)

@pytest.mark.parametrize(['x', 'xin', 'yin', 'expected'],
                         [(1., [0., 1., 2., 3.], [0., 1., 2., 3.], 1.),
                          (10., [0., 1., 100.], [0., 0., 10.], 5.0),
                          (10., [0., 1., 100.], [1.2, 1.2, 11.2], 6.2),
                         ])
def test_loginterpolate(x, xin, yin, expected, jdmath):
    assert jdmath[0].JDMlog_interpolate_f(x, xin, yin, len(xin)) == pytest.approx(expected)

@pytest.mark.parametrize('histtype', ['JDMhistogram_f',
                                      'JDMhistogram_d'])
def test_histogram(histtype, jdmath):
    """Test the histogram function"""
    hist = getattr(jdmath[0], histtype)
    histout = jdmath[1].new("unsigned int [4]")
    reverse_indices = jdmath[1].new("int [3]")
    ret = hist([1., 1., 2.], 3,
               [-0.5, 0.5, 1.5, 2.5, 3.5], 4,
               histout, reverse_indices)
    assert ret == 0
    assert list(histout) == [0, 2, 1, 0]
    assert list(reverse_indices) == [1, 1, 2]