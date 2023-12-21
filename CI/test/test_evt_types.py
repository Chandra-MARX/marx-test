# Licensed under GPL 2 - see LICENSE file
"""Test the output of marx2fits

"""
import subprocess
from glob import glob
import pycrates
from .utils import check_no_warnings

def run_marx_sim(args: dict[str, str | float]) -> None:
    marxcall = [f'{k}={v}' for k, v in args.items()]
    out = subprocess.run(['marx'] + marxcall, check=True, capture_output=True)
    check_no_warnings(out)
    out = subprocess.run(['marx2fits', '--pixadj=EDSER', args['OutputDir'],
                     args['OutputDir']/'evt.fits'],
                    check=True, capture_output=True)
    check_no_warnings(out)
    out = subprocess.run(['marxasp', f"MarxDir={args['OutputDir']}",
                     f'''OutputFile={args['OutputDir']/"asol1.fits"}'''],
                    check=True, capture_output=True)
    check_no_warnings(out)

def test_datatypes_of_evt_files(obsid11005, default_marx_sim):
    """The column types in an event file generated by marx2fits should match CIAO
    
    If they are different, it would be considerably more difficult to merge or compare event lists.
    However, the marx output contains extra columns with simulation output,
    such as the true x/y/z location of a photon, which we ignore.
    """
    evt2 = glob(str(obsid11005 / 'primary'/ '*evt2.fits.gz'))[0]
    evt2 = pycrates.read_file(evt2)
    sim = pycrates.read_file(str(default_marx_sim / 'point.fits'))
    for colname in evt2.get_colnames():
        if colname == 'pha_ro':
            # not written by marx
            continue

        # Need to split name for virtual columns such as chip(chipx, chipy)
        # which are accessed by column names "chip"
        colciao = evt2.get_column(colname.split('(')[0])
        colmarx = sim.get_column(colname.split('(')[0])

        # For ease of reading a failed test, we do not want to assert for
        # every column. If more than one column fails, we want to see all of them.
        # So, we collect all the failures and assert at the end.
        out = {}
        if not (colciao.unit, colciao.get_ndim(), colciao.values.dtype) == \
               (colmarx.unit, colmarx.get_ndim(), colmarx.values.dtype):
            out[colname] = (colciao.unit, colciao.get_ndim(), colciao.values.dtype,
                            colmarx.unit, colmarx.get_ndim(), colmarx.values.dtype)
    assert out == {}
    
    


    