import os
from warnings import warn
import gzip
import shutil

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from astropy.table import Table
import numpy as np


def detectorfromkeyword(keyword: str) -> str:
    '''Determine detector type from fits header keyword

    For ACIS, this uses a heuristic to decide which chips are on,
    since Chandra allows arbitrary combinations,
    but MARX only know ACIS-I or ACIS-S.

    Parameters
    ----------
    keyword : string
        Chandra fits file header keyword like 'ACIS-012356'.

    Returns
    -------
    detector : string
        one of 'ACIS-I', 'ACIS-S', 'HRC-I', 'HRC-S', or ''
    '''
    # DetectorType
    if keyword.startswith('HRC'):
        return keyword
    
    if keyword.startswith('ACIS'):
        # Use a heuristic to decide which chips are on, since Chandra
        # allows arbitrary combinations, but MARX only knows ACIS-I or ACIS-S
        acisi = len(set(keyword).intersection(set('0123')))
        aciss = len(set(keyword).intersection(set('456789')))
        if acisi > aciss:
            return 'ACIS-I'
        else:
            return 'ACIS-S'
    
    return ''


def target_coos_from_asol(asolfile: str) -> tuple[float, float]:
    '''Get coordinates for the target named in asol file.

    If the target name cannot be resolved, return the nominal aimpoint,
    assuming a on-axis pointing.

    Parameters
    ----------
    asolfile : str
        Path and name of an asol file

    Returns
    -------
    ra, dec : float
        RA and Dec in degrees
    '''
    asol = fits.getheader(asolfile, 1)
    try:
        objname = asol['OBJECT'].split('/')[0]
        skyco = SkyCoord.from_name(objname)
        return skyco.ra.value, skyco.dec.value
    except NameResolveError:
        warn('Target name {0} not resolved by SIMBAD. Assuming on-axis target.'.format(objname))
        return asol['RA_NOM'], asol['DEC_NOM']


def marxpars_from_asol(asolfile: str, evt2file: str) -> dict:
    '''Set MARX parameters from asol and evt file.

    This function parses the header of a fits file and uses the
    information in the header to set as many marx parameters as possible.

    Parameters
    ----------
    asolfile : str
        Path and name of an asol file

    evt2file : str
        Path and name of an evt2 file

    Returns
    -------
    marx_pars : dict
        Dictionary with marx parameters as far as they can be extracted from
        the data in the asol file.
    '''
    asol = fits.getheader(asolfile, 1)
    evt = fits.getheader(evt2file, 1)

    '''Unzip fits file.

    MARX cannot read zipped fits files, so we need to unzip the .fits.gz asol
    files that we downloaded from the archive. On the other hand, CIAO
    tools work on both zipped or unzipped files, so there is no need to
    unzip all of them, just the files that MARX reads as input.
    '''
    if not os.path.exists(asolfile.removesuffix('.gz')):
        with gzip.open(asolfile, 'rb') as f_in:
            with open(asolfile.removesuffix('.gz'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    marx_pars = {}
    marx_pars['RA_Nom'] = asol['RA_NOM']
    marx_pars['Dec_Nom'] = asol['DEC_NOM']
    marx_pars['Roll_Nom'] = asol['ROLL_NOM']
    marx_pars['GratingType'] = asol['GRATING']
    marx_pars['ExposureTime'] = asol['TSTOP'] - asol['TSTART']
    marx_pars['DitherModel'] = 'FILE'
    marx_pars['DitherFile'] = asolfile.removesuffix('.gz')
    marx_pars['TStart'] = asol['TSTART']

    # ACIS Exposure time (might vary for sub-array read-out)
    if evt['INSTRUME'][0:4] == 'ACIS':
        marx_pars['ACIS_Exposure_Time'] = evt['EXPTIME']
    # Target coordiantes
    # ra, dec = target_coos_from_asol(asolfile)
    marx_pars['SourceRA'] = evt['RA_TARG']
    marx_pars['SourceDEC'] = evt['DEC_TARG']

    # DetectorType
    det = detectorfromkeyword(evt['DETNAM'])
    if det != '':
        marx_pars['DetectorType'] = det
        # Find offsets from nominal pointing
        aimpts = Table.read(os.path.join(os.environ['CALDB'], 'data', 'chandra', 'default', 'aimpts', 'telD1999-07-23aimptsN0002.fits'))
        det2aimpts = {'ACIS-I': 'AI2', 'ACIS-S': 'AS1', 'HRC-I': 'HI1', 'HRC-S': 'HS1'}
        aimp = np.array([n.strip() for n in aimpts['AIMPOINT_NAME']])
        ind = aimp == det2aimpts[det]
        marx_pars['DetOffsetX'] = evt['SIM_X'] - aimpts[ind]['AIMPOINT'][0][0]
        # No Y offset possible
        marx_pars['DetOffsetZ'] = evt['SIM_Z'] - aimpts[ind]['AIMPOINT'][0][2]
    else:
        warn('detector {0} not understood. No DetectorType or aimpoint offset set.'.format(asol['DETNAM']))
    return marx_pars
