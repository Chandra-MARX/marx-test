from warnings import warn

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError


def detectorfromkeyword(keyword):
    '''Determine detector type form fits header keyword

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
    if keyword[:3] == 'HRC':
        return keyword
    elif keyword[:4] == 'ACIS':
        # Use a heuristic to decide which chips are on, since Chandra
        # allows arbitrary combinations, but MARX only know ACIS-I or ACIS-S
        acisi = len(set(keyword).intersection(set('0123')))
        aciss = len(set(keyword).intersection(set('456789')))
        if acisi >= aciss:
            return 'ACIS-I'
        else:
            return 'ACIS-S'
    else:
        return ''


def marxpars_from_asol(asolfile):
    '''Set MARX parameters from asol or evt file.

    This function parses the header of a fits file and uses the
    information in the header to set as many marx parameters as possible.

    Parameters
    ----------
    asolfile : string
        Path and name of an asol or evt file

    Returns
    -------
    marx_pars : dict
        Dictionary with marx parameters as far as they can be extracted from
        the data in the asol file.
    '''
    asol = fits.getheader(asolfile, 1)

    marx_pars = {'RA_Nom': asol['RA_NOM'],
                 'Dec_Nom': asol['DEC_NOM'],
                 'Roll_Nom': asol['ROLL_NOM'],
                 'GratingType': asol['GRATING'],
                 'ExposureTime': asol['TSTOP'] - asol['TSTART'],
                 'DitherModel': 'FILE',
                 'DitherFile': asolfile,
                 'TStart': asol['TSTART'],
    }
    # Target coordiantes
    try:
        objname = asol['OBJECT'].split('/')[0]
        skyco = SkyCoord.from_name(objname)
        marx_pars['SourceRA'] = skyco.ra.value
        marx_pars['SourceDEC'] = skyco.dec.value
    except NameResolveError:
        warn('Target name {0} not resolved. Assuming on-axis target.'.format(objname))
        marx_pars['SourceRA'] = asol['RA_NOM']
        marx_pars['SourceDEC'] = asol['DEC_NOM']

    # DetectorType
    det = detectorfromkeyword(asol['DETNAM'])
    if det != '':
        marx_pars['DetectorType'] = det
    else:
        warn('detector {0} not understood. No DetectorType set.'.format(asol['DETNAM']))

    return marx_pars


def spectrum_from_flux_correction(asolfile, evtfile, x, y, region, psffrac=1):
    '''Estimate a source spectrum through flux correction.

    We will simply extract the counts in a region containing the observed PSF
    and divide it by the effective area. This is known as "flux-correction",
    and involves no spectral fitting. Strictly speaking, the validity of this
    technique assumes that the spectrum does not vary much over the scale of
    the RMF.

    Parameters
    ----------
    asolfile : string
        Path and filename of asolfile
    evtfile : string
        Path and filename of evt2 file
    x, y : float
        x, y position of source in detector coordiantes
    region : string
        CIAO region for source extraction
    psffrac : float
        Correction factor for a PSF not fully included in ``region``

    Returns
    -------
    ciaocommands : list of strings
        CIAO commands ready to paste into a shell.
        When executed, the codewill write two files:

        - ``input_spec_marx.tbl``
        - ``input_spec_saotrace.rdb``
    '''
    evt = fits.getheader(evtfile, 1)

    det = detectorfromkeyword(evt['DETNAM'])
    if det[:4] != 'ACIS':
        raise ValueError('This works only for ACIS. Detector found: {0}'.format(det))

    return '''
asphist infile={asolfile} outfile=obs.asp evtfile={evtfile} clobber=yes
# It would be more accuarte to use mkwarf.
# However, that's slower and for our purposes an approximate solution is good enough.
mkarf detsubsys={det} grating={grating} outfile=obs.arf obsfile={evtfile} asphistfile=obs.asp sourcepixelx={x} sourcepixely={y} engrid="0.3:8.0:0.1" maskfile=NONE pbkfile=NONE dafile=NONE verbose=1 mode=h clobber=yes

# We have to use the same energy binning in energy space here that we used for the arf!
# So, first convert the PI to energy (to a precision that's good enough for this example.
dmtcalc "{evtfile}[EVENTS][sky={region}]" evt2_with_energy.fits  expr="energy=(float)pi*0.0149" clobber=yes
dmextract "evt2_with_energy.fits[bin energy=.3:7.999:0.1]" obs.spec clobber=yes opt=generic
dmcopy "obs.arf[cols energ_lo,energ_hi,specresp]" input_spec.fits clobber=yes
dmpaste input_spec.fits "obs.spec[cols counts]" input_spec_1.fits clobber=yes
dmtcalc input_spec_1.fits input_spec_2.fits  expr="flux=counts/((float)specresp * {psffrac}*{exptime})" clobber=yes

# devide by binwidth to turn the flux into a flux DENSITY for marx
dmtcalc input_spec_2.fits input_spec_3.fits  expr="fluxdens=flux/0.1" clobber=yes
dmcopy "input_spec_3.fits[cols energ_hi,fluxdens]" "input_spec_marx.tbl[opt kernel=text/simple]" clobber=yes
dmcopy "input_spec_3.fits[cols energ_lo,energ_hi,flux]" "input_spec_saotrace.tbl[opt kernel=text/simple]" clobber=yes

# We now use a combination of some of the most obscure UNIX tools to
# bring the SAOTRACE input spectrum into the right format
# see http://cxc.harvard.edu/cal/Hrma/RDB/FileFormat.html
# Remove the leading "#" from the line with the column names

awk '{{ sub(/\\# ENERG_LO/,"ENERG_LO"); print }}' < input_spec_saotrace.tbl > input_spec_saotrace.temp
# Then, replace spaces with tabs
tr ' ' \\\\t < input_spec_saotrace.temp > input_spec_saotrace.rdb

'''.format(det=det, grating=evt['GRATING'], evtfile=evtfile, x=x, y=y, asolfile=asolfile, region=region, psffrac=psffrac, exptime=evt['TSTOP'] - evt['TSTART']).split('\n')
