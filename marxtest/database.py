'''Connect to a database to store test results


'''
import sqlite3
from functools import wraps
import subprocess
import os
import re
from ConfigParser import NoOptionError
from warnings import warn

con = None
configdict = {}
configid = None


class InconsistentDBError(Exception):
    pass


def setup(f):
    '''Connect to database.

    Also, this looks up the id of the current configuration and creates a new
    row in the "config" table if necessary.'''
    @wraps(f)
    def wrapper(conf, *args, **kwds):

        global con, configdict, configid

        if len(configdict) == 0:
            set_configdict(conf)

        if con is None:
            con = sqlite3.connect(conf.get('Output', 'sqlitedb'))
            con.row_factory = sqlite3.Row

            con.execute("CREATE TABLE IF NOT EXISTS config (marx STRING, ciao STRING, caldb string,  saotrace STRING, compiler STRING, compilerversion STRING, compilerflags STRING, host STRING);")
            con.execute("CREATE TABLE IF NOT EXISTS expresult (testclass STRING NOT NULL, name STRING, version INTEGER, title STRING CHECK(LENGTH(title) <= 20), description STRING, value, unit STRING, acceptable, PRIMARY KEY (testclass, name, version));")
            con.execute("CREATE TABLE IF NOT EXISTS result (configid INTEGER, expid INTEGER, value, sigma, sigma2);")

        cur = con.execute("SELECT rowid FROM config WHERE marx=:marx AND ciao=:ciao AND caldb=:caldb AND saotrace=:saotrace AND compiler=:compiler AND compilerversion=:compilerversion AND compilerflags=:compilerflags and host=:host", configdict)

        rows = cur.fetchall()
        if len(rows) == 0:
            cur = con.execute('INSERT INTO config VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                              [configdict[n] for n in ['marx', 'ciao', 'caldb',
                                                       'saotrace', 'compiler',
                                                       'compilerversion',
                                                       'compilerflags', 'host']])
            configid = cur.lastrowid
        elif len(rows) == 1:
            configid = rows[0]['rowid']
        else:
            raise InconsistentDBError('Most than one row in table config matches current configuration.')

        return f(conf, *args, **kwds)
    return wrapper


def reset_cache():
    '''Reset database connection and cached version information'''
    global con, configdict, configid
    con = None
    configdict = {}
    configid = None


def set_configdict(conf):
    '''Collect configuration information from different places.'''
    configdict['marx'] = marx_version(conf)
    ciao, caldb = ciaocaldb_version(conf)
    configdict['ciao'] = ciao
    configdict['caldb'] = caldb
    configdict['saotrace'] = saotrace_version(conf)
    for n in ['compiler', 'compilerversion', 'compilerflags']:
        try:
            configdict[n] = conf.get('marx', n)
        except NoOptionError:
            configdict[n] = 'unknown'
    configdict['host'] = get_host()


def get_host():
    '''Return host name or "unkown".'''
    return os.environ.get('HOST', 'unknown')


def marx_version(conf):
    '''Return the MARX version

    Parameters
    ----------
    conf : `~ConfigParser.ConfigParser`

    Returns
    -------
    version : string
        MARX version
    '''
    out = subprocess.check_output([os.path.join(conf.get('marx', 'binpath'), 'marx'),
                                   '--version'], stderr=subprocess.STDOUT)
    ver = re.match(r'MARX version ([0-9.]+)', out)
    return ver.groups()[0]


def ciaocaldb_version(conf):
    '''Return the CIAO and CalDB version

    Parameters
    ----------
    conf : `~ConfigParser.ConfigParser`

    Returns
    -------
    ciao, caldb : string
        CIAO and CALDB version
    '''

    out = subprocess.check_output([conf.get('CIAO', 'setup')],
                                  stderr=subprocess.STDOUT, shell=True)
    ciao = re.search(r'CIAO ([0-9.]+)', out)
    caldb = re.search(r'CALDB[\s]+:[\s]+([0-9.]+)', out)

    return ciao.groups()[0], caldb.groups()[0]


def saotrace_version(conf):
    '''Guess SAOTrace version from installation path.

    Parameters
    ----------
    conf : `~ConfigParser.ConfigParser`

    Returns
    -------
    saotrace : string
        SAOTrace version
    '''
    out = subprocess.check_output([conf.get('SAOTrace', 'setup') +
                                   '\nwhich trace-nest'],
                                  shell=True)
    saotrace = re.search('saotrace-([0-9.]+)', out)
    return saotrace.groups()[0]


@setup
def insert_test_run(conf, testclass, testname, testversion,
                    value, sigma=None, sigma2=None):
    """Insert a test result into the results table.

    Parameters
    ----------
    conf : `~ConfigParser.ConfigParser`
        Test configuration (file paths etc.)
    testclass: string
        Name of the test class that defines a test
    testname: string
        A single test might output more than one number. So, in addition to the
        name of the test class, a "name" for the result can be set.
    testversion: string
        If the definition of a test changes so that results cannot be
        compared to previous runs any longer, the version needs to be
        incremented.
    value: string or float or int
        Results of the test
    sigma: float
        Uncertainty on ``value`` e.g. for fit results.
    sigma2: float
        In some cases a single value for ``sigma`` might not be sufficient,
        e.g. if error up and down a different. In this case, ``sigma2`` holds
        the second number.
    """

    cur = con.execute("SELECT rowid FROM expresult WHERE testclass=:testclass AND name=:name AND version=:version",
                      {'testclass': testclass, 'name': testname,
                       'version': testversion})
    expid = cur.fetchall()
    if len(expid) == 0:
        raise InconsistentDBError('No entry for {0}: {1}, version {2} in table expresult.'.format(testclass, testname, testversion))
    elif len(expid) == 1:
        expid = expid[0][0]
    else:
        raise InconsistentDBError('Several lines in expid match  {0}: {1}, version {2}'.format(testclass, testname, testversion))

    with con:
        con.execute("INSERT INTO result VALUES (?, ?, ?, ?, ?)",
                    (configid, expid, value, sigma, sigma2))


@setup
def insert_expected_result(conf, testclass, testname, testversion,
                           value, title, description=None, unit=None,
                           acceptable=None):
    """Insert a test description result into the expresults table.

    Parameters
    ----------
    conf : `~ConfigParser.ConfigParser`
        Test configuration (file paths etc.)
    testclass: string
        Name of the test class that defines a test
    testname: string
        A single test might output more than one number. So, in addition to the
        name of the test class, a "name" for the result can be set.
    testversion: string
        If the definition of a test changes so that results cannot be
        compared to previous runs any longer, the version needs to be
        incremented.
    value: string or float or int
        Results of the test
    title: string
        A title used when results are displayed on websites etc.
        Limited to 20 characters.
    description: string
        Describe the test if more details are needed.
    unit: string
        Unit of value. Used for display only.
    acceptable: float
        Values are rarely exact. This field defined the maximum deviation
        from ``value`` that still makes a test pass.
    """
    cur = con.execute("SELECT rowid FROM expresult WHERE testclass=:testclass AND name=:name AND version=:version",
                      {'testclass': testclass, 'name': testname,
                       'version': testversion})
    expid = cur.fetchone()
    if expid is None:
        with con:
            cur = con.execute("INSERT INTO expresult(testclass, name, version, title, description, value, unit, acceptable) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                              (testclass, testname, testversion, title,
                               description, value, unit, acceptable))
