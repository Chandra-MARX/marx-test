import importlib
from ConfigParser import ConfigParser
from os.path import join as pjoin
import os
import sys
import jinja2

from .utils import ChangeDir


__all__ = ['write_summary_rst', 'run_and_output']


def setup_env(conffile):
    '''
    Parameters
    ----------
    conffile : string
        Path and name of configuration file
    '''
    conf = ConfigParser()
    conf.read(conffile)
    outpath = conf.get('Output', 'outpath')

    jinjaenv = jinja2.Environment(loader=jinja2.FileSystemLoader(conf.get('Output', 'templates')))

    # Get tests fomr test module
    sys.path.append(conf.get('tests', 'path'))
    main_module = importlib.import_module(conf.get('tests', 'modulename'))

    return conf, outpath, jinjaenv, main_module.modulelist


def write_summary_rst(conffile):
    '''Write TOC pages.

    This function imports (but does not run) all test modules,
    looks for the tests included in those files
    and writes two rst pages with a TOC:

    - ``testindex.rst`` is a list of all modules (each of which could
      contain multiple tests).
    - ``codelist.rst`` is a list of all code pages (each test outputs
      an rst page with the code used to run this test). Those code pages
      are linked from the test results, but sphinx complains if they are not
      part of any TOC. That's what ``codelist.rst`` is for.

    If tests are run individually (e.g. when writing and debugging tests), then
    this function is useful to rebuild the cover pages without rerunning all
    tests (which may take a long time).

    Parameters
    ----------
    conffile : string
        Path and name of configuration file
    '''
    conf, outpath, jinjaenv, modulelist = setup_env(conffile)

    testindex = jinjaenv.get_template('testindex.rst')
    codeindex = jinjaenv.get_template('codelist.rst')

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # write module and code pages
    codelist = []
    for module in modulelist:
        imp_module = importlib.import_module(conf.get('tests', 'modulename') + '.' + module)
        for t in imp_module.tests:
            testclass = getattr(imp_module, t)
            testinst = testclass(conf)
            codelist.append(testinst.name)

    # write intro
    with open(pjoin(outpath, 'index.rst'), 'w') as f:
        f.write(testindex.render(modulelist=modulelist))

    # write codelist because every page must be in some TOC tree
    with open(pjoin(outpath, 'listofcode.rst'), 'w') as f:
        f.write(codeindex.render(codelist=codelist))


def run_and_output(conffile, run=True, modules=None, tests=None):
    '''Run tests and output results.

    Parameters
    ----------
    conffile : string
        Path and name of configuration file.
    run : boolean
        If ``True`` tests will be executed, otherwise test results
        already on disk from a previous run will be parsed to generate
        the rst output. This option is useful, if just the templates for
        the rst files have changed. It may not work for all tests.
    modules : list of strings or None
        If set, this defines the list of test modules to run and write.
        If ``None``, use ``tests.modulelist``.
    tests : list of strings or None
        If set, only tests with a name listed in ``tests`` will be run and
        generate rst output. Can be used in combination with ``modules``.
    '''
    conf, outpath, jinjaenv, modulelist = setup_env(conffile)
    testlisttemp = jinjaenv.get_template('testlist.rst')
    codetemp = jinjaenv.get_template('testcode.rst')

    if not os.path.exists(pjoin(outpath, 'code')):
        os.makedirs(pjoin(outpath, 'code'))

    # write module and code pages
    if modules is None:
        modules = modulelist
    for module in modules:
        imp_module = importlib.import_module('tests.' + module)
        t_list = []
        for t in imp_module.tests:
            # If tests is set, only run tests mentioned in that list
            if (tests is not None) and (t not in tests):
                continue
            testclass = getattr(imp_module, t)
            testinst = testclass(conf)
            if run:
                testinst.run()
            t_list.append(testinst)
            with open(pjoin(outpath, 'code', t + '.rst'), 'w') as f:
                with ChangeDir(testinst.basepath):
                    f.write(codetemp.render(testinst=testinst))
        with open(pjoin(outpath, module + '.rst'), 'w') as f:
            f.write(testlisttemp.render(module=imp_module,
                                        testinstances=t_list,
                                        figpath='figures')
                    )
