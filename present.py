import importlib
from os.path import join as pjoin
import os
import jinja2

from .tests import modulelist


def write_summary_rst(env):
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
    env : dict
        Dictionary with a ket ``outpath``.
    '''
    jinjaenv = jinja2.Environment(loader=jinja2.PackageLoader('marxtest', 'templates'))
    testindex = jinjaenv.get_template('testindex.rst')
    codeindex = jinjaenv.get_template('codelist.rst')

    if not os.path.exists(env['outpath']):
        os.makedirs(env['outpath'])

    # write module and code pages
    codelist = []
    for module in modulelist:
        imp_module = importlib.import_module('marxtest.tests.' + module)
        for t in imp_module.tests:
            testclass = getattr(imp_module, t)
            testinst = testclass(env)
            codelist.append(testinst.name)

    # write intro
    with open(pjoin(env['outpath'], 'index.rst'), 'w') as f:
        f.write(testindex.render(modulelist=modulelist))

    # write codelist because every page must be in some TOC tree
    with open(pjoin(env['outpath'], 'listofcode.rst'), 'w') as f:
        f.write(codeindex.render(codelist=codelist))


def run_and_output(env, run=True, modules=None, tests=None):
    '''

    Parameters
    ----------
    env : dict
        Dictionary with a key ``outpath``.
    run : boolean
        If ``True`` tests will be executed, otherwise test results
        already on disk from a previous run will be parsed to generate
        the rst output.
    modules : list of strings or None
        If set, this defines the list of test modules to run and write.
        If ``None``, use ``tests.modulelist``.
    tests : list of strings or None
        If set, only tests with a name listed in ``tests`` will be run and
        generate rst output. Can be used in combination with ``modules``.
    '''
    jinjaenv = jinja2.Environment(loader=jinja2.PackageLoader('marxtest', 'templates'))
    testlisttemp = jinjaenv.get_template('testlist.rst')
    codetemp = jinjaenv.get_template('testcode.rst')

    if not os.path.exists(pjoin(env['outpath'], 'code')):
        os.makedirs(pjoin(env['outpath'], 'code'))

    # write module and code pages
    if modules is None:
        modules = modulelist
    for module in modules:
        imp_module = importlib.import_module('marxtest.tests.' + module)
        t_list = []
        for t in imp_module.tests:
            # If tests is set, only run tests mentioned in that list
            if (tests is not None) and (t not in tests):
                continue
            testclass = getattr(imp_module, t)
            testinst = testclass(env)
            if run:
                testinst.run()
            t_list.append(testinst)
            with open(pjoin(env['outpath'], 'code', t + '.rst'), 'w') as f:
                f.write(codetemp.render(testinst=testinst))
        with open(pjoin(env['outpath'], module + '.rst'), 'w') as f:
            f.write(testlisttemp.render(module=imp_module,
                                        testinstances=t_list,
                                        figpath='figures')
                    )
