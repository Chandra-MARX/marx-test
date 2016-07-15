import importlib
from os.path import join as pjoin
import os
import jinja2

from . import tests


def run_and_writerst(path, env, run=True):
    jinjaenv = jinja2.Environment(loader=jinja2.PackageLoader('marxtest', 'templates'))
    testindex = jinjaenv.get_template('testindex.rst')
    codeindex = jinjaenv.get_template('codelist.rst')
    testlisttemp = jinjaenv.get_template('testlist.rst')
    codetemp = jinjaenv.get_template('testcode.rst')

    if not os.path.exists(pjoin(env['outpath'], 'code')):
        os.makedirs(pjoin(env['outpath'], 'code'))

    # write module and code pages
    codelist = []
    for module in tests.modulelist:
        imp_module = importlib.import_module('marxtest.tests.' + module)
        t_list = []
        for t in imp_module.tests:
            testclass = getattr(imp_module, t)
            testinst = testclass(env)
            if run:
                testinst.run()
            t_list.append(testinst)
            codelist.append(testinst.name)
            with open(pjoin(env['outpath'], 'code', t + '.rst'), 'w') as f:
                f.write(codetemp.render(testinst=testinst))
        with open(pjoin(env['outpath'], module + '.rst'), 'w') as f:
            f.write(testlisttemp.render(module=imp_module,
                                        testinstances=t_list,
                                        figpath='figures')
                    )
    # write intro
    with open(pjoin(env['outpath'], 'index.rst'), 'w') as f:
        f.write(testindex.render(modulelist=tests.modulelist))

    # write codelist because every page must be in some TOC tree
    with open(pjoin(env['outpath'], 'listofcode.rst'), 'w') as f:
        f.write(codeindex.render(codelist=codelist))
