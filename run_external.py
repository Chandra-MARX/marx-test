'''
Eventually, these string definitions should be moved to a configuration file,
but for now I just define them here.
'''
external_settings = {None : '',
                     'CIAO': '''
#!/bin/sh
# unset DISPLAY
# export CIAO_MAJOR_VER=-4.7

# Use -o option to overrride here. Otherwise, this script fails with an error if a CIAO
# environment is already set up and make stops.
# This is true, evne if the previous CIAO environmet is the same version.
# Due to the >/dev/null the error message is invisible, so it's better to make sure that
# something is set.
. /nfs/cxc/a1/setup/ciao-setup.sh -o > /dev/null
"$@"
''',
                     'marxparfile': '/nfs/melkor/d1/guenther/marx/installed/dev/share/marx/pfiles/marx.par',
                     'marxaspparfile': '/nfs/melkor/d1/guenther/marx/installed/dev/share/marx/pfiles/marxasp.par',
                     'marxpath': '/nfs/melkor/d1/guenther/marx/installed/dev/bin/',

                     'SAOTrace': '''
#!/bin/sh
# unset DISPLAY
. /nfs/mkx/a1/setup/saotrace-setup.sh > /dev/null
"$@"
''',

}
