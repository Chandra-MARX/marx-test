from __future__ import print_function
import os
from glob import glob
import subprocess
import argparse
import json
import sys


class ChangeDir(object):
    """ Context manager to step into a directory temporarily.

    See: https://pythonadventures.wordpress.com/2013/12/15/chdir-a-context-manager-for-switching-working-directories/
    """
    def __init__(self, directory):
        self.old_dir = os.getcwd()
        self.new_dir = directory

    def __enter__(self):
        os.chdir(self.new_dir)

    def __exit__(self, *args):
        os.chdir(self.old_dir)


def get_immidiate_subdirs():
    return next(os.walk('.'))[1]


def checkout_hash(directory, githash):
    with open(os.devnull, 'wb', 0) as DEVNULL:
        with ChangeDir(os.path.join(directory, 'dist')):
            print('{0}: git fetch'.format(directory), end=', ')
            sys.stdout.flush()
            subprocess.call(['git', 'fetch'],
                            stdout=DEVNULL, stderr=subprocess.STDOUT)
            print('git checkout', end=', ')
            sys.stdout.flush()
            subprocess.call(['git', 'checkout', githash],
                            stdout=DEVNULL, stderr=subprocess.STDOUT)
            print('make', end=', ')
            sys.stdout.flush()
            subprocess.call(['make', '-j{0}'.format(args.j)],
                            stdout=DEVNULL, stderr=subprocess.STDOUT)
            print('make install')
            sys.stdout.flush()
            subprocess.call(['make', 'install'],
                            stdout=DEVNULL, stderr=subprocess.STDOUT)


def time_command(cmd):
    '''time a command, e.g. >>> time_command(['ls', '-l'])'''
    with open(os.devnull, 'wb', 0) as DEVNULL:
        p = subprocess.Popen(cmd, stdout=DEVNULL, stderr=subprocess.STDOUT)
        ru = os.wait4(p.pid, 0)[2]
    return ru.ru_utime + ru.ru_stime


def execute_timing_tests(parpath, hostname, githash):
    times = []
    allpar = glob(os.path.join(parpath, '*.par'))
    allpar = [os.path.abspath(a) for a in allpar]
    for env in get_immidiate_subdirs():
        print('Running timing tests in {0}'.format(env), end=' ')
        with ChangeDir(env):
            for par in allpar:
                runtime = time_command(['install/bin/marx',
                                        '@@{0}'.format(par)])
                times.append([hostname, githash,
                              env, os.path.basename(par)[:-4], runtime])
                print('.',end="")
                sys.stdout.flush()
        print('')
    return times


parser = argparse.ArgumentParser(description='''Execute speed test for certain git commit.

This script needs to be run in a directory that contains one subdirectory
for each environment (combination of compiler and compilation options).
Each subdirectory in turn contains a folder "dist" which holds the git
repro of marx. This script assumes that "./configure" has been run on
those repros to install marx in a directory "install", so that the
structure will be:
env1/
env1/dist
env1/install
env2/
env2/dist
env2/install
...

Furthermore, the root directory must have a "Makefile" that contains a
target "subdirs" to run the speed tests.
''')
parser.add_argument('githash', type=str,
                    help='hash for the git commit')
parser.add_argument('file', type=str,
                    help='filename to save results')
parser.add_argument('parpath', type=str,
                    help='path to marx .par files')
parser.add_argument('-j', type=int, default=1,
                    help='number of jobs that make runs simultaneously')

args = parser.parse_args()

for d in get_immidiate_subdirs():
    checkout_hash(d, args.githash)

host = os.getenv('HOST', 'unknown host')
times = execute_timing_tests(args.parpath, args.githash, host)


# read jsonfile if it exists
try:
    with open(args.file, 'r') as jsonfile:
        timetabs = json.load(jsonfile)
except IOError:
    # file not present yet
    timetabs = []

# add new data
timetabs.extend(times)

# save file
with open(args.file, 'w') as jsonfile:
    json.dump(timetabs, jsonfile)
