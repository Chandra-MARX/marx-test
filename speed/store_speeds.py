
import os
from glob import glob
import subprocess
import argparse
import json


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



def summarize_runtimes():
    '''Read all *.time files and extract the user time'''
    timenames = glob("*/*.time")
    timenames.sort()

    host = os.getenv('HOST', 'unknown host')

    times = []
    for t in timenames:
        with open(t, 'r') as f:
            time = f.readline().split('u')[0]
            times.append([host,
                          os.path.basename(t),
                          os.path.basename(os.path.dirname(t)),
                          time])
    return times


def get_immidiate_subdirs():
    return next(os.walk('.'))[1]


def checkout_hash(directory, githash):
    with ChangeDir(os.path.join(directory, 'dist')):
        subprocess.call(['git', 'fetch'])
        subprocess.call(['git', 'checkout', githash])
        subprocess.call(['make', '-j{0}'.format(args.j)])
        subprocess.call(['make', 'install'])


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
parser.add_argument('-j', type=int, default=1,
                    help='number of jobs that make runs simultaneously')

args = parser.parse_args()

for d in get_immidiate_subdirs():
    checkout_hash(d, args.githash)

subprocess.call(['make', '-j{0}'.format(args.j), 'subdirs'])
times = summarize_runtimes()

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
