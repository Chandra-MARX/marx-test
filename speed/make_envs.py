# python of makefile or shell script?
import os
import subprocess
import argparse
import json

from ..utils import ChangeDir

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Create environments for MARX testing\n

Each environment is a combination of compiler and flags.
They are set up in separate directories.
''')
    parser.add_argument('envs', type=str,
                        help='filename with environment specifications\n' +
                        'The format is a json list of lists,' +
                        'where each list has at least' +
                        'three entries: name, CC, CFLAGS\n' +
                        'Example:\n' +
                        '[["gccdefault",  "gcc", "-g -O2"]]')
    parser.add_argument('gitpath', type=str,
                        help='git repository (local path or url)')

    args = parser.parse_args()

    with open(args.envs, 'r') as f:
        envs = json.load(f)

    for env in envs:
        os.mkdir(env[0])
        abspath = os.path.abspath(env[0])

        with ChangeDir(env[0]):
            subprocess.call(['git', 'clone', os.path.join('..', args.gitpath)])
        with ChangeDir(os.path.join(env[0], 'dist')):
            os.environ['CC'] = env[1]
            os.environ['CFLAGS'] = env[2]
            installpath = os.path.join(abspath, 'install')
            subprocess.call(['./configure',
                             '--prefix={0}'.format(installpath)])
