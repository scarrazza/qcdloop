from distutils.core import setup
from Cython.Build import cythonize
import subprocess
import sys

def call_command(command):
    l = command.split()
    try:
        result = subprocess.check_output(l)
    except OSError as e:
        print("Could not call %s: %s.\n"
              "Please make sure the relevant command is installed."
              % (l[0], e.strerror) )
        sys.exit(1)
    return result.decode().rstrip()

setup(
    name = 'qcdloop',
    author = 'Stefano Carrazza et al.',
    author_email = 'stefano.carrazza@cern.ch',
    version = call_command('qcdloop-config --version'),
    ext_modules = cythonize('*.pyx'),
)
