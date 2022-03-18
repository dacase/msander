# Note: always find header files in $MSANDERHOME/include

from distutils.core import setup, Extension
import os
from os.path import join
import sys


# Update/correct sys env vars to be correct for OS X.
if sys.platform == 'darwin':
    os.environ['CXX'] = 'clang++'
    os.environ['CC'] = 'clang'
# PGI does not work with Python, it appears. So force gcc, which is almost
# universally available
elif os.environ.get('CC', '').endswith('pgcc'):
    os.environ['CXX'] = 'g++'
    os.environ['CC'] = 'gcc'

sander_extension_source = 'sander/src/pysandermodule.c'
sanderles_extension_source = 'sanderles/src/pysandermodule.c'

# Replace any existing sanderles with sander if needed
if (not os.path.isfile(sanderles_extension_source)) or os.path.getmtime(sander_extension_source) > os.path.getmtime(sanderles_extension_source):
    print('Copying sander extension source to sanderles...')
    if sys.platform == 'win32':
        os.system('if exist sanderles rmdir /S /Q sanderles')
        os.system('xcopy /E /C /I sander sanderles')

    else:
        os.system('rm -fr sanderles')
        os.system('cp -r sander sanderles')

# Retrieve value of AMBERHOME env variable. (None if not defined).
amberhome = os.getenv('MSANDERHOME')
print('using amberhome as %s' % amberhome)

definitions = []
lesdefinitions = [('LES', None)]

# create arrays of dependencies
depends = ['sander/src/pysandermoduletypes.c']
lesdepends = ['sanderles/src/pysandermoduletypes.c']

# MSANDERHOME supplied.  Define new directory tree rooted at MSANDERHOME.
incdir = [join(amberhome, 'include')]
libdir = [join(amberhome, 'lib')]
depends.append(join(incdir[0], 'CompatibilityMacros.h'))
lesdepends.append(join(incdir[0], 'CompatibilityMacros.h'))


packages = ['sander', 'sanderles']

pysander = Extension('sander.pysander',
                     sources=[sander_extension_source],
                     include_dirs=incdir, library_dirs=libdir,
                     libraries=['sander'],
                     depends=depends,
                     define_macros=definitions
)

pysanderles = Extension('sanderles.pysander',
                        sources=[sanderles_extension_source],
                        include_dirs=incdir, library_dirs=libdir,
                        libraries=['sanderles'],
                        depends=lesdepends,
                        define_macros=lesdefinitions)

setup(name='sander',
      version="20.0",
      license='GPL v2 or later',
      author='Jason Swails',
      author_email='jason.swails -at- gmail.com',
      description='SANDER energy/force evaluation',
      ext_modules=[pysander, pysanderles],
      packages=packages,
)
