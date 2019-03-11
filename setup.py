from setuptools import find_packages, setup, Command
import stat
import os
from shutil import rmtree
import sys

# Package meta-data.
NAME = 'sdhc'
DESCRIPTION = 'Tools for computing spectral heat current distribution using LAMMPS simulations'
URL = 'https://github.com/ksaaskil/shc-python-tools'
EMAIL = 'ksaaskil@gmail.com'
AUTHOR = 'Kimmo S\"a\"askilahti'
REQUIRES_PYTHON = '>=2.7.8'
SRC_DIR = 'sdhc'  # Relative location wrt setup.py

# Required packages.
REQUIRED = ['numpy']

DEV = ['pytest', 'sphinx', 'sphinx_rtd_theme']

EXTRAS = {'dev': DEV }

# Entry point for CLI (relative to setup.py)
ENTRY_POINTS = ['sdhc=sdhc.__main__:main']

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
with open(os.path.join(here, 'README.md')) as f:
    long_description = '\n' + f.read()

# Load the package's __version__.py module as a dictionary.
about = dict()
with open(os.path.join(here, SRC_DIR, '__version__.py')) as f:
    exec(f.read(), about)

class SetupCommand(Command):
    """Base class for setup.py commands with no arguments"""
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def rmdir_if_exists(self, directory):
        self.status("Deleting {}".format(directory))
        rmtree(directory, ignore_errors=True)


class BuildDistCommand(SetupCommand):
    """Support setup.py upload."""
    description = "Build the package."

    def run(self):

        self.status("Removing previous builds...")
        self.rmdir_if_exists(os.path.join(here, 'dist'))

        self.status("Building Source and Wheel (universal) distribution...")
        os.system("{executable} setup.py sdist bdist_wheel --universal".format(executable=sys.executable))
        sys.exit()


def build_docs():
    os.chdir("docs")
    os.system("sphinx-apidoc -f -e -o source/ ../sdhc/")
    os.system("sphinx-build -M html -D version={version} source build".format(version=about['__version__']))


class BuildDocumentationCommand(SetupCommand):
    """Builds the sphinx documentation"""
    description = "Builds the sphinx documentation."

    def run(self):
        self.status("Removing previous builds...")

        build_dir = os.path.join(here, 'docs/build')
        self.rmdir_if_exists(build_dir)

        version_dir = os.path.join(here, 'docs', 'version={version}'.format(version=about['__version__']))
        self.rmdir_if_exists(version_dir)  # Need to delete this before building HTML docs

        self.status("Building documentation...")
        build_docs()

        self.status("Docs were built to `docs/build`.")
        sys.exit()


class TestCommand(SetupCommand):
    """Support setup.py test."""
    description = "Run local test if they exist"

    def run(self):
        os.system("pytest")
        sys.exit()


setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=('tests',)),
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='Apache 2.0',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Operating System :: Unix'
    ],
    entry_points={'console_scripts': ENTRY_POINTS},
    cmdclass={'dist': BuildDistCommand,
              'test': TestCommand,
              'doc': BuildDocumentationCommand}
)
