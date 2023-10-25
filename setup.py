from setuptools import setup

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding='utf-8') as f:
        return f.read()

setup(
    name = 'sign_vector_conditions',
    version = readfile("VERSION").strip(), # the VERSION file is shared with the documentation
    description = 'a SageMath package to work with sign vector conditions for chemical reaction networks',
    long_description = readfile("README.md"),
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/MarcusAichmayr/sign_vector_conditions',
    author = 'Marcus Aichmayr',
    author_email = 'aichmayr@mathematik.uni-kassel.de',
    license = 'GPLv3',
    classifiers = [
      # How mature is this project? Common values are
      #   3 - Alpha
      #   4 - Beta
      #   5 - Production/Stable
      'Development Status :: 5 - Stable',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: GNU General Public License v3',
      'Programming Language :: Python :: 3.8.5',
    ], # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords = ['chemical reaction networks', 'sign vector conditions, oriented matroids'],
    packages = ['sign_vector_conditions'],
    setup_requires   = ['elementary_vectors'],
    install_requires   = ['elementary_vectors', 'sphinx'],
)
