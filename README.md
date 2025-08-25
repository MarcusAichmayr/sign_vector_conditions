# Sign vector conditions for chemical reaction networks

## Description

a SageMath package to work with sign vector conditions for chemical reaction networks

## License

Distributed under the terms of the GNU General Public License (GPL, see the
LICENSE file), either version 3 or (at your option) any later version

- http://www.gnu.org/licenses/

## Requirements

Sage 10.0 or later is recommended.

The package [elementary_vectors](https://github.com/MarcusAichmayr/elementary_vectors) is necessary for this package to work.

## Installation

### Install from GitHub (recommended)

To download and install the latest development version on a system where Sage
was built from source or installed from official packages, run

    sage -pip install git+https://github.com/MarcusAichmayr/sign_vector_conditions.git

### Install from PyPI (recommended)

To install the package for Python directly, run:

    pip install sign-vector-conditions

### Local install from source

Download the source from the git repository:

    git clone https://github.com/MarcusAichmayr/sign_vector_conditions.git

Change to the root directory of the repository and run:

    make install

### Local install from source (no Sage installation required)

Download the source from the git repository:

    git clone https://github.com/MarcusAichmayr/sign_vector_conditions.git

Change to the root directory of the repository and run:

    python3 -m venv venv
    . venv/bin/activate
    pip install -v -e ".[passagemath]"

## Documentation

The documentation of this package is available on GitHub:

https://marcusaichmayr.github.io/sign_vector_conditions/index.html

To generate it, run

    make doc

or

    make doc-pdf

at the root directory of the repository.

## Testing

To run the test suite, install the package and run the command

    make test

at the root directory of the repository.
