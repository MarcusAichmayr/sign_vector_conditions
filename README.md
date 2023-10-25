# Sign vector conditions for chemical reaction networks

## Description

a SageMath package to work with sign vector conditions for chemical reaction networks

This package is still in development. Please inform me if something doesn't work as expected.

## License

Distributed under the terms of the GNU General Public License (GPL, see the
LICENSE file), either version 3 or (at your option) any later version

- http://www.gnu.org/licenses/

## Requirements

Sage 9.0 or later is recommended. Some features should work with older versions.

The package [elementary_vectors](https://github.com/MarcusAichmayr/elementary_vectors) is necessary for this package to work.

## Installation

### Install from GitHub (recommended)

To download and install the latest development version on a system where Sage
was built from source or installed from official packages, run

    $ sage -pip install git+https://github.com/MarcusAichmayr/sign_vector_conditions.git

or

    $ sage -pip install --user git+https://github.com/MarcusAichmayr/sign_vector_conditions.git

The optional `--user` flag causes the package to be installed in your `.sage` directory instead of the Sage installation tree.

### Local install from source

Download the source from the git repository:

    $ git clone https://github.com/MarcusAichmayr/sign_vector_conditions.git

Change to the root directory of the repository and run:

    $ sage -pip install --upgrade --no-index -v .

You can also run instead the shorthand:

    $ make install

### Documentation

The documentation of this package can be found on GitHub:
https://marcusaichmayr.github.io/sign_vector_conditions/index.html

To generate the documentation of this package, run

    $ make doc

or

    $ make doc-pdf

at the root directory of the repository.

### Testing

To run the test suite, install the package and run the command

    $ make test

at the root directory of the repository.
