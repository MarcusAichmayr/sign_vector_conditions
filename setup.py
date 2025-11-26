from setuptools import setup


# Get information from separate files (README)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


setup(
    name="sign_crn",
    version="2.2",
    description="SageMath package for (chemical) reaction networks using sign vector conditions",
    long_description=readfile("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/MarcusAichmayr/sign_crn",
    author="Marcus S. Aichmayr",
    author_email="aichmayr@mathematik.uni-kassel.de",
    license="GPL-3.0-or-later",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],  # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords=[
        "reaction networks",
        "crn",
        "equilibrium",
        "generalized mass-action",
        "robustness",
        "oriented matroids",
        "sign vector conditions",
    ],
    packages=["sign_crn"],
    install_requires=["elementary_vectors>=2.2", "sign_vectors>=1.2", "certlin>=1.2"],
    extras_require={
        "passagemath": [
            "passagemath-symbolics",
            "passagemath-flint",
            "passagemath-graphs",
            "passagemath-repl",
        ],
    },
)
